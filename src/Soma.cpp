#include "plugin.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>

// SOMA — a bursting/chaotic neuron oscillator on the Hindmarsh-Rose (HR) system.
//
// The sibling of Axon: HR adds a third, slow state variable `z` (adaptation) to
// the two-variable FitzHugh-Nagumo core, which is exactly what produces bursts
// (trains of spikes separated by quiescence) and, in a window of current, chaos.
// `x` (membrane potential) is the audio output; pitch is the *speed* the
// dimensionless simulation runs at, calibrated so the within-burst spike rate is
// C4 at 0 V (open-loop — the period is emergent, so CURRENT/BURST/ADAPT pull it).
//
//   dx/dt = y - a·x³ + b·x² - z + I      (fast: membrane potential, audio out)
//   dy/dt = c - d·x² - y                 (fast recovery / spiking)
//   dz/dt = r·( s·(x - x_R) - z )         (slow adaptation / bursting)
//
// a,b,c,d,x_R are fixed at the standard values; CURRENT=I, BURST=r, ADAPT=s are
// the exposed bifurcation controls. The integrator (RK4 + pitch-adaptive
// substepping) is the same engine as Axon, with one extra equation — the
// factoring the Axon plan called for.

static const int TRAIL = 512;   // length of the (x,z) phase-portrait trail

struct Soma : Module {

    enum ParamId {
        PITCH_PARAM,
        CURRENT_PARAM,
        BURST_PARAM,
        ADAPT_PARAM,
        CURRENT_ATT_PARAM,
        BURST_ATT_PARAM,
        PARAMS_LEN
    };

    enum InputId {
        VOCT_INPUT,
        CURRENT_INPUT,
        BURST_INPUT,
        TRIG_INPUT,
        INPUTS_LEN
    };

    enum OutputId {
        OUT_OUTPUT,
        SPIKE_OUTPUT,
        Z_OUTPUT,
        OUTPUTS_LEN
    };

    enum LightId {
        LIGHTS_LEN
    };

    // ─── Persistent state ───────────────────────────────────────────────────
    float xx = -1.6f, yy = -11.8f, zz = 2.f;   // HR state, seeded near the rest fixed point
    float trigPulse = 0.f;
    dsp::SchmittTrigger trigIn;
    dsp::SchmittTrigger spikeDet;
    dsp::PulseGenerator spikeGen;
    dsp::TRCFilter<float> dcBlock;
    float lastFs = 0.f;
    float trigDecay = 0.f;            // cached per-sample TRIG pulse decay (fn of fs)
    float antiDenorm = 1e-18f;        // alternating sub-LSB dither, anti-denormal on dcBlock

    // ─── Display trail (x,z phase portrait) ─────────────────────────────────
    // Double-buffered, lock-free snapshot (see Axon): the audio thread fills the
    // back buffer and flips dispBuf with a release store; the UI reads the front
    // buffer after an acquire load, with the head index carried alongside so the
    // trail and head dot stay coherent.
    float trailX[TRAIL] = {}, trailZ[TRAIL] = {};
    int   trailIdx = 0, trailDecim = 0;
    float dispX[2][TRAIL] = {}, dispZ[2][TRAIL] = {};
    int   dispHead[2] = {};
    std::atomic<int> dispBuf{0};
    int   dispClock = 0;

    // ─── Fixed HR parameters ────────────────────────────────────────────────
    static constexpr float A = 1.f, B = 3.f, C = 1.f, D = 5.f, XR = -1.6f;

    // ─── Tunable constants (RATE_CAL from tools/soma_stability_test.cpp) ─────
    static constexpr float RATE_CAL    = 14.925501f; // within-burst spike period at default → C4 at 0 V
    static constexpr float HSUB_MAX    = 0.05f;
    static constexpr int   MIN_SUB     = 4;
    static constexpr int   MAX_SUB     = 64;
    static constexpr float TRIG_AMP    = 1.0f;        // injected current pulse (kicks a burst from rest)
    static constexpr float TRIG_TAU_MS = 5.f;
    static constexpr float SPIKE_THRESH = 1.0f;
    static constexpr float OUT_GAIN    = 1.5f;
    static constexpr float Z_CENTER    = 2.0f;        // z midpoint, so Z_OUT is bipolar around a burst
    static constexpr float Z_GAIN      = 2.5f;
    static constexpr float CV_DEPTH    = 0.2f;        // ±5V CV × att → ±1.0 (linear I; log₂ for r)
    static constexpr float STATE_MAX   = 25.f;        // backstop; normal |y| peaks ~12, max observed ~21
    static constexpr float R_MIN       = 0.001f, R_MAX = 0.05f;

    Soma() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);

        configParam(PITCH_PARAM, -4.f, 4.f, 0.f, "Pitch", " Hz", 2.f, dsp::FREQ_C4);
        configParam(CURRENT_PARAM, 0.4f, 4.0f, 2.0f, "Current (excitability)");
        // BURST stores log₂(r); display base 2 shows the adaptation rate r directly.
        configParam(BURST_PARAM, std::log2(R_MIN), std::log2(R_MAX), std::log2(0.006f),
                    "Burst / adaptation rate r", "", 2.f);
        configParam(ADAPT_PARAM, 1.0f, 5.0f, 4.0f, "Adaptation strength s");
        configParam(CURRENT_ATT_PARAM, -1.f, 1.f, 0.f, "Current CV");
        configParam(BURST_ATT_PARAM, -1.f, 1.f, 0.f, "Burst CV");

        configInput(VOCT_INPUT, "1V/oct pitch");
        configInput(CURRENT_INPUT, "Current CV");
        configInput(BURST_INPUT, "Burst rate CV");
        configInput(TRIG_INPUT, "Trigger (fire a burst)");

        configOutput(OUT_OUTPUT, "Audio (membrane potential x)");
        configOutput(SPIKE_OUTPUT, "Spike (trigger on each spike)");
        configOutput(Z_OUTPUT, "Adaptation z (burst-envelope CV)");
    }

    void onReset() override {
        xx = -1.6f; yy = -11.8f; zz = 2.f;
        trigPulse = 0.f;
        trigIn.reset();
        spikeDet.reset();
        spikeGen.reset();
        dcBlock.reset();
    }

    // HR derivatives — the FHN f() with one extra (slow) line, as the Axon plan
    // anticipated. Same RK4 step below, now over three variables.
    static inline void f(float x, float y, float z, float Itot, float r, float s,
                         float& dx, float& dy, float& dz) {
        dx = y - A*x*x*x + B*x*x - z + Itot;
        dy = C - D*x*x - y;
        dz = r * (s * (x - XR) - z);
    }

    static inline void rk4(float& x, float& y, float& z, float h,
                           float Itot, float r, float s) {
        float ax,ay,az, bx,by,bz, cx,cy,cz, dx,dy,dz;
        f(x,           y,           z,           Itot, r, s, ax,ay,az);
        f(x+.5f*h*ax,  y+.5f*h*ay,  z+.5f*h*az,  Itot, r, s, bx,by,bz);
        f(x+.5f*h*bx,  y+.5f*h*by,  z+.5f*h*bz,  Itot, r, s, cx,cy,cz);
        f(x+h*cx,      y+h*cy,      z+h*cz,      Itot, r, s, dx,dy,dz);
        x += h/6.f*(ax + 2*bx + 2*cx + dx);
        y += h/6.f*(ay + 2*by + 2*cy + dy);
        z += h/6.f*(az + 2*bz + 2*cz + dz);
    }

    void process(const ProcessArgs& args) override {
        const float fs = args.sampleRate;

        // SR-derived coefficients, recomputed only on a sample-rate change (no
        // onSampleRateChange handler needed). Caching trigDecay keeps a
        // transcendental off the per-sample path.
        if (fs != lastFs) {
            dcBlock.setCutoffFreq(20.f / fs);
            trigDecay = std::exp(-1.f / (TRIG_TAU_MS * 1e-3f * fs));
            lastFs = fs;
        }

        // ── params → physics ──
        float I = params[CURRENT_PARAM].getValue()
                + inputs[CURRENT_INPUT].getVoltage() * params[CURRENT_ATT_PARAM].getValue() * CV_DEPTH;
        // BURST is log₂(r); CV adds in the log domain (a multiplicative nudge on r).
        float rLog = params[BURST_PARAM].getValue()
                   + inputs[BURST_INPUT].getVoltage() * params[BURST_ATT_PARAM].getValue() * CV_DEPTH;
        float r = clamp(std::exp2(rLog), R_MIN, R_MAX);
        float s = params[ADAPT_PARAM].getValue();

        // ── excitable trigger: rising edge injects a decaying current pulse ──
        if (trigIn.process(inputs[TRIG_INPUT].getVoltage(), 0.1f, 1.f))
            trigPulse = TRIG_AMP;
        float Itot = I + trigPulse;
        // cached decay coefficient; flush to zero once negligible (anti-denormal).
        trigPulse *= trigDecay;
        if (std::fabs(trigPulse) < 1e-30f) trigPulse = 0.f;

        // ── pitch = simulation speed (open-loop; tracks within-burst spike rate) ──
        float pitchHz = dsp::FREQ_C4 * std::exp2(
                            params[PITCH_PARAM].getValue() + inputs[VOCT_INPUT].getVoltage());
        float dtau = RATE_CAL * pitchHz / fs;
        int   K    = clamp((int) std::ceil(dtau / HSUB_MAX), MIN_SUB, MAX_SUB);
        float h    = dtau / K;

        // ── integrate K RK4 substeps ──
        for (int k = 0; k < K; k++)
            rk4(xx, yy, zz, h, Itot, r, s);

        // ── backstop ──
        if (!std::isfinite(xx) || !std::isfinite(yy) || !std::isfinite(zz)) {
            xx = -1.6f; yy = -11.8f; zz = 2.f;
        }
        xx = clamp(xx, -STATE_MAX, STATE_MAX);
        yy = clamp(yy, -STATE_MAX, STATE_MAX);
        zz = clamp(zz, -STATE_MAX, STATE_MAX);

        // ── outputs ──
        // Alternating sub-LSB dither keeps the DC-blocker's recursive state off
        // denormals when x is parked at rest (sub-threshold, no triggers);
        // ~1e-18 V is inaudible but well above the denormal floor.
        antiDenorm = -antiDenorm;
        dcBlock.process(xx + antiDenorm);
        outputs[OUT_OUTPUT].setVoltage(5.f * std::tanh(dcBlock.highpass() * OUT_GAIN));

        if (spikeDet.process(xx, 0.f, SPIKE_THRESH))
            spikeGen.trigger(1e-3f);
        outputs[SPIKE_OUTPUT].setVoltage(spikeGen.process(args.sampleTime) ? 10.f : 0.f);

        // Z is the slow adaptation variable — a burst-envelope CV (not high-passed).
        outputs[Z_OUTPUT].setVoltage(clamp((zz - Z_CENTER) * Z_GAIN, -5.f, 5.f));

        // ── feed the (x,z) phase-portrait trail (decimated; bursts are slow) ──
        if (++trailDecim >= 6) {
            trailDecim = 0;
            trailX[trailIdx] = xx;
            trailZ[trailIdx] = zz;
            trailIdx = (trailIdx + 1) % TRAIL;
        }
        if (++dispClock >= (int)(fs / 45.f)) {       // publish display snapshot ~45 Hz
            dispClock = 0;
            int next = 1 - dispBuf.load(std::memory_order_relaxed);
            std::copy(trailX, trailX + TRAIL, dispX[next]);
            std::copy(trailZ, trailZ + TRAIL, dispZ[next]);
            dispHead[next] = trailIdx;
            dispBuf.store(next, std::memory_order_release);
        }
    }
};


// ─── Phase-portrait display ────────────────────────────────────────────────
// The Hindmarsh-Rose attractor projected onto the (x, z) plane: fast spikes
// sweep horizontally while the slow adaptation z drifts up and down, so a burst
// shows as a cluster of spikes climbing along z and the quiescent gap as the
// slow return. In the chaotic window the trail never quite repeats. The faint
// diagonal is the z-nullcline z = s·(x − x_R) (where the slow drift reverses).
struct SomaDisplay : Widget {
    Soma* module = nullptr;
    std::shared_ptr<Font> font;

    static constexpr float XMIN = -2.0f, XMAX = 2.2f;
    static constexpr float ZMIN = 0.3f, ZMAX = 4.0f;

    void draw(const DrawArgs& args) override {
        nvgBeginPath(args.vg);
        nvgRoundedRect(args.vg, 0, 0, box.size.x, box.size.y, mm2px(2.f));
        nvgFillColor(args.vg, nvgRGB(0x07, 0x07, 0x12));
        nvgFill(args.vg);
        nvgStrokeColor(args.vg, nvgRGB(0x2b, 0x2b, 0x4d));
        nvgStrokeWidth(args.vg, 1.2f);
        nvgStroke(args.vg);
        Widget::draw(args);
    }

    void drawLayer(const DrawArgs& args, int layer) override {
        if (layer == 1) {
            nvgScissor(args.vg, 0.f, 0.f, box.size.x, box.size.y);

            auto X = [&](float x) { return (x - XMIN) / (XMAX - XMIN) * box.size.x; };
            auto Y = [&](float z) { return box.size.y - (z - ZMIN) / (ZMAX - ZMIN) * box.size.y; };

            float s = module ? module->params[Soma::ADAPT_PARAM].getValue() : 4.f;

            // z-nullcline: z = s·(x − x_R).
            nvgBeginPath(args.vg);
            nvgMoveTo(args.vg, X(XMIN), Y(s * (XMIN - Soma::XR)));
            nvgLineTo(args.vg, X(XMAX), Y(s * (XMAX - Soma::XR)));
            nvgStrokeColor(args.vg, nvgRGBA(0x90, 0x50, 0x70, 0x55));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);

            if (module) {
                int b = module->dispBuf.load(std::memory_order_acquire);
                const float* dx = module->dispX[b];
                const float* dz = module->dispZ[b];
                int idx = module->dispHead[b];   // coherent with dx/dz
                nvgLineCap(args.vg, NVG_ROUND);
                for (int k = 1; k < TRAIL; k++) {
                    int i0 = (idx + k - 1) % TRAIL;
                    int i1 = (idx + k) % TRAIL;
                    float alpha = (float) k / TRAIL;
                    nvgBeginPath(args.vg);
                    nvgMoveTo(args.vg, X(dx[i0]), Y(dz[i0]));
                    nvgLineTo(args.vg, X(dx[i1]), Y(dz[i1]));
                    nvgStrokeColor(args.vg, nvgRGBA(0xff, 0x9b, 0x3a, (int)(alpha * 0xcc)));
                    nvgStrokeWidth(args.vg, 1.6f);
                    nvgStroke(args.vg);
                }
                int newest = (idx + TRAIL - 1) % TRAIL;
                float hx = X(dx[newest]), hy = Y(dz[newest]);
                nvgBeginPath(args.vg); nvgCircle(args.vg, hx, hy, 4.f);
                nvgFillColor(args.vg, nvgRGBA(0xff, 0xc8, 0x88, 0x55)); nvgFill(args.vg);
                nvgBeginPath(args.vg); nvgCircle(args.vg, hx, hy, 2.f);
                nvgFillColor(args.vg, nvgRGB(0xff, 0xe0, 0xc0)); nvgFill(args.vg);
            } else {
                // Browser preview: an illustrative bursting loop.
                nvgBeginPath(args.vg);
                for (int i = 0; i <= 80; i++) {
                    float t = 2.f * (float)M_PI * i / 80.f;
                    float px = X(0.3f + 1.2f * std::cos(t) + 0.3f * std::sin(6*t));
                    float py = Y(2.0f + 1.1f * std::sin(t));
                    if (i == 0) nvgMoveTo(args.vg, px, py); else nvgLineTo(args.vg, px, py);
                }
                nvgStrokeColor(args.vg, nvgRGBA(0xff, 0x9b, 0x3a, 0x99));
                nvgStrokeWidth(args.vg, 1.6f);
                nvgStroke(args.vg);
            }

            if (!font)
                font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
            if (font) {
                nvgFontFaceId(args.vg, font->handle);
                nvgFontSize(args.vg, mm2px(3.2f));
                nvgTextLetterSpacing(args.vg, mm2px(0.5f));
                nvgFillColor(args.vg, nvgRGBA(0xff, 0xc0, 0x9a, 0xcc));
                nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
                nvgText(args.vg, mm2px(2.6f), mm2px(2.2f), "SOMA", NULL);
                nvgTextLetterSpacing(args.vg, 0.f);
                nvgFontSize(args.vg, mm2px(2.2f));
                nvgFillColor(args.vg, nvgRGBA(0xb0, 0x80, 0x60, 0xaa));
                nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_BOTTOM);
                nvgText(args.vg, box.size.x - mm2px(2.4f), box.size.y - mm2px(1.8f), "x–z", NULL);
            }

            nvgResetScissor(args.vg);
        }
        Widget::drawLayer(args, layer);
    }
};


// ─── Widget ─────────────────────────────────────────────────────────────────

struct SomaWidget : ModuleWidget {
    std::shared_ptr<Font> font;

    void draw(const DrawArgs& args) override {
        ModuleWidget::draw(args);
        if (!font)
            font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
        if (!font) return;

        nvgSave(args.vg);
        nvgFontFaceId(args.vg, font->handle);
        nvgTextAlign(args.vg, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);

        const NVGcolor dim    = nvgRGB(0x77, 0x77, 0x99);
        const NVGcolor outclr = nvgRGB(0xcc, 0xcc, 0xee);

        auto lbl = [&](float x, float y, float sz, NVGcolor col, const char* s) {
            nvgFontSize(args.vg, mm2px(sz));
            nvgFillColor(args.vg, col);
            nvgText(args.vg, mm2px(x), mm2px(y), s, nullptr);
        };

        lbl( 9.f,   61.f, 1.8f, dim, "PITCH");
        lbl(24.32f, 61.f, 1.8f, dim, "CURRENT");
        lbl(39.64f, 61.f, 1.8f, dim, "BURST");
        lbl(54.96f, 61.f, 1.8f, dim, "ADAPT");

        lbl(24.32f, 90.f, 1.7f, dim, "I CV");
        lbl(39.64f, 90.f, 1.7f, dim, "r CV");

        lbl( 7.5f, 118.5f, 1.7f, dim,    "V/OCT");
        lbl(19.f,  118.5f, 1.7f, dim,    "TRIG");
        lbl(34.f,  118.5f, 1.8f, outclr, "OUT");
        lbl(45.5f, 118.5f, 1.7f, outclr, "SPIKE");
        lbl(56.f,  118.5f, 1.8f, outclr, "Z");

        nvgRestore(args.vg);
    }

    SomaWidget(Soma* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Soma.svg")));

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 122.0f))));

        SomaDisplay* disp = new SomaDisplay();
        disp->module = module;
        disp->box.pos  = mm2px(Vec(5.5f, 8.f));
        disp->box.size = mm2px(Vec(50.f, 34.f));
        addChild(disp);

        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec( 9.f,   54.f)), module, Soma::PITCH_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(24.32f, 54.f)), module, Soma::CURRENT_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(39.64f, 54.f)), module, Soma::BURST_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(54.96f, 54.f)), module, Soma::ADAPT_PARAM));

        addParam(createParamCentered<Trimpot>(    mm2px(Vec(24.32f, 72.f)), module, Soma::CURRENT_ATT_PARAM));
        addParam(createParamCentered<Trimpot>(    mm2px(Vec(39.64f, 72.f)), module, Soma::BURST_ATT_PARAM));
        addInput(createInputCentered<PJ301MPort>( mm2px(Vec(24.32f, 84.f)), module, Soma::CURRENT_INPUT));
        addInput(createInputCentered<PJ301MPort>( mm2px(Vec(39.64f, 84.f)), module, Soma::BURST_INPUT));

        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec( 7.5f, 112.f)), module, Soma::VOCT_INPUT));
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec(19.f,  112.f)), module, Soma::TRIG_INPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(34.f,  112.f)), module, Soma::OUT_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(45.5f, 112.f)), module, Soma::SPIKE_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(56.f,  112.f)), module, Soma::Z_OUTPUT));
    }
};

Model* modelSoma = createModel<Soma, SomaWidget>("Soma");

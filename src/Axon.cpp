#include "plugin.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>

// AXON — a spiking-neuron oscillator on the FitzHugh-Nagumo (FHN) system.
//
// The membrane voltage v is the audio output; pitch is the *speed* the
// dimensionless simulation is run at (open-loop calibrated, not phase-locked —
// the limit-cycle period is emergent, so I/eps/shape pull pitch a little; that
// is part of the instrument). Above a current threshold it free-runs as a
// relaxation oscillator; below threshold it rests and fires one spike per
// trigger.
//
// The integrator is the whole engineering job: RK4 with pitch-adaptive
// substepping so the stiff relaxation spike stays accurate/stable as pitch
// rises. f() and the RK4 step are factored so a Hindmarsh-Rose sibling (v2)
// can reuse them with one extra state variable.

// Length of the (v,w) phase-portrait trail shown on the display.
static const int TRAIL = 512;

struct Axon : Module {

    enum ParamId {
        PITCH_PARAM,
        CURRENT_PARAM,
        EPS_PARAM,
        SHAPE_PARAM,
        CURRENT_ATT_PARAM,
        EPS_ATT_PARAM,
        PARAMS_LEN
    };

    enum InputId {
        VOCT_INPUT,
        CURRENT_INPUT,
        EPS_INPUT,
        TRIG_INPUT,
        INPUTS_LEN
    };

    enum OutputId {
        OUT_OUTPUT,
        SPIKE_OUTPUT,
        W_OUTPUT,
        OUTPUTS_LEN
    };

    enum LightId {
        LIGHTS_LEN
    };

    // ─── Persistent state ───────────────────────────────────────────────────
    float vv = -1.2f, ww = -0.6f;     // FHN state, seeded near the rest fixed point
    float trigPulse = 0.f;            // decaying injected current from TRIG
    dsp::SchmittTrigger trigIn;       // TRIG input edge detector (hysteresis window)
    dsp::SchmittTrigger spikeDet;     // detects v crossing SPIKE_THRESH upward
    dsp::PulseGenerator spikeGen;     // shapes the SPIKE output pulse
    dsp::TRCFilter<float> dcBlock;    // DC blocker on OUT (limit-cycle mean ≠ 0)
    float lastFs = 0.f;               // detects SR change to refresh coefficients
    float trigDecay = 0.f;            // cached per-sample TRIG pulse decay (fn of fs)
    float antiDenorm = 1e-18f;        // alternating sub-LSB dither, anti-denormal on dcBlock

    // ─── Display trail (phase portrait) ─────────────────────────────────────
    // The audio thread appends decimated (v,w) points to a ring; the UI thread
    // reads dispV/dispW (refreshed ~45 Hz). Float tearing during the copy is
    // harmless for a visualiser and the read no longer races a per-sample writer.
    float trailV[TRAIL] = {}, trailW[TRAIL] = {};
    int   trailIdx = 0, trailDecim = 0;
    float dispV[TRAIL] = {}, dispW[TRAIL] = {};
    int   dispClock = 0;

    // ─── Tunable constants (set by ear/scope at M7; RATE_CAL from tools/) ────
    static constexpr float RATE_CAL    = 37.899004f; // dimensionless period at default → C4 at 0 V
    static constexpr float B_FIXED     = 0.8f;       // recovery linear coefficient b
    static constexpr float HSUB_MAX    = 0.05f;      // max dimensionless substep size
    static constexpr int   MIN_SUB     = 4;
    static constexpr int   MAX_SUB     = 64;
    static constexpr float TRIG_AMP    = 0.6f;       // injected current pulse height
    static constexpr float TRIG_TAU_MS = 3.f;        // pulse decay time constant (ms)
    static constexpr float SPIKE_THRESH = 1.0f;      // v level for a SPIKE pulse
    static constexpr float OUT_GAIN    = 1.0f;
    static constexpr float W_GAIN      = 2.5f;
    static constexpr float CV_DEPTH    = 0.1f;       // ±5V CV × att → ±0.5
    static constexpr float STATE_MAX   = 10.f;       // safety clamp on v,w

    Axon() {
        config(PARAMS_LEN, INPUTS_LEN, OUTPUTS_LEN, LIGHTS_LEN);

        configParam(PITCH_PARAM, -4.f, 4.f, 0.f, "Pitch", " Hz", 2.f, dsp::FREQ_C4);
        configParam(CURRENT_PARAM, -0.2f, 1.6f, 0.6f, "Current (excitability)");
        configParam(EPS_PARAM, 0.01f, 0.30f, 0.08f, "Timescale ε (spike sharpness)");
        configParam(SHAPE_PARAM, 0.4f, 1.0f, 0.7f, "Shape (asymmetry a)");
        configParam(CURRENT_ATT_PARAM, -1.f, 1.f, 0.f, "Current CV");
        configParam(EPS_ATT_PARAM, -1.f, 1.f, 0.f, "ε CV");

        configInput(VOCT_INPUT, "1V/oct pitch");
        configInput(CURRENT_INPUT, "Current CV");
        configInput(EPS_INPUT, "ε CV");
        configInput(TRIG_INPUT, "Trigger (fire a spike)");

        configOutput(OUT_OUTPUT, "Audio (membrane voltage v)");
        configOutput(SPIKE_OUTPUT, "Spike (trigger on each spike)");
        configOutput(W_OUTPUT, "Recovery w (slow CV)");
    }

    void onReset() override {
        vv = -1.2f;
        ww = -0.6f;
        trigPulse = 0.f;
        trigIn.reset();
        spikeDet.reset();
        spikeGen.reset();
        dcBlock.reset();
    }

    // FHN derivatives in dimensionless time. Factored so a Hindmarsh-Rose
    // sibling can add a third (dz) line and reuse the RK4 step below.
    static inline void f(float v, float w, float Itot, float eps, float a,
                         float& dv, float& dw) {
        dv = v - v * v * v / 3.f - w + Itot;
        dw = eps * (v + a - B_FIXED * w);
    }

    // One RK4 substep of size h (updates v,w in place).
    static inline void rk4(float& v, float& w, float h, float Itot, float eps, float a) {
        float k1v, k1w, k2v, k2w, k3v, k3w, k4v, k4w;
        f(v,                  w,                  Itot, eps, a, k1v, k1w);
        f(v + 0.5f * h * k1v, w + 0.5f * h * k1w, Itot, eps, a, k2v, k2w);
        f(v + 0.5f * h * k2v, w + 0.5f * h * k2w, Itot, eps, a, k3v, k3w);
        f(v + h * k3v,        w + h * k3w,        Itot, eps, a, k4v, k4w);
        v += h / 6.f * (k1v + 2.f * k2v + 2.f * k3v + k4v);
        w += h / 6.f * (k1w + 2.f * k2w + 2.f * k3w + k4w);
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
        float eps = clamp(params[EPS_PARAM].getValue()
                        + inputs[EPS_INPUT].getVoltage() * params[EPS_ATT_PARAM].getValue() * CV_DEPTH,
                          0.01f, 0.30f);
        float a   = params[SHAPE_PARAM].getValue();   // no CV by design
        float I   = params[CURRENT_PARAM].getValue()
                  + inputs[CURRENT_INPUT].getVoltage() * params[CURRENT_ATT_PARAM].getValue() * CV_DEPTH;

        // ── excitable trigger: rising edge injects a decaying current pulse ──
        // Hysteresis window (0.1..1V) so a DC-coupled / offset trigger can't latch.
        if (trigIn.process(inputs[TRIG_INPUT].getVoltage(), 0.1f, 1.f))
            trigPulse = TRIG_AMP;
        float Itot = I + trigPulse;
        // decay the pulse once per sample (cached coefficient); flush to zero once
        // negligible so a long rest doesn't grind the multiply on denormals.
        trigPulse *= trigDecay;
        if (std::fabs(trigPulse) < 1e-30f) trigPulse = 0.f;

        // ── pitch = simulation speed (open-loop) ──
        float pitchHz = dsp::FREQ_C4 * std::exp2(
                            params[PITCH_PARAM].getValue() + inputs[VOCT_INPUT].getVoltage());
        float dtau = RATE_CAL * pitchHz / fs;                 // dimensionless advance / sample
        int   K    = clamp((int) std::ceil(dtau / HSUB_MAX), MIN_SUB, MAX_SUB);  // adaptive substeps
        float h    = dtau / K;

        // ── integrate K RK4 substeps ──
        for (int k = 0; k < K; k++)
            rk4(vv, ww, h, Itot, eps, a);

        // ── backstop: nonlinear systems can run away if pushed ──
        if (!std::isfinite(vv) || !std::isfinite(ww)) { vv = -1.2f; ww = -0.6f; }
        vv = clamp(vv, -STATE_MAX, STATE_MAX);
        ww = clamp(ww, -STATE_MAX, STATE_MAX);

        // ── outputs ──
        // Alternating sub-LSB dither keeps the DC-blocker's recursive state off
        // denormals when v is parked at the rest fixed point (sub-threshold, no
        // triggers); ~1e-18 V is inaudible but well above the denormal floor.
        antiDenorm = -antiDenorm;
        dcBlock.process(vv + antiDenorm);
        outputs[OUT_OUTPUT].setVoltage(5.f * std::tanh(dcBlock.highpass() * OUT_GAIN));

        // SchmittTrigger's two-arg form gives hysteresis so a noisy spike peak
        // emits exactly one pulse, not a burst.
        if (spikeDet.process(vv, 0.f, SPIKE_THRESH))
            spikeGen.trigger(1e-3f);
        outputs[SPIKE_OUTPUT].setVoltage(spikeGen.process(args.sampleTime) ? 10.f : 0.f);

        // W is a slow correlated CV — intentionally not high-passed.
        outputs[W_OUTPUT].setVoltage(clamp(ww * W_GAIN, -5.f, 5.f));

        // ── feed the phase-portrait trail (decimated) ──
        if (++trailDecim >= 4) {
            trailDecim = 0;
            trailV[trailIdx] = vv;
            trailW[trailIdx] = ww;
            trailIdx = (trailIdx + 1) % TRAIL;
        }
        if (++dispClock >= (int)(fs / 45.f)) {       // refresh display snapshot ~45 Hz
            dispClock = 0;
            std::copy(trailV, trailV + TRAIL, dispV);
            std::copy(trailW, trailW + TRAIL, dispW);
        }
    }
};


// ─── Phase-portrait display ────────────────────────────────────────────────
// The FHN limit cycle traced in the (v,w) plane: a glowing closed orbit when
// free-running, a point that jumps out and relaxes back on each trigger when
// excitable. The v-nullcline (the cubic) and w-nullcline (the line) are drawn
// as faint guides — their intersection is the fixed point whose stability the
// CURRENT knob controls (the Hopf bifurcation). State is read lock-free from
// the audio thread's snapshot — fine for a visualiser.
struct PhaseDisplay : Widget {
    Axon* module = nullptr;
    std::shared_ptr<Font> font;

    // Plot ranges for v (x axis) and w (y axis). Sized with margin around the
    // orbit envelope (v∈[-2.0,2.0], w∈[-0.6,2.3] across the oscillating param
    // space) so the limit cycle sits comfortably inside the screen.
    static constexpr float VMIN = -2.9f, VMAX = 2.9f;
    static constexpr float WMIN = -1.1f, WMAX = 2.7f;

    void draw(const DrawArgs& args) override {
        // Screen background + bezel (base, non-illuminated layer).
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
            // Clip everything to the screen rectangle. The cubic v-nullcline (and
            // an orbit pushed past its envelope by CV / triggers) would otherwise
            // draw outside the box and spill onto the panel.
            nvgScissor(args.vg, 0.f, 0.f, box.size.x, box.size.y);

            auto X = [&](float v) { return (v - VMIN) / (VMAX - VMIN) * box.size.x; };
            auto Y = [&](float w) { return box.size.y - (w - WMIN) / (WMAX - WMIN) * box.size.y; };

            float a = module ? module->params[Axon::SHAPE_PARAM].getValue() : 0.7f;
            float I = module ? module->params[Axon::CURRENT_PARAM].getValue() : 0.6f;

            // v-nullcline: w = v - v³/3 + I  (the cubic). Where dv/dt = 0.
            nvgBeginPath(args.vg);
            for (int i = 0; i <= 80; i++) {
                float v = VMIN + (VMAX - VMIN) * i / 80.f;
                float w = v - v * v * v / 3.f + I;
                float px = X(v), py = Y(w);
                if (i == 0) nvgMoveTo(args.vg, px, py); else nvgLineTo(args.vg, px, py);
            }
            nvgStrokeColor(args.vg, nvgRGBA(0x40, 0x70, 0x90, 0x55));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);

            // w-nullcline: w = (v + a)/b  (the line). Where dw/dt = 0.
            nvgBeginPath(args.vg);
            {
                float w0 = (VMIN + a) / Axon::B_FIXED, w1 = (VMAX + a) / Axon::B_FIXED;
                nvgMoveTo(args.vg, X(VMIN), Y(w0));
                nvgLineTo(args.vg, X(VMAX), Y(w1));
            }
            nvgStrokeColor(args.vg, nvgRGBA(0x90, 0x50, 0x70, 0x55));
            nvgStrokeWidth(args.vg, 1.f);
            nvgStroke(args.vg);

            // The trajectory trail, oldest→newest, brightening along its length.
            if (module) {
                int idx = module->trailIdx;   // newest just before idx
                nvgLineCap(args.vg, NVG_ROUND);
                for (int k = 1; k < TRAIL; k++) {
                    int i0 = (idx + k - 1) % TRAIL;
                    int i1 = (idx + k) % TRAIL;
                    float alpha = (float) k / TRAIL;   // older = dimmer
                    nvgBeginPath(args.vg);
                    nvgMoveTo(args.vg, X(module->dispV[i0]), Y(module->dispW[i0]));
                    nvgLineTo(args.vg, X(module->dispV[i1]), Y(module->dispW[i1]));
                    nvgStrokeColor(args.vg, nvgRGBA(0x55, 0xc8, 0xff, (int)(alpha * 0xcc)));
                    nvgStrokeWidth(args.vg, 1.6f);
                    nvgStroke(args.vg);
                }
                // Bright head dot at the newest point.
                int newest = (idx + TRAIL - 1) % TRAIL;
                float hx = X(module->dispV[newest]), hy = Y(module->dispW[newest]);
                nvgBeginPath(args.vg); nvgCircle(args.vg, hx, hy, 4.f);
                nvgFillColor(args.vg, nvgRGBA(0x9a, 0xe4, 0xff, 0x55)); nvgFill(args.vg);
                nvgBeginPath(args.vg); nvgCircle(args.vg, hx, hy, 2.f);
                nvgFillColor(args.vg, nvgRGB(0xc0, 0xee, 0xff)); nvgFill(args.vg);
            } else {
                // Browser preview: a static demo orbit.
                nvgBeginPath(args.vg);
                for (int i = 0; i <= 64; i++) {
                    float th = 2.f * (float)M_PI * i / 64.f;
                    float px = X(1.4f * std::cos(th)), py = Y(0.4f + 0.7f * std::sin(th));
                    if (i == 0) nvgMoveTo(args.vg, px, py); else nvgLineTo(args.vg, px, py);
                }
                nvgStrokeColor(args.vg, nvgRGBA(0x55, 0xc8, 0xff, 0x99));
                nvgStrokeWidth(args.vg, 1.6f);
                nvgStroke(args.vg);
            }

            // Screen title.
            if (!font)
                font = APP->window->loadFont(asset::system("res/fonts/DejaVuSans.ttf"));
            if (font) {
                nvgFontFaceId(args.vg, font->handle);
                nvgFontSize(args.vg, mm2px(3.2f));
                nvgTextLetterSpacing(args.vg, mm2px(0.5f));
                nvgFillColor(args.vg, nvgRGBA(0x9a, 0xb0, 0xff, 0xcc));
                nvgTextAlign(args.vg, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
                nvgText(args.vg, mm2px(2.6f), mm2px(2.2f), "AXON", NULL);
                nvgTextLetterSpacing(args.vg, 0.f);
                nvgFontSize(args.vg, mm2px(2.2f));
                nvgFillColor(args.vg, nvgRGBA(0x60, 0x80, 0xb0, 0xaa));
                nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_BOTTOM);
                nvgText(args.vg, box.size.x - mm2px(2.4f), box.size.y - mm2px(1.8f), "v–w", NULL);
            }

            nvgResetScissor(args.vg);
        }
        Widget::drawLayer(args, layer);
    }
};


// ─── Widget ─────────────────────────────────────────────────────────────────

struct AxonWidget : ModuleWidget {
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

        // Knob row labels (knobs at y=54)
        lbl( 9.f,   61.f, 1.8f, dim, "PITCH");
        lbl(24.32f, 61.f, 1.8f, dim, "CURRENT");
        lbl(39.64f, 61.f, 1.8f, dim, "EPS");
        lbl(54.96f, 61.f, 1.8f, dim, "SHAPE");

        // CV jack labels (jacks at y=84)
        lbl(24.32f, 90.f, 1.7f, dim, "I CV");
        lbl(39.64f, 90.f, 1.7f, dim, "ε CV");

        // I/O row labels (row at y=112)
        lbl( 7.5f, 118.5f, 1.7f, dim,    "V/OCT");
        lbl(19.f,  118.5f, 1.7f, dim,    "TRIG");
        lbl(34.f,  118.5f, 1.8f, outclr, "OUT");
        lbl(45.5f, 118.5f, 1.7f, outclr, "SPIKE");
        lbl(56.f,  118.5f, 1.8f, outclr, "W");

        nvgRestore(args.vg);
    }

    AxonWidget(Axon* module) {
        setModule(module);
        setPanel(createPanel(asset::plugin(pluginInstance, "res/Axon.svg")));

        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 1.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(1.0f,   122.0f))));
        addChild(createWidget<ScrewSilver>(mm2px(Vec(54.96f, 122.0f))));

        // Phase-portrait scope across the top.
        PhaseDisplay* disp = new PhaseDisplay();
        disp->module = module;
        disp->box.pos  = mm2px(Vec(5.5f, 8.f));
        disp->box.size = mm2px(Vec(50.f, 34.f));
        addChild(disp);

        // Knob row (y=54): PITCH | CURRENT | EPS | SHAPE
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec( 9.f,   54.f)), module, Axon::PITCH_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(24.32f, 54.f)), module, Axon::CURRENT_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(39.64f, 54.f)), module, Axon::EPS_PARAM));
        addParam(createParamCentered<RoundBlackKnob>(mm2px(Vec(54.96f, 54.f)), module, Axon::SHAPE_PARAM));

        // CV strips under CURRENT and EPS: attenuverter (y=72) + input jack (y=84)
        addParam(createParamCentered<Trimpot>(    mm2px(Vec(24.32f, 72.f)), module, Axon::CURRENT_ATT_PARAM));
        addParam(createParamCentered<Trimpot>(    mm2px(Vec(39.64f, 72.f)), module, Axon::EPS_ATT_PARAM));
        addInput(createInputCentered<PJ301MPort>( mm2px(Vec(24.32f, 84.f)), module, Axon::CURRENT_INPUT));
        addInput(createInputCentered<PJ301MPort>( mm2px(Vec(39.64f, 84.f)), module, Axon::EPS_INPUT));

        // I/O row (y=112): V/OCT, TRIG (in) | OUT, SPIKE, W (out)
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec( 7.5f, 112.f)), module, Axon::VOCT_INPUT));
        addInput(createInputCentered<PJ301MPort>(  mm2px(Vec(19.f,  112.f)), module, Axon::TRIG_INPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(34.f,  112.f)), module, Axon::OUT_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(45.5f, 112.f)), module, Axon::SPIKE_OUTPUT));
        addOutput(createOutputCentered<PJ301MPort>(mm2px(Vec(56.f,  112.f)), module, Axon::W_OUTPUT));
    }
};

Model* modelAxon = createModel<Axon, AxonWidget>("Axon");

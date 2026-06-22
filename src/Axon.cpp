#include "plugin.hpp"
#include <algorithm>
#include <atomic>
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
//
// Polyphonic: up to 16 independent FHN voices. The voice count follows V/OCT
// (falling back to TRIG, so a poly trigger drives poly percussion even with no
// pitch cable). The display traces every active voice, each on its own hue
// stepped across a narrow band around the module's accent colour.

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

    // ─── Per-voice persistent state (polyphonic, up to 16 channels) ─────────
    // Each voice is an independent FHN cell, so the integration state and the
    // stateful DSP helpers (TRIG edge detect, spike detect/shaper, DC blocker)
    // are one-per-voice. The voice count is driven by V/OCT (or TRIG).
    static const int MAX_POLY = 16;
    float vv[MAX_POLY], ww[MAX_POLY];          // FHN state per voice
    float trigPulse[MAX_POLY] = {};            // decaying injected current from TRIG
    dsp::SchmittTrigger trigIn[MAX_POLY];      // TRIG edge detector (hysteresis window)
    dsp::SchmittTrigger spikeDet[MAX_POLY];    // detects v crossing SPIKE_THRESH upward
    dsp::PulseGenerator spikeGen[MAX_POLY];    // shapes the SPIKE output pulse
    dsp::TRCFilter<float> dcBlock[MAX_POLY];   // DC blocker on OUT (limit-cycle mean ≠ 0)
    int   channels = 1;               // active voice count
    float lastFs = 0.f;               // detects SR change to refresh coefficients
    float trigDecay = 0.f;            // cached per-sample TRIG pulse decay (fn of fs)
    float antiDenorm = 1e-18f;        // alternating sub-LSB dither, anti-denormal on dcBlock

    // ─── Display trail (phase portrait) ─────────────────────────────────────
    // The audio thread appends decimated (v,w) points to a per-voice ring, then
    // ~45 Hz publishes a coherent snapshot into a double buffer: it fills the
    // back buffer and flips dispBuf with a release store; the UI reads the front
    // buffer after an acquire load. Lock-free and race-free, and the head index,
    // active-channel count and effective CURRENT travel *with* the arrays, so the
    // trails, head dots and nullcline stay mutually consistent.
    float trailV[MAX_POLY][TRAIL] = {}, trailW[MAX_POLY][TRAIL] = {};
    int   trailIdx = 0, trailDecim = 0;
    float dispV[2][MAX_POLY][TRAIL] = {}, dispW[2][MAX_POLY][TRAIL] = {};
    int   dispHead[2] = {};
    int   dispChannels[2] = {1, 1};
    float dispCurr[2] = {0.6f, 0.6f};   // voice-0 effective CURRENT for the v-nullcline
    std::atomic<int> dispBuf{0};
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

        for (int c = 0; c < MAX_POLY; c++) { vv[c] = -1.2f; ww[c] = -0.6f; }
    }

    void onReset() override {
        for (int c = 0; c < MAX_POLY; c++) {
            vv[c] = -1.2f; ww[c] = -0.6f; trigPulse[c] = 0.f;
            trigIn[c].reset(); spikeDet[c].reset(); spikeGen[c].reset(); dcBlock[c].reset();
        }
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
            for (int c = 0; c < MAX_POLY; c++) dcBlock[c].setCutoffFreq(20.f / fs);
            trigDecay = std::exp(-1.f / (TRIG_TAU_MS * 1e-3f * fs));
            lastFs = fs;
        }

        // Voice count follows V/OCT, but fall back to TRIG so a poly trigger
        // (percussion, no pitch cable) still spreads across voices.
        channels = std::max({inputs[VOCT_INPUT].getChannels(),
                             inputs[TRIG_INPUT].getChannels(), 1});
        outputs[OUT_OUTPUT].setChannels(channels);
        outputs[SPIKE_OUTPUT].setChannels(channels);
        outputs[W_OUTPUT].setChannels(channels);

        // Shared knob values; per-voice CV is added inside the loop.
        const float a         = params[SHAPE_PARAM].getValue();   // no CV by design
        const float epsBase   = params[EPS_PARAM].getValue();
        const float epsAtt    = params[EPS_ATT_PARAM].getValue();
        const float Ibase     = params[CURRENT_PARAM].getValue();
        const float Iatt      = params[CURRENT_ATT_PARAM].getValue();
        const float pitchKnob = params[PITCH_PARAM].getValue();

        // Alternating sub-LSB dither (anti-denormal), toggled once per sample;
        // applied to every voice's DC blocker so a parked rest can't denormalise.
        antiDenorm = -antiDenorm;

        float I0 = Ibase;   // voice-0 effective current, for the display nullcline

        for (int c = 0; c < channels; c++) {
            // ── params → physics (per voice) ──
            float eps = clamp(epsBase
                            + inputs[EPS_INPUT].getPolyVoltage(c) * epsAtt * CV_DEPTH,
                              0.01f, 0.30f);
            float I   = Ibase
                      + inputs[CURRENT_INPUT].getPolyVoltage(c) * Iatt * CV_DEPTH;
            if (c == 0) I0 = I;

            // ── excitable trigger: rising edge injects a decaying current pulse ──
            // Hysteresis window (0.1..1V) so a DC-coupled / offset trigger can't latch.
            if (trigIn[c].process(inputs[TRIG_INPUT].getPolyVoltage(c), 0.1f, 1.f))
                trigPulse[c] = TRIG_AMP;
            float Itot = I + trigPulse[c];
            trigPulse[c] *= trigDecay;
            if (std::fabs(trigPulse[c]) < 1e-30f) trigPulse[c] = 0.f;

            // ── pitch = simulation speed (open-loop) ──
            float pitchHz = dsp::FREQ_C4 * std::exp2(
                                pitchKnob + inputs[VOCT_INPUT].getPolyVoltage(c));
            float dtau = RATE_CAL * pitchHz / fs;                 // dimensionless advance / sample
            int   K    = clamp((int) std::ceil(dtau / HSUB_MAX), MIN_SUB, MAX_SUB);
            float h    = dtau / K;

            // ── integrate K RK4 substeps ──
            float v = vv[c], w = ww[c];
            for (int k = 0; k < K; k++)
                rk4(v, w, h, Itot, eps, a);

            // ── backstop: nonlinear systems can run away if pushed ──
            if (!std::isfinite(v) || !std::isfinite(w)) { v = -1.2f; w = -0.6f; }
            v = clamp(v, -STATE_MAX, STATE_MAX);
            w = clamp(w, -STATE_MAX, STATE_MAX);
            vv[c] = v; ww[c] = w;

            // ── outputs ──
            dcBlock[c].process(v + antiDenorm);
            outputs[OUT_OUTPUT].setVoltage(5.f * std::tanh(dcBlock[c].highpass() * OUT_GAIN), c);

            // Two-arg Schmitt gives hysteresis so a noisy peak emits one pulse.
            if (spikeDet[c].process(v, 0.f, SPIKE_THRESH))
                spikeGen[c].trigger(1e-3f);
            outputs[SPIKE_OUTPUT].setVoltage(spikeGen[c].process(args.sampleTime) ? 10.f : 0.f, c);

            // W is a slow correlated CV — intentionally not high-passed.
            outputs[W_OUTPUT].setVoltage(clamp(w * W_GAIN, -5.f, 5.f), c);
        }

        // ── feed the phase-portrait trails (all active voices, decimated) ──
        if (++trailDecim >= 4) {
            trailDecim = 0;
            for (int c = 0; c < channels; c++) {
                trailV[c][trailIdx] = vv[c];
                trailW[c][trailIdx] = ww[c];
            }
            trailIdx = (trailIdx + 1) % TRAIL;
        }
        if (++dispClock >= (int)(fs / 45.f)) {       // publish display snapshot ~45 Hz
            dispClock = 0;
            int next = 1 - dispBuf.load(std::memory_order_relaxed);
            for (int c = 0; c < channels; c++) {
                std::copy(trailV[c], trailV[c] + TRAIL, dispV[next][c]);
                std::copy(trailW[c], trailW[c] + TRAIL, dispW[next][c]);
            }
            dispHead[next] = trailIdx;
            dispChannels[next] = channels;
            dispCurr[next] = I0;
            dispBuf.store(next, std::memory_order_release);
        }
    }
};


// ─── Phase-portrait display ────────────────────────────────────────────────
// The FHN limit cycle traced in the (v,w) plane: a glowing closed orbit when
// free-running, a point that jumps out and relaxes back on each trigger when
// excitable. The v-nullcline (the cubic) and w-nullcline (the line) are drawn
// as faint guides — their intersection is the fixed point whose stability the
// CURRENT knob controls (the Hopf bifurcation). State is read lock-free from
// the audio thread's snapshot. With several voices patched, each orbit is drawn
// on its own hue stepped across a narrow cyan band so the chord is legible.
struct PhaseDisplay : Widget {
    Axon* module = nullptr;
    std::shared_ptr<Font> font;

    // Plot ranges for v (x axis) and w (y axis). Sized with margin around the
    // orbit envelope (v∈[-2.0,2.0], w∈[-0.6,2.3] across the oscillating param
    // space) so the limit cycle sits comfortably inside the screen.
    static constexpr float VMIN = -2.9f, VMAX = 2.9f;
    static constexpr float WMIN = -1.1f, WMAX = 2.7f;

    // Voice → hue (0..1) within a band centred on Axon's cyan accent (~0.554).
    static float voiceHue(int v, int nv) {
        const float center = 0.554f, halfBand = 0.072f;
        if (nv <= 1) return center;
        return center - halfBand + 2.f * halfBand * v / (nv - 1);
    }

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

            int b = module ? module->dispBuf.load(std::memory_order_acquire) : 0;
            float a = module ? module->params[Axon::SHAPE_PARAM].getValue() : 0.7f;
            float I = module ? module->dispCurr[b] : 0.6f;   // CV-included, coherent with the trail
            int nv = module ? module->dispChannels[b] : 1;

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

            // The trajectory trails, one per active voice, oldest→newest brightening
            // along their length. More voices → dimmer trails so they don't overwhelm.
            if (module) {
                int idx = module->dispHead[b];   // newest just before idx (coherent with arrays)
                float trailA = clamp(204.f - (nv - 1) * 6.f, 112.f, 204.f);  // 0x70..0xcc
                nvgLineCap(args.vg, NVG_ROUND);
                for (int v = 0; v < nv; v++) {
                    const float* dv = module->dispV[b][v];
                    const float* dw = module->dispW[b][v];
                    float hue = voiceHue(v, nv);
                    for (int k = 1; k < TRAIL; k++) {
                        int i0 = (idx + k - 1) % TRAIL;
                        int i1 = (idx + k) % TRAIL;
                        float alpha = (float) k / TRAIL;   // older = dimmer
                        nvgBeginPath(args.vg);
                        nvgMoveTo(args.vg, X(dv[i0]), Y(dw[i0]));
                        nvgLineTo(args.vg, X(dv[i1]), Y(dw[i1]));
                        nvgStrokeColor(args.vg, nvgHSLA(hue, 0.85f, 0.62f, (int)(alpha * trailA)));
                        nvgStrokeWidth(args.vg, 1.6f);
                        nvgStroke(args.vg);
                    }
                    // Bright head dot at the newest point.
                    int newest = (idx + TRAIL - 1) % TRAIL;
                    float hx = X(dv[newest]), hy = Y(dw[newest]);
                    nvgBeginPath(args.vg); nvgCircle(args.vg, hx, hy, 4.f);
                    nvgFillColor(args.vg, nvgHSLA(hue, 0.85f, 0.72f, 0x55)); nvgFill(args.vg);
                    nvgBeginPath(args.vg); nvgCircle(args.vg, hx, hy, 2.f);
                    nvgFillColor(args.vg, nvgHSLA(hue, 0.55f, 0.88f, 0xff)); nvgFill(args.vg);
                }
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

                // Voice-count badge when polyphonic.
                if (nv > 1) {
                    char buf[8];
                    snprintf(buf, sizeof(buf), "%dv", nv);
                    nvgFontSize(args.vg, mm2px(2.2f));
                    nvgFillColor(args.vg, nvgRGBA(0x60, 0x80, 0xb0, 0xcc));
                    nvgTextAlign(args.vg, NVG_ALIGN_RIGHT | NVG_ALIGN_TOP);
                    nvgText(args.vg, box.size.x - mm2px(2.4f), mm2px(2.2f), buf, NULL);
                }

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

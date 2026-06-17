# Axon

A spiking-neuron oscillator for VCV Rack 2.

![Axon](doc/screenshot.png)

Axon is a voice built on the **FitzHugh–Nagumo** (FHN) model — a two-variable
reduction of the Hodgkin–Huxley neuron. The membrane voltage `v` is the audio
output. Above a current threshold the neuron free-runs as a relaxation
oscillator (a slow charge-up followed by a fast spike); below threshold it sits
at rest and fires exactly one spike each time you trigger it. One knob —
**CURRENT** — moves you across that boundary (a Hopf bifurcation), so the same
module is a drone oscillator at one setting and a percussion/transient voice at
another.

## How it works

The state is two coupled variables in dimensionless time:

```
dv/dt = v − v³/3 − w + I        (fast: membrane voltage)
dw/dt = ε·(v + a − b·w)          (slow: recovery)
```

`v` rushes along the cubic nullcline and snaps back (the spike); `w` is the slow
recovery that drags it down and sets up the next spike. **ε** (EPS) is the ratio
of the two timescales — small ε gives a sharp relaxation spike, large ε a
smoother, near-sinusoidal swing. **a** (SHAPE) shifts the asymmetry / threshold.
`b` is fixed at 0.8.

The system is **stiff** (a fast variable riding a slow one), so the integrator is
the real work: each sample takes a number of **RK4 substeps**, and that number
adapts to pitch so the substep size stays bounded (`h ≤ 0.05`) no matter how fast
you play. A finiteness reset + state clamp are the backstop if forcing ever pushes
it to run away. The `f()` derivative and the RK4 step are factored so a planned
**Hindmarsh–Rose** sibling (bursting/chaos) can reuse them with one extra equation.

### Pitch is the simulation *speed*

The limit-cycle period is **emergent** — it depends on CURRENT, EPS and SHAPE —
so pitch is **open-loop calibrated, not phase-locked.** V/OCT maps to how fast
dimensionless time advances, calibrated (`RATE_CAL`) so the default voicing reads
C4 at 0 V. Tracking is within ~1 cent across the useful range at default params,
but **changing CURRENT / EPS / SHAPE detunes the pitch somewhat** — that coupling
is deliberate and part of the instrument's character, not a bug.

## Controls

| Control | Range | Purpose |
| --- | --- | --- |
| **PITCH** | ±4 oct | simulation speed (audio pitch); 0 = C4 |
| **CURRENT** | −0.2 … 1.6 | injected current `I`; the excitability / bifurcation control. Rests at both ends, oscillates in a middle band (~0.33–1.42 at default shape) |
| **EPS** | 0.01 … 0.30 | timescale ratio `ε`; small = sharp spike, large = smoother/near-sine |
| **SHAPE** | 0.4 … 1.0 | waveform asymmetry `a` (threshold position) |

**CURRENT** and **EPS** each have an attenuverter + CV input (±5 V). **V/OCT**
sums with PITCH (no attenuverter). **TRIG** injects a short decaying current pulse
on each rising edge — from rest that fires one spike (percussion); inside the
oscillating band it perturbs the phase. There is no mode switch: the regime is
set purely by where CURRENT sits, and triggers are honoured in both.

Outputs: **OUT** — the membrane voltage `v`, soft-clipped to ±5 V (`tanh`) and
internally DC-blocked at ~20 Hz (the limit cycle's mean is not zero). **SPIKE** —
a 10 V / ~1 ms pulse on each upward threshold crossing of `v` (one per spike, with
hysteresis so a noisy peak can't double-fire). **W** — the recovery variable as a
slow correlated ±5 V CV, intentionally *not* high-passed.

## Display

The screen traces the **phase portrait** in the `(v, w)` plane — the FHN limit
cycle is a glowing closed orbit, and a trigger from rest shows as a single loop
that jumps out and relaxes back. The two faint guides are the **nullclines**
(the cubic `v`-nullcline and the straight `w`-nullcline); their intersection is
the fixed point whose stability CURRENT controls, so you can watch the orbit grow
as CURRENT crosses into the oscillating band. The trail is read lock-free from a
~45 Hz snapshot — fine for a visualiser.

## Patches

`tools/make_patches.py` writes four smoke-test patches into `patches/` (and copies
them into the Windows Rack patches folder if present):

- **axon_1_freerun** — default voicing → audio; play V/OCT
- **axon_2_blips** — sub-threshold CURRENT, an LFO square clocking TRIG (one spike
  per clock); take OUT for the spike voice
- **axon_3_selfevolving** — W self-patched into CURRENT CV: a slow wandering
  texture that rides its own recovery variable
- **axon_4_crossmod** — VCO SAW into CURRENT CV for FM-like sidebands

## Notes / known limits (v1)

- **Pitch is emergent / approximate.** CURRENT, EPS and SHAPE pull the pitch a
  little (see above) — deliberate, not a bug. Calibration targets C4 at the
  default voicing.
- **Aliasing.** Spikes are sharp; the substep oversampling reduces but does not
  eliminate aliasing, so high notes alias — a known v1 limit (v2: more
  oversampling / BLEP on the spike).
- **State is not saved.** `v`, `w` are transient and re-seed at rest on load;
  params persist. Deliberate.
- **DC.** OUT is DC-blocked (the limit-cycle mean ≠ 0). W is intentionally not
  blocked — it's a slow correlated CV.
- **SR changes** need no handler: the only SR-derived cached state (the DC-blocker
  cutoff) is refreshed when the sample rate changes; everything else recomputes
  per sample.

`tools/stability_test.cpp` is a standalone replica of the kernel: it measures the
dimensionless period to set `RATE_CAL`, sweeps CURRENT × EPS × SHAPE × pitch and
asserts `v`,`w` stay finite/bounded, and checks V/OCT tracking.
`tools/render_wav.cpp` auditions voicings offline (writes WAVs).

```bash
g++ -O2 -o /tmp/t tools/stability_test.cpp && /tmp/t     # exit 0 = pass
python3 tools/panel_diagram.py                            # panel footprint check
```

## Deferred to v2

Hindmarsh–Rose sibling module (the integrator is factored for it); hard-sync /
reset input; closed-loop pitch tracking; per-spike velocity output; polyphony.

## References

- R. FitzHugh, *Impulses and physiological states in theoretical models of nerve
  membrane*, Biophysical Journal 1 (6), 1961.
- J. Nagumo, S. Arimoto, S. Yoshizawa, *An active pulse transmission line
  simulating nerve axon*, Proc. IRE 50 (10), 1962.
- [FitzHugh–Nagumo model — Wikipedia](https://en.wikipedia.org/wiki/FitzHugh%E2%80%93Nagumo_model)
- [Scholarpedia: FitzHugh–Nagumo model](http://www.scholarpedia.org/article/FitzHugh-Nagumo_model)

## Build

```bash
# Linux
make RACK_DIR=~/Rack2-SDK/Rack-SDK dist

# Windows cross-compile (from WSL)
RACK_DIR=~/Rack2-SDK-win/Rack-SDK \
  CC=x86_64-w64-mingw32-gcc-posix CXX=x86_64-w64-mingw32-g++-posix \
  STRIP=x86_64-w64-mingw32-strip MACHINE=x86_64-w64-mingw32 make dist
```

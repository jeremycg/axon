# Changelog

All notable changes to this plugin are documented here. Versions follow the
Rack-2 `MAJOR.MINOR.REVISION` convention.

## 2.1.0

- **New module: Soma** — a Hindmarsh-Rose bursting/chaotic neuron oscillator,
  the sibling of Axon. Adds the slow adaptation variable `z` for spike trains,
  bursts, and chaos; CURRENT / BURST (log-mapped `r`) / ADAPT (`s`) controls;
  OUT / SPIKE / Z (burst-envelope CV) outputs; an `(x, z)` phase-portrait display.
- **Display** is now a lock-free double-buffered snapshot: the trail, head dot,
  and nullcline are read from one coherent frame (no data race, no stale-index
  mismatch), and Axon's `v`-nullcline tracks the CV-modulated CURRENT.
- **Performance / robustness:** the per-sample TRIG-pulse decay coefficient is
  cached (was recomputing `exp()` every sample); `trigPulse` flushes to zero and
  the DC blocker gets a sub-LSB anti-denormal dither so a long sub-threshold rest
  can't stall on denormals.

## 2.0.0

- **Initial release: Axon** — a FitzHugh-Nagumo spiking-neuron oscillator.
  Membrane voltage `v` as audio out; pitch as the (open-loop calibrated)
  simulation speed; RK4 with pitch-adaptive substepping; free-running relaxation
  oscillation and excitable (trigger-fired) single-spike behaviour from one
  CURRENT control. OUT (soft-clipped + DC-blocked) / SPIKE / W outputs and a
  `(v, w)` phase-portrait display.

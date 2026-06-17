#!/usr/bin/env python3
"""Generate the §7 Axon smoke-test patches as .vcv files.

Axon positional ids (must match enum order in src/Axon.cpp):
  ParamId : 0 PITCH, 1 CURRENT, 2 EPS, 3 SHAPE, 4 CURRENT_ATT, 5 EPS_ATT
  InputId : 0 VOCT, 1 CURRENT, 2 EPS, 3 TRIG
  OutputId: 0 OUT, 1 SPIKE, 2 W
Axon is 12 HP wide → place downstream modules at x >= 12.

Third-party enum orders (verified from Fundamental v2 source):
  LFO : FREQ_PARAM=2; SQR_OUTPUT=3                (6 HP)
  VCO : FREQ_PARAM=2; PITCH_INPUT=0; SAW_OUTPUT=2  (10 HP)
"""

import json, os, io, glob, shutil, subprocess, sys, tarfile, tempfile, random

random.seed(11)
def uid():
    return random.randint(1_000_000_000_000_000, 9_007_199_254_740_991)

def axon(params, pos):
    return {"id": uid(), "plugin": "Axon", "model": "Axon",
            "version": "2.0.0", "params": params, "pos": pos}

def ap(pitch=0.0, current=0.6, eps=0.08, shape=0.7, current_att=0.0, eps_att=0.0):
    return [
        {"id": 0, "value": float(pitch)},
        {"id": 1, "value": float(current)},
        {"id": 2, "value": float(eps)},
        {"id": 3, "value": float(shape)},
        {"id": 4, "value": float(current_att)},
        {"id": 5, "value": float(eps_att)},
    ]

def audio(pos):
    return {"id": uid(), "plugin": "Core", "model": "AudioInterface",
            "version": "2.6.6", "params": [],
            "data": {"audio": {"driver": -1, "deviceName": "", "sampleRate": 44100.0,
                               "blockSize": 256, "inputOffset": 0, "outputOffset": 0},
                     "dcFilter": True},
            "pos": pos}

def cable(om, oid, im, iid, ci):
    colors = ["#f3374b", "#ffb437", "#00b56e", "#3695ef"]
    return {"id": uid(), "outputModuleId": om, "outputId": oid,
            "inputModuleId": im, "inputId": iid, "color": colors[ci % len(colors)],
            "inputPlugOrder": ci, "outputPlugOrder": ci}

def write_patch(name, modules, cables, master_id):
    patch = {"version": "2.6.6", "zoom": 0.5, "gridOffset": [0.0, 0.0],
             "modules": modules, "cables": cables, "masterModuleId": master_id}
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "patches")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, name)
    with tempfile.TemporaryDirectory() as tmp:
        jp = os.path.join(tmp, "patch.json")
        with open(jp, "w") as f:
            json.dump(patch, f, indent=2)
        tar_buf = io.BytesIO()
        with tarfile.open(fileobj=tar_buf, mode="w:") as tf:
            tf.add(jp, arcname="patch.json")
        r = subprocess.run(["zstd", "-19", "-o", out_file, "-f"],
                           input=tar_buf.getvalue(), capture_output=True)
        if r.returncode != 0:
            print("zstd error:", r.stderr.decode(), file=sys.stderr); sys.exit(1)
    print(f"  {name}: {len(modules)} modules, {len(cables)} cables, {os.path.getsize(out_file)} bytes")
    for win in glob.glob("/mnt/c/Users/*/AppData/Local/Rack2/patches"):
        shutil.copy2(out_file, os.path.join(win, name))
        print(f"    installed -> {win}/{name}")

# ── 1. Free-run tone: default voicing → audio ────────────────────────────────
def patch_freerun():
    x = axon(ap(current=0.6, eps=0.08, shape=0.7), [0, 0])
    a = audio([12, 0])   # Axon is 12 HP
    cs = [cable(x["id"], 0, a["id"], 0, 0), cable(x["id"], 0, a["id"], 1, 1)]
    write_patch("axon_1_freerun.vcv", [x, a], cs, a["id"])

# ── 2. Excitable blips: LFO square clocks TRIG, sub-threshold current ─────────
def patch_blips():
    lfo = {"id": uid(), "plugin": "Fundamental", "model": "LFO", "version": "2.6.4",
           "params": [{"id": 2, "value": 1.0}], "pos": [0, 0]}        # FREQ=1 → 2 Hz
    x = axon(ap(current=0.1, eps=0.08, shape=0.7), [6, 0])            # LFO is 6 HP
    a = audio([18, 0])                                                # Axon spans 6..18
    cs = [
        cable(lfo["id"], 3, x["id"], 3, 0),   # LFO SQR -> Axon TRIG
        cable(x["id"], 0, a["id"], 0, 1),     # Axon OUT -> L
        cable(x["id"], 0, a["id"], 1, 2),     # Axon OUT -> R
    ]
    write_patch("axon_2_blips.vcv", [lfo, x, a], cs, a["id"])

# ── 3. Self-evolving: W feeds CURRENT CV (rides its own recovery variable) ────
def patch_selfevolving():
    x = axon(ap(current=0.7, eps=0.08, shape=0.7, current_att=0.8), [0, 0])
    a = audio([12, 0])
    cs = [
        cable(x["id"], 2, x["id"], 1, 0),     # Axon W -> Axon CURRENT CV (self-patch)
        cable(x["id"], 0, a["id"], 0, 1),
        cable(x["id"], 0, a["id"], 1, 2),
    ]
    write_patch("axon_3_selfevolving.vcv", [x, a], cs, a["id"])

# ── 4. Cross-mod: VCO SAW into CURRENT CV for FM-like sidebands ───────────────
def patch_crossmod():
    vco = {"id": uid(), "plugin": "Fundamental", "model": "VCO", "version": "2.6.4",
           "params": [{"id": 2, "value": 0.0}], "pos": [0, 0]}        # VCO is 10 HP
    x = axon(ap(current=0.6, eps=0.08, shape=0.7, current_att=0.4), [10, 0])
    a = audio([22, 0])                                                # Axon spans 10..22
    cs = [
        cable(vco["id"], 2, x["id"], 1, 0),   # VCO SAW -> Axon CURRENT CV
        cable(x["id"], 0, a["id"], 0, 1),
        cable(x["id"], 0, a["id"], 1, 2),
    ]
    write_patch("axon_4_crossmod.vcv", [vco, x, a], cs, a["id"])

if __name__ == "__main__":
    print("Generating Axon smoke-test patches:")
    patch_freerun(); patch_blips(); patch_selfevolving(); patch_crossmod()
    print("Done.")

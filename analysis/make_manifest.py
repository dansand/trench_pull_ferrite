"""make_manifest.py — the central provenance record of every shipped model directory.

The driver stamps a provenance.txt into each model it produces, but the shipped models predate the stamp and carry
none.  This script records, per model directory under data/: the generating command and its key parameters (from the
driver's production block), the stress frame the stored sigma_* fields are in, the SHA-256 and size of every file,
and the model's origin.  Hashes are computed, never typed.  Two outputs, the same content:

    data/DATA_MANIFEST.json   machine-readable (gpe_analysis.Model.stress_frame reads it)
    data/DATA_MANIFEST.md     the same, as tables

    python analysis/make_manifest.py            # regenerate both from what is in data/
    python analysis/make_manifest.py --check    # verify every shipped file against the committed JSON; exit 1 on any
                                                # missing, extra or changed file (reproduce.sh runs this)
"""
import hashlib, json, os, re, sys
from datetime import date

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
FILES = ("gpe_model.vtu", "gpe_topo.csv", "tuned_V.txt")

# Fixed facts of the production runs (model/paper_models.jl, production block; model/idealized_beam.jl).
BASE = dict(E_Pa=70.0e9, nu=0.25, h_km=60, nx=800, nz=48, nsteps=24, L_km=1600, rho_a=3300, rho_w=1000, g=9.81)

# Origin — recorded facts, verified 2026-09-12 by SHA-256 against the private development archive
# (ferrite_plate_flexure, HEAD 6e7fb08 of 2026-09-11): every shipped file is byte-identical to its archive copy.
ARCHIVE = "the private development archive (ferrite_plate_flexure, HEAD 6e7fb08, 2026-09-11)"
ORIGIN = {
    "suite1_strength":  f"committed in {ARCHIVE} at finite_strain/out/paper/set1_rheology/<model> (commit 0c6fe38, 2026-07-08)",
    "suite4_thickness": f"committed in {ARCHIVE} at finite_strain/out/paper/set4_thickness/<model> (commit 26e6c9b, 2026-08-11)",
    "suite2_load":      f"on-disk output of `paper_models.jl v_sweep` in the working tree of {ARCHIVE} at finite_strain/out/paper/set3_v_sweep/<model>; never committed there",
    "suite3_background": f"on-disk output of `paper_models.jl nd_sweep` in the working tree of {ARCHIVE} at finite_strain/out/paper/set2_nd_sweep/<model>; never committed there",
    "idealized_beam":   "re-solved with the release code (model/idealized_beam.jl) on 2026-09-14; identical to the archive copy (finite_strain/out/idealized_beam, commit 10d0d39) on every shared field, plus the three Cauchy-gradient fields the archive export predated",
    "idealized_beam_elastic": "re-solved with the release code (model/idealized_beam.jl) on 2026-09-14; identical to the archive copy (finite_strain/out/idealized_beam_elastic, commit 758e4dc) on every shared field, plus the three Cauchy-gradient fields the archive export predated",
}


def describe(rel):
    """(command, stress_frame, parameters) for a model directory, from its suite and name."""
    suite, _, name = rel.partition("/")
    cmd = "julia --project=. model/paper_models.jl "
    if suite == "suite1_strength":
        p = dict(BASE, V_TN=4.0)
        if name == "tresca_deep_150_60km_V4":   p.update(strength="Tresca, uniform", sigma_Y_MPa=150)
        elif name == "elastic_deep_60km_V4":    p.update(strength="elastic (sigma_Y -> inf: 1e12 Pa)")
        elif name == "dd_vm_asym_60km_V4":      p.update(strength="DD-VM asymmetric (depth-dependent von Mises)", mu=0.15, cohesion_MPa=31.534, hardening_H_Pa=7.0e8, nsteps=22)
        elif name == "dd_vm_sym_60km_V4":       p.update(strength="DD-VM symmetric (depth-dependent von Mises)", mu=0.35, cohesion_MPa=36.695, hardening_H_Pa=7.0e8, nsteps=22)
        else: raise KeyError(rel)
        return cmd + "suite1", "massless", p
    if suite == "suite2_load":
        v = float(re.search(r"_V(\d+p?\d*)$", name).group(1).replace("p", "."))
        return cmd + "v_sweep", "massless", dict(BASE, strength="Tresca, uniform", sigma_Y_MPa=150, V_TN=v)
    if suite == "suite3_background":
        n = float(re.search(r"_mem(-?\d+)$", name).group(1))
        return cmd + "nd_sweep", "massless", dict(BASE, strength="Tresca, uniform", sigma_Y_MPa=150, V_TN=4.0, N_mem_TN=n)
    if suite == "suite4_thickness":
        hk = int(re.search(r"_(\d+)km$", name).group(1)); V0 = {30: 2.05, 40: 2.70, 50: 3.35}[hk]
        p = dict(BASE, h_km=hk, strength="Tresca, uniform", sigma_Y_MPa=150, w_target_m=3232.9, V0_TN=V0)
        tv = os.path.join(DATA, rel, "tuned_V.txt")
        if os.path.isfile(tv):
            for line in open(tv):
                k, _, v = line.strip().partition("=")
                if k == "V_TN": p["V_TN_tuned"] = float(v)
        return cmd + "thickness", "massless", p
    if suite in ("idealized_beam", "idealized_beam_elastic"):
        p = dict(nondimensional=True, L=10.0, h=1.0, nx=200, nz=40, E=60.0, nu=0.25, sigma_Y=(1.0e6 if suite.endswith("elastic") else 1.0),
                 delta_max=2.2, nsteps=44, element_order=2, foundation="none", gravity="none")
        return "julia --project=. model/idealized_beam.jl", "massless", p
    raise KeyError(rel)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def model_dirs():
    out = []
    for suite in sorted(os.listdir(DATA)):
        sd = os.path.join(DATA, suite)
        if not os.path.isdir(sd) or suite.startswith("_") or suite in ("reference_figures", "convergence"):
            continue                      # data/convergence is regeneration-only (Table S3) and never shipped: not a manifest member
        if os.path.isfile(os.path.join(sd, "gpe_model.vtu")):
            out.append(suite); continue
        for name in sorted(os.listdir(sd)):
            if os.path.isfile(os.path.join(sd, name, "gpe_model.vtu")):
                out.append(f"{suite}/{name}")
    return out


def build():
    models = {}
    for rel in model_dirs():
        cmd, frame, params = describe(rel)
        suite = rel.partition("/")[0]
        files = {f: {"sha256": sha256(os.path.join(DATA, rel, f)), "bytes": os.path.getsize(os.path.join(DATA, rel, f))}
                 for f in FILES if os.path.isfile(os.path.join(DATA, rel, f))}
        models[rel] = {"command": cmd, "stress_frame": frame, "parameters": params, "files": files,
                       "origin": ORIGIN[suite].replace("<model>", rel.partition("/")[2] or rel)}
    return {"generated": str(date.today()), "generator": "analysis/make_manifest.py", "models": models}


def write_md(man, path):
    L = ["# Data manifest — provenance of every shipped model", "",
         f"Generated {man['generated']} by `python analysis/make_manifest.py` (hashes computed, never typed); verify with "
         "`python analysis/make_manifest.py --check`. The machine-readable twin is `DATA_MANIFEST.json`, which "
         "`gpe_analysis.Model.stress_frame()` reads. The shipped models predate the driver's `provenance.txt` stamp; "
         "this file is their provenance record.", "",
         "**Stress frame** `massless` for every model: the stored `sigma_*` fields are the 2nd Piola–Kirchhoff stress "
         "of a run without gravity or lithostatic prestress (see README §4 for what the fields are).", ""]
    by_suite = {}
    for rel, e in man["models"].items():
        by_suite.setdefault(rel.partition("/")[0], []).append((rel, e))
    for suite, items in by_suite.items():
        L += [f"## `{suite}`", "", f"- command: `{items[0][1]['command']}`", f"- origin: {items[0][1]['origin'].replace(items[0][0].partition('/')[2] or items[0][0], '<model>')}", ""]
        L += ["| model | parameters | file | bytes | sha256 |", "|---|---|---|---|---|"]
        for rel, e in items:
            name = rel.partition("/")[2] or rel
            pr = ", ".join(f"{k}={v}" for k, v in e["parameters"].items())
            first = True
            for f, d in e["files"].items():
                L.append(f"| {'`'+name+'`' if first else ''} | {pr if first else ''} | `{f}` | {d['bytes']} | `{d['sha256']}` |")
                first = False
        L.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(L))


def check():
    with open(os.path.join(DATA, "DATA_MANIFEST.json")) as f:
        man = json.load(f)
    bad = 0
    present = set(model_dirs())
    for rel in sorted(set(man["models"]) | present):
        if rel not in man["models"]:
            print(f"EXTRA    {rel}: model directory not in the manifest"); bad += 1; continue
        if rel not in present:
            print(f"MISSING  {rel}: manifest entry has no model directory"); bad += 1; continue
        for f, d in man["models"][rel]["files"].items():
            p = os.path.join(DATA, rel, f)
            if not os.path.isfile(p):
                print(f"MISSING  {rel}/{f}"); bad += 1
            elif sha256(p) != d["sha256"]:
                print(f"CHANGED  {rel}/{f}: sha256 differs from the manifest"); bad += 1
        for f in FILES:
            if os.path.isfile(os.path.join(DATA, rel, f)) and f not in man["models"][rel]["files"]:
                print(f"EXTRA    {rel}/{f}: not in the manifest"); bad += 1
    n = sum(len(e["files"]) for e in man["models"].values())
    print(f"manifest check: {len(man['models'])} models, {n} files, {bad} problem(s)")
    return bad


if __name__ == "__main__":
    if "--check" in sys.argv[1:]:
        sys.exit(1 if check() else 0)
    man = build()
    with open(os.path.join(DATA, "DATA_MANIFEST.json"), "w") as f:
        json.dump(man, f, indent=1)
    write_md(man, os.path.join(DATA, "DATA_MANIFEST.md"))
    print(f"wrote data/DATA_MANIFEST.json and .md: {len(man['models'])} models, {sum(len(e['files']) for e in man['models'].values())} files")

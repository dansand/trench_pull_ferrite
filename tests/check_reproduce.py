"""check_reproduce.py — assert the headline numbers of a reproduce.sh run from its logs.   python tests/check_reproduce.py .reproduce_logs

Reads the per-step logs written by reproduce.sh and the executed notebook, and checks, with tolerances:
  ΔGPE* = 2.542 TN/m ± 0.5 %   identity |ΔN_D − ΔGPE*|/ΔGPE* < 0.05 %   arm ΔGPE*/(Δρ g w_T) = 34.9 ± 0.2 km   (notebook §4–5,
  and render_profiles' panel-(d) self-check)   S1/S2 benchmark assertions passed   the Suite-1 reconstruction table
  (render_gpe_compare: plate-top arm ≤ 3 %, sea-level arm ≤ 10 %)   S1 deflection offset 1.6 ± 0.3 %   the isostatic-column
  assumption (reference |change| ≤ 0.5 %, all models ≤ 5 %)   every figure file rewritten after the run started.
Exit 1 on any failure; prints every check either way."""
import glob, json, os, re, sys

LOG = sys.argv[1] if len(sys.argv) > 1 else ".reproduce_logs"
FIGURES = ["hero_tresca_deep60", "gpe_correlation", "gpe_compare_suite1", "profiles", "thickness_compare", "benchmark_boef",
           "benchmark_mp", "core_profiles_deep60", "corrected_density", "hero_tresca_30km", "hero_tresca_40km",
           "hero_dd_vm_asym", "hero_dd_vm_sym"]
bad = 0


def check(ok, msg):
    global bad
    print(("  ok    " if ok else "  FAIL  ") + msg)
    bad += (not ok)


def log(name):
    p = os.path.join(LOG, name + ".log")
    return open(p, encoding="utf-8", errors="replace").read() if os.path.isfile(p) else ""


print("check_reproduce:")
# --- benchmarks (S1, S2): the scripts assert their own misfits; here we only require that they said so
for name in ("benchmark", "mp_benchmark"):
    check("benchmark assertions passed" in log(name), f"{name}: in-script misfit assertions passed")
m = re.search(r"end shear, from FE\)\s*=\s*([\d.]+) TN/m", log("benchmark"))
check(m is not None and abs(float(m.group(1)) - 4.0) <= 0.04, f"S1 end resultant = {m and m.group(1)} TN/m (4.000 ± 1 %)")
m = re.search(r"offset\s*([+-][\d.]+)%", log("benchmark"))
check(m is not None and abs(float(m.group(1)) - 1.6) <= 0.3, f"S1 thick-plate deflection offset = {m and m.group(1)} % (1.6 ± 0.3)")
m = re.search(r"sections on the curve:\s*(\d+)\s+mean .*?=\s*([\d.]+)%\s+max\s+([\d.]+)%", log("mp_benchmark"))
check(m is not None and int(m.group(1)) == 148 and float(m.group(2)) <= 0.15 and float(m.group(3)) <= 1.0,
      f"S2 M–κ: {m and m.groups()} (148 sections, mean ≤ 0.15 %, max ≤ 1 %)")

# --- reconstruction table (render_gpe_compare): four Suite-1 rows
rows = re.findall(r"^\s*(.+?)\s*:\s*plate-top arm\s+([\d.]+)%,\s+sea-level arm\s+([\d.]+)%", log("gpe_compare"), re.M)
check(len(rows) == 4, f"gpe_compare reconstruction table has 4 rows (found {len(rows)})")
for title, etop, esea in rows:
    check(float(etop) <= 3.0 and float(esea) <= 10.0, f"reconstruction {title.strip()}: plate-top {etop} % (≤ 3), sea-level {esea} % (≤ 10)")
tr = [r for r in rows if "Tresca" in r[0]]
check(bool(tr) and float(tr[0][1]) <= 3.0 and 5.0 <= float(tr[0][2]) <= 11.0,
      f"reconstruction Tresca (uniform): plate-top {tr and tr[0][1]} % (~2), sea-level {tr and tr[0][2]} % (~8)")

# --- render_profiles panel-(d) self-check: the Tresca trench_pull value printed there
m = re.search(r"Tresca.*?trench_pull=\+?([\d.]+)", log("profiles"))
check(m is not None and abs(float(m.group(1)) / 2.542 - 1) <= 0.005, f"profiles: Tresca trench_pull = {m and m.group(1)} TN/m (2.542 ± 0.5 %)")

# --- the notebook: ΔGPE*, identity residual, arm, as printed by START_HERE §4–5
nbp = os.path.join(LOG, "START_HERE_executed.ipynb")
if os.path.isfile(nbp):
    out = ""
    for c in json.load(open(nbp, encoding="utf-8"))["cells"]:
        for o in c.get("outputs", []):
            out += "".join(o.get("text", [])) if "text" in o else ""
    m = re.search(r"ΔGPE\* = ([\d.]+) TN/m;\s+ΔN_D = ([\d.]+) TN/m;\s+residual[^=]*=\s*([\d.]+) %", out)
    check(m is not None and abs(float(m.group(1)) / 2.542 - 1) <= 0.005, f"notebook: ΔGPE* = {m and m.group(1)} TN/m (2.542 ± 0.5 %)")
    check(m is not None and float(m.group(3)) < 0.05, f"notebook: identity residual = {m and m.group(3)} % (< 0.05 %)")
    m = re.search(r"effective arm ΔGPE\*/\(Δρ g w_T\) = ([\d.]+) km", out)
    check(m is not None and abs(float(m.group(1)) - 34.9) <= 0.2, f"notebook: arm = {m and m.group(1)} km (34.9 ± 0.2)")
    check("Error" not in "".join(o.get("output_type", "") for c in json.load(open(nbp, encoding="utf-8"))["cells"] for o in c.get("outputs", [])),
          "notebook: no error outputs")
elif os.path.isfile(os.path.join(LOG, "NOTEBOOK_SKIPPED")):
    print("  skip  notebook not executed (--no-nb)")
else:
    check(False, "notebook: execution FAILED — no executed copy was produced (see notebook.log)")

# --- the isostatic-column assumption (analysis/isostatic_column_test.py): cost of taking x_I as lithostatic
m = re.search(r"Reference model \(Tresca, h = 60 km, V = 4 TN/m\): ([+-][\d.]+) %", log("isostatic_column"))
check(m is not None and abs(float(m.group(1))) <= 0.5, f"isostatic-column assumption, reference model: {m and m.group(1)} % (|change| ≤ 0.5 %)")
rows = re.findall(r"^\| suite\d\S* \| `[^`]+` \|(?:[^|]*\|){5}\s*([+-][\d.]+) \|", log("isostatic_column"), re.M)
check(len(rows) == 20 and all(abs(float(v)) <= 5.0 for v in rows),
      f"isostatic-column assumption, all {len(rows)} production models: max |change| {max((abs(float(v)) for v in rows), default=float('nan')):.2f} % (≤ 5 %)")

# --- every figure rewritten after the run started
t0 = float(open(os.path.join(LOG, "START")).read()) if os.path.isfile(os.path.join(LOG, "START")) else 0.0
for f in FIGURES:
    p = f"figures/{f}.png"
    check(os.path.isfile(p) and os.path.getmtime(p) >= t0 - 1, f"figure rewritten: {p}")
for f in ("equilibration_compare", "ridge_trench_overview_v2", "taux_cases"):
    p = f"schematic/{f}.pdf"
    if os.path.isfile(os.path.join(LOG, "schematic.log")):
        check(os.path.isfile(p) and os.path.getmtime(p) >= t0 - 1, f"schematic rebuilt: {p}")

print(f"check_reproduce: {bad} failure(s)")
sys.exit(1 if bad else 0)

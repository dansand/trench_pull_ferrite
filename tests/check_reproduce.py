"""check_reproduce.py — assert the numbers of a reproduce.sh run.   python tests/check_reproduce.py .reproduce_logs

Every number the paper quotes is read from the tables the scripts wrote in the same pass as their figures
(tables/*.csv, see gpe_analysis.write_table) — nothing is scraped from printed text.  Checks, with tolerances:
  model_summary: reference ΔGPE* = 2.542 TN/m ± 0.5 %, identity < 0.05 %, arm 34.9 ± 0.2 km, all 20 models identity < 0.1 %
  benchmark_boef: end resultant 4.000 ± 1 %, V misfit ≤ 1 %, deflection offset 1.6 ± 0.3 %, parabola ≤ 0.5 %
  benchmark_mp: 148 sections, mean ≤ 0.15 %, max ≤ 1 %        gpe_compare_reconstruction: plate-top ≤ 3 %, sea-level ≤ 10 %
  profiles_selfcheck: Tresca trench_pull = model_summary's (same function, same data) and area agrees to < 0.5 %
  isostatic_column: reference |change| ≤ 0.5 %, all models ≤ 5 %      notebook: ΔGPE*, identity, arm as printed, no errors
  every figure AND every table rewritten after the run started.
Exit 1 on any failure; prints every check either way."""
import glob, json, os, re, sys

sys.path.insert(0, "analysis")
from gpe_analysis import read_table

LOG = sys.argv[1] if len(sys.argv) > 1 else ".reproduce_logs"
FIGURES = ["hero_tresca_deep60", "gpe_correlation", "gpe_compare_suite1", "profiles", "thickness_compare", "benchmark_boef",
           "benchmark_mp", "core_profiles_deep60", "corrected_density", "hero_tresca_30km", "hero_tresca_40km",
           "hero_dd_vm_asym", "hero_dd_vm_sym"]
TABLES = ["model_summary", "isostatic_column", "frame_check", "edge_exclusion", "benchmark_boef", "benchmark_mp", "gpe_compare_reconstruction", "gpe_correlation",
          "thickness_compare", "profiles_selfcheck", "corrected_density", "hero_tresca_deep60", "hero_tresca_30km", "hero_tresca_40km",
          "hero_dd_vm_asym", "hero_dd_vm_sym"]
REF = "tresca_deep_150_60km_V4"
bad = 0


def check(ok, msg):
    global bad
    print(("  ok    " if ok else "  FAIL  ") + msg)
    bad += (not ok)


def log(name):
    p = os.path.join(LOG, name + ".log")
    return open(p, encoding="utf-8", errors="replace").read() if os.path.isfile(p) else ""


def table(name):
    p = f"tables/{name}.csv"
    return read_table(p) if os.path.isfile(p) else ({}, [])


def scalars(name):
    return {r["quantity"]: r["value"] for r in table(name)[1]}


print("check_reproduce:")
# --- model_summary: the headline numbers and the identity on every model
_, ms = table("model_summary")
ref = next((r for r in ms if r["model"] == REF and r["suite"] == "suite1_strength"), None)
check(len(ms) == 20, f"model_summary has 20 production models (found {len(ms)})")
check(ref is not None and abs(ref["dGPE_TN"] / 2.542 - 1) <= 0.005, f"reference ΔGPE* = {ref and ref['dGPE_TN']} TN/m (2.542 ± 0.5 %)")
check(ref is not None and ref["identity_pct"] < 0.05, f"reference identity residual = {ref and ref['identity_pct']} % (< 0.05 %)")
check(ref is not None and abs(ref["arm_km"] - 34.9) <= 0.2, f"reference arm = {ref and ref['arm_km']} km (34.9 ± 0.2)")
check(bool(ms) and all(r["identity_pct"] < 0.1 for r in ms), f"identity residual < 0.1 % on all models (max {max((r['identity_pct'] for r in ms), default=float('nan')):.3f} %)")

# --- benchmarks S1, S2 (the scripts also assert these in-script; the tables are what the SI quotes)
for name in ("benchmark", "mp_benchmark"):
    check("benchmark assertions passed" in log(name), f"{name}: in-script misfit assertions passed")
b = scalars("benchmark_boef")
check(abs(b.get("end_resultant_from_FE", 0) - 4.0) <= 0.04, f"S1 end resultant = {b.get('end_resultant_from_FE')} TN/m (4.000 ± 1 %)")
check(b.get("V_misfit_max_0_4alpha", 99) <= 1.0, f"S1 V misfit = {b.get('V_misfit_max_0_4alpha')} % (≤ 1)")
check(abs(b.get("w0_offset_vs_thin_beam", 99) - 1.6) <= 0.3, f"S1 thick-plate deflection offset = {b.get('w0_offset_vs_thin_beam')} % (1.6 ± 0.3)")
check(b.get("shear_parabola_misfit_max", 99) <= 0.5, f"S1 shear parabola misfit = {b.get('shear_parabola_misfit_max')} % (≤ 0.5)")
b = scalars("benchmark_mp")
check(b.get("sections_on_curve") == 148 and b.get("mean_misfit_M_over_Mp", 99) <= 0.15 and b.get("max_misfit_M_over_Mp", 99) <= 1.0,
      f"S2 M–κ: {b.get('sections_on_curve')} sections, mean {b.get('mean_misfit_M_over_Mp')} %, max {b.get('max_misfit_M_over_Mp')} % (148, ≤ 0.15, ≤ 1)")

# --- Fig 5 reconstruction table
_, rows = table("gpe_compare_reconstruction")
check(len(rows) == 4, f"gpe_compare reconstruction table has 4 rows (found {len(rows)})")
for r in rows:
    check(r["plate_top_arm_err_pct"] <= 3.0 and r["sea_level_arm_err_pct"] <= 10.0,
          f"reconstruction {r['label']}: plate-top {r['plate_top_arm_err_pct']} % (≤ 3), sea-level {r['sea_level_arm_err_pct']} % (≤ 10)")

# --- Fig 6 self-check, cross-checked against model_summary (same function, same data ⇒ must agree)
_, rows = table("profiles_selfcheck")
tr = next((r for r in rows if REF in str(r["model"])), None)
check(tr is not None and tr["diff_pct"] < 0.5, f"profiles: Tresca panel-(d) area vs trench_pull differ by {tr and tr['diff_pct']} % (< 0.5)")
check(tr is not None and ref is not None and abs(tr["trench_pull_TN"] / ref["dGPE_TN"] - 1) < 1e-3,
      f"cross-check: profiles trench_pull {tr and tr['trench_pull_TN']} = model_summary ΔGPE* {ref and ref['dGPE_TN']} (< 0.1 %)")

# --- the isostatic-column assumption
_, rows = table("isostatic_column")
ir = next((r for r in rows if r["model"] == REF and r["suite"] == "suite1_strength"), None)
check(ir is not None and abs(ir["change_pct"]) <= 0.5, f"isostatic-column assumption, reference model: {ir and ir['change_pct']} % (|change| ≤ 0.5 %)")
check(len(rows) == 20 and all(abs(r["change_pct"]) <= 5.0 for r in rows),
      f"isostatic-column assumption, all {len(rows)} models: max |change| {max((abs(r['change_pct']) for r in rows), default=float('nan')):.2f} % (≤ 5 %)")

# --- the notebook: as printed by START_HERE §4–5, cross-checked against model_summary
nbp = os.path.join(LOG, "START_HERE_executed.ipynb")
if os.path.isfile(nbp):
    cells = json.load(open(nbp, encoding="utf-8"))["cells"]
    out = "".join("".join(o.get("text", [])) for c in cells for o in c.get("outputs", []) if "text" in o)
    m = re.search(r"ΔGPE\* = ([\d.]+) TN/m;\s+ΔN_D = ([\d.]+) TN/m;\s+residual[^=]*=\s*([\d.]+) %", out)
    check(m is not None and ref is not None and abs(float(m.group(1)) / ref["dGPE_TN"] - 1) <= 1e-3, f"notebook: ΔGPE* = {m and m.group(1)} TN/m (= model_summary)")
    check(m is not None and float(m.group(3)) < 0.05, f"notebook: identity residual = {m and m.group(3)} % (< 0.05 %)")
    m = re.search(r"effective arm ΔGPE\*/\(Δρ g w_T\) = ([\d.]+) km", out)
    check(m is not None and abs(float(m.group(1)) - 34.9) <= 0.2, f"notebook: arm = {m and m.group(1)} km (34.9 ± 0.2)")
    check(not any(o.get("output_type") == "error" for c in cells for o in c.get("outputs", [])), "notebook: no error outputs")
elif os.path.isfile(os.path.join(LOG, "NOTEBOOK_SKIPPED")):
    print("  skip  notebook not executed (--no-nb)")
else:
    check(False, "notebook: execution FAILED — no executed copy was produced (see notebook.log)")

# --- every figure and every table rewritten after the run started
t0 = float(open(os.path.join(LOG, "START")).read()) if os.path.isfile(os.path.join(LOG, "START")) else 0.0
for f in FIGURES:
    p = f"figures/{f}.png"
    check(os.path.isfile(p) and os.path.getmtime(p) >= t0 - 1, f"figure rewritten: {p}")
present = sorted(os.path.splitext(os.path.basename(p))[0] for p in glob.glob("tables/*.csv"))
for t in TABLES:
    p = f"tables/{t}.csv"
    check(os.path.isfile(p) and os.path.getmtime(p) >= t0 - 1, f"table rewritten: {p}")
extra = [t for t in present if t not in TABLES and t != "convergence"]      # convergence is regeneration-only (shipped copy)
check(not extra, f"every table in tables/ is one the harness asserts (unlisted: {extra})")
for t in extra:
    check(os.path.getmtime(f"tables/{t}.csv") >= t0 - 1, f"table rewritten: tables/{t}.csv")
for p in ("tables/paper_numbers.tex", "tables/PAPER_NUMBERS.md"):
    check(os.path.isfile(p) and os.path.getmtime(p) >= t0 - 1, f"paper numbers rewritten: {p}")
for f in ("equilibration_compare", "ridge_trench_overview_v2", "taux_cases"):
    p = f"schematic/{f}.pdf"
    if os.path.isfile(os.path.join(LOG, "schematic.log")):
        check(os.path.isfile(p) and os.path.getmtime(p) >= t0 - 1, f"schematic rebuilt: {p}")

print(f"check_reproduce: {bad} failure(s)")
sys.exit(1 if bad else 0)

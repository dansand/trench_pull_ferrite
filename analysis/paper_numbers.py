"""paper_numbers.py — the verified numbers, ready to drop into the manuscript.

Reads tables/*.csv (the source of truth) and writes
    tables/paper_numbers.tex   one \\newcommand per quoted number, at the precision the paper uses
    tables/PAPER_NUMBERS.md    the placement guide: macro, value, source (table, row, column), where it goes in the
                               manuscript, and — with --compare — whether the text currently agrees
Run by reproduce.sh.  With `--compare DIR` (DIR = the manuscript root holding main.tex, sections/, si.tex) it reads the
LaTeX and reports, for every registered location, whether the literal the text carries matches the table value.

    python analysis/paper_numbers.py
    python analysis/paper_numbers.py --compare ~/projects/mypapers/trench_pull_force/2026_codex/full_manuscript

Registry entries (below) are the ONLY hand-maintained part: name, source table, selector, format, description, and the
manuscript locations.  Adding a quoted number = adding an entry; the value itself is never typed.
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gpe_analysis import read_table

REF = "tresca_deep_150_60km_V4"


def row(name, **key):
    _, rows = read_table(f"tables/{name}.csv")
    for r in rows:
        if all(str(r.get(k)) == str(v) or r.get(k) == v for k, v in key.items()):
            return r
    raise KeyError((name, key))


def meta(name):
    return read_table(f"tables/{name}.csv")[0]


def scalar(name, quantity):
    return row(name, quantity=quantity)["value"]


# (macro, value-thunk, format, description, [(file, approx line, literal currently in the text)])
REGISTRY = [
    ("RefTrenchDepthKm", lambda: row("model_summary", model=REF, suite="suite1_strength")["w_T_m"] / 1e3, ".2f",
     "reference model trench depth [km]", [("sections/results.tex", 13, "3.23")]),
    ("RefPull", lambda: row("model_summary", model=REF, suite="suite1_strength")["dGPE_TN"], ".2f",
     "reference trench pull ΔGPE* [TN/m], 2 dp", [("sections/results.tex", 14, "2.54"), ("si.tex", 886, "2.54")]),
    ("RefPullThree", lambda: row("model_summary", model=REF, suite="suite1_strength")["dGPE_TN"], ".3f",
     "reference trench pull ΔGPE* [TN/m], 3 dp", [("si.tex", 232, "2.542")]),
    ("RefIdentityPct", lambda: row("model_summary", model=REF, suite="suite1_strength")["identity_pct"], ".3f",
     "reference |ΔN_D − ΔGPE*|/ΔGPE* [%]", [("si.tex", 232, "0.019")]),
    ("RefArmKm", lambda: row("model_summary", model=REF, suite="suite1_strength")["arm_km"], ".1f",
     "reference effective arm ΔGPE*/(Δρ g w_T) [km]", []),
    ("RefXIkm", lambda: row("model_summary", model=REF, suite="suite1_strength")["x_I_km"], ".0f", "reference first isostatic column [km from trench]", []),
    ("RefXMkm", lambda: row("model_summary", model=REF, suite="suite1_strength")["x_M_km"], ".0f", "reference max-moment column [km from trench]", []),
    ("RefYieldThicknessPct", lambda: row("model_summary", model=REF, suite="suite1_strength")["yielded_thickness_xM_pct"], ".0f", "reference: yielded share of the thickness at max M [%]", []),
    ("PullElastic", lambda: row("model_summary", model="elastic_deep_60km_V4")["dGPE_TN"], ".2f", "Suite 1 elastic pull [TN/m]", [("sections/results.tex", 132, "2.23")]),
    ("PullDDVMasym", lambda: row("model_summary", model="dd_vm_asym_60km_V4")["dGPE_TN"], ".2f", "Suite 1 DD-VM asym pull [TN/m]", []),
    ("PullDDVMsym", lambda: row("model_summary", model="dd_vm_sym_60km_V4")["dGPE_TN"], ".2f", "Suite 1 DD-VM sym pull [TN/m]", [("sections/results.tex", 132, "2.84")]),
    ("DepthElasticKm", lambda: row("model_summary", model="elastic_deep_60km_V4")["w_T_m"] / 1e3, ".2f", "Suite 1 elastic trench depth [km]", [("sections/results.tex", 131, "2.90")]),
    ("DepthDDVMsymKm", lambda: row("model_summary", model="dd_vm_sym_60km_V4")["w_T_m"] / 1e3, ".2f", "Suite 1 DD-VM sym trench depth [km]", [("sections/results.tex", 130, "3.35")]),
    ("PullHthirty", lambda: row("thickness_compare", h_km=30.0)["pull_TN"], ".2f", "Suite 4 h=30 km pull [TN/m]", [("sections/results.tex", 167, "1.40"), ("si.tex", 900, "1.40")]),
    ("PullHforty", lambda: row("thickness_compare", h_km=40.0)["pull_TN"], ".2f", "Suite 4 h=40 km pull [TN/m]", []),
    ("PullHfifty", lambda: row("thickness_compare", h_km=50.0)["pull_TN"], ".2f", "Suite 4 h=50 km pull [TN/m]", []),
    ("ThicknessSlope", lambda: meta("thickness_compare")["through_origin_slope_dGPE_over_V"], ".2f", "Suite 4 through-origin ΔGPE*/V", [("sections/results.tex", 183, "0.65")]),
    ("BoefEndResultant", lambda: scalar("benchmark_boef", "end_resultant_from_FE"), ".3f", "S1 end resultant [TN/m]", [("si.tex", 130, "4.000"), ("si.tex", 148, "4.000"), ("si.tex", 259, "4.000")]),
    ("BoefVmisfitPct", lambda: scalar("benchmark_boef", "V_misfit_max_0_4alpha"), ".1f", "S1 max V misfit [%]", [("si.tex", 130, "0.8"), ("si.tex", 148, "0.8")]),
    ("BoefOffsetPct", lambda: scalar("benchmark_boef", "w0_offset_vs_thin_beam"), ".1f", "S1 thick-plate deflection offset [%]", [("si.tex", 106, "1.6"), ("si.tex", 127, "1.6")]),
    ("ParabolaMisfitPct", lambda: scalar("benchmark_boef", "shear_parabola_misfit_max"), ".2f", "S1 shear-parabola max deviation [%]", [("si.tex", 132, "0.39"), ("si.tex", 149, "0.39")]),
    ("MpSections", lambda: scalar("benchmark_mp", "sections_on_curve"), ".0f", "S2 sections on the curve", [("si.tex", 178, "148"), ("si.tex", 192, "148")]),
    ("MpMeanPct", lambda: scalar("benchmark_mp", "mean_misfit_M_over_Mp"), ".2f", "S2 mean |M/Mp − analytic| [%]", [("si.tex", 179, "0.09")]),
    ("MpMaxPct", lambda: scalar("benchmark_mp", "max_misfit_M_over_Mp"), ".2f", "S2 max |M/Mp − analytic| [%]", [("si.tex", 180, "0.85")]),
    ("ConvCoarsePull", lambda: row("convergence", model="bench_400x24")["dGPE_TN"], ".3f", "Table S3 400×24 pull", [("si.tex", 231, "2.550")]),
    ("ConvCoarseResid", lambda: row("convergence", model="bench_400x24")["identity_residual_pct"], ".3f", "Table S3 400×24 residual", [("si.tex", 231, "0.097")]),
    ("ConvFinePull", lambda: row("convergence", model="bench_1200x72")["dGPE_TN"], ".3f", "Table S3 1200×72 pull", [("si.tex", 233, "2.540")]),
    ("ConvFineResid", lambda: row("convergence", model="bench_1200x72")["identity_residual_pct"], ".3f", "Table S3 1200×72 residual", [("si.tex", 233, "0.010")]),
    ("ConvLoadIncPct", lambda: 100 * abs(row("convergence", model="bench_800x48_ns12")["dGPE_TN"] / row("convergence", configuration="24 (baseline)")["dGPE_TN"] - 1), ".3f",
     "Table S3: pull change 12 vs 24 load increments [%] (was 0.04 % from the 3-decimal transcription; full precision gives 0.005 %)", [("si.tex", 214, "0.04")]),
    ("ReconTopTrescaPct", lambda: row("gpe_compare_reconstruction", label="Tresca (uniform)")["plate_top_arm_err_pct"], ".1f", "Fig 5 Tresca plate-top-arm reconstruction error [%]", []),
    ("ReconSeaTrescaPct", lambda: row("gpe_compare_reconstruction", label="Tresca (uniform)")["sea_level_arm_err_pct"], ".1f", "Fig 5 Tresca sea-level-arm error [%]", []),
    ("ReconTopMaxPct", lambda: max(r["plate_top_arm_err_pct"] for r in read_table("tables/gpe_compare_reconstruction.csv")[1]), ".1f", "Fig 5 max plate-top-arm error, Suite 1 [%]", []),
    ("ReconSeaMaxPct", lambda: max(r["sea_level_arm_err_pct"] for r in read_table("tables/gpe_compare_reconstruction.csv")[1]), ".1f", "Fig 5 max sea-level-arm error, Suite 1 [%]", []),
    ("IsoRefPct", lambda: abs(row("isostatic_column", model=REF, suite="suite1_strength")["change_pct"]), ".2f",
     "isostatic-column assumption: |change| in the reference pull [%]", [("si.tex", 285, "0.04")]),
    ("IsoSuiteOneMaxPct", lambda: max(abs(r["change_pct"]) for r in read_table("tables/isostatic_column.csv")[1] if r["suite"] == "suite1_strength"), ".1f",
     "isostatic-column assumption: max |change| across Suite 1 [%]", []),
    ("IsoAllMaxPct", lambda: max(abs(r["change_pct"]) for r in read_table("tables/isostatic_column.csv")[1]), ".1f",
     "isostatic-column assumption: max |change| over all production models [%]", []),
    ("SfivePullTrench", lambda: row("corrected_density", column="trench")["dGPE_vs_isostatic_TN"], ".2f", "Fig S5 ΔGPE at the trench vs isostatic [TN/m]", [("si.tex", 886, "2.54")]),
    ("SfivePullMaxM", lambda: row("corrected_density", column="max M")["dGPE_vs_isostatic_TN"], ".2f", "Fig S5 ΔGPE at max M vs isostatic [TN/m]", [("si.tex", 887, "0.53")]),
    ("SfiveArmTrenchKm", lambda: row("corrected_density", column="trench")["dipole_arm_km"], ".0f", "Fig S5 dipole arm at the trench [km]", []),
    ("SfiveArmMaxMKm", lambda: row("corrected_density", column="max M")["dipole_arm_km"], ".0f", "Fig S5 dipole arm at max M [km]", []),
    ("SfiveMidPlateKm", lambda: meta("corrected_density")["deflected_mid_plate_km"], ".0f", "Fig S5 deflected mid-plate depth [km]", [("si.tex", 871, "33")]),
    ("SfiveBaseKm", lambda: meta("corrected_density")["deflected_base_km"], ".0f", "Fig S5 deflected base depth [km]", [("si.tex", 871, "63")]),
    ("HeroCentroidMeanKm", lambda: scalar("hero_tresca_deep60", "rhohat_centroid_depth_mean_over_window"), ".1f", "Fig 3 ρ̂ centroid depth, window mean [km]", []),
    # --- rounded / derived forms the text uses (from the same tables) ---
    ("RefPullRounded", lambda: row("model_summary", model=REF, suite="suite1_strength")["dGPE_TN"], ".1f", "reference pull, 1 dp (\"about 2.5\")",
     [("main.tex", 47, "2.5"), ("main.tex", 84, "2.5"), ("main.tex", 132, "2.5"), ("main.tex", 183, "2.5"), ("sections/discussion.tex", 89, "2.5"), ("sections/discussion.tex", 134, "2.5"), ("sections/discussion.tex", 145, "2.5")]),
    ("RefTrenchDepthRounded", lambda: row("model_summary", model=REF, suite="suite1_strength")["w_T_m"] / 1e3, ".1f", "reference trench depth, 1 dp (\"3.2-km trench\")",
     [("main.tex", 84, "3.2"), ("main.tex", 183, "3.2"), ("sections/model_setup.tex", 173, "3.2"), ("sections/discussion.tex", 89, "3.2")]),
    ("RefPullPerKm", lambda: row("model_summary", model=REF, suite="suite1_strength")["dGPE_TN"] / (row("model_summary", model=REF, suite="suite1_strength")["w_T_m"] / 1e3), ".1f",
     "reference pull per km of deflection [TN/m per km]", [("sections/discussion.tex", 88, "0.8")]),
    ("RefArmRoundedKm", lambda: row("model_summary", model=REF, suite="suite1_strength")["arm_km"], ".0f", "reference effective arm, rounded [km] (\"approximately 35\")",
     [("sections/model_setup.tex", 163, "35"), ("sections/results.tex", 57, "35"), ("sections/results.tex", 80, "35")]),
    ("RefArmOverH", lambda: row("model_summary", model=REF, suite="suite1_strength")["arm_over_h"], ".2f", "reference arm / h (\"approximately 0.6h\")", [("sections/results.tex", 81, "0.6"), ("sections/results.tex", 171, "0.6")]),
    ("ArmOverHMin", lambda: min(r["arm_over_h"] for r in read_table("tables/model_summary.csv")[1]), ".2f", "min arm/h over all production models (the 0.55 ± 0.10 coefficient)", [("sections/discussion.tex", 133, "0.55")]),
    ("ArmOverHMax", lambda: max(r["arm_over_h"] for r in read_table("tables/model_summary.csv")[1]), ".2f", "max arm/h over all production models", []),
    ("PullSuiteTwoMin", lambda: min(r["dGPE_TN"] for r in read_table("tables/model_summary.csv")[1] if r["suite"] == "suite2_load"), ".2f", "Suite 2 smallest pull [TN/m]", [("sections/results.tex", 44, "0.51")]),
    ("PullSuiteTwoMax", lambda: max(r["dGPE_TN"] for r in read_table("tables/model_summary.csv")[1] if r["suite"] == "suite2_load"), ".2f", "Suite 2 largest pull [TN/m]", [("sections/results.tex", 44, "3.20")]),
    ("PullMemMinusThree", lambda: row("model_summary", model="tresca_deep_150_60km_V4_mem-3")["dGPE_TN"], ".2f", "Suite 3 pull at N_D = −3 TN/m", [("sections/results.tex", 67, "3.31")]),
    ("TensionBackgroundMaxChangePct", lambda: max(abs(100 * (row("model_summary", model=f"tresca_deep_150_60km_V4_mem{k}")["dGPE_TN"] / row("model_summary", model=REF, suite="suite1_strength")["dGPE_TN"] - 1)) for k in (1, 2, 3)), ".1f",
     "largest |change| in pull for tension-like backgrounds +1…+3 TN/m [%] (\"about 1 %\")", [("sections/results.tex", 66, "1"), ("main.tex", 186, "1")]),
    ("SpreadFromTrendPct", lambda: max(abs(100 * (r["dGPE_TN"] / (row("model_summary", model=REF, suite="suite1_strength")["dGPE_TN"] / row("model_summary", model=REF, suite="suite1_strength")["w_T_m"] * r["w_T_m"]) - 1))
                                     for r in read_table("tables/model_summary.csv")[1] if r["suite"] in ("suite1_strength", "suite3_background")), ".1f",
     "largest departure of Suites 1 and 3 from the reference pull/deflection ratio [%] (\"within about 10 %\")",
     [("sections/results.tex", 62, "10"), ("sections/results.tex", 132, "10"), ("sections/discussion.tex", 91, "10"), ("main.tex", 82, "10"), ("main.tex", 185, "10")]),
    ("MidPlateApproxMinPct", lambda: min(r["mid_plate_approx_err_pct"] for r in read_table("tables/gpe_compare_reconstruction.csv")[1]), ".0f", "Fig 5: smallest h/2 under-estimate [%] (\"10–20 %\")", [("sections/results.tex", 101, "10")]),
    ("MidPlateApproxMaxPct", lambda: max(r["mid_plate_approx_err_pct"] for r in read_table("tables/gpe_compare_reconstruction.csv")[1]), ".0f", "Fig 5: largest h/2 under-estimate [%]", [("sections/results.tex", 101, "20")]),
    ("ReconTopMinPct", lambda: min(r["plate_top_arm_err_pct"] for r in read_table("tables/gpe_compare_reconstruction.csv")[1]), ".1f", "Fig 5 min plate-top-arm error, Suite 1 [%] (\"1.7–2.3 %\")", [("sections/results.tex", 98, "1.7"), ("sections/results.tex", 110, "1.7")]),
    ("ProdVsFinePct", lambda: 100 * abs(row("convergence", configuration="800 × 48 (baseline)")["dGPE_TN"] / row("convergence", model="bench_1200x72")["dGPE_TN"] - 1), ".2f", "production vs 1200×72 pull [%] (\"within 0.08 %\")", [("si.tex", 212, "0.08")]),
    ("HardeningIncrementMPa", lambda: max(r["hardening_H_Pa"] for r in read_table("tables/model_summary.csv")[1]) * 2e-3 / 1e6, ".1f", "hardening increment H·ε_p at the text's representative ε_p = 2×10⁻³ [MPa]", [("si/model_formulation.tex", 85, "1.4")]),
    ("HardeningIncrementPct", lambda: 100 * max(r["hardening_H_Pa"] for r in read_table("tables/model_summary.csv")[1]) * 2e-3 / 1e6 / 150, ".1f", "the same as % of 150 MPa (\"less than 1 %\")", [("si/model_formulation.tex", 86, "1")]),
    ("HardeningMaxIncrementMPa", lambda: max(r["hardening_increment_MPa"] for r in read_table("tables/model_summary.csv")[1]), ".1f", "hardening increment at the ACTUAL largest plastic strain in the DD-VM models [MPa]", []),
    ("MaxPlasticStrain", lambda: max(r["max_plastic_strain"] for r in read_table("tables/model_summary.csv")[1]), ".4f", "largest plastic strain in any production model", []),
    ("ElasticWZeroM", lambda: scalar("benchmark_boef", "w0_FE"), ".0f", "S1 elastic end deflection [m]", [("si.tex", 146, "2899")]),
    ("HoverAlpha", lambda: scalar("benchmark_boef", "h_over_alpha"), ".2f", "S1 slenderness h/α", [("si.tex", 109, "0.48")]),
    ("CumAvgDensityBase", lambda: meta("corrected_density")["cumulative_average_density_at_base_kg_m3"], ".0f", "Fig S5 common cumulative-average density at the base [kg/m³]", [("si.tex", 877, "3302")]),
    ("PullDDVMasymTwo", lambda: row("model_summary", model="dd_vm_asym_60km_V4")["dGPE_TN"], ".2f", "Suite 1 DD-VM asym pull [TN/m] (Fig S7 caption)", [("si.tex", 914, "2.29")]),
    ("IsoRatioDipolePct", lambda: 100 * abs(row("isostatic_column", model=REF, suite="suite1_strength")["Szz_xI_TN"]) / abs(row("model_summary", model=REF, suite="suite1_strength")["dGPE_TN"]), ".2f",
     "|Σzz(x_I)| as % of the trench dipole moment (\"0.1 % of the dipole moment\")", [("si.tex", 286, "0.1")]),
    # --- frame check (SI Text S2) ---
    ("FrameMaxSlopeDeg", lambda: max(r["max_slope_deg"] for r in read_table("tables/frame_check.csv")[1]), ".1f", "largest surface slope over all production models [deg] (\"below 2.5°\")", [("si.tex", 301, "2.5")]),
    ("FrameMaxCosDepPct", lambda: max(r["max_cos2theta_departure_pct"] for r in read_table("tables/frame_check.csv")[1]), ".1f", "largest 1 − cos 2θ over all models [%] (\"less than 0.3 %\")", [("si.tex", 302, "0.3")]),
    ("FrameMixTrenchRefTN", lambda: row("frame_check", model=REF, suite="suite1_strength")["mixing_term_at_trench_TN"], ".2f", "2V sin 2θ at the trench column, reference model [TN/m] (\"≈0.45\")", [("si.tex", 304, "0.45")]),
    ("FrameMixDecayRefKm", lambda: row("frame_check", model=REF, suite="suite1_strength")["mixing_first_below_0p02_km"], ".0f", "distance inboard at which 2V sin 2θ first falls below 0.02 TN/m, reference [km] (\"within about 50 km\")", [("si.tex", 305, "50")]),
    ("FrameMixStaysBelowRefKm", lambda: row("frame_check", model=REF, suite="suite1_strength")["mixing_below_0p02_beyond_km"], ".0f", "distance inboard beyond which 2V sin 2θ STAYS below 0.02 TN/m, reference [km]", []),
    ("FrameSlopeXIRefDeg", lambda: row("frame_check", model=REF, suite="suite1_strength")["slope_at_xI_deg"], ".1f", "surface slope at x_I, reference [deg] (\"about 0.2°\")", [("si.tex", 307, "0.2")]),
    ("FrameMixXIRefTN", lambda: abs(row("frame_check", model=REF, suite="suite1_strength")["mixing_term_at_xI_TN"]), ".2f", "|2V sin 2θ| at x_I, reference [TN/m] (\"≈0.01\")", [("si.tex", 307, "0.01")]),
    # --- loaded-edge exclusion (SI Text S1) ---
    ("EdgeChangeMinPct", lambda: min(abs(r["change_raw_to_plateau_pct"]) for r in read_table("tables/edge_exclusion.csv")[1]), ".0f", "smallest |arm change| when the edge window is excluded, Suite 1 [%] (\"4–8 %\")", [("si/model_formulation.tex", 49, "4")]),
    ("EdgeChangeMaxPct", lambda: max(abs(r["change_raw_to_plateau_pct"]) for r in read_table("tables/edge_exclusion.csv")[1]), ".0f", "largest |arm change| when the edge window is excluded, Suite 1 [%]", [("si/model_formulation.tex", 49, "8")]),
]


def fmt(v, f):
    return format(v, f)


def main(compare=None):
    os.makedirs("tables", exist_ok=True)
    tex = ["% tables/paper_numbers.tex — GENERATED by analysis/paper_numbers.py from tables/*.csv (reproduce.sh). Do not edit.",
           "% Usage: \\input{paper_numbers}  then e.g.  a trench pull of \\numRefPull~TN/m.  Regenerate: ./reproduce.sh, re-copy."]
    md = ["# Paper numbers — placement guide", "",
          "Generated by `python analysis/paper_numbers.py` from `tables/*.csv` (the source of truth), on every `reproduce.sh` run.",
          "`tables/paper_numbers.tex` defines one macro per number at the paper's precision: copy it next to the LaTeX sources,",
          "`\\input{paper_numbers}` in the preamble, and write the macro instead of the literal. A number that is not in this",
          "list is not a supported model result. 'Where' is the location in `2026_codex/full_manuscript` as of 2026-09-13 and the",
          "literal the text carried then; `--compare DIR` re-checks those locations against the tables.", "",
          "| macro | value | source | description | where (file:line, literal in text) | status |", "|---|---|---|---|---|---|"]
    report = []
    import inspect
    bad = [n for n, *_ in REGISTRY if not re.fullmatch(r"[A-Za-z]+", n)]
    assert not bad, f"LaTeX macro names must be letters only (no digits/underscores): {bad}"
    assert len({n for n, *_ in REGISTRY}) == len(REGISTRY), "duplicate macro name"
    for name, thunk, f, desc, where in REGISTRY:
        v = fmt(thunk(), f)
        tex.append(f"\\newcommand{{\\num{name}}}{{{v}}}")
        m = re.search(r'(?:row|scalar|meta)\("(\w+)"|read_table\("tables/(\w+)\.csv"', inspect.getsource(thunk))
        src = f"`tables/{(m.group(1) or m.group(2))}.csv`" if m else "—"
        if compare and where:
            parts = []
            for file, line, lit in where:
                path = os.path.join(compare, file)
                if not os.path.isfile(path):
                    parts.append("file missing"); continue
                lines = open(path, encoding="utf-8").read().splitlines()
                hit_lit = [i + 1 for i, l in enumerate(lines) if lit in l]
                hit_val = [i + 1 for i, l in enumerate(lines) if v in l]
                hit_mac = [i + 1 for i, l in enumerate(lines) if f"\\num{name}" in l]
                near = lambda hits: [h for h in hits if abs(h - line) <= 3]
                if hit_mac:
                    parts.append(f"integrated (macro at line {hit_mac[0]})")
                elif v == lit:
                    parts.append("match" if near(hit_lit) else (f"match (now at line {hit_lit[0]})" if hit_lit else "literal gone"))
                elif _rounded_match(v, lit):
                    parts.append(f"rounded match (text has {lit})" if hit_lit else "literal gone")
                else:
                    parts.append(f"DIFFERS: text has {lit}" if hit_lit else (f"updated (line {hit_val[0]})" if hit_val else "literal gone"))
            status = "; ".join(parts)
            report.append((name, v, where, status))
        else:
            status = "not compared (run with --compare DIR)" if where else "not yet quoted"
        wh = "; ".join(f"`{fl}:{ln}` ({lit})" for fl, ln, lit in where) or "—"
        md.append(f"| `\\num{name}` | {v} | {src} | {desc} | {wh} | {status} |")
    with open("tables/paper_numbers.tex", "w") as fh:
        fh.write("\n".join(tex) + "\n")
    with open("tables/PAPER_NUMBERS.md", "w") as fh:
        fh.write("\n".join(md) + "\n")
    print(f"wrote tables/paper_numbers.tex ({len(REGISTRY)} macros) and tables/PAPER_NUMBERS.md")
    if compare:
        return compare_report(compare, report)
    return 0


def _rounded_match(v, lit):
    """True when the manuscript literal is the table value quoted to fewer decimals (2.54 vs 2.542, 35 vs 34.9)."""
    try:
        a, b = float(v), float(lit)
    except ValueError:
        return False
    d = max(len(lit.split(".")[1]) if "." in lit else 0, 0)
    return round(a, d) == round(b, d) or abs(a - b) <= 0.5 * 10 ** (-d) + 1e-12


def compare_report(compare, report):
    """Print the comparison in three blocks (ok / rounded / DIFFERS) and return a nonzero exit code when anything
    needs a human: a missing manuscript file, a literal that is gone, or a genuine difference.  Also diffs the
    manuscript's copied paper_numbers.tex against the release file."""
    OK = ("match", "integrated"); SOFT = ("rounded match", "updated")
    ok, soft, bad = [], [], []
    for name, v, where, status in report:
        parts = status.split("; ")
        if all(p.startswith(OK) for p in parts):       ok.append((name, v, status))
        elif all(p.startswith(OK + SOFT) for p in parts): soft.append((name, v, status))
        else:                                            bad.append((name, v, status))
    print("\ncomparison against", compare)
    print(f"  ok      {len(ok):3d}  (literal matches the table, or the macro is already in the text)")
    print(f"  rounded {len(soft):3d}  (text quotes the table value to fewer decimals, or the literal moved)")
    for name, v, status in soft:
        print(f"           \\num{name} = {v:>8}   {status}")
    print(f"  DIFFERS {len(bad):3d}  (a genuine difference, a missing file, or a literal that is gone — needs a human)")
    for name, v, status in bad:
        print(f"           \\num{name} = {v:>8}   {status}")
    rc = 1 if bad else 0
    ms = os.path.join(compare, "paper_numbers.tex")
    if not os.path.isfile(ms):
        print(f"  manuscript copy of paper_numbers.tex: MISSING at {ms}"); rc = 1
    else:
        rel = {l.split("}")[0] for l in open("tables/paper_numbers.tex") if l.startswith("\\newcommand")}
        man = {l.split("}")[0] for l in open(ms) if l.startswith("\\newcommand")}
        relv = dict(re.findall(r"\\newcommand\{(\\num\w+)\}\{([^}]*)\}", open("tables/paper_numbers.tex").read()))
        manv = dict(re.findall(r"\\newcommand\{(\\num\w+)\}\{([^}]*)\}", open(ms).read()))
        missing = sorted(set(relv) - set(manv)); stale = sorted(k for k in relv if k in manv and relv[k] != manv[k])
        if not missing and not stale:
            print("  manuscript copy of paper_numbers.tex: current (every macro present with the release value)")
        else:
            rc = 1
            if missing: print(f"  manuscript copy of paper_numbers.tex: {len(missing)} macro(s) missing — re-copy tables/paper_numbers.tex: {', '.join(missing[:6])}{' …' if len(missing) > 6 else ''}")
            if stale:   print(f"  manuscript copy of paper_numbers.tex: {len(stale)} macro(s) STALE (value differs from the release): {', '.join(stale[:6])}")
    print("  exit", rc, "(0 = the manuscript is consistent with this release; 1 = something above needs attention)")
    return rc


if __name__ == "__main__":
    cmp = sys.argv[sys.argv.index("--compare") + 1] if "--compare" in sys.argv else None
    sys.exit(main(cmp) or 0)

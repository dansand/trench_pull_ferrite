"""paper_numbers.py — the headline numbers as LaTeX macros, read from tables/*.csv.

Writes tables/paper_numbers.tex (one \\newcommand per macro, at the precision the paper uses) and tables/PAPER_NUMBERS.md
(macro, value, source table, description).  The CSV tables are the authoritative numbers; the macros exist so the
manuscript's headline values come from the tables rather than being retyped.  Copy tables/paper_numbers.tex next to the
LaTeX sources and \\input it.  Every other number in the manuscript is checked by hand against the tables once before
submission.  Run by reproduce.sh.

    python analysis/paper_numbers.py

The registry below (macro, table selector, format, description) is the only hand-maintained part; a macro's value is
never typed.
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


# (macro, value-thunk, format, description) — the macros the manuscript uses
REGISTRY = [
    ("RefTrenchDepthKm", lambda: row("model_summary", model=REF, suite="suite1_strength")["w_T_m"] / 1e3, ".2f",
     "reference model trench depth [km]"),
    ("RefPull", lambda: row("model_summary", model=REF, suite="suite1_strength")["dGPE_TN"], ".2f",
     "reference trench pull ΔGPE* [TN/m], 2 dp"),
    ("RefPullThree", lambda: row("model_summary", model=REF, suite="suite1_strength")["dGPE_TN"], ".3f",
     "reference trench pull ΔGPE* [TN/m], 3 dp"),
    ("RefIdentityPct", lambda: row("model_summary", model=REF, suite="suite1_strength")["identity_pct"], ".3f",
     "reference |ΔN_D − ΔGPE*|/ΔGPE* [%]"),
    ("RefXIkm", lambda: row("model_summary", model=REF, suite="suite1_strength")["x_I_km"], ".0f",
     "reference first isostatic column [km from trench]"),
    ("WedgeAbovePlateTN", lambda: row("profiles_selfcheck", model="data/suite1_strength/tresca_deep_150_60km_V4")["above_plate_wedge_TN"], ".2f",
     "above-plate (water-wedge) term ½Δρ g w_T² omitted by the massless column, reference model [TN/m]"),
    ("WedgeAbovePlatePct", lambda: row("profiles_selfcheck", model="data/suite1_strength/tresca_deep_150_60km_V4")["above_plate_wedge_pct"], ".1f",
     "the same as % of the reference pull"),
    ("PullElastic", lambda: row("model_summary", model="elastic_deep_60km_V4")["dGPE_TN"], ".2f",
     "Suite 1 elastic pull [TN/m]"),
    ("PullDDVMasym", lambda: row("model_summary", model="dd_vm_asym_60km_V4")["dGPE_TN"], ".2f",
     "Suite 1 DD-VM asym pull [TN/m]"),
    ("PullDDVMsym", lambda: row("model_summary", model="dd_vm_sym_60km_V4")["dGPE_TN"], ".2f",
     "Suite 1 DD-VM sym pull [TN/m]"),
    ("DepthElasticKm", lambda: row("model_summary", model="elastic_deep_60km_V4")["w_T_m"] / 1e3, ".2f",
     "Suite 1 elastic trench depth [km]"),
    ("DepthDDVMsymKm", lambda: row("model_summary", model="dd_vm_sym_60km_V4")["w_T_m"] / 1e3, ".2f",
     "Suite 1 DD-VM sym trench depth [km]"),
    ("PullHthirty", lambda: row("thickness_compare", h_km=30.0)["pull_TN"], ".2f",
     "Suite 4 h=30 km pull [TN/m]"),
    ("ThicknessSlope", lambda: meta("thickness_compare")["through_origin_slope_dGPE_over_V"], ".2f",
     "Suite 4 through-origin ΔGPE*/V"),
    ("BoefEndResultant", lambda: scalar("benchmark_boef", "end_resultant_from_FE"), ".3f",
     "S1 end resultant [TN/m]"),
    ("BoefVmisfitPct", lambda: scalar("benchmark_boef", "V_misfit_max_0_4alpha"), ".1f",
     "S1 max V misfit [%]"),
    ("BoefOffsetPct", lambda: scalar("benchmark_boef", "w0_offset_vs_thin_beam"), ".1f",
     "S1 thick-plate deflection offset [%]"),
    ("ParabolaMisfitPct", lambda: scalar("benchmark_boef", "shear_parabola_misfit_max"), ".2f",
     "S1 shear-parabola max deviation [%]"),
    ("MpSections", lambda: scalar("benchmark_mp", "sections_on_curve"), ".0f",
     "S2 sections on the curve"),
    ("MpMeanPct", lambda: scalar("benchmark_mp", "mean_misfit_M_over_Mp"), ".2f",
     "S2 mean |M/Mp − analytic| [%]"),
    ("MpMaxPct", lambda: scalar("benchmark_mp", "max_misfit_M_over_Mp"), ".2f",
     "S2 max |M/Mp − analytic| [%]"),
    ("ConvCoarsePull", lambda: row("convergence", model="bench_400x24")["dGPE_TN"], ".3f",
     "Table S3 400×24 pull"),
    ("ConvCoarseResid", lambda: row("convergence", model="bench_400x24")["identity_residual_pct"], ".3f",
     "Table S3 400×24 residual"),
    ("ConvFinePull", lambda: row("convergence", model="bench_1200x72")["dGPE_TN"], ".3f",
     "Table S3 1200×72 pull"),
    ("ConvFineResid", lambda: row("convergence", model="bench_1200x72")["identity_residual_pct"], ".3f",
     "Table S3 1200×72 residual"),
    ("ConvLoadIncPct", lambda: 100 * abs(row("convergence", model="bench_800x48_ns12")["dGPE_TN"] / row("convergence", configuration="24 (baseline)")["dGPE_TN"] - 1), ".3f",
     "Table S3: pull change 12 vs 24 load increments [%] (was 0.04 % from the 3-decimal transcription; full precision gives 0.005 %)"),
    ("SfivePullTrench", lambda: row("corrected_density", column="trench")["dGPE_vs_isostatic_TN"], ".2f",
     "Fig S5 ΔGPE at the trench vs isostatic [TN/m]"),
    ("SfivePullMaxM", lambda: row("corrected_density", column="max M")["dGPE_vs_isostatic_TN"], ".2f",
     "Fig S5 ΔGPE at max M vs isostatic [TN/m]"),
    ("SfiveMidPlateKm", lambda: meta("corrected_density")["deflected_mid_plate_km"], ".0f",
     "Fig S5 deflected mid-plate depth [km]"),
    ("SfiveBaseKm", lambda: meta("corrected_density")["deflected_base_km"], ".0f",
     "Fig S5 deflected base depth [km]"),
    ("RefPullRounded", lambda: row("model_summary", model=REF, suite="suite1_strength")["dGPE_TN"], ".1f",
     "reference pull, 1 dp (\"about 2.5\")"),
    ("RefArmRoundedKm", lambda: row("model_summary", model=REF, suite="suite1_strength")["arm_km"], ".0f",
     "reference effective arm, rounded [km] (\"approximately 35\")"),
    ("RefArmOverH", lambda: row("model_summary", model=REF, suite="suite1_strength")["arm_over_h"], ".2f",
     "reference arm / h (\"approximately 0.6h\")"),
    ("ArmOverHMin", lambda: min(r["arm_over_h"] for r in read_table("tables/model_summary.csv")[1]), ".2f",
     "min arm/h over all production models (the 0.55 ± 0.10 coefficient)"),
    ("ArmOverHMax", lambda: max(r["arm_over_h"] for r in read_table("tables/model_summary.csv")[1]), ".2f",
     "max arm/h over all production models"),
    ("MidPlateApproxMinPct", lambda: min(r["mid_plate_approx_err_pct"] for r in read_table("tables/gpe_compare_reconstruction.csv")[1]), ".0f",
     "Fig 5: smallest h/2 under-estimate [%] (\"10–20 %\")"),
    ("MidPlateApproxMaxPct", lambda: max(r["mid_plate_approx_err_pct"] for r in read_table("tables/gpe_compare_reconstruction.csv")[1]), ".0f",
     "Fig 5: largest h/2 under-estimate [%]"),
    ("FrameMaxSlopeDeg", lambda: max(r["max_slope_deg"] for r in read_table("tables/frame_check.csv")[1]), ".1f",
     "largest surface slope over all production models [deg] (\"below 2.5°\")"),
    ("FrameMaxCosDepPct", lambda: max(r["max_cos2theta_departure_pct"] for r in read_table("tables/frame_check.csv")[1]), ".1f",
     "largest 1 − cos 2θ over all models [%] (\"less than 0.3 %\")"),
    ("FrameMixTrenchRefTN", lambda: row("frame_check", model=REF, suite="suite1_strength")["mixing_term_at_trench_TN"], ".2f",
     "2V sin 2θ at the trench column, reference model [TN/m] (\"≈0.45\")"),
    ("FrameMixDecayRefKm", lambda: row("frame_check", model=REF, suite="suite1_strength")["mixing_first_below_0p02_km"], ".0f",
     "distance inboard at which 2V sin 2θ first falls below 0.02 TN/m, reference [km] (\"within about 50 km\")"),
    ("FrameMixStaysBelowRefKm", lambda: row("frame_check", model=REF, suite="suite1_strength")["mixing_below_0p02_beyond_km"], ".0f",
     "distance inboard beyond which 2V sin 2θ STAYS below 0.02 TN/m, reference [km]"),
    ("FrameSlopeXIRefDeg", lambda: row("frame_check", model=REF, suite="suite1_strength")["slope_at_xI_deg"], ".1f",
     "surface slope at x_I, reference [deg] (\"about 0.2°\")"),
    ("FrameMixXIRefTN", lambda: abs(row("frame_check", model=REF, suite="suite1_strength")["mixing_term_at_xI_TN"]), ".2f",
     "|2V sin 2θ| at x_I, reference [TN/m] (\"≈0.01\")"),
]


def main():
    os.makedirs("tables", exist_ok=True)
    bad = [n for n, *_ in REGISTRY if not re.fullmatch(r"[A-Za-z]+", n)]
    assert not bad, f"LaTeX macro names must be letters only: {bad}"
    assert len({n for n, *_ in REGISTRY}) == len(REGISTRY), "duplicate macro name"
    import inspect
    tex = ["% tables/paper_numbers.tex — GENERATED by analysis/paper_numbers.py from tables/*.csv (reproduce.sh). Do not edit.",
           "% Usage: \\input{paper_numbers}  then e.g.  a trench pull of \\numRefPull~TN/m."]
    md = ["# Paper numbers", "",
          "Generated by `python analysis/paper_numbers.py` from `tables/*.csv` (the authoritative numbers) on every",
          "`reproduce.sh` run. One LaTeX macro per headline value, at the paper's precision; copy `tables/paper_numbers.tex`",
          "next to the LaTeX sources and `\\input{paper_numbers}`. Numbers not listed here are read from the tables directly.", "",
          "| macro | value | source table | description |", "|---|---|---|---|"]
    for name, thunk, f, desc in REGISTRY:
        v = format(thunk(), f)
        tex.append(f"\\newcommand{{\\num{name}}}{{{v}}}")
        m = re.search(r'(?:row|scalar|meta)\("(\w+)"|read_table\("tables/(\w+)\.csv"', inspect.getsource(thunk))
        src = f"`tables/{(m.group(1) or m.group(2))}.csv`" if m else "—"
        md.append(f"| `\\num{name}` | {v} | {src} | {desc} |")
    with open("tables/paper_numbers.tex", "w") as fh:
        fh.write("\n".join(tex) + "\n")
    with open("tables/PAPER_NUMBERS.md", "w") as fh:
        fh.write("\n".join(md) + "\n")
    print(f"wrote tables/paper_numbers.tex ({len(REGISTRY)} macros) and tables/PAPER_NUMBERS.md")


if __name__ == "__main__":
    main()

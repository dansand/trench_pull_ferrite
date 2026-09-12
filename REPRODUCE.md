# REPRODUCE — every manuscript figure: script → data → command

The frozen figure set of *The 'trench pull' force: constraints from elasto-plastic bending models*. Every
figure regenerates from the data in `data/` without re-running the solver, except the S2 convergence table
(noted below). **Run every command from the repository root** — Python in the `ferrite-figs` env, Julia as
`julia --project=.` (see `README.md`). Each script writes its figure into `data/reference_figures/`, which ships
with the reference renders to compare against (output-equivalent, not byte-identical).

Data folders are named after the manuscript's suites (`data/suite1_strength` … `data/suite4_thickness`, `data/appendixB_esweep`);
the driver's *commands* stay descriptive: `v_sweep` writes `suite2_load`, `nd_sweep` writes `suite3_background`.

## Main text

| Fig | Manuscript file | Command (from repo root) | Data | Regenerate data |
|---|---|---|---|---|
| 1 | `ridge_trench_overview_v2.pdf` | `cd schematic && tectonic ridge_trench_overview_v2.tex` (hand TikZ) | — | — |
| 2 | `equilibration_compare.pdf` | `cd schematic && python gen_equilibration_compare.py && tectonic equilibration_compare.tex` | — (computed curves) | — |
| 3 | `hero_tresca_deep60.png` | `python figures/render_hero.py` | `data/suite1_strength/tresca_deep_150_60km_V4` | `model/paper_models.jl suite1` |
| 4 | `gpe_correlation.png` | `python figures/render_gpe_correlation.py` | `suite1_strength`, `suite3_background`, `suite2_load` | `suite1 nd_sweep v_sweep` |
| 5 | `gpe_compare_suite1.png` | `python figures/render_gpe_compare.py` | `suite1_strength/*` (4) | `suite1` |
| 6 | `profiles.png` | `python figures/render_profiles.py` | `suite1_strength/*` (4) | `suite1` |
| 7 | `ferrite_thickness_compare.png` | `python figures/render_thickness_compare.py` → `thickness_compare.png` (renamed on copy) | `suite1_strength/tresca_*`, `suite4_thickness/*` (3) | `suite1 thickness` |

## Supporting information

| SI | Manuscript file | Command (from repo root) | Data | Regenerate data |
|---|---|---|---|---|
| S1 | `benchmark_boef.png` | `python figures/render_benchmark.py` | `suite1_strength/elastic_deep_60km_V4`, `data/idealized_beam_elastic` | `suite1`; `model/idealized_beam.jl` |
| S2 | `benchmark_mp.png` | `python figures/render_mp_benchmark.py` | `data/idealized_beam` | `model/idealized_beam.jl` |
| S2 | convergence table | `python figures/render_convergence.py` | `data/convergence/` — **not included**; values in `data/CONVERGENCE.md` | **run first:** `julia --project=. model/paper_models.jl convergence` |
| S3 | `fig_taux_cases_grid.pdf` | `cd schematic && tectonic taux_cases.tex` → `taux_cases.pdf` (renamed on copy) | — | — |
| S4 | `core_profiles_deep60.png` | `python figures/render_core_profiles.py` | `suite1_strength/tresca_deep_150_60km_V4` | `suite1` |
| S5 | `corrected_density.png` | `python figures/render_corrected_density.py` | `suite1_strength/tresca_deep_150_60km_V4` | `suite1` |
| S6, S7 | `ferrite_hero_h40.png`, `ferrite_hero_h30.png` | `HERO_WINDOW_KM=260 HERO_MFIX=0 python figures/render_hero.py data/suite4_thickness/tresca_150_40km data/reference_figures/hero_tresca_40km.png` · `HERO_WINDOW_KM=200 python figures/render_hero.py data/suite4_thickness/tresca_150_30km data/reference_figures/hero_tresca_30km.png` (renamed on copy) | `suite4_thickness/tresca_150_{40,30}km` | `thickness` |
| S8, S9 | `hero_dd_vm_asym.png`, `hero_dd_vm_sym.png` | `python figures/render_hero.py data/suite1_strength/dd_vm_asym_60km_V4 data/reference_figures/hero_dd_vm_asym.png` (and `…sym…`) | `suite1_strength/dd_vm_*` | `suite1` |
| parked | `edge_effect_SI.png` — loaded-edge disclosure; in `si_parked.tex` only, not the live SI (kept: it backs the loaded-edge numbers the SI text quotes) | `python figures/render_edge_si.py` | `suite1_strength/*` (4) | `suite1` |
| App. B | (number only: E-sweep arm ≈ 35 km, §4.1) | `python figures/render_esweep_test.py` → `esweep.png` | `appendixB_esweep/*` + the three suites | `all` |

`figures/render_hero.py` is one parameterised renderer — the same figure layout for any model directory
(`render_hero.py MODEL_DIR OUT [DIFF_DIR]`; `HERO_WINDOW_KM`, `HERO_MFIX` set the plotted window for thin plates).
The hero-family SI figures are *model* variants, not different figure designs. From Python/Jupyter:
`import render_hero; fig = render_hero.build(model_dir, save=False)`.

## Records backing numbers quoted in the SI (not figures)
- `data/BOUNDARY_ARTIFACT.md` — loaded-edge sensitivity (29 % force- vs displacement-controlled split; raw
  trench 12–14 % above the interior; convergence within 40–60 km of the loaded edge).
- `data/CONVERGENCE.md` — the S2 convergence / consistency table values.
- `figures/render_gpe_compare.py` prints the plate-top-arm (~2 %) and sea-level-arm (~8 %) reconstruction errors.

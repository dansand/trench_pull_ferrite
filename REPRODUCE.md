# REPRODUCE — every manuscript figure: script → data → command

The frozen figure set of *The 'trench pull' force: constraints from elasto-plastic bending models*. Every
figure regenerates from the data in `data/` without re-running the solver. **Run every command from the repository root** — Python in the `ferrite-figs` env, Julia as
`julia --project=.` (see `README.md`). Each script writes its figure into `figures/`, which ships with the
reference renders to compare against (output-equivalent, not byte-identical).

**One command does all of it: `./reproduce.sh`** — every command below in order, the headline numbers asserted
(ΔGPE\* 2.542 TN/m, identity < 0.05 %, arm 34.9 km, the S1/S2 benchmark misfits, the reconstruction table),
`START_HERE.ipynb` executed top to bottom, `pytest tests/` — nonzero exit on any failure.

Data folders are named after the manuscript's suites (`data/suite1_strength` … `data/suite4_thickness`);
the driver's *commands* stay descriptive: `v_sweep` writes `suite2_load`, `nd_sweep` writes `suite3_background`. The
provenance of every shipped model (command, parameters, SHA-256 of every file, origin) is in `data/DATA_MANIFEST.md`.

## Main text

| Fig | Manuscript file | Command (from repo root) | Table written in the same pass | Data | Regenerate data |
|---|---|---|---|---|---|
| 1 | `ridge_trench_overview_v2.pdf` | `cd schematic && tectonic ridge_trench_overview_v2.tex` (hand TikZ) | — | — | — |
| 2 | `equilibration_compare.pdf` | `cd schematic && python gen_equilibration_compare.py && tectonic equilibration_compare.tex` | — | — (computed curves) | — |
| 3 | `hero_tresca_deep60.png` (four panels: topography, σxx−σzz with crosses, ρ̂, resultants) | `python scripts/render_hero.py` | `hero_tresca_deep60.csv` | `data/suite1_strength/tresca_deep_150_60km_V4` | `model/paper_models.jl suite1` |
| 4 | `gpe_correlation.png` | `python scripts/render_gpe_correlation.py` | `gpe_correlation.csv` | `suite1_strength`, `suite3_background`, `suite2_load` | `suite1 nd_sweep v_sweep` |
| 5 | `gpe_compare_suite1.png` | `python scripts/render_gpe_compare.py` | `gpe_compare_reconstruction.csv` | `suite1_strength/*` (4) | `suite1` |
| 6 | `profiles.png` | `python scripts/render_profiles.py` | `profiles_selfcheck.csv` | `suite1_strength/*` (4) | `suite1` |
| 7 | `ferrite_thickness_compare.png` | `python scripts/render_thickness_compare.py` → `thickness_compare.png` (renamed on copy) | `thickness_compare.csv` | `suite1_strength/tresca_*`, `suite4_thickness/*` (3) | `suite1 thickness` |

The three TikZ schematics (Figs 1, 2, S3) are drawings, not model output, and `reproduce.sh` does not build them:
`cd schematic && python gen_equilibration_compare.py && tectonic equilibration_compare.tex && tectonic ridge_trench_overview_v2.tex && tectonic taux_cases.tex`.

## Supporting information

Numbering follows the canonical SI source (`2026_codex/full_manuscript/si.tex`): Figures S1–S9; Tables S1 (symbols),
S2 (verification), S3 (convergence).

| SI | Manuscript file | Command (from repo root) | Table written in the same pass | Data | Regenerate data |
|---|---|---|---|---|---|
| Fig. S1 | `benchmark_boef.png` | `python scripts/render_benchmark.py` (asserts: V misfit ≤ 1 %, end resultant 4.000 ± 1 % TN/m, parabola ≤ 0.5 %) | `benchmark_boef.csv` | `suite1_strength/elastic_deep_60km_V4`, `data/idealized_beam_elastic` | `suite1`; `model/idealized_beam.jl` |
| Fig. S2 | `benchmark_mp.png` | `python scripts/render_mp_benchmark.py` (asserts: 148 sections, mean ≤ 0.15 %, max ≤ 1 %) | `benchmark_mp.csv` | `data/idealized_beam` | `model/idealized_beam.jl` |
| Table S3 | convergence table | `python scripts/render_convergence.py` | `convergence.csv` (+ `.md`) | `data/convergence/*` (5 models) | `convergence` |
| Fig. S3 | `fig_taux_cases_grid.pdf` | `cd schematic && tectonic taux_cases.tex` → `taux_cases.pdf` (renamed on copy) | — | — | — |
| Fig. S4 | `core_profiles_deep60.png` | `python scripts/render_core_profiles.py` | — | `suite1_strength/tresca_deep_150_60km_V4` | `suite1` |
| Fig. S5 | `corrected_density.png` | `python scripts/render_corrected_density.py` | — | `suite1_strength/tresca_deep_150_60km_V4` | `suite1` |
| Fig. S6 | `hero_tresca_deep60_full.png` — the reference model in the full five-panel layout (adds the τzx panel) | `HERO_PANELS=dsr python scripts/render_hero.py data/suite1_strength/tresca_deep_150_60km_V4 figures/hero_tresca_deep60_full.png` | `hero_tresca_deep60_full.csv` | `suite1_strength/tresca_deep_150_60km_V4` | `suite1` |
| Fig. S7 | `ferrite_hero_h30.png` | `HERO_PANELS=dsr HERO_WINDOW_KM=200 HERO_MFIX=0 python scripts/render_hero.py data/suite4_thickness/tresca_150_30km figures/hero_tresca_30km.png` (renamed on copy) | `hero_tresca_30km.csv` | `suite4_thickness/tresca_150_30km` | `thickness` |
| Fig. S8, S9 | `hero_dd_vm_asym.png`, `hero_dd_vm_sym.png` | `HERO_PANELS=dsr python scripts/render_hero.py data/suite1_strength/dd_vm_asym_60km_V4 figures/hero_dd_vm_asym.png` (and `…sym…`) | `hero_dd_vm_{asym,sym}.csv` | `suite1_strength/dd_vm_*` | `suite1` |
| supporting (not in the live SI) | `hero_tresca_40km.png` — the 40 km member of Suite 4, same layout as S6; referenced only by the SI source snapshot | `HERO_PANELS=dsr HERO_WINDOW_KM=260 HERO_MFIX=0 python scripts/render_hero.py data/suite4_thickness/tresca_150_40km figures/hero_tresca_40km.png` | `hero_tresca_40km.csv` | `suite4_thickness/tresca_150_40km` | `thickness` |

`scripts/render_hero.py` is one parameterised renderer — the same figure layout for any model directory
(`render_hero.py MODEL_DIR OUT [DIFF_DIR]`; `HERO_WINDOW_KM`, `HERO_MFIX` set the plotted window for thin plates;
`HERO_PANELS` selects the map panels — `dr` (default, Fig. 3: σxx−σzz and ρ̂) or `dsr` (the SI figures, adding τzx)).
The hero-family SI figures are *model* variants, not different figure designs. From Python/Jupyter:
`import render_hero; fig = render_hero.build(model_dir, save=False)`.

The two benchmark scripts read the exported material-frame stress S directly (the idealized beam is a material-section
benchmark); everything else goes through the deformed-Cauchy extractor — see `CLAUDE.md` and the scripts' docstrings.

## Tables — every number the paper quotes, as CSV
`tables/*.csv`, one per figure or record, written by the same script in the same pass as its figure from the same
arrays (`gpe_analysis.write_table`); a sidecar `tables/<name>.json` names each one's script, figure and models. `reproduce.sh`
asserts every table was rewritten and reads its checks from them. Standalone tables: `model_summary.csv` (one row per
production model: h, V, N_D, w_T, columns, ΔGPE*, ΔN_D, identity, arm — `analysis/model_summary.py`),
`isostatic_column.csv` (the cost of assuming x_I lithostatic — `analysis/isostatic_column_test.py`), `convergence.csv`.

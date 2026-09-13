# REPRODUCE — every manuscript figure: script → data → command

The frozen figure set of *The 'trench pull' force: constraints from elasto-plastic bending models*. Every
figure regenerates from the data in `data/` without re-running the solver, except the convergence table
(noted below). **Run every command from the repository root** — Python in the `ferrite-figs` env, Julia as
`julia --project=.` (see `README.md`). Each script writes its figure into `figures/`, which ships with the
reference renders to compare against (output-equivalent, not byte-identical).

**One command does all of it: `./reproduce.sh`** — every command below in order, the headline numbers asserted
(ΔGPE\* 2.542 TN/m, identity < 0.05 %, arm 34.9 km, the S1/S2 benchmark misfits, the reconstruction table),
`START_HERE.ipynb` executed top to bottom, `pytest tests/` — nonzero exit on any failure.

Data folders are named after the manuscript's suites (`data/suite1_strength` … `data/suite4_thickness`);
the driver's *commands* stay descriptive: `v_sweep` writes `suite2_load`, `nd_sweep` writes `suite3_background`. The
provenance of every shipped model (command, parameters, SHA-256 of every file, origin) is in `data/DATA_MANIFEST.md`.

## Main text

| Fig | Manuscript file | Command (from repo root) | Data | Regenerate data |
|---|---|---|---|---|
| 1 | `ridge_trench_overview_v2.pdf` | `cd schematic && tectonic ridge_trench_overview_v2.tex` (hand TikZ) | — | — |
| 2 | `equilibration_compare.pdf` | `cd schematic && python gen_equilibration_compare.py && tectonic equilibration_compare.tex` | — (computed curves) | — |
| 3 | `hero_tresca_deep60.png` | `python scripts/render_hero.py` | `data/suite1_strength/tresca_deep_150_60km_V4` | `model/paper_models.jl suite1` |
| 4 | `gpe_correlation.png` | `python scripts/render_gpe_correlation.py` | `suite1_strength`, `suite3_background`, `suite2_load` | `suite1 nd_sweep v_sweep` |
| 5 | `gpe_compare_suite1.png` | `python scripts/render_gpe_compare.py` | `suite1_strength/*` (4) | `suite1` |
| 6 | `profiles.png` | `python scripts/render_profiles.py` | `suite1_strength/*` (4) | `suite1` |
| 7 | `ferrite_thickness_compare.png` | `python scripts/render_thickness_compare.py` → `thickness_compare.png` (renamed on copy) | `suite1_strength/tresca_*`, `suite4_thickness/*` (3) | `suite1 thickness` |

## Supporting information

Numbering follows the canonical SI source (`2026_codex/full_manuscript/si.tex`): Figures S1–S8; Tables S1 (symbols),
S2 (verification), S3 (convergence).

| SI | Manuscript file | Command (from repo root) | Data | Regenerate data |
|---|---|---|---|---|
| Fig. S1 | `benchmark_boef.png` | `python scripts/render_benchmark.py` (asserts: V misfit ≤ 1 %, end resultant 4.000 ± 1 % TN/m, parabola ≤ 0.5 %) | `suite1_strength/elastic_deep_60km_V4`, `data/idealized_beam_elastic` | `suite1`; `model/idealized_beam.jl` |
| Fig. S2 | `benchmark_mp.png` | `python scripts/render_mp_benchmark.py` (asserts: 148 sections, mean ≤ 0.15 %, max ≤ 1 %) | `data/idealized_beam` | `model/idealized_beam.jl` |
| Table S3 | convergence table | `python scripts/render_convergence.py` | `data/convergence/` — **not included**; values in `data/CONVERGENCE.md` | **run first:** `julia --project=. model/paper_models.jl convergence` |
| Fig. S3 | `fig_taux_cases_grid.pdf` | `cd schematic && tectonic taux_cases.tex` → `taux_cases.pdf` (renamed on copy) | — | — |
| Fig. S4 | `core_profiles_deep60.png` | `python scripts/render_core_profiles.py` | `suite1_strength/tresca_deep_150_60km_V4` | `suite1` |
| Fig. S5 | `corrected_density.png` | `python scripts/render_corrected_density.py` | `suite1_strength/tresca_deep_150_60km_V4` | `suite1` |
| Fig. S6 | `ferrite_hero_h30.png` | `HERO_WINDOW_KM=200 python scripts/render_hero.py data/suite4_thickness/tresca_150_30km figures/hero_tresca_30km.png` (renamed on copy) | `suite4_thickness/tresca_150_30km` | `thickness` |
| Fig. S7, S8 | `hero_dd_vm_asym.png`, `hero_dd_vm_sym.png` | `python scripts/render_hero.py data/suite1_strength/dd_vm_asym_60km_V4 figures/hero_dd_vm_asym.png` (and `…sym…`) | `suite1_strength/dd_vm_*` | `suite1` |
| supporting (not in the live SI) | `hero_tresca_40km.png` — the 40 km member of Suite 4, same layout as S6; referenced only by the SI source snapshot | `HERO_WINDOW_KM=260 HERO_MFIX=0 python scripts/render_hero.py data/suite4_thickness/tresca_150_40km figures/hero_tresca_40km.png` | `suite4_thickness/tresca_150_40km` | `thickness` |

`scripts/render_hero.py` is one parameterised renderer — the same figure layout for any model directory
(`render_hero.py MODEL_DIR OUT [DIFF_DIR]`; `HERO_WINDOW_KM`, `HERO_MFIX` set the plotted window for thin plates).
The hero-family SI figures are *model* variants, not different figure designs. From Python/Jupyter:
`import render_hero; fig = render_hero.build(model_dir, save=False)`.

The two benchmark scripts read the exported material-frame stress S directly (the idealized beam is a material-section
benchmark); everything else goes through the deformed-Cauchy extractor — see `CLAUDE.md` and the scripts' docstrings.

## Records backing numbers quoted in the SI (not figures)
- `data/BOUNDARY_ARTIFACT.md` — loaded-edge sensitivity (29 % force- vs displacement-controlled split; raw
  trench 12–14 % above the interior; convergence within 40–60 km of the loaded edge).
- `data/CONVERGENCE.md` — the Table S3 convergence / consistency values.
- `data/ISOSTATIC_COLUMN.md` — the cost of assuming the first isostatic column lithostatic (Σzz(x_I) = 0, as the analytic
  estimate does) instead of measuring it: −0.13 % on the reference model, up to ~2 % in Suite 1 and ~4 % with a background
  N_D (`python analysis/isostatic_column_test.py`, run by `reproduce.sh`). The SI quotes this.
- `scripts/render_gpe_compare.py` prints the plate-top-arm (~2 %) and sea-level-arm (~8 %) reconstruction errors.

# MANIFEST — the frozen manuscript figure set: figure → script → data → command

Every figure in *The 'trench pull' force: constraints from elasto-plastic bending models*, the script that
draws it, the model data it reads, and the driver command that regenerates that data. All data listed is
**included in this repository**, so every figure reproduces without re-running the solver — except the S2
convergence table (noted below). Run all commands from the repository root; Python via the
`ferrite-figs` env, Julia via `julia --project=.` (see `README.md`).

Manuscript ↔ driver suite naming: the manuscript's **Suite 2 (load sweep)** is the driver's `v_sweep` →
`out/paper/set3_v_sweep`; the manuscript's **Suite 3 (background N_D)** is `nd_sweep` → `out/paper/set2_nd_sweep`.

## Main text

| Fig | Manuscript file | Script → output | Model data | Regenerate data |
|---|---|---|---|---|
| 1 | `ridge_trench_overview_v2.pdf` | `schematic/ridge_trench_overview_v2.tex` (hand TikZ) → same | — | `cd schematic && tectonic ridge_trench_overview_v2.tex` |
| 2 | `equilibration_compare.pdf` | `schematic/gen_equilibration_compare.py` → `.tex` → same | — (computed curves) | `cd schematic && python gen_equilibration_compare.py && tectonic equilibration_compare.tex` |
| 3 | `hero_tresca_deep60.png` | `render_hero.py` → same | `set1_rheology/tresca_deep_150_60km_V4` | `paper_models.jl suite1` |
| 4 | `gpe_correlation.png` | `render_gpe_correlation.py` → same | `set1_rheology`, `set2_nd_sweep`, `set3_v_sweep` | `suite1 nd_sweep v_sweep` |
| 5 | `gpe_compare_suite1.png` | `render_gpe_compare.py` → same (its only output) | `set1_rheology/*` (4) | `suite1` |
| 6 | `profiles.png` | `render_profiles.py` → same | `set1_rheology/*` (4) | `suite1` |
| 7 | `ferrite_thickness_compare.png` | `render_thickness_compare.py` → `thickness_compare.png` (renamed on copy) | `set1_rheology/tresca_*`, `set4_thickness/*` (3) | `suite1 thickness` |

## Supporting information

| SI | Manuscript file | Script → output | Model data | Regenerate data |
|---|---|---|---|---|
| S1 | `benchmark_boef.png` | `render_benchmark.py` → same | `set1_rheology/elastic_deep_60km_V4`, `out/idealized_beam_elastic` | `suite1`; `idealized_beam.jl` |
| S2 | `benchmark_mp.png` | `render_mp_benchmark.py` → same | `out/idealized_beam` | `idealized_beam.jl` |
| S2 | convergence table | `render_convergence.py` | `out/paper/convergence/` — **not included** (values in `out/paper/CONVERGENCE.md`) | `paper_models.jl convergence` (**required first**) |
| S4 | `core_profiles_deep60.png` | `render_core_profiles.py` → same | `set1_rheology/tresca_deep_150_60km_V4` | `suite1` |
| S5 | `corrected_density.png` | `render_corrected_density.py out/paper/figures/corrected_density.png` (OUT arg required) | `set1_rheology/tresca_deep_150_60km_V4` | `suite1` |
| S6 | `edge_effect_SI.png` | `render_edge_si.py` → same | `set1_rheology/*` (4) | `suite1` |
| SI | `fig_taux_cases_grid.pdf` | `schematic/taux_cases.tex` (hand TikZ) → `taux_cases.pdf` (renamed on copy) | — | `cd schematic && tectonic taux_cases.tex` |
| SI | `hero_dd_vm_asym.png`, `hero_dd_vm_sym.png` | `render_hero.py MODEL_DIR OUT` → same | `set1_rheology/dd_vm_*` | `suite1` |
| SI | `ferrite_hero_h40.png`, `ferrite_hero_h30.png` | `HERO_WINDOW_KM=260 HERO_MFIX=0 render_hero.py …40km`, `HERO_WINDOW_KM=200 render_hero.py …30km` → `hero_tresca_40km.png`, `hero_tresca_30km.png` (renamed on copy) | `set4_thickness/tresca_150_{40,30}km` | `thickness` |
| App. B | (number only: E-sweep arm ≈ 35 km, §4.1) | `render_esweep_test.py` → `esweep.png` | `esweep/*` + `set1_rheology`, `set2_nd_sweep`, `set3_v_sweep` | `all` |

`render_hero.py` is one parameterised renderer: the same figure layout for any model directory it is given
(`render_hero.py MODEL_DIR OUT [DIFF_DIR]`; `HERO_WINDOW_KM`, `HERO_MFIX` adjust the plotted window for the
thin plates). The hero-family SI figures are therefore *model* variants, not different figure designs.

## Records backing numbers quoted in the SI (not figures)
- `out/paper/BOUNDARY_ARTIFACT.md` — loaded-edge sensitivity (29 % force- vs displacement-controlled split;
  raw trench 12–14 % above the interior; convergence within 40–60 km).
- `out/paper/CONVERGENCE.md` — the S2 convergence / consistency table values.
- `render_gpe_compare.py` prints the plate-top-arm (~2 %) and sea-level-arm (~8 %) reconstruction errors.

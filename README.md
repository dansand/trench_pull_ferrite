# trench_pull_ferrite

Models, analysis and figure code for

> **The 'trench pull' force: constraints from elasto-plastic bending models** — Sandiford (JGR: Solid Earth, submitted 2026).

Plane-strain, total-Lagrangian finite-element models of a bending lithospheric plate (Ferrite.jl), and the
Python post-processing that turns them into the paper's figures. **All model output data is included**
(~113 MB; largest file 4 MB), so every figure reproduces from this repository alone, without re-running
the finite-element solves. Re-running the solves is also supported (see §4).

This is the lean, curated release. The full development history lives in a separate private archive.

## 1. Layout

| | |
|---|---|
| `gpe_mvm.jl`, `gpe_plastic.jl` | the solver: massless plane-strain flexure, plane-strain J2 / in-plane Tresca plasticity |
| `paper_models.jl` | the **production driver** — every model in the paper is a command of this script |
| `idealized_beam.jl` | the beam benchmarks (SI Text S1) |
| `gpe_analysis.py` | the analysis library: deformed-mesh Cauchy integration, trench pull ΔGPE\*, equivalent density ρ̂ |
| `render_*.py` | one script per manuscript figure (see `MANIFEST.md`) |
| `schematic/` | TikZ sources for Figures 1, 2 and the SI stress-regime grid |
| `out/paper/` | **all model data** (`set1_rheology`, `set2_nd_sweep`, `set3_v_sweep`, `set4_thickness`, `esweep`) and the reference figure renders (`figures/`) |
| `out/idealized_beam*` | benchmark beam output |
| `FIGURE_LEDGER.md`, `MANIFEST.md` | figure → script → data → command; the frozen manuscript figure set |
| `demo_trench_pull.ipynb` | a short walk-through of the analysis on the reference model |

Naming note: the production driver's output directories are `set2_nd_sweep` (background-stress sweep) and
`set3_v_sweep` (load sweep); in the **manuscript** these are **Suite 3** and **Suite 2** respectively. The
mapping is in `MANIFEST.md`.

## 2. Environments

**Julia** (solver) — Julia 1.10.5, Ferrite.jl 1.4.1, Tensors.jl 1.17.0, pinned in `Project.toml`/`Manifest.toml`:
```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
**Python** (analysis + figures) — Python 3.12, pyvista 0.46.3, numpy 1.26.4, matplotlib 3.9.1, scipy 1.13.1:
```bash
conda env create -f environment.yml && conda activate ferrite-figs
```
**TikZ** (Figures 1, 2, SI grid) — compile with [`tectonic`](https://tectonic-typesetting.github.io/) (no
pdflatex needed).

Run every command below **from the repository root**.

## 3. Reproduce the figures (from the included data)

`MANIFEST.md` lists every figure with its exact command. In short:

```bash
# main text
python render_hero.py                                                    # Fig 3  hero_tresca_deep60.png
python render_gpe_correlation.py                                         # Fig 4  gpe_correlation.png
python render_gpe_compare.py                                             # Fig 5  gpe_compare_suite1.png
python render_profiles.py                                                # Fig 6  profiles.png
python render_thickness_compare.py                                       # Fig 7  thickness_compare.png
# supporting information
python render_benchmark.py            # S1 benchmark_boef.png       python render_mp_benchmark.py   # S2 benchmark_mp.png
python render_core_profiles.py        # S4 core_profiles_deep60.png python render_edge_si.py        # S6 edge_effect_SI.png
python render_corrected_density.py out/paper/figures/corrected_density.png                          # S5
python render_esweep_test.py          # App. B  esweep.png (the E-sweep arm quoted in §4.1)
# hero variants (same script, different model): SI hero suite + thin-plate heroes
python render_hero.py out/paper/set1_rheology/dd_vm_asym_60km_V4 out/paper/figures/hero_dd_vm_asym.png
python render_hero.py out/paper/set1_rheology/dd_vm_sym_60km_V4  out/paper/figures/hero_dd_vm_sym.png
HERO_WINDOW_KM=260 HERO_MFIX=0 python render_hero.py out/paper/set4_thickness/tresca_150_40km out/paper/figures/hero_tresca_40km.png
HERO_WINDOW_KM=200             python render_hero.py out/paper/set4_thickness/tresca_150_30km out/paper/figures/hero_tresca_30km.png
# schematics
( cd schematic && python gen_equilibration_compare.py && tectonic equilibration_compare.tex \
                && tectonic ridge_trench_overview_v2.tex && tectonic taux_cases.tex )
```
Each script overwrites its file in `out/paper/figures/`; the shipped copies are the reference renders to
compare against. Expect output-equivalent figures (identical data, layout and numbers); byte-identical PNGs
are not guaranteed across matplotlib/font builds. Several scripts print the numbers quoted in the paper
(e.g. `render_gpe_compare.py` prints the ~2% plate-top-arm and ~8% sea-level-arm reconstruction errors).

Some figures are **renamed when copied into the manuscript**: `thickness_compare.png` →
`ferrite_thickness_compare.png`; `hero_tresca_30km.png`/`hero_tresca_40km.png` → `ferrite_hero_h30.png`/`ferrite_hero_h40.png`;
`schematic/taux_cases.pdf` → `fig_taux_cases_grid.pdf`.

**One exception — the S2 convergence table.** `render_convergence.py` reads `out/paper/convergence/`, which
is **not included** (regeneration-only). The table's values are recorded in `out/paper/CONVERGENCE.md`; to
regenerate the data run `julia --project=. paper_models.jl convergence` first (several solves).

## 4. Re-run the models

The driver is skip-existing (it will not recompute a model whose output exists — move/delete a directory to
force a rerun) and never overwrites:
```bash
julia --project=. paper_models.jl suite1          # Suite 1: four strength models at V = 4 TN/m  (the reference set)
julia --project=. paper_models.jl v_sweep         # manuscript Suite 2: load sweep V = 1 … 4.5   → out/paper/set3_v_sweep
julia --project=. paper_models.jl nd_sweep        # manuscript Suite 3: background N_D = −3 … +3 → out/paper/set2_nd_sweep
julia --project=. paper_models.jl thickness       # Suite 4: h = 30/40/50 km at matched deflection (secant-tuned V)
julia --project=. paper_models.jl esweep          # Appendix B: Young's-modulus sweep
julia --project=. paper_models.jl convergence     # SI Table S2 (not part of `all`)
julia --project=. paper_models.jl all             # suite1 + nd_sweep + v_sweep + esweep + thickness
julia --project=. idealized_beam.jl               # the benchmark beams (S1, S2)
```
Each model directory carries a `provenance.txt` stamp (driver, parameters, `stress_frame="massless"`).

Reference model, for orientation: uniform Tresca σ_Y = 150 MPa, h = 60 km, V = 4 TN/m → trench deflection
3233 m, ΔGPE\* = 2.54 TN/m, with the equilibrium identity ΔN_D = ΔGPE\* satisfied to 0.02%.

## 5. Citation

If you use this code, cite the paper above and the solver:
Carlsson, K., Ekre, F. & Ferrite.jl contributors, *Ferrite.jl* v1.4.1, doi:10.5281/zenodo.13862652.

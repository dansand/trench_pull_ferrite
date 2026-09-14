# trench_pull_ferrite

<img src="animation/trench_pull_loading.gif" width="820">

*A 60 km elasto-plastic plate loaded at its trench edge: the deflection grows with the applied shear, and with it a
pressure deficit beneath the trench — a GPE-like resultant, ΔGPE\* = −Δσ̄zz. That resultant is the trench pull. Static
equilibrium requires it to be balanced, and in these models it is balanced by an exchange with the other horizontal
component: an equal and opposite normal-stress-difference resultant, ΔN_D, arises. (24 load steps; `animation/`.)*

*The models start from a plate carrying no horizontal force at all — true even of the horizontal normal stress, since there
are no body forces. Loading the left edge vertically then couples the vertical load to the horizontal resultants: an
exchange between the vertical and the horizontal force. The nature and degree of that coupling — the efficiency of the
exchange, how far it departs from the simple picture of Section 2 of the associated manuscript, and how well it agrees
with the analytical predictions there — is what these models measure.*

Models, analysis and figure code for

> **The 'trench pull' force: constraints from elasto-plastic bending models** — Sandiford (JGR: Solid Earth, submitted 2026).

<img src="figures/hero_tresca_deep60.png" width="820">

*The reference model (uniform Tresca, σ_Y = 150 MPa, V = 4 TN/m): topography, stress, the shear stress τzx, the equivalent density
ρ̂ = τzx,x/g that supports the pressure deficit, and the resultants — the trench pull ΔGPE\* = 2.54 TN/m, with the equilibrium
identity ΔN_D = ΔGPE\* holding to 0.02 %.*

Plane-strain, total-Lagrangian finite-element models of a bending lithospheric plate (Ferrite.jl), and the Python
post-processing that turns them into the paper's figures. **All model output data is included** (~100 MB; largest
file 4 MB), so every figure reproduces from this repository alone, without re-running the solver. Re-running the
solver is also supported (§4).

**New here? Open `START_HERE.ipynb`.**

## 1. Layout

```
model/        the solver and the production driver (Julia)
                gpe_plastic.jl               the solver: total-Lagrangian massless plate on a Δρg foundation; in-plane Tresca / depth-dependent von Mises plasticity
                paper_models.jl              every model in the paper is a command of this script
                idealized_beam.jl            the beam benchmarks (SI S1, S2)
analysis/     gpe_analysis.py — deformed-mesh Cauchy integration, trench pull ΔGPE*, equivalent density ρ̂
scripts/      one render_*.py per manuscript figure (see REPRODUCE.md)
figures/      the reference renders — what the scripts write; the manuscript's copies are taken from here
tables/       every number the paper quotes, as CSV — written by the same script, in the same pass, as its figure
                model_summary.csv (one row per model), isostatic_column.csv, convergence.csv (Table S3), the
                benchmark misfits, the Fig. 5 reconstruction errors, the Fig. 7 values, the hero reference lines
schematic/    TikZ sources for Figures 1, 2 and the SI stress-regime grid
data/         all model output, one folder per manuscript suite: suite1_strength, suite2_load, suite3_background,
              suite4_thickness; idealized_beam* (benchmarks)
              DATA_MANIFEST.md / .json — per model: generating command, parameters, SHA-256 of every file, origin
              BOUNDARY_ARTIFACT.md — the loaded-edge investigation record the SI quotes
START_HERE.ipynb   the analysis step by step, then on any model, then the paper's figures from their scripts
REPRODUCE.md       figure → script → data → command, for every figure in the paper
reproduce.sh       one command that regenerates and checks every figure (and runs the notebook and the tests)
tests/test_quick.py  fast checks (seconds): push-forward, centroid, resultants, benchmark misfits, import hygiene
```
Data folders are named after the manuscript's suites; the driver's *commands* stay descriptive (`v_sweep` → `data/suite2_load`, `nd_sweep` → `data/suite3_background`).

## 2. Environments

**Julia** (solver) — Julia 1.10.5, Ferrite.jl 1.4.1, Tensors.jl 1.17.0, pinned in `Project.toml`/`Manifest.toml`:
```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
**Python** (analysis + figures) — Python 3.12, pyvista 0.46.3, numpy 1.26.4, matplotlib 3.9.1, scipy 1.13.1:
```bash
conda env create -f environment.yml && conda activate ferrite-figs
python -m ipykernel install --user --name ferrite-figs   # so START_HERE.ipynb can select this env as its kernel
```
Headless machines: the figure scripts draw with matplotlib (Agg, no display needed) — except `scripts/render_hero.py`,
whose field panels are rendered off-screen by PyVista/VTK and need an OpenGL context. On a display-less Ubuntu 22.04 VM with
the Mesa libraries installed (`libegl1-mesa`, `libgl1-mesa-dri`) VTK falls back to Mesa's software EGL renderer by itself and
the full `./reproduce.sh` passes (verified 2026-09-13); without Mesa, run it under Xvfb (`xvfb-run -a python scripts/render_hero.py`).
**TikZ** (Figures 1, 2, SI grid) — [`tectonic`](https://tectonic-typesetting.github.io/); no pdflatex needed.

**Run every command from the repository root.** Scripts locate `analysis/` and `data/` relative to it.

## 3. Reproduce the figures (from the included data)

`REPRODUCE.md` has every figure's exact command. The short version:
```bash
python scripts/render_hero.py                 # Fig 3     python scripts/render_benchmark.py       # S1
python scripts/render_gpe_correlation.py      # Fig 4     python scripts/render_mp_benchmark.py    # S2
python scripts/render_gpe_compare.py          # Fig 5     python scripts/render_core_profiles.py   # S4
python scripts/render_profiles.py             # Fig 6     python scripts/render_corrected_density.py  # S5
python scripts/render_thickness_compare.py    # Fig 7
( cd schematic && python gen_equilibration_compare.py && tectonic equilibration_compare.tex \
                && tectonic ridge_trench_overview_v2.tex && tectonic taux_cases.tex )       # Figs 1, 2, SI grid
```
**`./reproduce.sh` regenerates and checks every figure**: it runs every command above in order, asserts the headline
numbers (ΔGPE\* = 2.542 TN/m, identity < 0.05 %, arm 34.9 km, the S1/S2 benchmark misfits, the reconstruction
table), executes `START_HERE.ipynb` top to bottom and runs `pytest tests/`, exiting nonzero on any failure. The numbers it checks are
read from `tables/*.csv`, which every figure script writes alongside its figure from the same arrays — a table and
its figure cannot disagree, and the harness asserts both were rewritten.

Each script writes its figure into `figures/`; the shipped copies are the reference renders.
Expect output-equivalent figures (identical data, layout and numbers) — byte-identical PNGs are not guaranteed
across matplotlib/font builds. Several scripts print the numbers quoted in the paper (e.g. `render_gpe_compare.py`
prints the ~2 % plate-top-arm and ~8 % sea-level-arm reconstruction errors).

Some figures are **renamed when copied into the manuscript**: `thickness_compare.png` → `ferrite_thickness_compare.png`;
`hero_tresca_{30,40}km.png` → `ferrite_hero_h{30,40}.png`; `schematic/taux_cases.pdf` → `fig_taux_cases_grid.pdf`.

The one figure-side exception: the **Table S3 convergence table** reads `data/convergence/`, which is not included
(regeneration-only, several solves). Its values are in `tables/convergence.csv` (and `.md`); to regenerate, run
`julia --project=. model/paper_models.jl convergence` first.

## 3b. Tables — the numbers the paper quotes

`tables/` holds every model-derived number the manuscript and SI quote, as CSV, one file per figure or record. Each is
written by the **same script, in the same pass, from the same arrays** as its figure, so a table and its figure cannot
disagree; `./reproduce.sh` rewrites all of them and reads its numerical checks from them. Each file's `#` header names the
script, the figure and the models it came from. `tables/README.md` is the handoff for the manuscript: which file backs
which number, columns and units, and how to validate. Standalone tables: `model_summary.csv` (one row per production
model: thickness, load, background N_D, w_T, the three columns, ΔGPE\*, ΔN_D, identity residual, effective arm),
`isostatic_column.csv`, `convergence.csv` (Table S3).

## 4. Re-run the models

The driver is skip-existing and never overwrites: it checks only for `gpe_model.vtu` in the target directory, so a
finished model is skipped and a directory with partial output (no VTU) is re-solved — reruns go into an empty or
absent directory; move a finished directory aside to force a rerun. Output goes to `data/<suite>/<model>/`.
Every model the driver produces gets a `provenance.txt` stamp (commit, parameters, completion flag). **The shipped
models predate that stamp and carry none**; their provenance — generating command, parameters, SHA-256 of every
file, origin — is recorded centrally in `data/DATA_MANIFEST.md` (`python analysis/make_manifest.py --check`
verifies the shipped files against it). The shipped models were produced with Julia 1.10.5 and Ferrite.jl 1.4.1 (the
pinned `Manifest.toml`); re-solving the baseline with that environment on macOS arm64 reproduced the shipped files byte for
byte (w_T 3233 m, ΔGPE\* 2.542 TN/m, identity residual 0.019 %). On other platforms expect the numbers, not necessarily
the bytes.
```bash
julia --project=. model/paper_models.jl suite1        # Suite 1: four strength models at V = 4 TN/m (the reference set)
julia --project=. model/paper_models.jl v_sweep       # manuscript Suite 2: load sweep V = 1 … 4.5   → data/suite2_load
julia --project=. model/paper_models.jl nd_sweep      # manuscript Suite 3: background N_D = −3 … +3 → data/suite3_background
julia --project=. model/paper_models.jl thickness     # Suite 4: h = 30/40/50 km at matched deflection (secant-tuned V)
julia --project=. model/paper_models.jl convergence   # Table S3 (not part of `all`)
julia --project=. model/paper_models.jl all           # suite1 + nd_sweep + v_sweep + thickness
julia --project=. model/idealized_beam.jl             # the benchmark beams (S1, S2)
```
Reference model, for orientation: uniform Tresca σ_Y = 150 MPa, h = 60 km, V = 4 TN/m → trench deflection 3233 m,
ΔGPE\* = 2.54 TN/m, with the equilibrium identity ΔN_D = ΔGPE\* satisfied to 0.02 %.

### What is in a model's `gpe_model.vtu` — read this before using the fields directly
The VTU is a VTK XML unstructured grid (ParaView/PyVista-readable) of the **reference (undeformed) mesh**, with the
displacement `u` and these point fields:

| field | what it is |
|---|---|
| `sigma_xx`, `sigma_zz`, `sigma_xz` [Pa] | the **second Piola–Kirchhoff stress S** on the reference configuration — **not Cauchy stress**. The plate tilts by up to a few percent, so S and Cauchy differ by O(slope) mixing — measured on the reference model: up to 5 % of the field scale for τzx, 8 % for σxx−σzz, 9 % for σzz (95th percentile ≈ 4 %). |
| `sigma_xz_cauchy` [Pa] | the true Cauchy shear σxz = J⁻¹F·S·Fᵀ, projected to nodes |
| `dsxz_dx`, `dszz_dz` [Pa/m] | FE spatial gradients **of the Cauchy** σxz and σzz on the deformed geometry; `dsxz_dx` = τzx,x is the equivalent density ρ̂·g |
| `pressure`, `von Mises` [Pa], `plastic_strain`, `yielded` | derived scalars / the yield indicator (0/1) |

**Deformed coordinates are not stored**: the mesh in the file is the reference rectangle; the bent plate is reference + `u`,
formed in post-processing (`gpe_analysis.deformed_line`, `START_HERE.ipynb`'s `plot_field`) or in ParaView with *Warp By Vector* on `u`.
Likewise the stress: `sigma_*` is the 2nd Piola–Kirchhoff stress S (work-conjugate to the Green–Lagrange strain, expressed in material
axes, unchanged by rigid rotation); the Cauchy stress σ = J⁻¹ F S Fᵀ carries it into the lab axes on the deformed body — at these
strains (~0.2 %) essentially a rotation by the local plate tilt.

To work with the **full Cauchy tensor**, use `analysis/gpe_analysis.py`: `Model(dir).cauchy_fields()` pushes all three
components forward (`σ = J⁻¹ F S Fᵀ`, `F = I + ∂u/∂X`) — that is what every number in the paper is computed from.
(Only the shear was exported as Cauchy from Julia because ρ̂ needs its *gradient*, which is clean on the FE side;
the naming is a historical artefact, not a principled scheme.)

## 5. Animations (not in the manuscript)

`animation/` shows the reference model being loaded incrementally — the trench deepening, the equivalent density ρ̂
and the shear stress τzx growing with load — as three short movies, shipped as `.gif` (renders inline on GitHub)
and `.mp4`: `trench_pull_loading`, `trench_pull_rhohat`, `trench_pull_sxz`. They are illustrative and no manuscript
figure depends on them. To regenerate: `julia --project=. animation/gen_frames.jl` solves the reference Tresca
model at 24 load fractions (a model run, minutes; writes `animation/frames/`, not shipped), then
`python animation/render_anim.py`, `render_anim_rho.py`, `render_anim_sxz.py` assemble the movies.

## 6. Citation

Cite the paper above and the solver: Carlsson, K., Ekre, F. & Ferrite.jl contributors, *Ferrite.jl* v1.4.1,
doi:10.5281/zenodo.13862652.

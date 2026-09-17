# CLAUDE.md

Working notes for anyone (or any coding assistant) editing this repository. Written for a clone on any machine: nothing
below depends on the author's other directories. If you only want to reproduce the figures, `README.md` and
`REPRODUCE.md` are enough; this file is about how to change things without breaking them.

## What this is
The reproducibility package for *The 'trench pull' force: constraints from elasto-plastic bending models* (Sandiford;
preprint doi:10.22541/essoar.174413825.53806221/v1, a revised version is under review). It contains the finite-element
models (Julia, Ferrite.jl), all model output, and the scripts that turn that output into every manuscript figure and
the model-derived numbers the paper quotes (the manuscript-side ledger lists the few literals that are not backed by a
table here). It was curated from a larger private development history; that history is not needed
and not public. **This repository is the single source of truth for the paper's models and numbers.**

## Layout (run everything from the repository root)
`model/` Julia solver (`gpe_plastic.jl`) and the production driver `paper_models.jl` (one command per suite) ·
`analysis/gpe_analysis.py` the extraction library (+ `make_manifest.py`, `model_summary.py`, `paper_numbers.py` and
the assumption tests) · `scripts/render_*.py` one per manuscript figure · `figures/` the shipped renders (what the
scripts write) · `tables/` the model-derived numbers the paper quotes, as CSV with a JSON provenance sidecar, written by the same
script in the same pass as its figure (`tables/README.md` explains the contract) · `schematic/` TikZ ·
`data/` all model output, one directory per model, with `DATA_MANIFEST.md/.json` (command, parameters, SHA-256, origin) ·
`START_HERE.ipynb` a guided tour · `REPRODUCE.md` figure → script → data → command · `reproduce.sh` regenerates and
checks everything · `tests/` (`test_quick.py`, `test_quick.jl`, `check_reproduce.py`) · `animation/` loading movies
(illustrative; `gen_frames.jl` is a model run of 24 solves; `frames/` is not shipped).

## Rules that keep the package trustworthy
- **One script → one figure (+ its table).** A script that produces a number the paper quotes writes it to
  `tables/<figure>.csv` through `gpe_analysis.write_table`, from the same arrays as the figure, never by a separate
  computation. `reproduce.sh` asserts that every table was rewritten and reads all its checks from the tables.
  `render_hero.py` is the one parameterised renderer (the same layout for any model directory); several designs of a
  figure in one script is not allowed.
- **Figure scripts stay importable:** `matplotlib.use("Agg")` only under `if __name__ == "__main__"`, so the notebook can
  import them without losing its backend. Helpers that draw on a given `ax` are preferred.
- **Model runs only through the driver.** `julia --project=. model/paper_models.jl <command>` skips existing model
  directories and never overwrites; move a directory aside to force a rerun. Solves take minutes to hours; there are no
  ad-hoc model runs. An assistant must ask before starting any solve.
- **Data are provenance-documented.** Every model directory has an entry in the manifest (command, parameters, SHA-256
  of the primary output files, origin); models solved with the release code also carry a `provenance.txt` stamp
  (parameters, solver version, completion flag) — the five convergence models do, the older shipped models predate it. After adding or regenerating a model, run `python analysis/make_manifest.py` and
  check with `--check`. Provenance is the pinned solver versions plus the manifest, not commit hashes.
- **Don't drop figure elements** (curves, lines, panels) when restyling without saying so.
- **Commits and pushes are made by the author**, not by an assistant. Commit messages are one short line.

## Conventions
- Python from `environment.yml` (`ferrite-figs`: numpy 1.26, so `np.trapz`, not `np.trapezoid`). Julia 1.10.5 with the
  shipped `Project.toml`/`Manifest.toml` (Ferrite.jl 1.4.1); both environments have been verified from clean installs.
  TikZ compiles with `tectonic` (no pdflatex needed).
- Terminology follows the manuscript: **equivalent density** ρ̂ = τzx,x/g (never "pseudo-density"); the shear-stress
  gradient is written `τzx,x`; **centroids are always signed** (∫z·τ dz / ∫τ dz), never magnitude-weighted.
- Two distinct "arm" quantities, kept apart in comments and labels: **A**, the force-based effective arm
  ΔGPE\*/(Δρ g w_T), the one the paper reports (≈ 35 km, ≈ 0.58 h); **B**, the signed centre of mass of the τzx,x
  distribution, drawn for illustration and undefined where the net charge of a column vanishes.
- Trench pull is always trench-referenced: between the trench column (the leftmost complete deformed column) and the
  **first isostatic column** x_I (first w = 0), via `gpe_analysis.trench_pull`. The identity ΔN_D = ΔGPE\* holds to
  about 0.02 % and is asserted.
- Data folders are named after the manuscript's suites (`suite1_strength`, `suite2_load`, `suite3_background`,
  `suite4_thickness`, plus `convergence` and the two `idealized_beam*` benchmarks); the driver commands stay descriptive
  (`v_sweep` → suite2, `nd_sweep` → suite3).
- **The trench-edge yield-strength ramp** (+100 MPa, e-folding 10 km, in `paper_models.jl`, recorded in the manifest)
  exists because the parabolic face traction (peak 1.5V/h = 100 MPa at V = 4 TN/m) exceeds the in-plane Tresca shear
  capacity σY/2 = 75 MPa; without it the loaded face yields in shear, a deformation mode unrelated to the bending the
  paper studies. It is not what makes the trench column elastic (that is the bending moment vanishing at the free end).
  This was settled by a ramp-free comparison (author's records, 2026-09-15); do not reopen a reference-model redesign
  on this ground.
- **VTU stress fields are not all Cauchy.** `sigma_xx/zz/xz` in `gpe_model.vtu` are the 2nd Piola–Kirchhoff stress on
  the reference mesh; only `sigma_xz_cauchy` and the gradients `dsxz_dx`, `dszz_dz` are Cauchy. Always go through
  `Model.cauchy_fields()` (σ = J⁻¹ F S Fᵀ) for stresses: the rotation mixes components at the order of the surface
  slope, tens of per cent on the shear. **One deliberate exception:** the two idealized-beam benchmarks
  (`render_benchmark.py`, `render_mp_benchmark.py`) read the exported material-frame stress directly, because their
  analytic comparisons (Hetényi deflection, shear parabola, M–κ) are material-section quantities on the reference
  thickness. Both scripts assert the SI numbers, so a wrong frame fails loudly.

## Manuscript handshake
The manuscript reads figures from its own directory under the names listed in `REPRODUCE.md` (some are renamed on
copy). Regenerated figures are copied there by hand, never written by a script. The numbers the paper quotes are consumed
as LaTeX macros from `tables/paper_numbers.tex`; `python analysis/paper_numbers.py --compare <manuscript dir>` reports
which manuscript literals match, differ from, or already use the macros (`tables/README.md`, "how the numbers get into
the manuscript").

## Acceptance test
`./reproduce.sh` passes: `pytest tests/` (Python and, if Julia is present, `test_quick.jl`), every figure script from
the root against `data/`, the notebook top to bottom, `make_manifest.py --check`, and the headline numbers asserted by
`tests/check_reproduce.py` (ΔGPE\* = 2.542 TN/m ± 0.5 %, identity < 0.05 %, arm 34.9 ± 0.2 km, the S1/S2 benchmark
misfits, the reconstruction and convergence tables, the isostatic-column bounds, every figure and table rewritten).
Run it after any change to scripts, analysis, data or paths; a passing run is the definition of done. A printed number
nobody re-reads is not a check.

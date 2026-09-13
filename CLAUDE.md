# CLAUDE.md

Guidance for Claude Code in this repository.

## What this is
The **lean, public reproducibility package** for *The 'trench pull' force: constraints from elasto-plastic bending
models* (Sandiford, JGR: Solid Earth). Everything needed to reproduce every manuscript figure from the included
data, and to re-run the finite-element models. It was curated out of the private development archive
`~/projects/ferrite` (`ferrite_plate_flexure`), which keeps the full history, diagnostics and dev tooling — look
there for history, never develop there. **This repo is the single source of truth.**

## Layout (role-based; run everything from the repo root)
`model/` Julia solver + `paper_models.jl` production driver · `analysis/gpe_analysis.py` (+ `make_manifest.py`) ·
`scripts/render_*.py` (one per manuscript figure) · `figures/` the shipped reference renders (what the scripts write) ·
`tables/` every number the paper quotes as CSV, written by the same script in the same pass as its figure ·
`schematic/` TikZ · `data/` all model output + `DATA_MANIFEST.md/.json` (provenance: command, parameters, SHA-256, origin) ·
`START_HERE.ipynb` · `REPRODUCE.md` (figure → script → data → command) · `reproduce.sh` (regenerate + check everything) ·
`tests/test_quick.py` (+ `test_quick.jl`; seconds) · `animation/` (loading movies — illustrative,
no manuscript figure depends on them; `gen_frames.jl` is a **model run** (24 solves) → ask before running; `frames/` is gitignored).

## Hard rules
- **The user makes every commit and push.** Never `git commit` or `git push`. Stage if useful, then show the command.
- **Never run a finite-element solve unprompted.** Solves take minutes–hours and only ever go through
  `model/paper_models.jl <command>` (skip-existing, never clobbers) — no ad-hoc model runs. Ask first.
- **One script → one figure (+ its table).** No ghost outputs. A figure script that produces a number the paper quotes
  writes it to `tables/<figure>.csv` via `gpe_analysis.write_table` in the same pass, from the same arrays — never a
  separate computation; `reproduce.sh` asserts every table was rewritten and reads its checks from the tables. `render_hero.py` is the one parameterised renderer (same layout for
  any model directory) — that is fine; several *designs* of a figure in one script is not.
- Figure scripts must stay **importable**: `matplotlib.use("Agg")` lives under `if __name__ == "__main__"`, never at
  module level (it would hijack Jupyter's backend). Plotting helpers that take an `ax` are preferred.
- Don't drop figure elements (curves, lines, panels) when restyling without saying so.
- After rendering an image, open it: `open -a Preview <path>`.

## Conventions
- Python from the pinned `environment.yml` env (`ferrite-figs`, numpy 1.26 → use `np.trapz`, not `np.trapezoid`).
  Julia 1.10.5, `julia --project=.`, Ferrite.jl 1.4.1. TikZ compiles with **`tectonic`** (no pdflatex).
- Terminology: **equivalent density** ρ̂ = τzx,x/g (never "pseudo-density"); the shear-gradient is written
  `τzx,x`; **centroids are always signed** (∫z·τ/∫τ), never magnitude-weighted.
- Two distinct "arm" quantities — keep them apart in code comments and labels: **A**, the effective force-based
  arm ΔGPE\*/(Δρ g w) (the reported one, ≈35 km); **B**, the centre of mass of the τzx,x distribution (drawn,
  illustrative, undefined where the net dipole charge vanishes).
- Trench pull is always trench-referenced, between the trench column and the **first isostatic column** x_I
  (`gpe_analysis.trench_pull`); the identity ΔN_D = ΔGPE\* holds to ~0.02 %.
- Data folders are named after the manuscript's suites (`suite1_strength`, `suite2_load`, `suite3_background`,
  `suite4_thickness`); the driver commands stay descriptive (`v_sweep` → suite2, `nd_sweep` → suite3).
- **VTU stress fields are NOT all Cauchy.** `sigma_xx/zz/xz` in `gpe_model.vtu` are the 2nd Piola–Kirchhoff stress S on the
  reference mesh; only `sigma_xz_cauchy` (+ the gradients `dsxz_dx`, `dszz_dz`) are Cauchy. Always go through
  `Model.cauchy_fields()` (σ = J⁻¹F S Fᵀ) for stresses — rotation mixes components at O(slope), tens of percent on the shear.
  **One exception, by design:** the two idealized-beam benchmarks (`render_benchmark.py`, `render_mp_benchmark.py`)
  read the exported S directly — their analytic comparisons (Hetényi, shear parabola, M–κ) are material-section
  quantities on the reference thickness. Switching them to Cauchy broke them once (2026-09-11 → Codex audit); both
  scripts now assert the SI numbers, so a wrong frame fails loudly.

## Manuscript handshake
The manuscript is `~/projects/mypapers/trench_pull_force/2026_codex/full_manuscript/` (being finalised with Codex; its LaTeX
reads figures from its own `figures/` directory — some are renamed on copy, see `REPRODUCE.md`). `2026_version/` is a
static relic — never read it as current. When figures are regenerated for the paper: copy them there. Provenance is the
pinned solver versions (Julia 1.10.5, Ferrite.jl 1.4.1 in `Manifest.toml`) and `data/DATA_MANIFEST.md`; do not chase
commit hashes through the README or the paper's records (Dan, 2026-09-13 — that is more detail than the package needs).

## Acceptance test
`./reproduce.sh` passes: `pytest tests/`, every figure script from the root against `data/`, the notebook top to bottom in
the `ferrite-figs` env, `make_manifest.py --check`, and the headline numbers ASSERTED (not printed) by
`tests/check_reproduce.py` — ΔGPE\* = 2.542 TN/m ± 0.5 %, identity < 0.05 % (0.019 %), arm 34.9 km, the S1/S2 benchmark
misfits, the Suite-1 reconstruction table. Run it after any change to scripts, analysis or paths; a passing run is the
definition of done. (Lesson of the 2026-09-12 audit: a printed number nobody re-reads is not a check.)

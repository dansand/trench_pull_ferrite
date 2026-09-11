# Figure ledger — the frozen manuscript figure set (living tracker)

Manuscript is FROZEN (finishing week of 2026-09-11); this is the authoritative figure set.
Working rule: **every script generates only real manuscript figure(s)** — ghost outputs are removed
(superseded sources → local-only `scratch/`). The lean release is exported by `make_release.sh` (repo root)
into the new repo **`trench_pull_ferrite`**, built from an explicit allowlist; `MANIFEST.md` there carries
the exact figure → script → data → command table.

Clutter key: **clean** = one script → one figure · **param** = one script → several *real* figures by
argument (acceptable).

## Main text (2026_version/main.tex)

| Fig | File (in manuscript) | Script | Clutter | Styled this week |
|----|----|----|----|----|
| 1 | `ridge_trench_overview_v2.pdf` | `schematic/ridge_trench_overview_v2.tex` (hand TikZ) | clean | ✓ sidewall sign |
| 2 | `equilibration_compare.pdf` | `schematic/gen_equilibration_compare.py` | clean | ✓ full redesign |
| 3 | `hero_tresca_deep60.png` | `render_hero.py` | param (also SI hero variants) | — |
| 4 | `gpe_correlation.png` | `render_gpe_correlation.py` | clean | — |
| 5 | `gpe_compare_suite1.png` | `render_gpe_compare.py` | **clean** (decluttered 2026-09-11: ghost `gpe_compare.png` + `main_tresca` removed; SI-cited errors now printed) | ✓ B&W restyle |
| 6 | `profiles.png` | `render_profiles.py` | clean | ✓ x_I relabel |
| 7 | `ferrite_thickness_compare.png` | `render_thickness_compare.py` (emits `thickness_compare.png`, renamed on copy) | clean | — |

## Supporting information

| File (in manuscript) | Script | Clutter |
|----|----|----|
| `benchmark_boef.png` | `render_benchmark.py` | clean |
| `benchmark_mp.png` | `render_mp_benchmark.py` | clean |
| convergence table (S2) | `render_convergence.py` | clean — **data `out/paper/convergence/` NOT on disk (regen-only)** |
| `core_profiles_deep60.png` | `render_core_profiles.py` | clean |
| `corrected_density.png` | `render_corrected_density.py` (argv OUT) | clean |
| `edge_effect_SI.png` | `render_edge_si.py` | clean |
| `fig_taux_cases_grid.pdf` | `schematic/taux_cases.tex` (emits `taux_cases.pdf`, renamed) | clean |
| `ferrite_hero_h30.png`, `ferrite_hero_h40.png` | `render_hero.py` (argv variants) | param |
| `hero_dd_vm_asym.png`, `hero_dd_vm_sym.png` | `render_hero.py` (argv variants) | param |
| (App. B number, §4.1) | `render_esweep_test.py` → `esweep.png` | clean — kept: computes the cited E-sweep arm |

## Separation convention (applied automatically)
- **local-only junk** → `scratch/` (gitignored): superseded/duplicate sources, never pushed.
- **release set** → `trench_pull_ferrite` (via `make_release.sh` allowlist): pushed to GitHub/Zenodo.
- **everything else** stays in this repo (`ferrite_plate_flexure`), the private development archive.

## Resolved / done (2026-09-11)
- **`ridge_trench_overview` (non-`_v2`)** — superseded duplicate of Fig 1 → `scratch/retired_figures/`; old tracked copies staged for deletion.
- **`trailing_plate_budget.{tex,pdf}`** — parked TikZ schematic, not in the paper → `scratch/retired_figures/`; staged for deletion.
- **`render_gpe_compare.py`** — decluttered to its one figure; ghost `out/paper/figures/gpe_compare.png` staged for deletion; mypapers copy restored.
- **`render_hero.py`** — confirmed a **parameterised renderer** (same layout for any model dir). Kept as-is.
- **`render_esweep_test.py`** — kept (its number is cited even though its figure is not shown).

## Not in the release (remain in this archive only)
Ghost render scripts: `render_mid_decomp.py`, `render_np_diag.py` (RETIRED from main text 2026-07-08),
`render_suite_overlay.py`, `render_ratio_test.py`, `render_gpe.py` (superseded `fig_2b/2c`),
`render_paper_figures.py`; `_sandbox` dev scripts `render_core_thinning / gradient_compare / levels / overlay`;
stale notebooks/build scripts; `diagnostics/`, `animation/`, `attic/`, `small_strain/`, `legacy/`.

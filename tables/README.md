# tables/ — every number the paper quotes, as CSV

**For the paper agent.** This directory is the handoff between the model analysis and the manuscript. Every number the
main text or the Supporting Information quotes from the models is in one of these files. Read the numbers from here,
never from a figure or from printed output.

## The rule that makes the tables trustworthy

A table is written by the **same script**, in the **same pass**, from the **same in-memory arrays** as the figure it
accompanies. Nothing is computed twice. The three standalone tables (`model_summary`, `isostatic_column`, `convergence`)
are written by their own scripts on the same models with the same functions. `./reproduce.sh` regenerates every table,
asserts that every table and every figure were rewritten in that run, and reads all of its numerical checks from these
files (`tests/check_reproduce.py`). So a table can only be stale if the harness has not been run, and a stale table fails
the harness.

Each CSV is plain (header row, then rows) so GitHub and spreadsheets render it. Its provenance — the script that wrote
it, the figure it accompanies, the model directories it read, and any scalar results — is in a sidecar of the same name,
`tables/<name>.json`. Floats are written with six significant figures.

Read them with `gpe_analysis.read_table(path)` → `(meta, rows)`, or with any CSV reader; the scalars are in the sidecar.

## The tables

| file | figure / record | what each row is | columns (units) | backs in the manuscript |
|---|---|---|---|---|
| `model_summary.csv` | — (`analysis/model_summary.py`) | one production model (20) | suite, model, strength, h_km, V_TN, N_mem_TN, w_T_m, x_T_km, x_I_km, x_M_km, dGPE_TN, dND_TN, identity_pct, arm_km, arm_over_h, max_plastic_strain, hardening_H_Pa, hardening_increment_MPa, yielded_thickness_xM_pct, yielded_face_zone_pct | the reference values ΔGPE* = 2.542 TN/m, identity 0.019 %, arm 34.9 km; any per-model value; the yielded share of the thickness at max M and within 10 km of the face |
| `isostatic_column.csv` | — (`analysis/isostatic_column_test.py`) | one production model (20) | suite, model, h_km, x_I_km, Szz_xI_TN, dGPE_measured_TN, dGPE_if_xI_lithostatic_TN, change_pct, net_charge_at_xI_as_deflection_m, peak_abs_szz_on_xI_MPa | SI: the cost of assuming the first isostatic column lithostatic (reference −0.13 %; up to ~2 % in Suite 1, ~4 % with a background N_D) |
| `frame_check.csv` | — (`analysis/frame_check.py`) | one production model (20) | suite, model, h_km, max_slope_deg, max_cos2theta_departure_pct, mixing_term_at_trench_TN, mixing_first_below_0p02_km, mixing_below_0p02_beyond_km, slope_at_xI_deg, mixing_term_at_xI_TN, x_I_km | SI Text S2 'frame of the resultants': surface slope, cos 2θ departure, the 2V sin 2θ mixing term at the trench, its inboard decay, and at x_I |
| `edge_exclusion.csv` | — (`analysis/edge_exclusion.py`) | one Suite-1 model (4) | label, model, h_km, raw_trench_arm_km, edge_window_km, plateau_arm_km, plateau_min/max_km, change_raw_to_plateau_pct, plateau_arm_over_h, arm_at_40km_km, arm_at_60km_km | SI Text S1 loaded-edge sensitivity: the arm change when the edge window is excluded; the inboard arm values |
| `convergence.csv` | Table S3 (`scripts/render_convergence.py`) | one convergence run (6 rows, 5 models; the 800 × 48 model appears in both the mesh and the load-increment series) | model, configuration, w_T_m, dGPE_TN, identity_residual_pct | Table S3 |
| `benchmark_boef.csv` | Fig. S1 | one quantity | quantity, value, unit | end resultant 4.000 TN/m, V misfit 0.8 %, deflection offset 1.6 %, parabola 0.39 % |
| `benchmark_mp.csv` | Fig. S2 | one quantity | quantity, value, unit | 148 sections, mean 0.09 %, max 0.85 % |
| `gpe_compare_reconstruction.csv` | Fig. 5 | one Suite-1 model (4) | model, label, plate_top_arm_err_pct, sea_level_arm_err_pct, mid_plate_approx_err_pct | the ~2 % plate-top-arm and ~8 % sea-level-arm reconstruction errors; the h/2 misses |
| `gpe_correlation.csv` | Fig. 4 | one plotted model (17) | suite, label, model, V_TN, N_mem_TN, w_T_km, dGPE_TN; meta: uniform_plate_slope_TN_per_km, reference_w_T_km, reference_dGPE_TN | every point in Fig. 4; the uniform-plate slope Δρ g h/2 |
| `thickness_compare.csv` | Fig. 7 | one thickness member (4) | h_km, V_TN, w_T_km, x_I_km, x_M_km, pull_TN, pull_over_V; meta: through_origin_slope_dGPE_over_V | Fig. 7 panel (a) values; ΔGPE* ≈ 0.65 V |
| `profiles_selfcheck.csv` | Fig. 6 | one Suite-1 model (4) | model, label, trench_curve_area_TN, trench_pull_TN, diff_pct | the panel-(d) self-check: trench-curve area equals the pull to < 0.2 % |
| `hero_tresca_deep60.csv`, `hero_tresca_30km.csv`, `hero_tresca_40km.csv`, `hero_dd_vm_asym.csv`, `hero_dd_vm_sym.csv` | Fig. 3, S6, S7, S8 and the 40 km supporting render | one quantity | quantity, value, unit | reference-line positions (trench, max M, first isostatic, forebulge) and the ρ̂ centroid depth over the window |

Units: TN = 10¹² N per metre of strike (TN/m); km, m, MPa as named; percentages are relative unless the column says otherwise.
`x_*_km` are distances from the trench. `N_mem_TN` is the background in-plane force (positive = tension). `arm_km` is the
force-based effective arm ΔGPE*/(Δρ g w_T), not a centre of mass.

## The contract (agreed with the paper side, 2026-09-13)

    model analysis → tables/*.csv → tables/paper_numbers.tex → manuscript / SI

Any number presented as a result of these numerical models — including verification, convergence, reconstruction,
sensitivity, fitted scalings and post-processing — is quoted through a macro, whose value comes from a table. No
model-derived value is typed, rounded, updated or inferred in the LaTeX. A required value with no macro is a request to
the model side (a registry entry in `analysis/paper_numbers.py`; the value comes from a table, never from a keyboard).
Rounding is a registry choice, not an editorial one: where the text wants "about 2.5" the registry carries a rounded
macro (`\numRefPullRounded`) beside the precise one (`\numRefPull`). Out of scope by nature: prescribed model inputs,
purely analytical estimates, literature-derived values, and results of other models.

## How the numbers get into the manuscript — `paper_numbers.tex`

The tables are the source; the manuscript consumes them through **one generated LaTeX file**, so no model number is ever
typed into the text:

1. `analysis/paper_numbers.py` (run by `reproduce.sh`) reads the tables and writes `tables/paper_numbers.tex`: one
   `\newcommand{\num<Name>}{<value>}` per quoted number, already at the precision the paper uses (e.g. `\numRefPull` =
   2.54, `\numRefPullThree` = 2.542, `\numMpSections` = 148). `tables/PAPER_NUMBERS.md` lists every macro with its value,
   description, and the place in the manuscript it belongs.
2. Copy `tables/paper_numbers.tex` next to the LaTeX sources (as the figures are copied) and add `\input{paper_numbers}`
   to the preamble of `main.tex` and `si.tex`.
3. Replace each literal with its macro — `2.54~\si{\tera\newton\per\meter}` becomes `\numRefPull~\si{\tera\newton\per\meter}`.
   The placement guide names the file and line for every literal that existed on 2026-09-13.
4. When the package changes, `./reproduce.sh` regenerates the file; re-copy it and recompile. A change to a quoted number
   then reaches the paper without anyone retyping it, and a number the paper needs that has no macro is a request to the
   model side (add it to the registry in `paper_numbers.py`; the value comes from a table, never from a keyboard).

`python analysis/paper_numbers.py --compare <manuscript dir>` re-reads the LaTeX at every registered location and reports
`match` or `DIFFERS` for each; run it before any submission.

## How to validate a manuscript number against the tables

1. Run `./reproduce.sh` (or trust the committed tables — they are committed only after a passing run).
2. Find the number's row and column above; the sidecar `.json` names the figure and models it came from.
3. Quote the value at the paper's precision. If the paper's number is not in any table, it is not a model result this
   package supports — say so rather than reading it off a figure.
4. `tests/check_reproduce.py` lists every tolerance the harness applies; a number the paper quotes at tighter precision
   than the harness checks should be checked by eye against the table.

There are no prose records: a number with no table has no source in this package (the loaded-edge displacement-controlled
figures once quoted from an archive note were removed on 2026-09-14 for that reason).

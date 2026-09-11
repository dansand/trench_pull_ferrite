# Demo: trench pull as a GPE dipole

`demo_trench_pull.ipynb` is a short, runnable tutorial that loads one flexure model with `gpe_analysis.py`
and reproduces the two headline results of the study:

1. **The force-balance identity** ΔN_D = ΔGPE\* — the horizontal driving force equals the deficit in the
   depth-integrated vertical normal stress (they agree to <0.05%).
2. **The corrected-density dipole** — every column carries the same integrated mass, so the trench pull is a
   dipole of the corrected-density *difference* between the trench and an isostatic reference.

Build/refresh the notebook from its source:

```bash
python build_demo_nb.py
jupyter nbconvert --to notebook --execute --inplace demo_trench_pull.ipynb
```

## Committed model subset

Only the **Suite-1 rheology models** are committed (the `.vtu` mesh + `gpe_topo.csv` per model, ~14 MB total),
enough to run the notebook and the single-model / four-rheology figures:

- `out/paper/set1_rheology/tresca_deep_150_60km_V4/` — the locked baseline (Tresca, σY=150 MPa, V=4 TN/m)
- `out/paper/set1_rheology/elastic_deep_60km_V4/`
- `out/paper/set1_rheology/dd_vm_asym_60km_V4/`
- `out/paper/set1_rheology/dd_vm_sym_60km_V4/`

These are tracked despite the `out/` entry in `.gitignore` (force-added). The rest of the study — the N_D
sweep (`set2_nd_sweep`), V sweep (`set3_v_sweep`), and E sweep (`esweep`) — is **not** committed; the
multi-suite figures (`render_gpe_correlation.py`, `render_suite_overlay.py`, `render_esweep_test.py`) need
those. Regenerate the full study from the driver:

```bash
julia --project=. paper_models.jl suite1   # rheology suite + N_D sweep + E sweep
julia --project=. paper_models.jl v_sweep     # V sweep
```

See `out/paper/MODELS.md` for the full model list, values, and figure↔model mapping.

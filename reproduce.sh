#!/usr/bin/env bash
# reproduce.sh — regenerate and CHECK every manuscript figure from the shipped data, in one command, from the repo root.
#
#   ./reproduce.sh              # everything: tests, every REPRODUCE.md command in order, the notebook, the number checks
#   ./reproduce.sh --no-nb      # skip executing START_HERE.ipynb (the slowest step)
#   The three TikZ schematics (Figs 1, 2, S3) are not model output and are built separately — see REPRODUCE.md.
#   Table S3 needs data/convergence/ (shipped; regenerate with `julia --project=. model/paper_models.jl convergence`).
#
# Runs in the ferrite-figs env; PYTHON=/path/to/python overrides the interpreter for the scripts AND the notebook
# kernel (which must have ipykernel installed, as the pinned env does).  Every command's output is logged to
# .reproduce_logs/<step>.log; tests/check_reproduce.py then asserts the headline numbers against their tolerances
# (ΔGPE* = 2.542 TN/m ± 0.5 %, identity < 0.05 %, arm 34.9 km, the S1/S2 benchmark misfits, the reconstruction table)
# and that every figure was rewritten.  Exit status is nonzero on ANY failure.
set -uo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python}
LOG=.reproduce_logs; rm -rf "$LOG"; mkdir -p "$LOG"
date +%s > "$LOG/START"
fail=0
run() {                       # run <name> <command...>  — logs, reports, never stops the sequence
    local name=$1; shift
    printf '%-20s ' "$name"
    if "$@" > "$LOG/$name.log" 2>&1; then echo ok; else echo "FAILED (see $LOG/$name.log)"; fail=1; fi
}

run tests             "$PY" -m pytest tests/ -q -p no:cacheprovider
run hero              "$PY" scripts/render_hero.py
run hero_full         env HERO_PANELS=dsr "$PY" scripts/render_hero.py data/suite1_strength/tresca_deep_150_60km_V4 figures/hero_tresca_deep60_full.png
run gpe_correlation   "$PY" scripts/render_gpe_correlation.py
run gpe_compare       "$PY" scripts/render_gpe_compare.py
run profiles          "$PY" scripts/render_profiles.py
run thickness_compare "$PY" scripts/render_thickness_compare.py
run benchmark         "$PY" scripts/render_benchmark.py
run mp_benchmark      "$PY" scripts/render_mp_benchmark.py
run convergence       "$PY" scripts/render_convergence.py
run core_profiles     "$PY" scripts/render_core_profiles.py
run corrected_density "$PY" scripts/render_corrected_density.py
run hero_h30          env HERO_PANELS=dsr HERO_WINDOW_KM=200 HERO_MFIX=0 "$PY" scripts/render_hero.py data/suite4_thickness/tresca_150_30km figures/hero_tresca_30km.png
run hero_h40          env HERO_PANELS=dsr HERO_WINDOW_KM=260 HERO_MFIX=0 "$PY" scripts/render_hero.py data/suite4_thickness/tresca_150_40km figures/hero_tresca_40km.png
run hero_dd_vm_asym   env HERO_PANELS=dsr "$PY" scripts/render_hero.py data/suite1_strength/dd_vm_asym_60km_V4 figures/hero_dd_vm_asym.png
run hero_dd_vm_sym    env HERO_PANELS=dsr "$PY" scripts/render_hero.py data/suite1_strength/dd_vm_sym_60km_V4 figures/hero_dd_vm_sym.png
run isostatic_column  "$PY" analysis/isostatic_column_test.py
run model_summary     "$PY" analysis/model_summary.py
run frame_check       "$PY" analysis/frame_check.py
run edge_exclusion    "$PY" analysis/edge_exclusion.py
run paper_numbers     "$PY" analysis/paper_numbers.py
run manifest          "$PY" analysis/make_manifest.py --check
no_nb=0
for a in "$@"; do case "$a" in --no-nb) no_nb=1;; *) echo "unknown option $a"; exit 2;; esac; done
if [[ $no_nb -eq 0 ]]; then
    # The notebook runs in a kernel built on THIS interpreter ($PY), not whichever "python3" kernelspec happens to be
    # registered: a throwaway kernelspec under $LOG is put first on JUPYTER_PATH and selected by name.
    PYABS=$(command -v "$PY"); mkdir -p "$LOG/kernels/ferrite-reproduce"
    printf '{"argv": ["%s", "-m", "ipykernel_launcher", "-f", "{connection_file}"], "display_name": "ferrite-reproduce", "language": "python"}\n' \
        "$PYABS" > "$LOG/kernels/ferrite-reproduce/kernel.json"
    run notebook env JUPYTER_PATH="$(pwd)/$LOG" "$PY" -m jupyter nbconvert --to notebook --execute START_HERE.ipynb \
                 --ExecutePreprocessor.kernel_name=ferrite-reproduce --ExecutePreprocessor.timeout=1800 \
                 --output START_HERE_executed.ipynb --output-dir "$LOG"
else
    touch "$LOG/NOTEBOOK_SKIPPED"
fi

echo; "$PY" tests/check_reproduce.py "$LOG" || fail=1
if [[ $fail -eq 0 ]]; then echo; echo "reproduce.sh: ALL CHECKS PASSED"; else echo; echo "reproduce.sh: FAILURES (see above and $LOG/)"; fi
exit $fail

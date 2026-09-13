#!/usr/bin/env bash
# reproduce.sh — regenerate and CHECK every manuscript figure from the shipped data, in one command, from the repo root.
#
#   ./reproduce.sh              # everything: tests, every REPRODUCE.md command in order, the notebook, the number checks
#   ./reproduce.sh --no-nb      # skip executing START_HERE.ipynb (the slowest step)
#
# Runs in the ferrite-figs env; PYTHON=/path/to/python overrides the interpreter for the scripts AND the notebook
# kernel (which must have ipykernel installed, as the pinned env does).  Every command's output is logged to
# .reproduce_logs/<step>.log; tests/check_reproduce.py then asserts the headline numbers against their tolerances
# (ΔGPE* = 2.542 TN/m ± 0.5 %, identity < 0.05 %, arm 34.9 km, the S1/S2 benchmark misfits, the reconstruction table)
# and that every figure was rewritten.  Exit status is nonzero on ANY failure.  The convergence table (Table S3) needs
# data/convergence/, which is regeneration-only, so it is not part of this run (data/CONVERGENCE.md is the record).
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
run gpe_correlation   "$PY" scripts/render_gpe_correlation.py
run gpe_compare       "$PY" scripts/render_gpe_compare.py
run profiles          "$PY" scripts/render_profiles.py
run thickness_compare "$PY" scripts/render_thickness_compare.py
run benchmark         "$PY" scripts/render_benchmark.py
run mp_benchmark      "$PY" scripts/render_mp_benchmark.py
run core_profiles     "$PY" scripts/render_core_profiles.py
run corrected_density "$PY" scripts/render_corrected_density.py
run hero_h30          env HERO_WINDOW_KM=200 "$PY" scripts/render_hero.py data/suite4_thickness/tresca_150_30km figures/hero_tresca_30km.png
run hero_h40          env HERO_WINDOW_KM=260 HERO_MFIX=0 "$PY" scripts/render_hero.py data/suite4_thickness/tresca_150_40km figures/hero_tresca_40km.png
run hero_dd_vm_asym   "$PY" scripts/render_hero.py data/suite1_strength/dd_vm_asym_60km_V4 figures/hero_dd_vm_asym.png
run hero_dd_vm_sym    "$PY" scripts/render_hero.py data/suite1_strength/dd_vm_sym_60km_V4 figures/hero_dd_vm_sym.png
run isostatic_column  "$PY" analysis/isostatic_column_test.py
run manifest          "$PY" analysis/make_manifest.py --check
if command -v tectonic > /dev/null; then
    run schematic bash -c "cd schematic && $PY gen_equilibration_compare.py && tectonic equilibration_compare.tex && tectonic ridge_trench_overview_v2.tex && tectonic taux_cases.tex"
else
    echo "schematic            skipped (tectonic not on PATH)"
fi
if [[ "${1:-}" != "--no-nb" ]]; then
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

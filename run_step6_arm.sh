#!/usr/bin/env bash
# STAGE 1.5.3 step 6 (G-c step 2) -- ONE ARM: build, train, snapshot.
#   ./run_step6_arm.sh A    # control  : AGE_MASKED_DATASETS = frozenset()
#   ./run_step6_arm.sh B    # treatment: AGE_MASKED_DATASETS = {"hff_sc"}
#   ./run_step6_arm.sh C 0  # control  : HFF labels SHUFFLED, seed 0 (arg 2)
#
# The arms write SEPARATE fold roots (cellfate_loocv_<donor>_armA/_armB) via
# CELLFATE_FOLD_SUFFIX, which scorecard.py honours too. The first step-6 run let arm B
# overwrite arm A, which cost arm A's scalers.json -- its deconfounder coefficient had to
# be reported from a proxy build. Both arms now survive on disk (~1.6 GB per arm).
set -euo pipefail
ARM="${1:?usage: run_step6_arm.sh A|B}"
SEED="${2:-0}"
case "$ARM" in
  A) TAG=gc2_A_keep_hff ;;
  B) TAG=gc2_B_mask_hff ;;
  C) TAG="gc2_C_shuffle_hff_s${SEED}" ;;   # label-permutation control; seed goes in the tag
  *) echo "arm must be A, B or C" >&2; exit 2 ;;
esac
PY=/d/.venv-cellfate/Scripts/python.exe
export PYTHONUTF8=1
export CELLFATE_FOLD_SUFFIX="_arm${ARM}"
echo "=== step 6 arm $ARM | fold roots cellfate_loocv_<donor>${CELLFATE_FOLD_SUFFIX} -> snapshot '$TAG' ==="
"$PY" local_runners/run_loocv.py "D:\GSE242423" "D:\Gill" --arm "$ARM" --age-window-k 4       --shuffle-seed "$SEED"
"$PY" scorecard.py snapshot --tag "$TAG"
echo "=== arm $ARM done and snapshotted as $TAG ==="

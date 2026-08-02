#!/usr/bin/env bash
# STAGE 1.5.3 step 6 (G-c step 2) -- ONE ARM, then its snapshot.
#   ./run_step6_arm.sh A    # control  : AGE_MASKED_DATASETS = frozenset()
#   ./run_step6_arm.sh B    # treatment: AGE_MASKED_DATASETS = {"hff_sc"}
#
# Arm B OVERWRITES arm A's builds (scorecard.py:132 resolves cellfate_loocv_<donor> exactly),
# so the snapshot is chained onto the run rather than left as a separate step a human must
# remember. Losing arm A's snapshot means repeating hours of compute.
set -euo pipefail
ARM="${1:?usage: run_step6_arm.sh A|B}"
case "$ARM" in
  A) TAG=gc2_A_keep_hff ;;
  B) TAG=gc2_B_mask_hff ;;
  *) echo "arm must be A or B" >&2; exit 2 ;;
esac
PY=/d/.venv-cellfate/Scripts/python.exe
export PYTHONUTF8=1
echo "=== step 6 arm $ARM -> snapshot '$TAG' ==="
"$PY" local_runners/run_loocv.py "D:\GSE242423" "D:\Gill" --arm "$ARM" --age-window-k 4
"$PY" scorecard.py snapshot --tag "$TAG"
echo "=== arm $ARM done and snapshotted as $TAG ==="

# HOW TO RUN THE SCRIPTS IN THIS FOLDER

## What is in here, and what is not

`plan_tests/` holds the **per-stage verification gates** — the scripts a plan in `plans/` names as
the thing that decides whether that stage passed. They answer one question each, print a verdict,
and write a JSON report into `results/`.

| script | the plan it gates | the question it answers |
|---|---|---|
| `verify_1a.py` | `plans/STAGE_1_CALIBRATION.md` (S1a) | did the donor label column land, and is inner-LODO possible? |
| `verify_stage1_5.py` | `plans/STAGE_1_5_HARMONIZATION_AUDIT.md` | did the silent ΔAge zero-point fallback fire on the real build? |
| `smoke_stage1.py` | `plans/STAGE_1_CALIBRATION.md` | end-to-end smoke test on synthetic data shaped like the real failure |

**Not in here, on purpose:**

| | where | why |
|---|---|---|
| `pytest` unit tests | `tests/` | they run in CI, need no data, and are a different thing entirely |
| exploratory / numbered tests | `experiments/` | `test18_forward_gate.py`, `test5_ridge_gap.py`, the `diag_*` scripts |
| operational tools | repo root | `scorecard.py`, `retrain_stage1.py`, `audit_metrics.py` — you run these constantly and they are imported by other code |

## The one rule: run from the REPO ROOT

Same as `experiments/` and `local_runners/`. These scripts use paths relative to the repo root —
`src/`, `configs/`, `runs/`, `results/` — so the working directory must be `D:\cellfate-rx`.

```powershell
cd D:\cellfate-rx
D:\.venv-cellfate\Scripts\Activate.ps1       # prompt must show (.venv-cellfate)

python plan_tests\verify_stage1_5.py "D:\GSE242423" "D:\Gill"
python plan_tests\verify_1a.py
python plan_tests\smoke_stage1.py
```

Launching from *inside* `plan_tests\` makes them look for `src\` and `configs\` in here, not find
them, and fail with a `ModuleNotFoundError` or a "not found". Nothing is broken; it is looking in
the wrong place.

## Where the output goes

All three write their JSON report to **`results/`** at the repo root — `verify_1a_results.json`,
`verify_stage1_5_results.json`. Those paths are `__file__`-relative, not
working-directory-relative, so the report lands in the same place regardless of where you launched
from. Only the *reading* of `src/` and `configs/` depends on the working directory.

## A note for anyone editing these

`verify_stage1_5.py` is **imported by `tests/test_baseline_census.py`**, which puts this directory
on `sys.path`. Everything above `main()` in that file is deliberately data-free and safe to import;
keep it that way. If you move or rename it, `tests/test_baseline_census.py` breaks.

"""STAGE 3a — THE GATE, run on C-7-clean labels (and on the contaminated ones, for contrast).

    python experiments/stage3a_forward_gate.py                 # both arms if both exist
    python experiments/stage3a_forward_gate.py _c7             # one arm

READ-ONLY. Writes `results/stage3a_forward_gate_results.json`. `src/` untouched, no build
touched. `experiments/test18_forward_gate.py` is imported UNMODIFIED; only its `resolve_root`
is redirected, exactly as the Test 7 reproduction did.

WHY THIS WRAPPER EXISTS RATHER THAN JUST RUNNING TEST 18
---------------------------------------------------------
`STAGE_3_TOOL.md` §3a's STOP branch is **terminal**: *"do not write tool code. Ship the scoring
model; go to Stage 5."* Its target is `y_age` for the held-out Gill donor, pooled across folds
-- and STAGE_1_5_6 §5.14 proved that quantity is inflated ~3x in **five of six** arm-A folds by
one degenerate GEO column, with a 16.67 yr between-fold spread that is not biology.

**A terminal decision must not be taken on labels known to be contaminated.** So 3a runs on the
`_c7` folds, where the degenerate control is rejected and donor N2's ΔAge is masked by rule 4.

The contaminated arm is run too, and reported beside it, for one reason only: **if the two arms
disagree, the disagreement is itself the finding** -- it would mean the label defect was strong
enough to flip a terminal gate, which is worth recording whichever way 3a lands. The `_c7` arm
is the operative verdict; `_armA` is context, never the decision.

An earlier proposal -- run 3a twice, once with donor N2 excluded -- is **withdrawn and must not
be revived**. Excluding N2 removes N2's *rows* but leaves N2's *control* inside the harmonizer
of every fold that does not hold N2 out, so it drops the one fold whose harmonizer is clean and
keeps the five contaminated ones. Exactly backwards. Recorded in §5.16.

THE BAR IS §3a's OWN, UNCHANGED
--------------------------------
    GO       Δt improves prediction beyond noise AND the sweep moves > 2 yr  -> build 3b-3d
    WEAK GO  the sweep moves but Δt does not beat state-only                 -> build, tempered
    STOP     neither                                                        -> do not write
                                                                               tool code
Nothing here re-derives or relaxes it.
"""
from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(REPO / "experiments"))
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

OUT = _RESULTS / "stage3a_forward_gate_results.json"
ARMS_DEFAULT = ("_c7", "_armA")
OPERATIVE = "_c7"        # the clean arm; the one the verdict is taken from


def _run_arm(suffix: str) -> dict | None:
    """Import test 18 UNMODIFIED, point it at one arm's fold roots, run its own `main`.

    `main` is called rather than its parts re-assembled here: test 18 owns the LODO ridge, the
    Δt sweep and the GO/WEAK GO/STOP rule, and re-implementing any of that in a wrapper is how
    a wrapper silently becomes a different test. Its printed report is captured verbatim so the
    verdict recorded is the one test 18 actually produced.
    """
    import contextlib
    import io as _io

    if "test18_forward_gate" in sys.modules:
        del sys.modules["test18_forward_gate"]
    t18 = importlib.import_module("test18_forward_gate")
    t18.resolve_root = lambda name, _s=suffix: str(REPO / f"{name}{_s}")

    missing = [d for d in t18.DONORS if not (REPO / f"cellfate_loocv_{d}{suffix}").exists()]
    if missing:
        print(f"   [{suffix}] missing fold roots: {missing} -- skipped")
        return None

    buf = _io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            t18.main()
        err = None
    except Exception as exc:                                        # noqa: BLE001
        err = f"{type(exc).__name__}: {exc}"
    report = buf.getvalue()
    print(report)
    if err:
        print(f"   [{suffix}] FAILED: {err}")

    verdict = None
    for token in ("STOP", "WEAK GO", "GO"):
        if token in report:
            verdict = token
            break
    return {"suffix": suffix, "verdict": verdict, "error": err, "report": report}


def main() -> int:
    from cellfate.common.console import install_pretty_console
    install_pretty_console()

    arms = tuple(sys.argv[1:]) or ARMS_DEFAULT
    print("\n" + "=" * 78)
    print("STAGE 3a — the forward-Δt gate")
    print("=" * 78)
    print(f"operative arm: {OPERATIVE} (C-7 clean).  arms requested: {list(arms)}")
    print("A STOP here is TERMINAL, so it is taken on clean labels only; the contaminated")
    print("arm is context. See STAGE_1_5_6 §5.14 for why arm A's y_age cannot decide this.")

    out: dict = {"script": "stage3a_forward_gate", "operative_arm": OPERATIVE, "arms": {}}
    for suffix in arms:
        print(f"\n--- arm {suffix} ---")
        res = _run_arm(suffix)
        if res is not None:
            out["arms"][suffix] = res
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(f"\n   wrote {OUT.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

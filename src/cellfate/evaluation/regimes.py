"""The three held-out generalisation regimes (Document 5, S1).

Each names a distinct way of holding data out; the evaluator reports every regime
separately, never a single pooled number (success criterion 3).
"""

from __future__ import annotations

# regime key (as written in splits/) -> human description
REGIMES = {
    "scaffold": "leave-drug-out (unseen chemical scaffolds)",
    "cell_line": "leave-cell-line-out (unseen cell lines)",
    "both": "both-unseen (unseen scaffold x unseen cell line)",
    "random": "random cell-level split (fit validation, NOT generalization)",
    "holdout": "leave-named-cell-line(s)-out for test; rest cell-level train/val/calib",
}


def describe(regime: str) -> str:
    return REGIMES.get(regime, regime)


def iter_regimes(selected):
    for r in selected:
        if r not in REGIMES:
            raise ValueError(f"unknown regime {r!r}; choose from {tuple(REGIMES)}")
        yield r

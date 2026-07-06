"""Leaderboards, reliability diagrams, and JSON/Markdown reports (Document 5, S5)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .metrics import mean_finite

# key metrics shown on the leaderboard (flat keys produced by the evaluator)
_LEADER_COLS = ("mean_prauc", "mean_auroc", "brier", "ece", "reg_mae", "reg_rmse", "reg_pearson")


def summarise_estimator(flat: dict) -> dict:
    """Collapse an estimator's flat metric dict into leaderboard columns."""
    return {
        "mean_prauc": mean_finite(flat.get(f"prauc_{c}") for c in range(3)),
        "mean_auroc": mean_finite(flat.get(f"auroc_{c}") for c in range(3)),
        "brier": flat.get("brier", float("nan")),
        "ece": flat.get("ece", float("nan")),
        "reg_mae": flat.get("reg_mae", float("nan")),
        "reg_rmse": flat.get("reg_rmse", float("nan")),
        "reg_pearson": flat.get("reg_pearson", float("nan")),
    }


def leaderboard(estimator_metrics: dict[str, dict]) -> list[tuple[str, dict]]:
    """(name, summary) rows sorted by mean PR-AUC descending (model included)."""
    rows = [(name, summarise_estimator(m)) for name, m in estimator_metrics.items()]
    rows.sort(key=lambda kv: (-(kv[1]["mean_prauc"] if np.isfinite(kv[1]["mean_prauc"]) else -1)))
    return rows


def _fmt(v) -> str:
    return f"{v:.4f}" if isinstance(v, (int, float)) and np.isfinite(v) else "n/a"


def _markdown_table(rows: list[tuple[str, dict]]) -> str:
    head = "| estimator | " + " | ".join(_LEADER_COLS) + " |\n"
    head += "|" + "---|" * (len(_LEADER_COLS) + 1) + "\n"
    body = "".join(
        "| " + name + " | " + " | ".join(_fmt(summ[c]) for c in _LEADER_COLS) + " |\n"
        for name, summ in rows
    )
    return head + body


def reliability_diagram(y_true, p, path, n_bins: int = 15) -> bool:
    """Save a reliability diagram PNG. Returns False (without raising) if matplotlib
    is unavailable, so reporting degrades gracefully to JSON/MD only."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - optional dependency
        return False
    y_true = np.asarray(y_true)
    conf, pred = p.max(1), p.argmax(1)
    acc = (pred == y_true).astype(float)
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    xs, ys = [], []
    for lo, hi in zip(bins[:-1], bins[1:], strict=True):
        m = (conf > lo) & (conf <= hi)
        if m.any():
            xs.append(conf[m].mean())
            ys.append(acc[m].mean())
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot([0, 1], [0, 1], "--", color="gray", label="perfect")
    ax.plot(xs, ys, "o-", label="model")
    ax.set_xlabel("confidence")
    ax.set_ylabel("accuracy")
    ax.set_title("Reliability")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return True


def write_report(regime: str, R: dict, out_dir) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{regime}.json").write_text(json.dumps(R, indent=2, default=float))

    estimator_metrics = {k: v for k, v in R.items() if isinstance(v, dict) and "ece" in v}
    rows = leaderboard(estimator_metrics)
    md = [f"# Evaluation report - regime: {regime}\n",
          "## Baseline leaderboard (sorted by mean PR-AUC)\n",
          _markdown_table(rows)]
    if "coverage" in R:
        md.append(f"\n**Conformal coverage (model):** {_fmt(R['coverage'])}\n")
    if "ranking" in R:
        md.append(f"**RES ranking (model):** Spearman {_fmt(R['ranking']['spearman'])}, "
                  f"precision@k {_fmt(R['ranking']['precision_at_k'])}\n")
    (out / f"{regime}.md").write_text("".join(md))

    if "model" in estimator_metrics and "_y_true" in R and "_p_model" in R:
        reliability_diagram(R["_y_true"], np.asarray(R["_p_model"]), out / f"{regime}_reliability.png")
    return out / f"{regime}.json"


def write_summary(results: dict, gates: dict, out_dir) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    payload = {"gates": gates,
               "regimes": {r: {k: v for k, v in R.items() if not k.startswith("_")}
                           for r, R in results.items() if r != "external"},
               "external": results.get("external", {})}
    (out / "summary.json").write_text(json.dumps(payload, indent=2, default=float))

    lines = ["# CellFate-Rx - acceptance summary\n"]
    for regime, g in gates.items():
        passed = all(g.values())
        lines.append(f"\n## {regime} - {'PASS' if passed else 'FAIL'}\n")
        for crit, ok in g.items():
            lines.append(f"- {'[x]' if ok else '[ ]'} {crit}\n")
    ext = results.get("external", {})
    lines.append(f"\n## external validation\n- status: {ext.get('status', 'n/a')}\n")
    (out / "summary.md").write_text("".join(lines))
    return out / "summary.json"

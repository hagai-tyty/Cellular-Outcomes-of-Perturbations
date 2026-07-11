"""Opt-in human-readable console output for the run scripts.

The library logs structured JSON lines (machine-parseable). That is great for
parsing but noisy to watch -- e.g. one JSON line per shard. Calling
``install_pretty_console()`` swaps the JSON stream handler for one that:

  * collapses high-frequency ``chunk.done`` events into a single in-place
    progress line (``building shards... 32 done``) instead of 51 JSON lines, and
  * prints every other event as a clean one-liner instead of raw JSON.

This changes ONLY presentation; the structured events and all pipeline logic are
untouched. Library default (JSON) is unchanged, so tests are unaffected.
"""
from __future__ import annotations

import logging
import sys

# friendly one-line templates for the common events (fallback is generic)
_PRETTY = {
    "panel.fit":        lambda f: f"gene panel ready ({f.get('n')} genes)",
    "panel.loaded":     lambda f: f"gene panel loaded ({f.get('n')} genes)",
    "harmonizer.fit":   lambda f: (f"harmonizer fit on {f.get('datasets')} "
                                   f"({f.get('n_genes')} genes; held out {f.get('heldout')}, "
                                   f"{f.get('heldout_controls_excluded')} control(s) excluded)"),
    "dataset.done":     lambda f: f"dataset built: {f.get('n_samples'):,} cells in {f.get('n_shards')} shards",
    "data.loaded":      lambda f: (f"data loaded: {f.get('n_train'):,} train / "
                                   f"{f.get('n_val'):,} val / {f.get('n_calib'):,} calib"),
    "member.trained":   lambda f: f"model {int(f.get('idx', 0)) + 1} trained (val_loss {f.get('val_loss'):.3f})",
    "bundle.done":      lambda f: f"ensemble ready ({f.get('n_members')} models, temp {f.get('temperature'):.3f})",
}


class _PrettyHandler(logging.Handler):
    def __init__(self) -> None:
        super().__init__()
        self._built = 0
        self._in_progress = False

    def emit(self, record: logging.LogRecord) -> None:
        fields = getattr(record, "extra_fields", {}) or {}
        msg = record.getMessage()
        out = sys.stdout
        try:
            if msg == "chunk.done":                       # collapse into one in-place line
                self._built += 1
                n = fields.get("n", 0)
                out.write(f"\r   building shards... {self._built} done "
                          f"({n} cells in last)      ")
                out.flush()
                self._in_progress = True
                return
            # any non-progress event: end the progress line first
            if self._in_progress:
                out.write("\n")
                self._in_progress = False
                self._built = 0
            if msg == "chunk.failed":
                out.write(f"   ! chunk failed: {fields.get('chunk')} ({fields.get('err')})\n")
            else:
                render = _PRETTY.get(msg)
                if render is not None:
                    out.write(f"   \u2713 {render(fields)}\n")
                else:                                      # generic clean fallback (no raw JSON)
                    extra = "  ".join(f"{k}={v}" for k, v in fields.items())
                    out.write(f"   \u00b7 {msg}{('  ' + extra) if extra else ''}\n")
            out.flush()
        except Exception:  # never let logging crash the run
            self.handleError(record)


def render_table(headers: list[str], rows: list[list[str]], title: str = "",
                 aligns: list[str] | None = None) -> str:
    """Render a clean fixed-width box-drawing table (copy-paste safe in any terminal).

    ``aligns`` is a per-column 'l'/'r' list (default: first column left, rest right).
    Column widths are computed from the actual content so borders always line up.
    """
    ncol = len(headers)
    aligns = aligns or (["l"] + ["r"] * (ncol - 1))
    cells = [[str(c) for c in row] for row in rows]
    widths = [len(headers[i]) for i in range(ncol)]
    for row in cells:
        for i in range(ncol):
            widths[i] = max(widths[i], len(row[i]))

    def fmt(row: list[str]) -> str:
        out = []
        for i in range(ncol):
            out.append(row[i].rjust(widths[i]) if aligns[i] == "r" else row[i].ljust(widths[i]))
        return "│ " + " │ ".join(out) + " │"

    top = "┌─" + "─┬─".join("─" * w for w in widths) + "─┐"
    sep = "├─" + "─┼─".join("─" * w for w in widths) + "─┤"
    bot = "└─" + "─┴─".join("─" * w for w in widths) + "─┘"
    lines = []
    if title:
        lines.append(title)
    lines += [top, fmt(headers), sep]
    lines += [fmt(r) for r in cells]
    lines.append(bot)
    return "\n".join(lines)


def install_pretty_console(logger_name: str = "cellfate") -> None:
    """Replace the JSON stream handler on the cellfate logger tree with the
    pretty handler. Safe to call once at the top of a run script."""
    names = (logger_name, "cellfate.data", "cellfate.training",
             "cellfate.evaluation", "cellfate.models", "cellfate.inference")
    shared = _PrettyHandler()          # one handler so the progress counter is shared
    for name in names:
        lg = logging.getLogger(name)
        lg.handlers = [h for h in lg.handlers if not isinstance(h, logging.StreamHandler)]
        lg.addHandler(shared)
        lg.setLevel(logging.INFO)
        lg.propagate = False


def _demo() -> None:
    """Print a sample of the pretty console UI using representative numbers, so
    the output can be previewed without running the real (~20 min) pipeline.

    Run:  python show_ui.py        (or)   python -m cellfate.common.console --show
    """
    import time

    from cellfate.common.logging import get_logger, log_event

    install_pretty_console()
    log = get_logger("cellfate.data")
    tlog = get_logger("cellfate.training")

    print("  (this is a PREVIEW of the console UI — no real run, numbers are illustrative)\n")
    print("[env] torch 2.5.1+cu121 | CUDA: True | GPU: NVIDIA GeForce GTX 1080")
    print("[data] Gill donors: ['N2','N3','O1','O2','Y1','Y2']  ->  holding out 'O1' as TEST")
    print("\n[1/5] BUILD combined dataset (regime=holdout, holdout=O1, harmonize=ON) — streaming ...")
    log_event(log, "panel.fit", n=2000, hash="783f269a214aa972")
    log_event(log, "harmonizer.fit", datasets=["gill_bulk", "hff_sc"], n_genes=5328,
              heldout=["O1"], heldout_controls_excluded=1)
    for i in range(51):                                   # animated in-place progress bar
        log_event(log, "chunk.done", chunk=f"reprogramming:HFF:b{i}", n=940 + (i % 15))
        time.sleep(0.03)
    log_event(log, "dataset.done", n_samples=42605, n_shards=51)

    # ---- sample [sanity] table (real O1-fold shape) ------------------------- #
    print("\n[sanity] where each cell line's cells landed, and their average predictions")
    srows = [["HFF", "42,481", "calib=4478 train=33613 val=4390", "0.436", "-8.27"],
             ["N2", "21", "calib=4 train=14 val=3", "0.007", "-14.44"],
             ["N3", "21", "calib=1 train=16 val=4", "0.270", "+30.15"],
             ["O1", "21", "test=21", "0.249", "+14.74  <- TEST"],
             ["O2", "21", "calib=2 train=18 val=1", "0.291", "+14.22"],
             ["Y1", "19", "calib=3 train=13 val=3", "0.222", "+14.32"],
             ["Y2", "21", "calib=2 train=14 val=5", "0.240", "+32.53"]]
    print(render_table(
        ["cell line", "cells", "where they went (train/val/calib/test)", "P(loss)", "ΔAge (yrs)"],
        srows, aligns=["l", "r", "l", "r", "r"]))
    print("   P(loss) = avg predicted chance the cell loses its identity (0-1).  "
          "ΔAge = avg years younger(−)/older(+).")
    print("\n   >>> HELD-OUT (test) donor: 'O1'  <<<  (model never trains on it)")

    print("\n[2/5] TRAIN 5-member ensemble on cuda (leave-cell-line-out) ...")
    log_event(tlog, "data.loaded", g=2000, n_train=33688, n_val=4406, n_calib=4490)
    for idx, vl in enumerate([7.976, 7.936, 7.924, 7.920, 7.973]):
        log_event(tlog, "member.trained", idx=idx, seed=idx, val_loss=vl)
        time.sleep(0.15)
    log_event(tlog, "bundle.done", n_members=5, temperature=0.5454)
    print("\n[check] trained on 33,688 cells | validated on 4,406 | calibrated on 4,490")
    print("        (temperature 0.55 = confidence rescaling; "
          "conformal_q 9.68 = ±interval width in years)")

    # ---- sample evaluation output ------------------------------------------ #
    print("\n[3/5] EVALUATE  (held-out donor — model never saw it) ...")
    print("\n   gates:  ranking PASS   calibration ----   coverage PASS   beats-all ----")
    trows = [["model", "n/a", "0.274", "5.39", "← our model"],
             ["mean", "n/a", "0.280", "28.38", ""],
             ["ridge", "n/a", "0.108", "8.25", "model wins MAE"],
             ["x_only", "n/a", "0.105", "8.03", "model wins MAE"],
             ["knn", "n/a", "0.254", "25.47", "model wins MAE"]]
    print("\n" + render_table(
        ["estimator", "PR-AUC", "ECE", "ΔAge MAE", "vs model"], trows,
        aligns=["l", "r", "r", "r", "l"],
        title="   MODEL vs BASELINES on the held-out donor:"))
    print("   PR-AUC: fate accuracy, higher better (0-1).   "
          "ECE: calibration error, lower better (<0.05 good).")
    print("   ΔAge MAE: avg error in years, lower better.   "
          "'model' is our tool; the rest are simple baselines to beat.")
    print("\n   WHAT THIS MEANS:")
    print("     • Ranking (the tool's main job): STRONG  (Spearman 0.68) — "
          "can it order perturbations correctly?")
    print("     • ΔAge magnitude: MAE 5.4  — beats the linear baseline (ridge) on this donor")
    print("     • Calibration on this unseen donor: good  (coverage 0.81, ECE 0.27)")
    print("     Note: one held-out donor is noisy; the leave-one-donor-out run "
          "(run_loocv.py) gives the honest mean±std.")
    print("\n[5/5] DONE.  (this was a preview — run 'python run_multi_local.py ...' for the real thing)")


if __name__ == "__main__":
    import sys
    if "--show" in sys.argv or "-show" in sys.argv or len(sys.argv) == 1:
        _demo()
    else:
        print("usage: python -m cellfate.common.console --show")


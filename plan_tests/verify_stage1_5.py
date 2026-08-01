"""
STAGE 1.5 VERIFICATION — did the silent ΔAge zero-point fallback fire on the real build?

    python plan_tests/verify_stage1_5.py                       # defaults: D:\\GSE242423  D:\\Gill
    python plan_tests/verify_stage1_5.py "D:\\GSE242423" "D:\\Gill"

WHY THIS EXISTS (STAGE_1_5_HARMONIZATION_AUDIT.md §0.2, §2 Group E).
ΔAge is control-relative: `ΔAge = age - mean(age over that line's vehicle controls)`. But
`aging.py:88` has a silent fallback — if a line lands in a chunk with **no** controls, its
zero-point flips from control-relative to self-centred, forcing that line's mean ΔAge toward 0
with no warning. Stage 2 then spends wet-lab money on the premise that the ±12.7 yr per-donor
offset is real biology. An offset that survives is EITHER that biology OR this fallback firing,
and nothing else distinguishes them. This gate settles it, cheaply, before the money is spent.

WHAT IT DOES. Replays `plan_all(sources)` (the same chunking the build used — deterministic),
fetches each chunk, and counts vehicle controls per chunk. A chunk that has perturbed cells but
**zero** controls fired the fallback.

    every chunk with perturbed cells also has >=1 control  ->  PASS  (offset is not an artefact)
    any such chunk has zero controls                        ->  FAIL  (a FINDING, not a bug to
                                                                        patch here — see §3)

Per §1, a FAIL is recorded in STAGE_1_DEVIATIONS.md and fixed under its OWN pre-registered Change;
"fixing" the zero-point here would move every target and invalidate the guards.

DESIGN. Following the verify_1a lesson (a decision function whose only exercised path is the one
that says PASS is not a gate), `decide_verdict()` is a PURE function separated from all I/O and is
unit-tested on every branch in `tests/test_harmonize.py`. This script only wires real data into it.

NOTE ON THE REAL-DATA WIRING (review this). Source construction mirrors `run_multi_local.py`
(GSE242423 + Gill via its `discover_*` helpers). If the runner's source setup changes, update
`build_sources()` below. Everything above `main()` is data-free and safe to import.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- #
# Pure logic (data-free, fully unit-tested) — nothing below imports repo data  #
# --------------------------------------------------------------------------- #


MIN_BASELINE_N = 2      # gate G-a: below this a zero-point has no error bar


@dataclass
class ChunkControlStat:
    """Control census for one built chunk (one source × cell_line)."""
    chunk_id: str
    cell_line: str
    n_cells: int
    n_control: int
    error: str | None = None
    # Stage 1.5.2 gate G-a. Purely additive and defaulted, so every existing construction
    # site and test keeps working and the PASS/FAIL rule below is untouched.
    control_batches: tuple[str, ...] = ()
    cell_batches: tuple[str, ...] = ()

    @property
    def n_perturbed(self) -> int:
        return self.n_cells - self.n_control

    @property
    def fallback_fired(self) -> bool:
        # self-centring only bites a line that HAS perturbed cells but no control to
        # anchor them; an all-control chunk has nothing whose zero-point could flip.
        return self.error is None and self.n_perturbed > 0 and self.n_control == 0

    @property
    def unreplicated_baseline(self) -> bool:
        """A zero-point built from a single sample — no error bar, and until G-a, silent.

        Stage 1.5 made `n = 0` visible. This is the `n = 1` case: the per-donor offset Stage 2
        exists to correct is ±12.7 yr, and ONE clock measurement's own CV error is 12.27 yr.
        """
        return (self.error is None and self.n_perturbed > 0
                and 0 < self.n_control < MIN_BASELINE_N)

    @property
    def cross_batch_baseline(self) -> bool:
        """Every control from one batch while the line's cells span several — finding D1.

        `ΔAge = age(Exp1 sample) − age(Exp2 baseline)` for ~half of every Gill donor's samples.
        That sits INSIDE the definition of `y_age`, not downstream of it.
        """
        return (self.error is None and self.n_control > 0
                and len(set(self.control_batches)) == 1
                and len(set(self.cell_batches)) > 1)


def decide_verdict(stats: list[ChunkControlStat]) -> dict:
    """PASS iff no chunk fired the no-control fallback. Pure; no I/O.

    - CANNOT_VERIFY — nothing was scanned (no chunks, or all errored on fetch).
    - FAIL          — >=1 chunk has perturbed cells and zero controls (fallback fired).
    - PASS          — every chunk with perturbed cells also carries >=1 control.
    """
    scanned = [s for s in stats if s.error is None]
    errored = [s for s in stats if s.error is not None]
    if not scanned:
        return {
            "status": "CANNOT_VERIFY",
            "reason": "no chunks could be scanned (none planned, or all failed to fetch)",
            "n_chunks": len(stats), "n_scanned": 0,
            "n_errored": len(errored), "fallback_chunks": [],
        }
    fired = [s for s in scanned if s.fallback_fired]
    status = "FAIL" if fired else "PASS"
    if fired:
        reason = (f"{len(fired)} chunk(s) had perturbed cells but NO vehicle controls — the "
                  f"aging.py:88 fallback fired, so their ΔAge zero-point is self-centred, not "
                  f"control-relative. Part of the per-donor offset may be an artefact.")
    else:
        reason = (f"every one of {len(scanned)} scanned chunks with perturbed cells also carries "
                  f">=1 control; the zero-point is control-relative throughout.")
    # Gate G-a additions. These are WARNINGS, deliberately kept OUT of `status`: the Stage 1.5
    # gate means one specific thing ("the no-control fallback did not fire") and four runs are
    # recorded against it. Folding a new condition into the same PASS would silently redefine
    # what those four PASSes meant. They are reported beside it instead.
    unrep = [s for s in scanned if s.unreplicated_baseline]
    xbatch = [s for s in scanned if s.cross_batch_baseline]
    warnings = []
    if unrep:
        warnings.append(
            f"{len(unrep)} chunk(s) rest on an UNREPLICATED baseline (n < {MIN_BASELINE_N}): "
            + ", ".join(f"{s.cell_line} (n={s.n_control})" for s in unrep)
            + " — a zero-point with no error bar")
    if xbatch:
        warnings.append(
            f"{len(xbatch)} chunk(s) have a CROSS-BATCH baseline: "
            + ", ".join(f"{s.cell_line} (controls {sorted(set(s.control_batches))}, cells "
                        f"{sorted(set(s.cell_batches))})" for s in xbatch)
            + " — ΔAge is a cross-batch difference for the non-matching cells (finding D1)")
    return {
        "status": status,
        "reason": reason,
        "n_chunks": len(stats), "n_scanned": len(scanned), "n_errored": len(errored),
        "fallback_chunks": [
            {"chunk_id": s.chunk_id, "cell_line": s.cell_line,
             "n_cells": s.n_cells, "n_control": s.n_control} for s in fired
        ],
        "baseline_warnings": warnings,
        "unreplicated_chunks": [
            {"chunk_id": s.chunk_id, "cell_line": s.cell_line, "n_control": s.n_control}
            for s in unrep],
        "cross_batch_chunks": [
            {"chunk_id": s.chunk_id, "cell_line": s.cell_line,
             "control_batches": sorted(set(s.control_batches)),
             "cell_batches": sorted(set(s.cell_batches))} for s in xbatch],
    }


# --------------------------------------------------------------------------- #
# Real-data replay (imports repo data machinery only when actually run)        #
# --------------------------------------------------------------------------- #
def scan_build(work) -> list[ChunkControlStat]:
    """Fetch every planned chunk and count its vehicle controls. Read-only: fetches
    raw cells, writes nothing. `work` is the output of `plan_all(sources)`."""
    import numpy as np

    stats: list[ChunkControlStat] = []
    for src, chunk in work:
        cid = chunk["id"]
        line = chunk.get("cell_line", "?")
        try:
            raw = src.fetch(chunk)
            is_ctrl = np.asarray(raw.obs["is_control"].to_numpy(), dtype=bool)
            # G-b stamps `batch` into obs; sources that do not are simply reported as single-batch,
            # which is honest — an unrecorded batch is not the same as a verified one.
            b = ([str(x) for x in raw.obs["batch"].to_numpy()] if "batch" in raw.obs.columns
                 else [""] * int(is_ctrl.size))
            stats.append(ChunkControlStat(
                cid, line, int(is_ctrl.size), int(is_ctrl.sum()),
                control_batches=tuple(x for x, c in zip(b, is_ctrl, strict=True) if c),
                cell_batches=tuple(b)))
        except Exception as exc:  # noqa: BLE001 — recorded per chunk, never aborts the scan
            stats.append(ChunkControlStat(cid, line, 0, 0, error=repr(exc)[:160]))
    return stats


def build_sources(gse_dir: str, gill_dir: str):
    """Reconstruct the exact GSE242423 + Gill sources the build used.

    Mirrors run_multi_local.py:100-114. Imported lazily so this module stays
    data-free at import time (the unit tests import only the pure logic above).
    """
    runners = Path(__file__).resolve().parents[1] / "local_runners"
    if str(runners) not in sys.path:
        sys.path.insert(0, str(runners))
    try:
        # discover_* + the two sampling constants live in the runner, next to the
        # source construction they feed — reuse them rather than duplicating.
        from run_multi_local import (  # type: ignore
            CELLS_PER_RUN,
            MAX_CELLS,
            discover_gill,
            discover_gse,
        )
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            "could not import the source-discovery helpers from local_runners/run_multi_local.py "
            f"({exc!r}); the real-data replay depends on the runner's setup — see the module "
            "docstring."
        ) from exc

    from cellfate.data.sources import GillReprogrammingSource, GSE242423SingleCellSource

    gse_samples, gse_genes = discover_gse(gse_dir)
    gill_expr, gill_series = discover_gill(gill_dir)
    gse = GSE242423SingleCellSource(gse_samples, gse_genes, cell_line="HFF", min_genes=500,
                                    max_cells_per_sample=MAX_CELLS, cells_per_run=CELLS_PER_RUN,
                                    seed=0)
    gill = GillReprogrammingSource(gill_expr, gill_series)
    return [gse, gill]


def _render(stats: list[ChunkControlStat], verdict: dict) -> None:
    from cellfate.common.console import install_pretty_console, render_table
    install_pretty_console()

    rows = []
    for s in sorted(stats, key=lambda x: x.chunk_id):
        if s.error is not None:
            # Six cells, matching the header. `render_table` indexes `row[i] for i in
            # range(len(headers))`, so a short row raises IndexError -- which would crash the
            # renderer on exactly the errored-chunk path `scan_build` goes out of its way to
            # survive ("recorded per chunk, never aborts the scan").
            rows.append([s.chunk_id, s.cell_line, "ERR", "ERR", s.error[:36], "-"])
            continue
        flag = "** FALLBACK FIRED **" if s.fallback_fired else "ok"
        notes = []
        if s.unreplicated_baseline:
            notes.append(f"n={s.n_control} UNREPLICATED")
        if s.cross_batch_baseline:
            notes.append(f"CROSS-BATCH (ctrl {sorted(set(s.control_batches))})")
        rows.append([s.chunk_id, s.cell_line, str(s.n_cells),
                     f"{s.n_control} ctrl / {s.n_perturbed} pert", flag, "; ".join(notes) or "-"])
    print("\n  PER-CHUNK CONTROL CENSUS")
    print(render_table(["chunk", "cell_line", "cells", "composition", "zero-point", "baseline (G-a)"],
                       rows, aligns=["l", "l", "r", "l", "l", "l"]))

    print(f"\n   VERDICT: {verdict['status']}")
    print(f"     => {verdict['reason']}")
    if verdict.get("baseline_warnings"):
        print("\n   BASELINE WARNINGS (gate G-a — reported BESIDE the verdict, never folded into")
        print("   it: the Stage 1.5 PASS means one specific thing and four runs are recorded")
        print("   against it, so redefining it retroactively would invalidate that record):")
        for w in verdict["baseline_warnings"]:
            print(f"       ! {w}")
    if verdict["fallback_chunks"]:
        print("     chunks that fired the fallback (their ΔAge is self-centred, mean forced to 0):")
        for c in verdict["fallback_chunks"]:
            print(f"       - {c['chunk_id']}  ({c['cell_line']}): "
                  f"{c['n_cells']} cells, {c['n_control']} controls")
        print("\n   This is a FINDING, not a bug to patch here (STAGE_1_5 §1/§3): record it in")
        print("   STAGE_1_DEVIATIONS.md and fix the zero-point under its OWN pre-registered Change.")


def main() -> None:
    # cp1255 (Hebrew) consoles cannot encode the box-drawing table; emit UTF-8 so a
    # print can never abort the gate. The JSON report is written regardless.
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass

    gse_dir = sys.argv[1] if len(sys.argv) > 1 else r"D:\GSE242423"
    gill_dir = sys.argv[2] if len(sys.argv) > 2 else r"D:\Gill"

    print("\nSTAGE 1.5 VERIFICATION — did the no-control ΔAge fallback fire on the real build?")
    print(f"  GSE242423 dir : {gse_dir}")
    print(f"  Gill dir      : {gill_dir}")

    if not (Path(gse_dir).exists() and Path(gill_dir).exists()):
        print("\n   !! source data not found at those paths. Pass them as arguments:")
        print('      python plan_tests/verify_stage1_5.py "D:\\GSE242423" "D:\\Gill"')
        print("   (Group A–D of tests/test_harmonize.py cover everything that does NOT need the data.)")
        return

    from cellfate.data.chunking import plan_all

    sources = build_sources(gse_dir, gill_dir)
    work = plan_all(sources)
    print(f"  planned {len(work)} chunks across {len(sources)} sources; fetching to census controls...")
    stats = scan_build(work)
    verdict = decide_verdict(stats)

    report = {
        "script": "verify_stage1_5",
        "utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "gse_dir": gse_dir, "gill_dir": gill_dir,
        "verdict": verdict,
        "chunks": [
            {"chunk_id": s.chunk_id, "cell_line": s.cell_line, "n_cells": s.n_cells,
             "n_control": s.n_control, "n_perturbed": s.n_perturbed,
             "fallback_fired": s.fallback_fired, "error": s.error}
            for s in stats
        ],
    }
    out = _RESULTS / "verify_stage1_5_results.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\n  saved -> {out}   |   VERDICT: {verdict['status']}")

    try:
        _render(stats, verdict)
    except Exception as exc:  # noqa: BLE001 — the JSON already holds the verdict
        print(f"  (table render skipped: {exc!r}; the verdict is in {out})")


if __name__ == "__main__":
    main()

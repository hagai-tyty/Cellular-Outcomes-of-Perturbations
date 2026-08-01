"""Transcriptomic aging clock + ΔAge (Document 2, S10).

The clock predicts a transcriptomic age from expression; ΔAge is the predicted
age relative to the matched vehicle-control baseline of the same cell line
(negative = rejuvenated). On cancer / transformed lines the clock is
out-of-distribution, so ΔAge is masked (``age_mask = False``) -- the safety head
still trains on those cells, the age head does not.

``LinearClock`` (age = w.x + b over panel genes) is the dependency-free default
and the interface real clocks (Buckley et al.; scAgeClock) plug into.
"""

from __future__ import annotations

import json
from abc import ABC, abstractmethod
from pathlib import Path

import numpy as np
import pandas as pd

from cellfate.common import constants as C
from cellfate.common.panel import GenePanel


class AgingClock(ABC):
    """Predicts a transcriptomic age (years) from expression.

    The clock consumes the **full normalised profile** and its own gene panel --
    it aligns its weights against the gene symbols it is handed, NOT the model's
    2000-HVG input. So it is decoupled from the model input: the model sees the
    HVG panel; the clock sees every gene it can match.
    """

    @abstractmethod
    def predict_age(self, expr: np.ndarray, genes: list[str]) -> np.ndarray:
        """Return (N,) predicted ages for an (N, len(genes)) matrix in ``genes`` order."""


class LinearClock(AgingClock):
    """age = sum_g w_g * x_g + b, with weights keyed by gene symbol."""

    def __init__(self, weights: dict[str, float], intercept: float = 0.0) -> None:
        self.weights = weights
        self.intercept = float(intercept)

    def predict_age(self, expr: np.ndarray, genes: list[str]) -> np.ndarray:
        w = np.array([self.weights.get(g, 0.0) for g in genes], dtype=np.float64)
        return np.asarray(expr, dtype=np.float64) @ w + self.intercept

    @classmethod
    def random(cls, panel: GenePanel, seed: int = 0, scale: float = 1.0) -> LinearClock:
        """A deterministic **random** clock over a panel.

        For synthetic/smoke runs and tests ONLY -- its ages are meaningless. Real
        runs must load fitted weights via :meth:`from_json` (see scripts/fit_clock.py).
        """
        rng = np.random.default_rng(seed)
        w = rng.normal(0.0, scale, size=len(panel)) / np.sqrt(len(panel))
        return cls({g: float(wi) for g, wi in zip(panel.genes, w, strict=True)}, intercept=40.0)

    @classmethod
    def from_json(cls, path: str | Path) -> LinearClock:
        """Load a fitted clock: ``{"weights": {gene: w}, "intercept": b, ...}``."""
        d = json.loads(Path(path).read_text(encoding="utf-8"))
        if "weights" not in d:
            raise ValueError(f"clock file {path} has no 'weights' key")
        weights = {str(k): float(v) for k, v in d["weights"].items()}
        if not weights:
            raise ValueError(f"clock file {path} has empty weights")
        return cls(weights, intercept=float(d.get("intercept", 0.0)))

    def to_json(self, path: str | Path, meta: dict | None = None) -> None:
        """Serialise fitted weights (+ optional provenance ``meta``) to JSON."""
        payload = {"weights": self.weights, "intercept": self.intercept}
        if meta:
            payload["meta"] = meta
        Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _control_baseline(values: np.ndarray, lines: np.ndarray, is_ctrl: np.ndarray,
                      census: dict | None = None,
                      composition: dict[str, np.ndarray] | None = None) -> np.ndarray:
    """Per-line mean over vehicle controls. Falls back to the line's own mean when
    a line has no controls in this chunk (values then centred within the line).

    When ``census`` is supplied it is filled in with, per line, how many controls the
    baseline actually rests on and what they were made of. **This is recording only --
    the arithmetic above is untouched and ΔAge is bit-identical with and without it**
    (``tests/test_baseline_census.py`` asserts that rather than assuming it).

    Why it exists (Stage 1.5.2 §0, gate G-a). Stage 1.5 made the ``n = 0`` case visible;
    ``n = 1`` was still silent. A zero-point built from one unreplicated sample has no error
    bar, and nothing in the pipeline said so -- while the per-donor offset Stage 2 exists to
    correct is ±12.7 yr, the same magnitude as ONE clock measurement's error. You cannot
    audit what is not recorded.

    ``composition`` maps a label (e.g. ``"batch"``) to a per-cell array; the distinct values
    among each line's controls are recorded. That is what makes finding D1 -- every baseline
    drawn from ``Exp2`` while ~50% of treated samples are ``Exp1`` -- visible in the output
    instead of reconstructible only by hand.
    """
    baseline = np.empty_like(values, dtype=np.float64)
    for line in np.unique(lines):
        in_line = lines == line
        ctrl = in_line & is_ctrl
        ref = values[ctrl] if ctrl.any() else values[in_line]
        baseline[in_line] = ref.mean()
        if census is not None:
            used = ctrl if ctrl.any() else in_line
            rec = {"n_control": int(ctrl.sum()), "n_cells": int(in_line.sum()),
                   "source": "controls" if ctrl.any() else "self_fallback",
                   "unreplicated": bool(ctrl.sum() == 1)}
            for name, arr in (composition or {}).items():
                a = np.asarray(arr)
                # Both are needed to say anything useful. The baseline's values alone cannot
                # distinguish "the controls are all Exp2 while the line spans Exp1 and Exp2"
                # (finding D1, a real defect) from "this column is constant within a line"
                # (donor_age, where a single value is the only possible answer).
                rec[name] = sorted({str(v) for v in a[used]})
                rec[f"{name}_in_line"] = sorted({str(v) for v in a[in_line]})
            census[str(line)] = rec
    return baseline


def age_label_policy(
    n: int,
    source: str,
    obs: pd.DataFrame,
    *,
    masked_datasets: frozenset[str] = frozenset(),
    clock_age_range: tuple[float, float] | None = None,
) -> tuple[np.ndarray, list[str | None]]:
    """Which cells get a usable ΔAge label, and -- when they do not -- WHY. Pure.

    Stage 1.5.3 C-1. Before this, the rule was one expression keyed on ``source`` alone, and
    both reprogramming sources report ``source = "reprogramming"`` -- so **no policy could
    mask HFF and keep Gill**, which is precisely what Stage 1.5.2's gate G-c step 2 requires.

    Three independent reasons, checked in decreasing order of certainty. The first is the
    pre-existing rule and is never weakened by the other two; a cell excluded for more than one
    reason reports the FIRST that applied, so the string is stable under reordering of the
    later rules.

      1. ``cancer_source``   -- the clock is out of distribution on transformed lines.
                                Today's only rule (``CANCER_SOURCES``), unchanged.
      2. ``dataset_policy``  -- gate G-c: this dataset's labels are withheld by decision, not
                                by cell type. Empty by default.
      3. ``donor_out_of_clock_range`` -- the donor's chronological age is outside the range the
                                clock was FITTED on (``configs/clocks/*.json`` ->
                                ``meta.age_range``). An UNKNOWN age never masks.

    Returns ``(age_mask, reasons)`` with ``reasons[i] is None`` exactly where ``age_mask[i]``
    is True -- the invariant :class:`~cellfate.common.schemas.Sample` validation enforces.
    """
    mask = np.full(n, True, dtype=bool)
    reasons: list[str | None] = [None] * n

    def _exclude(bad: np.ndarray, why: str) -> None:
        newly = np.asarray(bad, dtype=bool) & mask
        mask[newly] = False
        for i in np.flatnonzero(newly):
            reasons[i] = why

    if source in C.CANCER_SOURCES:
        _exclude(np.ones(n, dtype=bool), "cancer_source")

    # FAIL LOUD, NEVER OPEN. If a withholding policy is switched on and the column it needs is
    # absent, the silent outcome is to KEEP labels that were meant to be withheld -- the unsafe
    # direction, and invisible. Not hypothetical: G-b reached Gill and not HFF, so `donor_age`
    # was missing on HFF until C-3. The step order protects us; correctness must not depend on it.
    if masked_datasets:
        if "dataset_id" not in obs.columns:
            raise KeyError(
                "age_label_policy: masked_datasets was requested but obs has no 'dataset_id' "
                "column, so the policy cannot be applied. Refusing to silently keep labels "
                "that were meant to be withheld.")
        _exclude(obs["dataset_id"].isin(masked_datasets).to_numpy(), "dataset_policy")

    if clock_age_range is not None:
        if "donor_age" not in obs.columns:
            raise KeyError(
                "age_label_policy: clock_age_range was requested but obs has no 'donor_age' "
                "column (see C-3). Refusing to silently treat every cell as in-range.")
        a = pd.to_numeric(obs["donor_age"], errors="coerce").to_numpy(dtype=float)
        lo, hi = clock_age_range
        # NaN comparisons are False, so an UNKNOWN donor age cannot mask. Deliberate, and
        # distinct from the raise above: a missing COLUMN means the policy is inapplicable and
        # is an error; a missing VALUE in a present column is recorded absence, and absence of
        # evidence is not acted on.
        _exclude((a < lo) | (a > hi), "donor_out_of_clock_range")

    return mask, reasons


def census_warnings(census: dict, min_controls: int = 2) -> list[str]:
    """Human-readable problems in a baseline census. Pure, so it is testable on its own.

    Three things are worth saying out loud, in descending severity:
      * no controls at all -- the ``aging.py`` self-centring fallback fired (Stage 1.5's gate);
      * a single unreplicated control -- a zero-point with no error bar (G-a's reason to exist);
      * a baseline drawn from strictly fewer batches than the line spans -- finding D1's
        cross-batch structure, which sits *inside* the definition of ``y_age``.

    The third check compares the baseline's values against the **whole line's**, so it fires
    only on a genuine mismatch. A column that is constant within a line by construction
    (``donor_age``) can never trigger it -- an earlier version warned on every donor, which is
    noise that would have trained the reader to ignore the warnings that matter.
    """
    out = []
    for line, rec in sorted(census.items()):
        if rec["source"] == "self_fallback":
            out.append(f"{line}: NO controls in this chunk; ΔAge centred within the line "
                       f"({rec['n_cells']} cells)")
        elif rec["n_control"] < min_controls:
            out.append(f"{line}: baseline rests on n={rec['n_control']} control cell(s) -- "
                       f"an unreplicated zero-point, no error bar")
        for name, vals in rec.items():
            in_line = rec.get(f"{name}_in_line")
            if not isinstance(vals, list) or in_line is None:
                continue
            if len(vals) == 1 and len(in_line) > 1:
                out.append(f"{line}: every baseline cell has {name}={vals[0]}, but the line "
                           f"spans {in_line} -- ΔAge is a cross-{name} difference for the rest")
    return out


def recenter_on_control_arrays(
    values: np.ndarray, lines: np.ndarray, is_ctrl: np.ndarray
) -> np.ndarray:
    """Array form of :func:`recenter_on_controls` (no DataFrame needed).

    Subtracts the per-line vehicle-control baseline from ``values`` given the
    per-cell ``lines`` and boolean ``is_ctrl`` arrays.
    """
    return np.asarray(values, dtype=np.float64) - _control_baseline(values, lines, is_ctrl)


def recenter_on_controls(values: np.ndarray, obs: pd.DataFrame) -> np.ndarray:
    """Subtract the per-line vehicle-control baseline from ``values``.

    Used to re-anchor ΔAge after cell-cycle deconfounding: ``deconfound_age``
    removes the regression intercept and so re-centres the whole population to
    mean 0, which shifts the controls off zero. Re-applying the control baseline
    restores the invariant that ΔAge is *control-relative* (controls ~ 0), which
    is the zero-point the rejuvenation score depends on -- without reintroducing
    the cell-cycle slope that was just removed.
    """
    lines = obs["cell_line"].to_numpy()
    is_ctrl = obs["is_control"].to_numpy().astype(bool)
    return recenter_on_control_arrays(values, lines, is_ctrl)


def delta_age(
    clock: AgingClock,
    expr: np.ndarray,
    genes: list[str],
    obs: pd.DataFrame,
    source: str,
    census: dict | None = None,
    composition_cols: tuple[str, ...] = ("batch", "donor_age"),
    clock_age_range: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, list[str | None]]:
    """Compute (ΔAge, age_mask, age_mask_reason) from the **full normalised profile**.

    The clock is handed ``expr`` (N, len(genes)) with its gene symbols ``genes``
    -- the full profile, not the 2000-HVG model input -- so it can dot its
    weights against every gene it matches.

    ΔAge[i] = age[i] - mean(age over vehicle controls of the same cell line).
    If a cell line has no controls in this chunk, its own mean age is used as the
    baseline (ΔAge centred within the line). ``age_mask`` is all-False for
    sources in ``CANCER_SOURCES``.

    Note: if the caller subsequently deconfounds ΔAge for cell cycle, it must
    re-anchor with ``recenter_on_controls`` to preserve the control-relative
    zero-point (deconfounding otherwise re-centres the whole population).

    Pass ``census`` (a dict) to have the per-line baseline count and composition
    recorded into it -- Stage 1.5.2's gate G-a. Purely additive: the returned ΔAge
    does not depend on whether ``census`` was supplied. ``composition_cols`` names
    the ``obs`` columns to summarise; missing ones are skipped, so a source that
    does not stamp them is not an error.

    ``clock_age_range`` enables :func:`age_label_policy`'s third rule (Stage 1.5.3 C-2).
    ``None`` -- the default -- leaves it off, so nothing changes unless a caller opts in.
    """
    age = clock.predict_age(expr, genes)
    lines = obs["cell_line"].to_numpy()
    is_ctrl = obs["is_control"].to_numpy().astype(bool)
    comp = {c: obs[c].to_numpy() for c in composition_cols if c in obs.columns} or None
    d = age - _control_baseline(age, lines, is_ctrl, census=census, composition=comp)
    age_mask, age_mask_reason = age_label_policy(
        age.shape[0], source, obs,
        masked_datasets=C.AGE_MASKED_DATASETS, clock_age_range=clock_age_range)
    return d, age_mask, age_mask_reason

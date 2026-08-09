"""Bulk sample integrity — reject columns that are not transcriptomes (Change C-7).

WHY THIS EXISTS
---------------
`GSE165176_Log2_RPM_Sendai_reprogramming` contains columns that are not expression profiles.
`N2_Fib_Sendai_Exp2` -- donor N2's day-0 control, and therefore N2's entire ΔAge zero-point --
is **nearly a constant vector**: `min = median = mean = 11.490`, `max = 13.227`, a dynamic
range of 1.74 log2 units where every sound control spans 13-15, and a linear library 68x the
cohort. The clock reads it as **98.65 yr** for a donor of age **0** -- the highest of all six,
and +35.12 yr above the other five.

It does not stay in its own donor. The day-0 `_Fib_` sample is `is_control`
(`sources.py:417`), which makes it **two things at once**: that donor's zero-point, and one of
the five or six controls `sigma_gill` is fitted on. `sigma_gill / sigma_hff` is the gain applied
to **HFF's** labels -- 99.7% of the age-labelled corpus -- in every fold that does not hold N2
out. Step 3c measured the consequence: removing this one column moves HFF's day-14 ΔAge from
**-26.755 to -8.196**, `+18.558 yr`, against healthy control drops of -0.32 to -2.98 yr.

THE TWO CONDITIONS, AND WHY THEY ARE UNIT-BASED RATHER THAN QUANTILE-BASED
--------------------------------------------------------------------------
An earlier proposal thresholded `mean - min` at 1/5 of the cohort median. It was rejected
(STAGE_1_5_6 §5.10): it cuts a continuous distribution 8% from its neighbour, so on a new
cohort it would flag or miss arbitrarily. Both conditions here come from what the numbers
*mean*; this cohort only confirms that they separate.

  G1  LIBRARY SIZE. The matrix is Reads Per Million, so a sound column's linear values sum to
      ~1e6 BY DEFINITION of the units. Accept a decade either side.
  G2  DYNAMIC RANGE. Any real transcriptome spans orders of magnitude between its least- and
      most-expressed gene. Require >= 256-fold.

Verified over all 124 Gill columns before this file was written: G1 AND G2 reject **exactly 5**
with **0** false positives of 119. The other 119 span library 2.859e+05 - 3.880e+06 and range
9.00 - 15.26. G1 margins 2.58x below the ceiling and 1.69x above it; G2 gap 7.26 vs 9.00, no
overlap. **Each condition independently rejects all five**, so they are kept as two only
because they fail differently -- G1 catches a mis-scaled library, G2 a collapsed distribution.

SCOPE -- BULK ONLY, AND THAT IS NOT AN ARBITRARY RESTRICTION
-------------------------------------------------------------
G1 is defined on RPM. Single-cell sources yield **raw UMI counts** (`normalize_counts` runs
later, at `build_dataset.py:174`), which sum to ~1e3-1e4 per cell, so G1 would reject every
cell **by construction of the units** rather than by any property of the data. This module is
therefore called from bulk sources only, and single-cell sources never reach it.

PURE. No I/O, no config, no logging. Testable on synthetic arrays.
"""
from __future__ import annotations

import numpy as np

# A sound RPM column sums to ~1e6 by definition of the units. A decade either side is the
# tolerance; the observed cohort sits within 3.9x of 1e6 and the rejects are 17x to 2148x it.
G1_LIBRARY_BAND: tuple[float, float] = (1e5, 1e7)

# >= 256-fold between the least- and most-expressed gene. The observed cohort's sound columns
# span 9.00-15.26 log2 units; the rejects span 0.15-7.26.
G2_MIN_LOG2_RANGE: float = 8.0

REASON_LIBRARY = "library_out_of_band"
REASON_RANGE = "dynamic_range_collapsed"


def linear_library_size(log2_values: np.ndarray) -> float:
    """Sum of the column in LINEAR space, inverting log2 RPM exactly as the pipeline does.

    Mirrors `sources.py`'s own inversion (``2**x - 1``, negatives clipped) so the number this
    gate judges is the number the pipeline goes on to use -- not a differently-derived proxy.
    """
    v = np.asarray(log2_values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 0.0
    lin = np.power(2.0, v) - 1.0
    lin[lin < 0] = 0.0
    return float(lin.sum())


def log2_dynamic_range(log2_values: np.ndarray) -> float:
    """`max - min` over the column's log2 values. Zero for a constant column."""
    v = np.asarray(log2_values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 0.0
    return float(v.max() - v.min())


def bulk_column_verdict(log2_values: np.ndarray) -> tuple[bool, str | None]:
    """Admit or reject ONE bulk sample from its raw log2-RPM column. Pure.

    Returns ``(admitted, reason)`` with ``reason is None`` **exactly when** ``admitted`` is
    True -- the same invariant `age_label_policy` holds for ``(age_mask, reasons)``, so the two
    read the same way at their call sites.

    A column must satisfy BOTH conditions to be admitted. When both fail, the library reason is
    reported first: it is the more mechanical of the two, and a stable choice keeps the string
    reproducible rather than dependent on evaluation order.
    """
    lib = linear_library_size(log2_values)
    lo, hi = G1_LIBRARY_BAND
    if not (lo <= lib <= hi):
        return False, REASON_LIBRARY
    if log2_dynamic_range(log2_values) < G2_MIN_LOG2_RANGE:
        return False, REASON_RANGE
    return True, None


def screen_bulk_matrix(log2_matrix: np.ndarray, sample_names: list[str]) -> dict[str, str]:
    """Verdicts for every column of a (genes x samples) log2-RPM matrix. Pure.

    Returns ``{sample_name: reason}`` for REJECTED samples only, so an empty dict means the
    whole matrix is admissible and the caller needs no special case for the healthy path.
    """
    m = np.asarray(log2_matrix, dtype=np.float64)
    if m.ndim != 2:
        raise ValueError(f"expected a 2-D (genes x samples) matrix, got shape {m.shape}")
    if m.shape[1] != len(sample_names):
        raise ValueError(
            f"{m.shape[1]} columns but {len(sample_names)} sample names -- refusing to guess "
            "the pairing, because a mis-paired verdict would reject the wrong sample")
    out: dict[str, str] = {}
    for j, name in enumerate(sample_names):
        ok, why = bulk_column_verdict(m[:, j])
        if not ok:
            out[str(name)] = str(why)
    return out

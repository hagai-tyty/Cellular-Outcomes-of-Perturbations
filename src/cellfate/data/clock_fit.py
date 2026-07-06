"""Fit a real linear transcriptomic aging clock from an age-labelled expression
matrix (e.g. GSE113957 human dermal fibroblasts, ages 1-96).

The clock consumes the **full normalised profile** at label time (the pipeline
hands ``predict_age`` every gene, not the 2000-HVG model input), so fit it on ALL
available genes -- do NOT truncate to the model panel. Weights are keyed by gene
SYMBOL; at label time the clock dot-products them against whatever symbols it
matches in the data's full profile. ``panel_genes`` is an optional restriction and
is normally left ``None``.

The output is a weights JSON loadable by ``LinearClock.from_json``, so ``clock:``
in a data config can point straight at it. See scripts/fit_clock.py for the CLI.
"""

from __future__ import annotations

import numpy as np

from .aging import LinearClock
from .normalize import normalize_counts


def fit_linear_clock(
    counts: np.ndarray,
    genes: list[str],
    ages: np.ndarray,
    *,
    panel_genes: list[str] | None = None,
    alphas: np.ndarray | None = None,
    normalized: bool = False,
    seed: int = 0,
) -> tuple[LinearClock, dict]:
    """Fit ``age ~ expression`` (ridge, cross-validated) and return (clock, metrics).

    Parameters
    ----------
    counts : (n_samples, n_genes) raw counts (or already-normalized if ``normalized``).
    genes : gene SYMBOLS, length n_genes, matching the columns of ``counts``.
    ages : (n_samples,) chronological ages in years.
    panel_genes : if given, restrict/align the clock to these symbols (e.g. the
        978-gene L1000 panel the model uses); genes absent from the data are dropped.
    normalized : set True if ``counts`` is already log1p CP10k (skip normalization).

    The returned metrics carry an honest k-fold CV MAE / Pearson so you can report
    the clock's accuracy in the paper.
    """
    from sklearn.linear_model import RidgeCV
    from sklearn.metrics import mean_absolute_error
    from sklearn.model_selection import KFold, cross_val_predict

    X = np.asarray(counts, dtype=np.float64)
    ages = np.asarray(ages, dtype=np.float64)
    genes = [str(g) for g in genes]
    if X.shape[1] != len(genes):
        raise ValueError(f"counts has {X.shape[1]} gene columns but {len(genes)} names")
    if X.shape[0] != ages.shape[0]:
        raise ValueError(f"counts has {X.shape[0]} samples but {ages.shape[0]} ages")
    if X.shape[0] < 5:
        raise ValueError("need at least 5 samples to fit a clock")

    Xn = X if normalized else normalize_counts(X)   # log1p CP10k == predict_age input

    if panel_genes is not None:
        panel_genes = [str(g) for g in panel_genes]
        pos = {g: i for i, g in enumerate(genes)}
        used = [g for g in panel_genes if g in pos]
        if not used:
            raise ValueError("no overlap between data genes and the panel")
        Xn = Xn[:, [pos[g] for g in used]]
    else:
        used = genes

    alphas = alphas if alphas is not None else np.logspace(-1.0, 4.0, 24)
    model = RidgeCV(alphas=alphas).fit(Xn, ages)

    n_splits = int(min(5, len(ages)))
    cv = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    pred = cross_val_predict(RidgeCV(alphas=alphas), Xn, ages, cv=cv)
    mae = float(mean_absolute_error(ages, pred))
    pearson = float(np.corrcoef(pred, ages)[0, 1]) if len(ages) > 2 else float("nan")

    weights = {g: float(w) for g, w in zip(used, model.coef_, strict=True)}
    clock = LinearClock(weights, intercept=float(model.intercept_))
    metrics = {
        "n_samples": int(len(ages)),
        "n_genes": int(len(used)),
        "cv_mae_years": mae,
        "cv_pearson": pearson,
        "alpha": float(model.alpha_),
        "age_range": [float(ages.min()), float(ages.max())],
        "normalization": "log1p_cp10k",
    }
    return clock, metrics

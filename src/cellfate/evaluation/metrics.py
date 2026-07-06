"""Calibration-, imbalance-, and ranking-aware metrics (Document 5, S3).

Classification uses per-class AUROC/PR-AUC (the rare dangerous classes need
precision-recall, not accuracy), Brier, and ECE with reliability bins. Regression
reports MAE/RMSE/Pearson/R2 on age-valid cells only. Coverage checks the conformal
interval against its nominal level. Ranking quality (Spearman + precision@k) scores
RES against the measured rejuvenation effect -- the triage-utility metric. DES/PDS
are provided for a future full-state predictor (this model predicts class + ΔAge,
not the whole transcriptome, so they are not part of the main loop).
"""

from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

NAN = float("nan")


def per_class_auroc(y_true, p) -> dict[int, float]:
    y_true = np.asarray(y_true)
    out = {}
    for c in range(p.shape[1]):
        yc = (y_true == c).astype(int)
        out[c] = float(roc_auc_score(yc, p[:, c])) if 0 < yc.sum() < len(yc) else NAN
    return out


def per_class_prauc(y_true, p) -> dict[int, float]:
    y_true = np.asarray(y_true)
    out = {}
    for c in range(p.shape[1]):
        yc = (y_true == c).astype(int)
        out[c] = float(average_precision_score(yc, p[:, c])) if yc.sum() > 0 else NAN
    return out


def brier(y_onehot, p) -> float:
    return float(((np.asarray(p) - np.asarray(y_onehot)) ** 2).sum(1).mean())


def ece(y_true, p, n_bins: int = 15) -> float:
    y_true = np.asarray(y_true)
    conf, pred = p.max(1), p.argmax(1)
    acc = (pred == y_true).astype(float)
    bins = np.linspace(0.0, 1.0, n_bins + 1)
    e = 0.0
    for lo, hi in zip(bins[:-1], bins[1:], strict=True):
        m = (conf > lo) & (conf <= hi)
        if m.any():
            e += m.mean() * abs(acc[m].mean() - conf[m].mean())
    return float(e)


def regression_metrics(y, yhat, mask) -> dict[str, float]:
    mask = np.asarray(mask, dtype=bool)
    y = np.asarray(y, dtype=np.float64)[mask]
    yhat = np.asarray(yhat, dtype=np.float64)[mask]
    if len(y) == 0:
        return {"mae": NAN, "rmse": NAN, "pearson": NAN, "r2": NAN, "n": 0}
    mae = float(np.abs(y - yhat).mean())
    rmse = float(np.sqrt(((y - yhat) ** 2).mean()))
    pear = float(np.corrcoef(y, yhat)[0, 1]) if y.std() > 0 and yhat.std() > 0 else NAN
    denom = ((y - y.mean()) ** 2).sum()
    r2 = float(1.0 - ((y - yhat) ** 2).sum() / denom) if denom > 0 else NAN
    return {"mae": mae, "rmse": rmse, "pearson": pear, "r2": r2, "n": int(len(y))}


def coverage(y, lo, hi, mask) -> float:
    mask = np.asarray(mask, dtype=bool)
    y = np.asarray(y)[mask]
    lo = np.asarray(lo)[mask]
    hi = np.asarray(hi)[mask]
    if len(y) == 0:
        return NAN
    return float(((y >= lo) & (y <= hi)).mean())


def precision_at_k(res_scores, measured_quality, k: int) -> float:
    res_scores = np.asarray(res_scores)
    measured_quality = np.asarray(measured_quality)
    k = min(k, len(res_scores))
    if k == 0:
        return NAN
    top = np.argsort(-res_scores)[:k]
    thr = np.median(measured_quality)
    return float((measured_quality[top] >= thr).mean())


def ranking_metrics(res_scores, measured_age_shift, k: int = 10) -> dict[str, float]:
    res = np.asarray(res_scores, dtype=np.float64)
    quality = -np.asarray(measured_age_shift, dtype=np.float64)  # more rejuvenation -> higher quality
    if len(res) < 3 or np.std(res) == 0 or np.std(quality) == 0:
        rho = NAN
    else:
        rho = float(spearmanr(res, quality).correlation)
    return {"spearman": rho, "precision_at_k": precision_at_k(res, quality, k)}


def des_pds(pred_state, true_state, control, top_frac: float = 0.1) -> dict[str, float]:
    """Virtual-Cell-Challenge-style expression fidelity for a full-state predictor:
    DES = Jaccard overlap of the top differential-expression gene sets (vs control)
    between prediction and truth; PDS = L1 nearest-match discrimination rate."""
    ctrl = np.asarray(control, dtype=np.float64).mean(0)

    def de_set(state):
        g = np.abs(np.asarray(state, dtype=np.float64) - ctrl).mean(0)
        k = max(1, int(len(g) * top_frac))
        return set(np.argsort(-g)[:k].tolist())

    dp, dt = de_set(pred_state), de_set(true_state)
    des = len(dp & dt) / len(dp | dt) if (dp | dt) else NAN

    P = np.asarray(pred_state, dtype=np.float64)
    T = np.asarray(true_state, dtype=np.float64)
    hits = sum(int(np.argmin(np.abs(P - T[i]).sum(1)) == i) for i in range(len(T)))
    pds = hits / len(T) if len(T) else NAN
    return {"des": float(des), "pds": float(pds)}


def energy_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Energy (E-) distance between two point clouds a (n,d) and b (m,d).

    E = 2*mean||a_i - b_j|| - mean||a_i - a_i'|| - mean||b_j - b_j'||  (Euclidean).
    It is >= 0, and 0 iff the two empirical distributions coincide; larger values
    mean a bigger distributional shift. On scRNA-seq this is the standard measure of
    how far a perturbed cell population has moved from control -- used here as the
    *measured effect size* the RES ranking is scored against (Goal 4.5).
    """
    from scipy.spatial.distance import cdist

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != b.shape[1]:
        raise ValueError("a and b must be 2-D with matching feature dimension")
    d_ab = float(cdist(a, b).mean())
    d_aa = float(cdist(a, a).mean())
    d_bb = float(cdist(b, b).mean())
    return max(0.0, 2.0 * d_ab - d_aa - d_bb)   # clamp tiny negatives from sampling


def edistance_to_control(
    X: np.ndarray,
    groups,
    control,
    *,
    n_pcs: int = 50,
    max_cells: int = 1000,
    seed: int = 0,
) -> dict[str, float]:
    """E-distance from the control population to each perturbation group, after PCA.

    Returns ``{group: e_distance}`` (control excluded). Cells are PCA-reduced (as is
    standard before E-distance) and each group is subsampled to ``max_cells`` so the
    pairwise-distance computation stays tractable. This is the ground-truth effect
    magnitude to correlate the model's RES ranking against.
    """
    from sklearn.decomposition import PCA

    X = np.asarray(X, dtype=np.float64)
    groups = np.asarray(groups)
    rng = np.random.default_rng(seed)

    def _sub(mask: np.ndarray) -> np.ndarray:
        idx = np.where(mask)[0]
        if len(idx) > max_cells:
            idx = rng.choice(idx, size=max_cells, replace=False)
        return idx

    n_comp = int(min(n_pcs, X.shape[1], max(1, X.shape[0] - 1)))
    Z = PCA(n_components=n_comp, random_state=seed).fit_transform(X) if X.shape[1] > n_comp else X
    ctrl = _sub(groups == control)
    if len(ctrl) == 0:
        raise ValueError(f"no control cells labelled {control!r}")
    out: dict[str, float] = {}
    for g in np.unique(groups):
        if g == control:
            continue
        out[str(g)] = energy_distance(Z[ctrl], Z[_sub(groups == g)])
    return out


def mean_finite(values) -> float:
    vals = [float(v) for v in values if v is not None and np.isfinite(v)]
    return float(np.mean(vals)) if vals else NAN

"""STAGE 1.5.6 step 1c — MEASURE the harmonization gain on HFF. Predicted 2.26.

    python experiments/diag_harmonization_gain.py "D:\\Gill" "D:\\GSE242423\\GSE242423"

READ-ONLY. Writes `results/diag_harmonization_gain_results.json`. `src/` untouched, no labels move.

THE PREDICTION, MADE BEFORE THIS RAN
------------------------------------
Step 1b eliminated the deconfounder: cell-cycle deconfounding and control re-centring contribute
**+2.11 yr** at HFF day-14 and move ΔAge TOWARD zero. So they cannot explain the gap between the
clock applied directly (**-10.62**) and the pipeline's recorded `y_age` (**-24.02**).

The remaining candidate is **harmonization**, which the real build switches on
(`local_runners/run_multi_local.py:161`, ref = `gill_bulk`). `STAGE_1_5_HARMONIZATION_AUDIT.md` §2
Group B derived its effect in closed form:

    dAge = sum_g (x_pert,g - x_ctrl,g) . sigma_ref,g / (sigma_d,g + EPS) . w_g

`sigma_d` does not cancel; it survives as a per-dataset multiplicative GAIN, and HFF carries
`sigma_gill / sigma_hff`.

    PREDICTED   gain = 24.02 / 10.62 = 2.26

**If the measured gain is far from 2.26, the attribution is wrong and step 1b's conclusion stands
only as "not the deconfounder", with the source still unknown.**

HOW THE GAIN IS COMPUTED — THE PIPELINE'S OWN ARITHMETIC, NOT AN APPROXIMATION
------------------------------------------------------------------------------
`Harmonizer.transform` is `(x - mu_d) / (sigma_d + EPS)` and `project_to_clock` is
`x_scaled * sigma_ref + mu_ref`, so a harmonized HFF cell reaching the clock is

    x_clock = (x_hff - mu_hff) / (sigma_hff + EPS) * sigma_gill + mu_gill

and because ΔAge is a DIFFERENCE, both `mu` terms cancel exactly:

    dAge_harmonized = sum_g w_g . (x_pert,g - x_ctrl,g) . sigma_gill,g / (sigma_hff,g + EPS)
    dAge_direct     = sum_g w_g . (x_pert,g - x_ctrl,g)

So the gain is not assumed — it is the ratio of two quantities this script computes from the same
cells. **The variance floor is applied exactly as `harmonize.py:112` does it** (`sigma` floored at
the MEDIAN sigma over admissible genes), because that floor is what stops near-constant genes
producing an exploding ratio, and omitting it would inflate the answer.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
_RESULTS = Path(__file__).resolve().parents[1] / "results"
_RESULTS.mkdir(exist_ok=True)

IPSC_DAY = 21.0
PREDICTED_GAIN = 2.26
DIRECT_DAY14 = -10.62      # results/diag_pipeline_decompose_results.json (S2)
SHARD_DAY14 = -24.02       # results/diag_gc_hff_signature_results.json


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def floored_sigma(M: np.ndarray) -> np.ndarray:
    """Per-gene SD over control observations, floored at the MEDIAN — `harmonize.py:111-113`.

    The floor is not cosmetic: without it a near-constant gene gives sigma ~ 0 and the
    sigma_ref/sigma_d ratio explodes, which would manufacture the very gain this script is trying
    to measure.
    """
    sigma = M.std(axis=0)
    return np.maximum(sigma, float(np.median(sigma)))


def main() -> int:
    for _s in (sys.stdout, sys.stderr):
        try:
            _s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
    if len(sys.argv) < 3:
        print(__doc__)
        return 2

    from cellfate.data.normalize import normalize_counts
    from cellfate.data.qc import QCConfig, apply_qc
    from cellfate.data.sources import GSE242423SingleCellSource
    dv = _load("dv", ROOT / "experiments" / "diag_dage_variants.py")
    hli = _load("hli", ROOT / "experiments" / "diag_hff_label_identity.py")

    clock = json.loads((ROOT / "configs" / "clocks" / "fleischer_clock.json").read_text("utf-8"))
    W = {k: float(v) for k, v in clock["weights"].items()}

    # ---- Gill: the REFERENCE dataset. Its controls are the day-0 samples ---------------- #
    g_samples, g_genes, g_lin = dv.load_gill(Path(sys.argv[1]))
    g_norm = normalize_counts(g_lin)
    import re as _re
    g_day = np.array([float(m.group(1)) if (m := _re.search(r"_d(\d+)_", s)) else 0.0
                      for s in g_samples])
    g_ctrl = g_norm[g_day == 0.0]
    print(f"\n[gill]  {len(g_samples)} samples, {int((g_day == 0.0).sum())} controls (day 0)")

    # ---- HFF: the dataset carrying the gain --------------------------------------------- #
    hff = Path(sys.argv[2])
    genes_file = next(hff.glob("*genes.tsv.gz"))
    samples = []
    for mtx in sorted(hff.glob("*.matrix.mtx.gz")):
        bc = mtx.with_name(mtx.name.replace(".matrix.mtx.gz", ".barcodes.tsv.gz"))
        if bc.exists():
            samples.append({"matrix": str(mtx), "barcodes": str(bc),
                            "label": mtx.name.split(".")[0].split("_")[-1]})
    src = GSE242423SingleCellSource(samples, str(genes_file), min_genes=hli.MIN_GENES,
                                    max_cells_per_sample=hli.MAX_CELLS,
                                    cells_per_run=hli.CELLS_PER_RUN)
    qc = QCConfig(max_mito_frac=0.20, min_genes=hli.MIN_GENES)
    mats, days, h_genes = [], [], None
    for chunk in src.plan():
        raw = apply_qc(src.fetch(chunk), qc)
        if len(raw.obs) == 0:
            continue
        h_genes = raw.genes
        mats.append(normalize_counts(raw.counts))
        days.append(raw.obs["time_h"].to_numpy(dtype=float) / 24.0)
    h_day = np.concatenate(days)
    print(f"[hff]   {len(h_day)} cells, {int((h_day == 0.0).sum())} controls (day 0)")

    # ---- shared gene space, clock genes only --------------------------------------------- #
    gi_g = {g: i for i, g in enumerate(g_genes)}
    gi_h = {g: i for i, g in enumerate(h_genes)}
    shared = sorted(set(gi_g) & set(gi_h) & set(W))
    print(f"[genes] {len(shared)} clock genes present in BOTH datasets")

    ig = np.array([gi_g[g] for g in shared])
    ih = np.array([gi_h[g] for g in shared])
    w = np.array([W[g] for g in shared])

    sig_gill = floored_sigma(g_ctrl[:, ig])
    h_ctrl = np.concatenate([m[:, ih] for m, d in zip(mats, days, strict=True)])[h_day == 0.0]
    sig_hff = floored_sigma(h_ctrl)
    ratio = sig_gill / (sig_hff + 1e-8)

    # ---- ΔAge both ways, on the SAME cells ------------------------------------------------ #
    Xh = np.concatenate([m[:, ih] for m in mats])
    keep = h_day != IPSC_DAY
    Xh, hd = Xh[keep], h_day[keep]
    base = Xh[hd == 0.0].mean(axis=0)
    delta = Xh - base                                     # per-gene deviation from the control mean

    d_direct = delta @ w
    d_harm = delta @ (w * ratio)
    d14_direct = float(d_direct[hd == 14.0].mean())
    d14_harm = float(d_harm[hd == 14.0].mean())
    gain = d14_harm / d14_direct if d14_direct else float("nan")

    print(f"\n  sigma ratio (gill/hff) over {len(shared)} clock genes: "
          f"median {np.median(ratio):.3f}, mean {ratio.mean():.3f}")
    print(f"\n  HFF day-14 ΔAge, direct       = {d14_direct:+8.2f} yr")
    print(f"  HFF day-14 ΔAge, harmonized   = {d14_harm:+8.2f} yr")
    print(f"  => MEASURED GAIN = {gain:.3f}   (predicted {PREDICTED_GAIN})")
    print(f"  recorded shard value {SHARD_DAY14:+.2f}; harmonized route gives {d14_harm:+.2f}")

    ok = abs(gain - PREDICTED_GAIN) <= 0.5
    out = {"script": "diag_harmonization_gain", "utc": datetime.now(UTC).isoformat(),
           "n_shared_clock_genes": len(shared),
           "sigma_ratio": {"median": float(np.median(ratio)), "mean": float(ratio.mean()),
                           "p10": float(np.percentile(ratio, 10)),
                           "p90": float(np.percentile(ratio, 90))},
           "day14_direct": d14_direct, "day14_harmonized": d14_harm,
           "measured_gain": float(gain), "predicted_gain": PREDICTED_GAIN,
           "shard_reference": SHARD_DAY14,
           "verdict": "CONFIRMED" if ok else "ATTRIBUTION_WRONG"}
    (_RESULTS / "diag_harmonization_gain_results.json").write_text(json.dumps(out, indent=2), "utf-8")
    print(f"\n  => {out['verdict']}")
    print("  wrote results/diag_harmonization_gain_results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

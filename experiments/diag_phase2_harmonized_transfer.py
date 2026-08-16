"""PHASE 2 -- out-of-cohort age transfer using the project's OWN Harmonizer.  (read-only)

    python experiments/diag_phase2_harmonized_transfer.py

Pre-registered in `plans/THREE_TESTS_PREREG.md` Phase 2.

WHY
---
`diag_age_transfer` compared raw features and a crude per-cohort z-score and skipped
`cellfate.data.harmonize.Harmonizer`, which exists for exactly this. Raw transfer gave MAE 70-118
with NEGATIVE correlation -- a scale failure, not evidence that age is unreadable across cohorts.
This runs the same test through the machinery the pipeline actually uses.

The Harmonizer holds per-dataset CONTROL mu/sigma on a shared gene space and applies a Z-transform
plus the Gill Projection. Both sides here are control populations -- GSE113957 is entirely
untreated fibroblasts, and the held-out cohorts are day-0 fibroblasts -- so it is being used for
the case it was designed for.

"TRANSDUCTIVE", RESTATED HONESTLY
--------------------------------
I earlier called the z-score correction non-deployable because it needs the whole test cohort.
That is true of ANY batch alignment, harmonizer included: a batch effect cannot be estimated from
one sample. Samples are processed in batches in practice, so per-batch alignment is a DEPLOYMENT
CONSTRAINT, not a defect. What would be a defect is needing the test LABELS -- none of these do.

PRE-REGISTERED READING
  TRANSFER WORKS  harmonized MAE <= 20 yr AND Spearman >= 0.6 on the Gill cohort, majority of alphas
  DOES NOT        otherwise
GSE165177 is reported but EXCLUDED FROM THE VERDICT: its donors span 38-53 yr against a ~12 yr
instrument error, so it has no dynamic range to rank. Declared here, before the numbers.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
_RESULTS = Path(__file__).resolve().parents[1] / "results"
sys.path.insert(0, str(ROOT / "src"))


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


dre = _load("dre", "experiments/diag_residual_expression.py")
dac = _load("dac", "experiments/diag_age_capacity.py")
dat = _load("dat", "experiments/diag_age_transfer.py")

from cellfate.data.build_dataset import GenePanel  # noqa: E402
from cellfate.data.harmonize import Harmonizer  # noqa: E402

REF = "gse113957"


def main() -> None:
    meta = dac.load_meta()
    tr_expr, tr_genes = dac.load_expression(meta.gsm.tolist())
    keep = meta.disease.eq("Normal").to_numpy() & meta.age.notna().to_numpy()
    tr_expr, tr_age = tr_expr[keep], meta.age.to_numpy(float)[keep]

    print("=" * 100)
    print(f"PHASE 2 -- transfer through the project Harmonizer.  train n={len(tr_age)}")
    print(f"  pre-registered: WORKS iff MAE <= {dat.TRANSFER_MAE_BAR} AND spearman >= "
          f"{dat.TRANSFER_RHO_BAR} on Gill, majority of alphas")
    print("  GSE165177 reported but EXCLUDED from the verdict (38-53 yr range, no dynamic range)")
    print("=" * 100)

    panel = set(GenePanel.load(str(dre.PANEL_PATH)).genes)
    res: dict = {"n_train": int(len(tr_age)), "cohorts": {}}

    for cname, (ex, se) in (("GSE165176 Gill", dat.GILL), ("GSE165177", dat.TRANS)):
        d = dat.load_bulk_day0(ex, se, 12)
        if d.empty:
            continue
        te_age = d["__age__"].to_numpy(float)
        te = d.drop(columns="__age__")
        te_genes = list(te.columns)
        te_expr = te.to_numpy(dtype=np.float64)

        # Fit the harmonizer on both cohorts' CONTROL populations, reference = the training cohort.
        h = Harmonizer.fit({REF: [(tr_expr, tr_genes)], "held_out": [(te_expr, te_genes)]},
                           ref_dataset=REF)
        Htr = h.transform(tr_expr, tr_genes, REF)
        Hte = h.transform(te_expr, te_genes, "held_out")
        shared = [g for g in h.genes if g in panel]
        hi = {g: i for i, g in enumerate(h.genes)}
        cols = np.array([hi[g] for g in shared])
        base = float(np.abs(te_age - np.median(tr_age)).mean())

        print(f"\n[{cname}]  n={len(te_age)}  ages {sorted(te_age.astype(int))}  "
              f"harmonised gene space {len(h.genes)}, panel-shared {len(shared)}")
        print(f"  predicting the TRAINING median ({np.median(tr_age):.0f} yr) gives MAE {base:.2f}")
        print(f"  {'alpha':>10}{'MAE':>9}{'spearman':>10}{'pearson':>9}   pass")
        n_pass, per = 0, {}
        for a in dat.ALPHAS:
            p = dre.ridge_fit_predict(Htr[:, cols], tr_age, Hte[:, cols], a)
            mae = float(np.abs(p - te_age).mean())
            rho = (float(pd.Series(p).corr(pd.Series(te_age), method="spearman"))
                   if np.std(p) > 0 and np.std(te_age) > 0 else float("nan"))
            pe = (float(np.corrcoef(p, te_age)[0, 1])
                  if np.std(p) > 0 and np.std(te_age) > 0 else float("nan"))
            ok = bool(mae <= dat.TRANSFER_MAE_BAR and np.isfinite(rho)
                      and rho >= dat.TRANSFER_RHO_BAR)
            n_pass += ok
            per[str(a)] = {"mae": mae, "spearman": rho, "pearson": pe, "pass": ok}
            print(f"  {a:>10.0f}{mae:>9.2f}{rho:>10.3f}{pe:>9.3f}   {'YES' if ok else 'no'}")
        verdict = "TRANSFER WORKS" if n_pass > len(dat.ALPHAS) / 2 else "does not"
        counts = " (EXCLUDED from verdict)" if "165177" in cname else ""
        print(f"  -> harmonized: {verdict}  ({n_pass}/{len(dat.ALPHAS)} alphas){counts}")
        res["cohorts"][cname] = {"n": len(te_age), "baseline_mae": base, "alphas": per,
                                 "n_pass": n_pass, "verdict": verdict,
                                 "counts_toward_verdict": "165177" not in cname}

    _RESULTS.mkdir(exist_ok=True)
    p = _RESULTS / "diag_phase2_harmonized_transfer_results.json"
    p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    print(f"\nSaved: {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()

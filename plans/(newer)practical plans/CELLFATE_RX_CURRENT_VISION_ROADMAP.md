# CellFate-Rx — Current Vision Roadmap

**Purpose:** High-level roadmap only.  
This document is intentionally less specific than the individual stage plans. Each detailed stage must have its own pre-registration, tests, and PASS/FAIL rules before execution.

---

## STAGE 21 — Prospective Fate Feasibility

Use the datasets we already have to determine whether a legitimate **early molecular state → later fate** prediction task can be constructed without leakage, and whether early transcriptomic information adds anything beyond time, treatment, and donor metadata. This stage does **not** change CellFate-Rx. Its only job is to decide whether the prospective direction is scientifically testable with the available data.

**PASS → STAGE 22**

---

## STAGE 22 — Prospective CellFate-Rx

If Stage 21 passes, adapt the fate side of CellFate-Rx to a true prospective target: molecular state measured before the outcome, combined with the proposed perturbation, predicting an independently observed later fate. Compare the full model against simple `X-only`, `treatment/metadata-only`, and `X + treatment` baselines. The goal is to prove that the model is predicting the future rather than recognizing the current state.

**PASS → STAGE 23**

---

## STAGE 23 — Independent Replication

Test the frozen prospective result on an independent biological dataset, donor set, cell line, or experimental system that was not used to develop the Stage-22 model. No tuning should be performed against the replication set. This stage decides whether the prospective signal is a real generalizable result rather than something specific to one dataset.

**PASS → EXECUTE**

---

## STAGE 24 — Final Execution / Evidence Lock

Freeze the model, data definitions, evaluation code, baselines, calibration, and all manuscript-facing metrics. Run the complete final evaluation exactly once under the locked protocol and generate the final tables, figures, confidence intervals, ablations, negative controls, and reproducibility outputs. After this stage, the scientific results are treated as fixed; no result-driven model tuning is allowed.

**PASS → NOTE**

---

## STAGE 25 — Scientific Note / Claim Lock

Write the final scientific interpretation before drafting the full manuscript. Record exactly what the project proves, what it does not prove, which old claims were withdrawn, what the prospective result adds, and which limitations must travel with every headline number. This becomes the single source of truth for the paper so the manuscript cannot gradually overstate the evidence.

**PASS → PUBLISH STAGES**

---

# PUBLISH STAGES

## STAGE 26 — Manuscript Build

Turn the locked evidence into the paper. The prospective treatment-conditioned fate result becomes the central contribution, while the existing CellFate-Rx work supports it: timepoint confounding, ΔAge circularity, calibration, uncertainty, reproducibility, and the negative results that motivated the prospective evaluation design. Build the figures, methods, results, discussion, supplementary material, and reproducibility package from the frozen Stage-24 outputs.

**PASS → STAGE 27**

---

## STAGE 27 — Internal Review / Submission Lock

Audit the manuscript as if reviewing someone else's work. Check every headline against the frozen result files, verify that limitations are explicit, remove claims not directly supported by evidence, confirm that all analyses are reproducible, and perform a final prior-art/novelty check. Only corrections are allowed here; no new fishing analyses or model optimization.

**PASS → SUBMIT**

---

## STAGE 28 — Submission / Preprint

Submit the manuscript to the selected journal and, if desired, release the preprint and reproducibility materials. At this point the first CellFate-Rx research arc is complete. Reviewer-requested analyses should be handled as clearly labelled follow-up work rather than silently changing the original evidence base.

**PASS → PUBLICATION / REVISION CYCLE**

---

# AFTER PAPER 1 — PRODUCT / GENERATION 2

After submission, return to the larger CellFate-Rx vision: multiple perturbations, stronger held-out-treatment generalization, richer fate classes, better uncertainty, rejuvenation endpoints, and eventually an integrated decision score such as RES. These are **not requirements for Paper 1**. Paper 1 should establish one strong prospective capability first; the broader experimental tool can then be perfected on top of that validated foundation.

---

# One-line roadmap

```text
STAGE 21
Prospective feasibility
        ↓ PASS

STAGE 22
Prospective CellFate-Rx
        ↓ PASS

STAGE 23
Independent replication
        ↓ PASS

STAGE 24
FINAL EXECUTION / EVIDENCE LOCK
        ↓ PASS

STAGE 25
SCIENTIFIC NOTE / CLAIM LOCK
        ↓ PASS

STAGE 26
MANUSCRIPT BUILD
        ↓ PASS

STAGE 27
INTERNAL REVIEW / SUBMISSION LOCK
        ↓ PASS

STAGE 28
SUBMIT / PREPRINT / PUBLISH
```

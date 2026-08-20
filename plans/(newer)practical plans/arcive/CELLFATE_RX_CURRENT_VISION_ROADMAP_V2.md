# CellFate-Rx — Current Vision Roadmap v2

**Purpose:** High-level roadmap.  
The current priority is no longer to squeeze another predictive result from the existing local datasets. The critical path is to obtain a valid prospective dataset, prove a simple prospective signal, then adapt CellFate-Rx and replicate it. The existing audited CellFate-Rx work remains preserved and may support a separate methods/audit paper, but it is not allowed to substitute for independently grounded prospective evidence.

---

## STAGE 21 — Prospective Data Qualification

Audit the datasets already on disk to document exactly why they can or cannot support a true early-state → later-outcome experiment, then move directly to public lineage-resolved datasets. Qualify candidates such as reprogramming lineage datasets and treatment-response lineage datasets by checking pre-outcome RNA, clone/lineage linkage, an independently measured later outcome, treatment metadata, sufficient independent units, and reconstructable public files. Stage 21 passes only when one dataset has been downloaded, reconstructed, sanity-checked against its source study, and frozen into a valid `X + U -> Y_future` table.

**PASS → STAGE 22**

---

## STAGE 22 — Prospective Baseline Experiment

Before changing CellFate-Rx, test whether the qualified dataset actually contains useful prospective information. Compare metadata/treatment-only, transcriptome-only, and transcriptome + treatment models under grouped biological holdouts, with log loss as the primary metric and a structure-preserving molecular-shuffle null. The key question is whether early molecular state adds information beyond treatment and metadata; if it does not, deep-model development stops rather than trying architectures until something looks positive.

**PASS → STAGE 23**

---

## STAGE 23 — Prospective CellFate-Rx

Adapt the fate side of CellFate-Rx to the frozen prospective target from Stage 21. The model now receives molecular state measured before the outcome plus the proposed perturbation and predicts the independently observed future outcome. Compare it against the strongest simple Stage-22 baseline, preserve grouped holdouts and calibration, and require that the neural model adds real value rather than merely reproducing a linear predictor.

**PASS → STAGE 24**

---

## STAGE 24 — Independent Replication

Freeze the Stage-23 model and test the prospective result on a second biological dataset, cell system, donor set, or perturbation setting not used to discover the result. The purpose is to show that the prospective signal is not peculiar to one lineage experiment. No tuning against the replication data is allowed; any adaptation must be defined before the replication outcome is examined.

**PASS → STAGE 25**

---

## STAGE 25 — Final Evidence Lock

Freeze the model, target definitions, datasets, baselines, calibration, split rules, negative controls, and manuscript-facing metrics. Run the complete final evaluation under the locked protocol and generate the final tables, figures, confidence intervals, ablations, OOD results, and reproducibility artifacts. After this stage, no result-driven tuning is allowed; changes must be treated as a new experiment rather than silently improving the final numbers.

**PASS → STAGE 26**

---

## STAGE 26 — Scientific Claim Lock

Write the authoritative scientific interpretation before drafting the paper. State exactly what the prospective experiment proves, what remains a surrogate or limitation, which older CellFate-Rx claims were withdrawn, and how the earlier circularity, time-confounding, calibration, uncertainty, and scorecard findings motivate the final evaluation design. This document becomes the claim firewall for the manuscript.

**PASS → STAGE 27**

---

## STAGE 27 — Main Manuscript Build

Build the primary CellFate-Rx paper around the prospective result: `X_before + U -> independently observed Y_future`. The existing project work becomes supporting evidence showing why same-state labels, time-only performance, weak uncertainty and badly chosen decision metrics can overstate predictive capability. The manuscript should present the prospective benchmark, baselines, model comparison, replication, calibration, failure modes, and limitations from the frozen evidence only.

**PASS → STAGE 28**

---

## STAGE 28 — Internal Review / Submission Lock

Audit the manuscript like a hostile external reviewer. Check every headline against frozen result files, verify that no expression-derived surrogate is described as independent biological fate, confirm that all units of replication are biologically correct, rerun reproducibility checks, and perform the final novelty/prior-art review. Only corrections are allowed here; no new model fishing or metric selection.

**PASS → SUBMIT**

---

## STAGE 29 — Submission / Preprint

Submit the main prospective CellFate-Rx manuscript to the selected journal and release a preprint/reproducibility package if appropriate. Reviewer-requested analyses are handled as explicitly labelled follow-up work. Paper 1's evidence base remains frozen rather than being retroactively rewritten around reviewer requests.

**PASS → PUBLICATION / REVISION CYCLE**

---

# PARALLEL TRACK — Existing Methods / Audit Paper

The current CellFate-Rx work already contains a coherent methodological audit: same-state ΔAge circularity, timepoint confounding in fate, the methylation instrument floor, RES failing because uncertainty exceeds signal, calibration-target mismatch, scorecard decision bugs, and a heavily tested reproducibility framework. This may become a separate methods/audit manuscript if a focused novelty and venue review says the contribution is sufficiently distinct. It should proceed in parallel only if it does not delay the prospective data path; the prospective paper remains the main scientific target.

---

# AFTER THE MAIN PROSPECTIVE PAPER — PRODUCT / GENERATION 2

After submission, expand from one validated prospective task toward the original full product vision: multiple perturbations, held-out-treatment generalization, richer fate classes, stronger uncertainty, better OOD handling, rejuvenation endpoints, and eventually a validated integrated decision score such as RES. These are not prerequisites for the first prospective paper. The first paper needs one strong, independently grounded future-outcome capability that survives simple baselines and independent replication.

---

# One-line roadmap

```text
STAGE 21
QUALIFY + RECONSTRUCT PROSPECTIVE DATA
        ↓ PASS

STAGE 22
PROSPECTIVE BASELINES
        ↓ PASS

STAGE 23
PROSPECTIVE CELLFATE-RX
        ↓ PASS

STAGE 24
INDEPENDENT REPLICATION
        ↓ PASS

STAGE 25
FINAL EVIDENCE LOCK
        ↓ PASS

STAGE 26
SCIENTIFIC CLAIM LOCK
        ↓ PASS

STAGE 27
MAIN MANUSCRIPT
        ↓ PASS

STAGE 28
INTERNAL REVIEW / SUBMISSION LOCK
        ↓ PASS

STAGE 29
SUBMIT / PREPRINT / PUBLISH
```

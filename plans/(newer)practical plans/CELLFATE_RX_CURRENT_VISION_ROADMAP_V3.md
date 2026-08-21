# CellFate-Rx — Current Vision Roadmap v3

**Primary scientific target:** Build a genuinely new prospective treatment-conditioned fate model, not merely validate the existing classifier on another dataset.

The intended paper is built around:

```text
pre-intervention molecular state X
+
candidate perturbation U
        ↓
independently measured future outcome Y
```

and, where the experimental design permits:

```text
for the same / comparable starting state:
rank which candidate perturbation is most likely to succeed
```

The old CellFate-Rx audit remains valuable supporting evidence and may become a separate methods paper, but the main research arc now creates a new benchmark, a new prospective model, and a new state-conditioned treatment-ranking capability.

---

## STAGE 21 — Prospective Data Qualification

Audit the local datasets once, then qualify and reconstruct public lineage-resolved prospective data. The stage seeks **two complementary capabilities**: a reprogramming dataset with pre-outcome RNA and an independently measured later outcome (**Role A — required**), plus a multi-perturbation lineage dataset where treatment identity genuinely varies (**Role B — high value, non-blocking**). Stage 21 ends when the Role-A data are reconstructable, statistically usable, leakage-checked, and frozen for benchmark construction; a qualified Role B is folded in on the same terms when the search finds one.

**PASS → STAGE 22**, on either path:

```text
FULL_DATA_PATH   Role A + Role B qualified
CORE_DATA_PATH   Role A qualified, Role B not found inside the frozen search budget
                 -> proceed; the treatment-ranking contribution is SCOPED DOWN
DATA_BLOCKED     Role A not found -> do not build the prospective architecture yet
```

*Corrected 2026-08-21: this paragraph originally read "The stage deliberately seeks two
complementary capabilities … Stage 21 ends only when those data are reconstructable, statistically
usable, leakage-checked, and frozen for benchmark construction", i.e. both roles required. That was
V2's PASS rule. It is superseded by `STAGE_21_PROSPECTIVE_DATA_QUALIFICATION_V3.md` §9.8, which
makes Role B non-blocking so a missing multi-perturbation dataset scopes the treatment-ranking
contribution down instead of holding the paper. Executed 21C returned `FULL_DATA_PATH`, so both
roles are in fact qualified — the rule is corrected here regardless, because it governs what a PASS
means, not what happened to be found.*

---

## STAGE 22 — Prospective Benchmark Build

Create a standardized lineage-grounded benchmark from the qualified datasets. Map each experiment into a common `X_before + U -> Y_future` schema, define biologically correct split groups, freeze the outcome ontology, metrics, treatment metadata, leakage rules, and external-test partitions, and reproduce the source-study bookkeeping. The benchmark itself is one contribution: it separates true prospective prediction from post-treatment state classification and prevents random-cell or same-state shortcuts.

**PASS → STAGE 23**

---

## STAGE 23 — Learnability & Interaction Gate

Before building a neural model, prove that the benchmark contains the information needed for the CellFate-Rx thesis. Compare metadata/treatment-only, molecular-state-only, additive `X+U`, and a simple explicit `X×U` interaction baseline under grouped biological holdouts. Require that early state adds information beyond treatment and that treatment adds information beyond state; on the multi-treatment dataset, test whether state-treatment interaction improves prediction beyond purely additive effects. If the interaction thesis is unsupported, stop rather than inventing architecture around noise.

**PASS → STAGE 24**

---

## STAGE 24 — Build CellFate-Rx Prospective

This is the main invention stage. Build a new prospective architecture designed specifically for `X_before + U -> Y_future`, with separate molecular-state and perturbation encoders, an explicit state×perturbation interaction module, probabilistic future-outcome prediction, calibrated uncertainty, and OOD-aware abstention. The model must be compared against the strongest Stage-23 simple baseline and must earn its complexity through better predictive value, treatment ranking, calibration, or generalization rather than simply being deeper.

**PASS → STAGE 25**

---

## STAGE 25 — State-Conditioned Treatment Ranking Challenge

Test the capability that makes the new system useful: for the same or comparable starting molecular state, can it distinguish which candidate treatment is more likely to produce the desired future outcome? Use multi-treatment lineage data to evaluate pairwise treatment ranking, within-state treatment selection, and decision utility against treatment-only and state-only baselines. This should become one of the paper's central figures because it demonstrates that the model has learned `cell state × perturbation`, not merely "good cells" or "good treatments."

**PASS → STAGE 26**

---

## STAGE 26 — Held-Out Perturbation / Generalization Challenge

Where the data make it scientifically resolvable, hold out entire perturbations, treatment classes, experimental batches, or biological contexts and test whether the frozen model transfers. This is a harder claim than prediction under already-seen treatments, so it is separated from Stage 25. A failure here does not erase the known-treatment prospective result, but it determines whether the paper may claim generalization to unseen interventions.

**PASS OR SCOPED LIMIT → STAGE 27**

---

## STAGE 27 — Independent Biological Replication

Freeze the model and core evaluation rules, then test the main prospective result on a separate biological dataset/system not used to choose the architecture. Prefer a second lineage-resolved system with an independent future outcome. The goal is to show that the prospective signal and treatment-conditioning result are not peculiar to one experiment. No model fishing is allowed against the replication set.

**PASS → STAGE 28**

---

## STAGE 28 — Calibration, OOD & Decision Utility

Turn the predictor into a trustworthy decision instrument. Evaluate probability calibration against the actual hard future outcome, uncertainty/risk-coverage, abstention behavior, OOD detection, and whether rejecting uncertain cases improves decision quality. Reuse the lessons from the old Stage-16/17 calibration work, but calibrate only against the prospective ground truth. This stage is where the engineering history becomes a methodological advantage rather than baggage.

**PASS → STAGE 29**

---

## STAGE 29 — Final Evidence Lock

Freeze benchmark versions, model weights, splits, baselines, calibration, OOD rules, metrics, negative controls and manuscript-facing figures. Run the complete final evaluation once under the locked protocol. After this point, result-driven tuning is forbidden; anything new is a separately labelled experiment.

**PASS → STAGE 30**

---

## STAGE 30 — Scientific Claim Lock

Write the authoritative claim document before drafting the manuscript. Separate what is established under known treatments from what generalizes to unseen treatments, distinguish biological outcomes from transcriptomic surrogates, record every limitation, and state exactly which old CellFate-Rx claims were withdrawn. This becomes the claim firewall for every abstract, figure caption and conclusion.

**PASS → STAGE 31**

---

## STAGE 31 — Main Manuscript Build

Build the main paper around four contributions: **(1)** a prospective lineage-grounded benchmark, **(2)** a new treatment-conditioned CellFate-Rx architecture, **(3)** state-conditioned treatment ranking from pre-intervention molecular state, and **(4)** calibrated uncertainty / abstention under biologically correct holdouts. The old project's circularity, time-confounding, scorecard and calibration failures become supporting evidence explaining why this stricter evaluation framework is necessary.

**PASS → STAGE 32**

---

## STAGE 32 — Internal Review / Submission Lock

Review the manuscript as a hostile external reviewer. Verify every headline against frozen artifacts, confirm that clone/donor/replicate counts are not pseudo-replication, ensure no surrogate is called independent fate, rerun reproducibility checks, and perform the final novelty/prior-art audit. Only corrections are allowed; no performance fishing.

**PASS → SUBMIT**

---

## STAGE 33 — Submission / Preprint

Submit the prospective CellFate-Rx manuscript and release the preprint/reproducibility package if appropriate. Reviewer-requested analyses are recorded as explicit follow-up work rather than silently changing the original evidence base.

**PASS → PUBLICATION / REVISION CYCLE**

---

# PARALLEL TRACK — Methods / Audit Paper

The existing project contains a coherent methodological audit: same-state ΔAge circularity, timepoint confounding, instrument disagreement, RES collapsing because uncertainty exceeds signal, calibration-target mismatch, scorecard bugs, and strong reproducibility discipline. That may support a separate methods/audit manuscript after a focused novelty/venue review. It should run in parallel only if it does not delay Stages 21–27.

---

# AFTER THE PROSPECTIVE PAPER — GENERATION 2 PRODUCT

After the first prospective paper, return to the larger experimental-triage product:

```text
multiple perturbation modalities
richer fate classes
dose optimisation
stronger held-out-treatment transfer
rejuvenation endpoints
future ΔAge with independent anchors
validated RES / integrated utility
```

Paper 1 does not need all of that.

It needs one strong new capability:

> **given a pre-intervention molecular state and candidate treatment, predict and rank independently observed future outcomes under biologically valid prospective evaluation.**

---

# One-line roadmap

```text
STAGE 21
PROSPECTIVE DATA QUALIFICATION
        ↓ PASS

STAGE 22
PROSPECTIVE BENCHMARK BUILD
        ↓ PASS

STAGE 23
LEARNABILITY + X×U INTERACTION GATE
        ↓ PASS

STAGE 24
BUILD CELLFATE-RX PROSPECTIVE
        ↓ PASS

STAGE 25
STATE-CONDITIONED TREATMENT RANKING
        ↓ PASS

STAGE 26
HELD-OUT PERTURBATION / GENERALIZATION
        ↓ PASS OR SCOPED LIMIT

STAGE 27
INDEPENDENT BIOLOGICAL REPLICATION
        ↓ PASS

STAGE 28
CALIBRATION + OOD + DECISION UTILITY
        ↓ PASS

STAGE 29
FINAL EVIDENCE LOCK
        ↓ PASS

STAGE 30
SCIENTIFIC CLAIM LOCK
        ↓ PASS

STAGE 31
MAIN MANUSCRIPT
        ↓ PASS

STAGE 32
INTERNAL REVIEW / SUBMISSION LOCK
        ↓ PASS

STAGE 33
SUBMIT / PREPRINT / PUBLISH
```

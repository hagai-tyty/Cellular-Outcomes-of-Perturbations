# WORK ORDER — 2026-08-14. Five steps, in order, with decision rules fixed in advance.

**Revision 2**, after a second external review. Changes from revision 1, all adopted:

| # | what was wrong | fix |
|---|---|---|
| 1 | "95 % CI, n = 3" could read as a large-sample robust cluster interval | stated explicitly as a **t-based interval on three donor-level estimates**; with 3 clusters the method matters and asymptotic cluster-robust results do not apply |
| 2 | CpG coverage reported only as a fraction | also report **which** CpGs — two matrices at 60 % weight coverage can differ biologically |
| 3 | P1.2 claimed concordance ruled out "a transcriptomic-clock artefact" | softened: it rules out an artefact **unique to one modality**, not a shared one |
| 4 | the ratio `ΔAge_RNA / ΔAge_Meth` was the primary magnitude test | **demoted to secondary/descriptive.** A near-zero denominator makes it pathological. **Regression is now primary** |
| 5 | P1.4 regressed **sample-level** points while declaring donor the independent unit | **contradiction removed** — donor-level or mixed-effects with donor as a random effect; no sample-level R² quoted as though n = samples |
| 6 | "slope stable across donors" implied a test n = 3 cannot power | reworded descriptively |
| 7 | P3's progress coordinate could be learned from all cells, leaking the held-out day | **fit on training timepoints only, then project held-out cells into that fixed coordinate** |
| 8 | "overall safety defined afterwards" invited choosing the rule after seeing results | **pre-specified now: no composite is built** |

**Order: ~~P2~~ → P1 → P3 → P4 → P5.** P2 is complete.

**Out of scope, deliberately:** no retrain, no acquisition, no harmonizer change (that one would
confound P1), 3b/3c/3d unwritten, Stage 5 unentered.

---

# P2 — ✅ DONE. The outward-facing material is retired, not corrected.

The expert consultation was dropped, so there is nothing to correct — the question becomes what to
do with material that carries known-wrong numbers and will not be used.

- `plans/EXPERT_MESSAGE.md` — **deleted.** It was a draft of an action no longer being taken. Git
  history retains it; nothing that records a *finding* lives only there.
- `plans/EXPERT_BRIEF.md` — **kept, bannered RETIRED.** Deleting it would lose the provenance index
  mapping each claim to the artefact that produced it, which is still accurate. The banner states
  that every ΔAge figure inside is superseded twice — by the pseudoreplication fix and by the
  double-`log1p` fix — and gives the current values.
- `plans/DATA_REQUIREMENT_SECOND_TIMECOURSE.md` — **unchanged.** Already marked ON HOLD, and its
  technical spec is independent of the ΔAge numbers.

*Records are annotated, not erased; a draft of an abandoned action is deleted.*

---

# P1 — Does the transcriptomic ΔAge survive comparison with an independent ageing modality?

**Not** *"is −42 real"*. A second imperfect clock cannot answer that.

`GSE165179` is the methylation twin of `GSE165177`: same three donors, same arms, same days, same
sample names. **Primary clock: Horvath skin & blood 2018** (fitted for fibroblasts); multi-tissue
2013 reported alongside, never substituted after the fact.

## P1.0 — Step zero, before anything is computed

1. **Re-read what the recorded ρ = 0.568 methylation ceiling actually measured** before citing it
   as a bar. Do not cite it until re-read.
2. **CpG coverage**: report the fraction of each clock's CpGs present, the fraction of total
   |weight| they carry, **and the identity of the covered and missing CpGs**. Fractional coverage
   alone can hide two very different biological pictures.

## P1.1 — The adjudication floor. **Read before any comparison.**

Within each (donor, day) there are 2–3 untreated replicates.

> **`SD_meth` = pooled SD of methylation age among control replicates within (donor, day).**

| condition | consequence |
|---|---|
| `SD_meth` ≥ ½ × mean \|ΔAge_Meth\| | **THE YARDSTICK CANNOT ADJUDICATE.** Report it, stop, record the question as open. **Do not read P1.2–P1.5.** |
| otherwise | proceed, quoting `SD_meth` beside every number below |

This is the branch most likely to fire.

## P1.2 — DIRECTION

For each donor, the sign of mean ΔAge in each modality.

| result | reading |
|---|---|
| both modalities negative in 3/3 donors | **CONCORDANT IN DIRECTION** — evidence the shift is **not unique to the transcriptomic modality**. It does **not** exclude an artefact shared by both |
| signs disagree in any donor | **DISCORDANT** — the modalities do not provide convergent evidence for rejuvenation on this trajectory. Investigate before interpreting magnitude. Disagreement can arise from noise, calibration or modality-specific biology, and this branch does not distinguish them |

## P1.3 — MAGNITUDE. **Regression primary, ratio secondary.**

**Primary:** `ΔAge_RNA = α + β · ΔAge_Meth`, fitted per P1.4's clustering rules.

| result | reading |
|---|---|
| β ≈ 1 | the modalities are on a **similar scale** |
| β > 1 | the transcriptomic effect is **larger relative to** the methylation modality |
| β < 1 | **smaller relative to** it |

**Secondary, descriptive only:** the ratio `ΔAge_RNA / ΔAge_Meth`. **It must not drive any
conclusion** — when the denominator approaches zero the ratio and its interval become pathological,
which is entirely plausible here. Report it with the denominator beside it, always.

**No branch concludes −42 is correct or incorrect in absolute terms.** Both clocks are fitted on
chronological age in populations; both could mis-scale acute reprogramming the same way.

## P1.4 — TRAJECTORY. Do the clocks track the same movement?

Per donor, `day → ΔAge` in both modalities, then the relationship between them **respecting the
clustering**:

- **either** regress on the three **donor-level** trajectories,
- **or** a repeated-measures / mixed-effects model with **donor as a random effect**.

**An ordinary sample-level regression is not reported, and no R² is quoted with n equal to the
number of samples.** Samples within a donor are not independent; this project has already made that
error twice and it is not repeated here.

| result | reading |
|---|---|
| positive relationship, and the day-ordering of effect size agrees between modalities | the clocks follow **the same movement** |
| flat despite both being negative on average | they agree the cells move younger but **disagree about when** — two similar endpoints, not a shared trajectory |

## P1.5 — SCALE

From P1.4: slope, intercept, and a fit statistic appropriate to the clustered model.

| result | reading |
|---|---|
| donor-specific slopes **directionally consistent**, and a pooled model approximately linear | consistent with a **systematic scale factor** between modalities |
| curved or saturating | the modalities diverge at large effects; no single factor describes them |
| slopes inconsistent in direction, or a wide pooled interval | shape **not established** at n = 3; report and claim nothing |

*Three donors cannot power a formal test of slope heterogeneity, so no such test is claimed.*

## What P1 cannot do

Establish the true number of years. Exclude a bias shared by both clocks. Or separate "compression
is general" from "compression is specific to the chronological axis" — that needs an anchor neither
modality provides.

---

## P1 RESULT — 2026-08-14. Floor passed. Direction concordant. **RNA reads ~2.5–2.9× larger than methylation on the reprogramming axis.**

*Artefacts: `experiments/dage_meth_concordance.py`, `results/dage_meth_concordance_results.json`.
Methylation ages recomputed through `diag_methylation_anchor`'s own verified loader, linear
predictor, Horvath anti-transform and implied-intercept derivation — nothing re-derived.*

### P1.1 — the adjudication floor **PASSED**, which is not what I expected

| clock | control-replicate SD | mean \|ΔAge_meth\| | blocked if SD ≥ | verdict |
|---|---|---|---|---|
| skin & blood | **3.36 yr** | 10.47 | 5.23 | ✅ passed |
| multi-tissue | **3.41 yr** | 11.70 | 5.85 | ✅ passed |

I flagged this as "the branch most likely to fire." It did not. **The methylation zero-point is
reproducible to ~3.4 yr on contemporaneous replicates — twice as tight as the RNA zero-point
(9.82 yr) on the same samples.** Methylation can adjudicate here.

**But the two methylation clocks disagree with each other by RMS 9.07 yr** across 45 matched cells.
This is a **precision / measurement ceiling on the modality** — the recorded ρ = 0.568 limit
expressed in years — and it should be read as a scale for how much of an RNA-vs-methylation gap
could be instrumental. It is **not** the uncertainty of any particular ΔAge comparison: it is an
aggregate across cells, and individual contrasts have their own precision, which is reported
separately with each one.

### P1.2–P1.5, donor-clustered (n = 3, `t = 4.303`)

| stratum | ΔAge RNA | ΔAge meth (sb) | ΔAge meth (mt) | both < 0 | RNA/meth |
|---|---|---|---|---|---|
| **transient_int** (in the reprogramming phase) | **−64.67** [−77.8, −51.5] | **−23.50** [−41.2, −5.8] | **−26.14** [−38.5, −13.7] | **3/3** both clocks | **2.92 / 2.52** |
| **transient_fib** (returned to fibroblast identity) | −24.53 [−50.4, **+1.4**] | −7.57 [−31.0, **+15.9**] | −8.91 [−16.9, −0.9] | 2/3 · 3/3 | 0.22 / 3.30 |

**DIRECTION — concordant where it counts.** On the intermediates every donor is negative in every
modality, 3/3 across both clocks. The apparent conflict between our earlier −42 and Gill's ~30 was
never a modality disagreement; it was my own arm pooling.

**MAGNITUDE — on the intermediate trajectory, RNA shows a ~2.5–2.9× larger response than
methylation** — the one stratum where both modalities have a solid effect and both denominators
are safely away from zero.

**The ratio's instability is now demonstrated, not just anticipated.** For `transient_fib` it reads
**0.22** against skin & blood and **3.30** against multi-tissue — a 15× swing driven entirely by a
near-zero, sign-uncertain denominator. The review's insistence that the ratio be demoted to
descriptive was correct, and this is the evidence.

### 🔑 What P1 establishes

1. **The ~3.4× chronological-axis compression does not transfer to the reprogramming axis; on the
   intermediate trajectory, RNA shows a ~2.5–2.9× larger response than methylation.** Chronological
   sensitivity says nothing about reprogramming sensitivity — exactly the concern raised in review.
2. **Gill's actual claim — the RETURNED fibroblasts — is not established by our RNA clock**
   (CI includes zero), **and the two methylation clocks disagree about it** (sb includes zero, mt
   excludes it marginally). **Methylation cannot settle that stratum either.**
3. **The unambiguous molecular-age-associated effect in this dataset is confined to the
   intermediate stratum** — both modalities, every donor.

### What P1 does NOT establish

Any absolute number of years. That both clocks are not wrong together — both are fitted on
chronological age in populations and could mis-scale acute reprogramming the same way. Or a
trajectory claim: with 3–4 days per donor and n = 3, P1.4's day-ordering comparison is descriptive
and no slope-heterogeneity test is claimed.

**P1 is closed. Next: P3.**

---

# P3 — Does molecular progress predict risk better than calendar day? *(within HFF only)*

Two lines can run the same course at different speeds, so calendar day cannot transfer while
trajectory position might. This also **reinterprets an existing result**: earlier tests found time
redundant with state along one trajectory and filed it as failure; under this model that is the
expected observation.

**Leakage-proofed design.** Hold out one timepoint. **Fit the progress coordinate on the training
timepoints only** — including any dimensionality reduction or trajectory construction — then
**project the held-out cells into that fixed coordinate** and predict their risk. Nothing used to
build the coordinate may derive from the held-out day, or the "held-out" claim is void.

Compare three predictors: day only, progress only, both.

| result | action |
|---|---|
| progress beats day, held out | supported within one line → **P5 becomes worth doing** |
| progress ties day | no advantage where it should be easiest → **do not spend on P5** |
| progress loses to day | reconsider how progress is inferred before concluding anything |

**Stopping rule:** if it fails within the single line it was designed on, stop.

---

# P4 — Split identity-loss and apoptosis into separate heads

Currently collapsed into one flag. Our own data separates them: the transient arm ran P(loss) 0.613
/ P(death) 0.328 while the failed arms were near-pure loss.

**Pre-specified now, before any result:** the two risks are reported and modelled **separately**,
and **no composite endpoint is constructed.** Not `1 − (1 − P_loss)(1 − P_death)`, not a weighted
sum, not a union. They are biologically distinct failure modes with no reason to peak together, and
choosing a combination rule after seeing which one looks better is a forking path. **If a downstream
decision later requires a single number, that combination rule gets its own pre-registration.**

---

# P5 — Cross-line test on `GSE221739` — **gated on P3**

8 sampled days (D0–D15), 10x, replicate pools per timepoint.

**Stated limit:** BJ is a neonatal foreskin fibroblast line, like HFF. This tests
**protocol-and-line transfer, not donor-background transfer.**

**Do not start unless P3 returns "progress beats day."** The gate is not relaxed.

---

## The standing rule

> **No clock is treated as truth. Each clock is a noisy, potentially mis-scaled measurement whose
> generalisability must be demonstrated independently.**

Decision rules are written before the numbers exist; the donor is the unit for any claim that
generalises; and the instrument brought in to check the first instrument gets the same scrutiny as
the one under test.

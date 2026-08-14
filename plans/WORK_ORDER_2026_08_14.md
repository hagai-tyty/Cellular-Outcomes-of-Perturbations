# WORK ORDER — 2026-08-14. Five steps, in order, with decision rules fixed in advance.

**Supersedes** the informal plan of the same day. Revised after external review, which corrected
one real error: the first version treated the methylation clock as **ground truth** while the
project's own records say its agreement ceiling is limited. It is an **independent, imperfect
yardstick**, and every rule below is written accordingly.

**Order: P2 → P1 → P3 → P4 → P5.** P2 first because it is cheap, externally visible, and currently
wrong. P1 is the first *scientific* priority.

**Not in scope, deliberately:** no retrain, no data acquisition, no harmonizer change, 3b/3c/3d
stay unwritten, Stage 5 unentered. The harmonizer especially — changing it now would confound P1.

---

# P2 — Correct the outward-facing message *(first, because it leaves the building)*

`plans/EXPERT_MESSAGE.md` still quotes ΔAge −17.9 with a CI computed on 12 pseudoreplicated
donor-day cells. Both numbers are superseded.

**Replace with:** ΔAge(transient) **−42.45**, donor-clustered 95 % CI **[−67.39, −17.51]** (n = 3
donors), and the paired transient-vs-failed contrast reported as
*"direction consistent in 3/3 donors, interval includes zero at n = 3"* rather than as an
established effect. Add the two findings that did not exist when it was drafted: the clock
**compresses** the age axis ~3.4× out of domain, and the supplier hypothesis is dead.

**Done when:** no superseded number remains, and the compression finding is stated as an open
question rather than a result.

---

# P1 — Does the transcriptomic ΔAge survive comparison with an independent ageing modality?

**Not** *"is −42 real"*. That question cannot be answered by a second imperfect clock.

`GSE165179` is the **methylation twin** of `GSE165177`: same three donors, same arms, same days
(0/10/13/15/17), same sample names. Two Horvath clocks are already in the repo.
**Primary: skin & blood 2018** (fitted for fibroblasts). Multi-tissue 2013 reported alongside,
never substituted after the fact.

## P1.0 — Step zero, before anything is computed

Two things must be settled first, and both are checks on the *instrument*, not the biology.

1. **Re-read what the recorded ρ = 0.568 methylation ceiling actually measured.** It is carried in
   this project's notes as a limit on methylation-to-methylation agreement, but the contrast it was
   computed on must be confirmed before it is used as a bar. **Do not cite it until re-read.**
2. **CpG coverage.** Report how many of each clock's CpGs are present in the matrix and what
   fraction of total |weight| they carry — the same check that showed 57 % of genes carried 89 % of
   the transcriptomic clock's weight. Low coverage would confound everything downstream.

## P1.1 — The adjudication floor. **Read before any comparison.**

A yardstick can only settle a question if it is finer than the thing being measured. The
contemporaneous controls give this directly: within each (donor, day) there are 2–3 untreated
replicates.

> **Compute `SD_meth` = pooled SD of methylation age among control replicates within (donor, day).**

| condition | consequence |
|---|---|
| `SD_meth` ≥ ½ × the mean \|ΔAge_Meth\| being measured | **THE YARDSTICK CANNOT ADJUDICATE.** Report `SD_meth`, stop, and record that the question stays open. **Do not read P1.2–P1.5.** |
| `SD_meth` < ½ × mean \|ΔAge_Meth\| | proceed, quoting `SD_meth` alongside every number below |

This is the branch the review asked for, and it is the one most likely to fire.

## P1.2–P1.5 — Four separate questions, not one

For every non-control sample, in **both** modalities, against its **contemporaneous** control:

```
ΔAge_X(s) = clock_X(s) − mean( clock_X(controls: same donor, same day) )
```

Donor is the independent unit throughout — n = 3, `t(0.975, df=2) = 4.303`. Sample-level figures
are reported for shape only and never carry a generalisation claim.

### P1.2 — DIRECTION. Do both modalities agree the cells move younger?

| result | reading |
|---|---|
| both modalities negative in **3/3 donors** | **CONCORDANT IN DIRECTION** — evidence the effect is not a transcriptomic-clock artefact |
| modalities disagree in sign in any donor | **DISCORDANT** — at least one molecular-age measure is responding to something other than rejuvenation. **Escalate; magnitude questions become meaningless** |

### P1.3 — MAGNITUDE. How large is the transcriptomic effect *relative to* the other modality?

Report the donor-level ratio `r = ΔAge_RNA / ΔAge_Meth` with its n = 3 interval, **and** `SD_meth`
beside it.

| result | reading |
|---|---|
| interval on `r` contains 1 | the two modalities are **not distinguishable in magnitude** at this n |
| interval on `r` lies **above** 1 | the transcriptomic magnitude is **inflated relative to** the methylation modality — *not* proof that it is inflated relative to truth |
| interval on `r` lies **below** 1 | the transcriptomic magnitude is **attenuated relative to** methylation |

**No branch of this concludes that −42 is correct or incorrect in absolute terms.** Both clocks can
be wrong together.

### P1.4 — TRAJECTORY. Do the two clocks track the same biological movement?

Not endpoints — the whole course. Plot `day → ΔAge_RNA` and `day → ΔAge_Meth` per donor, then
regress paired sample-level ΔAge values.

| result | reading |
|---|---|
| positive slope, and the day-ordering of effect size agrees between modalities | the clocks are following **the same movement** |
| near-zero slope despite both being negative on average | they agree the cells move younger but **disagree about when** — two vaguely similar endpoint numbers, not a shared trajectory |

This is the question a single endpoint comparison cannot answer, and it is the most informative one.

### P1.5 — SCALE. Is the relationship linear, and is there a systematic factor?

From the P1.4 regression: report slope, intercept and R², donor-clustered.

| result | reading |
|---|---|
| approximately linear, slope stable across donors | a **systematic scale factor** between modalities — quantifiable and correctable in principle |
| curved or saturating | the modalities diverge at large effect sizes; no single factor describes the relationship |
| low R² with a large slope interval | shape is **not established** at this n; report and claim nothing |

## What P1 cannot do, stated in advance

It cannot establish the true number of years. It cannot rule out that both clocks share a bias —
they are both fitted on chronological age in populations, and both could mis-scale acute
reprogramming the same way. And it cannot separate "compression is general" from "compression is
specific to the chronological axis" — that needs an anchor neither modality provides.

---

# P3 — Does molecular progress predict risk better than calendar day? *(within HFF only)*

The reframing: the coordinate for safety should be **where a cell is in the reprogramming
trajectory**, not how many days have elapsed. Two lines can run the same course at different
speeds, so day cannot transfer while progress might.

This also **reinterprets an existing result rather than adding one**: earlier tests found time is
redundant with state along a single trajectory and that was filed as a failure. Under this model it
is the *expected* observation — the state already encodes the progress.

**Design.** Infer a per-cell progress coordinate on `GSE242423` (9 timepoints, ~42k cells).
**Leave-one-timepoint-out**: hold out a day, compute progress from the remaining days *without the
held-out day's labels*, and compare three predictors of risk — day only, progress only, both.

| result | reading | action |
|---|---|---|
| progress beats day, held out | the reframing is **supported within one line** | P5 becomes worth doing |
| progress ties day | no advantage where it should be easiest | the reframing buys nothing; **do not spend on P5** |
| progress loses to day | the coordinate is not capturing trajectory position | reconsider how progress is inferred before concluding anything |

**Stopping rule:** if it fails within the single line it was designed on, stop. Do not test transfer
of a construct that does not work at home.

---

# P4 — Split identity-loss and apoptosis into separate heads

Currently collapsed into one "unsafe" flag. Our own data says they separate: the transient arm ran
P(loss) 0.613 / P(death) 0.328 while the failed arms were near-pure loss. Two biologically distinct
failure modes that need not peak together, and collapsing them discards that.

**Deliverable:** P(identity loss) and P(apoptosis) reported and modelled separately across HFF's
nine timepoints, with overall safety defined *afterwards* from the two, rather than the biology
being forced into one binary up front.

---

# P5 — Cross-line test on `GSE221739` **— gated on P3 succeeding**

8 sampled days (D0–D15), 10x, replicate pools per timepoint, downloadable. Tests whether the
molecular ordering of unsafe states reproduces on another line and another induction system.

**Honest limit, stated now:** BJ is a neonatal foreskin fibroblast line, like HFF. This tests
**protocol-and-line transfer, not donor-background transfer.** It cannot answer the question a
genuinely different adult donor would.

**Do not start P5 unless P3 returns "progress beats day."**

---

## The standing rule this order enforces

Every step's decision rules are written **before** the numbers exist, the donor is the unit for any
claim that generalises, and no instrument is treated as truth — including the one brought in to
check the first.

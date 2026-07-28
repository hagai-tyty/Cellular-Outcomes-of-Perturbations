# REVIEW of `STAGE_1_5_1_REVISED.md` — verified independently, 2026-07-26

**Companion file. `STAGE_1_5_1_REVISED.md` is left byte-unmodified.**

Every load-bearing claim was re-derived from the raw data with independently written code (not by
running `diag_gill_replication.py`), and the plan's literature anchor was checked against the source
paper.

> **Verdict in one line: the DIAGNOSIS is correct and valuable; the PROPOSED FIX is not validated
> and should not be implemented as written.** Step 1 (harden) should proceed. Step 2–3 (redefine
> `is_control` to the non-responder arm) must not, until the confound below is resolved.

---

## 1. What I verified as CORRECT — this part is a real advance

| Claim | How checked | Result |
|---|---|---|
| Gill labels: 65 responders / 47 non-responders / 6 iPSC / 6 day-0 | parsed `GSE165176_series_matrix` directly | ✅ **exact** |
| responders day0→peak **+8.2 yr** | re-derived from scratch | ✅ **+8.2** |
| non-responders day0→peak **+36.5 yr** | re-derived | ✅ **+36.5** |
| paired separation **−28.3 yr, 5/6 negative** | re-derived | ✅ **−28.3, 5/6** (my CI [−49.9, −6.7] vs their [−51.6, −5.0] — trivial method difference) |
| **A1** pooling error — `diag_e1_trajectory.py` excludes only iPSC | read the code | ✅ **TRUE.** 47 non-responders *were* averaged in with 65 responders |
| **A2** control = day 0 — `sources.py:471` | read the code | ✅ **TRUE.** `is_ctrl = m["day"] == 0.0 or ctype == "Dermal fibroblast"` |
| **A3** monotonic Spearman on a dip-shaped effect | read the code | ✅ **TRUE** |

**A1 and A3 are genuine methodological errors and their discovery is a real contribution.** Pooling
47 cells-that-cannot-rejuvenate into the treatment arm, and testing a dip with a monotonic rank
statistic, are both mistakes that would independently produce the null results the Stage 1.5 arc got.
That part of the plan stands and should be kept.

---

## 2. 🔴 The proposed fix — four problems, all measured

### 2.1 The stated mechanism is factually wrong

> Plan §1(a): *"it is the clock reading loss of fibroblast identity as ageing."*

Measured identity state (mean marker-set expression, log1p-CP10k):

| group | pluripotency | fibroblast identity | predicted age |
|---|---|---|---|
| day-0 fibroblasts | 0.05 | 2.47 | 69.4 |
| **non-responders (d10–13)** | **0.10** | **2.82** | **102.2** |
| responders (d10–13) | 0.80 | 1.18 | 76.3 |
| iPSC | 1.32 | 0.32 | 49.1 |

**The non-responders never lost fibroblast identity.** They sit essentially where day-0 cells sit
(pluripotency 0.10 vs 0.05; fibroblast 2.82 vs 2.47) — if anything *more* fibroblast-like — yet read
**33 years older**. Their +36.5 yr therefore **cannot** be "loss of fibroblast identity". The plan's
explanation for its own central observation does not survive measurement.

### 2.2 The control is not identity-matched, so the confound is not cancelled

> Plan §1(b): *"Both arms undergo the same OSKM exposure and **the same identity change**; only one
> rejuvenates. Differencing them removes the artefact."*

The two arms differ in identity by **8×** on pluripotency (0.80 vs 0.10). Switching control from
day-0 to non-responders changes the responder-minus-control identity gap like this:

| contrast | Δ pluripotency | Δ fibroblast identity |
|---|---|---|
| responder − **day 0** (current control) | +0.75 | −1.29 |
| responder − **non-responder** (proposed) | +0.70 | **−1.64** |

**The proposed control does not reduce the identity gap. On the fibroblast axis it makes it 27%
larger.** The premise the fix rests on is measurably false.

And the clock is demonstrably reading that axis: across all 112 OSKM-exposed fibroblasts,
`corr(predicted_age, pluripotency) = −0.62`, `corr(predicted_age, fibroblast identity) = +0.62`.

### 2.3 61% of the effect is explained by identity state

Regressing predicted age on donor + arm, then adding the two identity scores (peak window, n=36):

| model | arm effect |
|---|---|
| donor-adjusted only | **−28.3 yr** ← the plan's headline |
| donor + **identity**-adjusted | **−11.1 yr** |

**61% of the effect is identity, not ageing.**

*Fair caveat, stated because it cuts against me:* pluripotency could be a **mediator** (OSKM →
pluripotency → genuine rejuvenation), in which case adjusting for it removes real effect. That
objection is legitimate for this regression — **but not for §2.1**, where non-responders have day-0
identity and still read +33 yr. That comparison is mediator-free and is the damning one.

### 2.4 The "replication of Gill" does not hold — and the sign is opposite

The plan's strongest rhetorical claim is *"Gill reported ~−30 yr. We measure −28.3 yr. That is a
replication."* Checked against [Gill et al. 2022, eLife 71624](https://elifesciences.org/articles/71624):

| | Gill | This plan |
|---|---|---|
| clock used | **their own custom random-forest predictor**, built *because* "existing transcription clocks failed to accurately predict the age of our negative control samples" | **Fleischer** — an existing transcription clock, i.e. the class Gill reports as failing |
| control compared against | **negative controls** (untreated cells) | **failing-to-reprogram cells** |
| what non-responders did | *"a moderate **reduction** in transcription age in cells that failed to transiently reprogram"* | **+36.5 yr older** |

Three separate mismatches: different instrument, different control arm, and **opposite sign on the
very samples the plan proposes to use as its control.** Matching Gill's number while disagreeing
with Gill's non-responder direction is not replication — the agreement is coincidental.

**Worse, Gill's actual negative control is not in the dataset.** Verified: `GSE165176` contains **no
untreated sample at any day > 0** — the only non-OSKM samples are the six day-0 fibroblasts. A
Gill-style negative control cannot be constructed from this data.

---

## 3. What the non-responders actually are

Not senescent-arrested, not identity-changed — **inflammatory**:

| group | senescence markers | SASP | proliferation |
|---|---|---|---|
| day-0 fibroblasts | 1.04 | 0.46 | 0.31 |
| **non-responders** | 0.91 | **1.06 (2.3×)** | 0.57 |
| responders | 0.52 | 0.18 | 1.05 |

Non-responders are an **OSKM-stressed inflammatory state**; responders are proliferative and
non-inflammatory. That is a large, real biological difference — but it is a difference in
*cell state*, and calling it "28 years of age" is the interpretive leap the whole Stage 1.5 arc has
been trying to avoid. **This is a real biological condition, not an artefact to be cancelled**, so
subtracting it removes signal as well as noise.

---

## 4. Verdict and what to change

**Not ready to execute Steps 2–4.** Redefining `is_control` to the non-responder arm would bake an
identity + inflammation confound permanently into `y_age` — and unlike the current day-0 control, it
is a **treatment-dependent, moving baseline**, which is harder to reason about and harder to deploy
(the plan's own §5.5 concedes deployment needs internal controls).

### Keep

- **A1 and A3 fixes — unconditionally.** Stop pooling non-responders into the treatment arm; stop
  using a monotonic statistic for a dip. Re-run E1/E1b/D2 **split by responder status with a
  window contrast**. This alone may resolve the escalation, costs nothing, and touches no labels.
- **Step 1 (harden the finding)** — the LOO and window checks are worth having either way.
- The two new ground rules (§10 negative controls, §11 state the shape first). Both are correct and
  well-earned, independent of this plan's fate.

### Change

- **Do not redefine `is_control` yet.** The premise it rests on ("same identity change") is measured
  false.
- **Add a pre-registered identity-confound test to Step 2's decision**, before any label change:
  does the separation survive identity adjustment, and is the surviving part large enough to matter?
  Measured today: **−11.1 yr**, which against label noise 17.9 yr is **SNR 0.62** — *below* the
  current labels' 0.64. On this evidence the fix does not improve SNR at all.
- **Drop the "replication of Gill" framing** or restate it precisely: we reproduce a similar
  *magnitude* with a different clock against a different control, while disagreeing with Gill on the
  non-responder direction.
- **Treat the downstream annotations (MASTER_PLAN, STAGE_2–6) as provisional** — they inherit this
  premise and should not be relied on until it is settled.

### The honest open question

We now have three incompatible readings of the same cells — day-0 control (+11.4 yr effect),
non-responder control (−28.3), identity-adjusted (−11.1). **Deciding between them needs an
identity-independent anchor**, which transcriptomics alone cannot supply. Gill solved it with
multi-omic evidence (methylation). That is the real shape of the problem, and it is worth stating
plainly rather than picking whichever control gives the most publishable number.

---

# 5. THE A1/A3 RE-RUN — EXECUTED (2026-07-26)

`experiments/diag_e1_corrected.py` (+21 tests). Fixes **only** the two undisputed errors — stop
pooling non-responders (A1), use a window contrast instead of a monotonic Spearman (A3) — and
**keeps the current day-0 control**, so it is independent of the disputed control swap.

**Question:** once the pooling and the statistic are fixed, do the *current* labels show
rejuvenation in responders?

## 5.1 Answer: no — at any window

Responders only, ΔAge vs their own day-0 baseline:

| window | n | mean ΔAge | 95% CI | verdict |
|---|---|---|---|---|
| 7–9 d | 6 | **+20.6** | [−10.6, +51.8] | NO_EFFECT |
| **10–13 d** (pre-registered) | 6 | **+8.2** | [−20.1, +36.5] | **NO_EFFECT** |
| 13–15 d | 6 | +5.7 | [−18.4, +29.8] | NO_EFFECT |
| 15–21 d | 6 | +10.6 | [−15.6, +36.8] | NO_EFFECT |
| 21–29 d | 6 | +22.6 | [−4.1, +49.4] | NO_EFFECT |

Leave-one-donor-out: **STABLE** — a real null, not one donor carrying it. **Every point estimate is
positive** (mildly ageing), none significant.

**So the escalation was NOT purely a method artefact.** A1 and A3 are real errors, but fixing them
does not rescue the current labels.

## 5.2 The decomposition — this is the finding

```
separation  =  responder arm  −  non-responder arm
   −28.3     =     +8.2       −      +36.5
```

> **100% of the "−28.3 yr rejuvenation" comes from the CONTROL arm rising. 0% comes from the
> treatment arm falling.**

The responders never get younger — they drift slightly *older* (+8.2). The entire effect is
non-responders reading +36.5 yr. This is arithmetic, not a hypothesis test, so it does not depend on
sample size or on the identity adjustment in §2.3.

Non-responders, for contrast, are the only arm with a significant signal — and it is **AGEING** at
every window from day 10 on (+36.5, +44.6, +42.4, +43.7).

**This is the strongest single argument against the proposed fix.** Redefining `is_control` to the
non-responder arm would define ΔAge as "how much *less* the control arm inflates" — a quantity whose
entire dynamic range is supplied by an artefact in the reference, not by the treatment.

## 5.3 Power — stated honestly, because it limits the claim

| true effect | power at n=6 with the observed spread (sd 27.0) |
|---|---|
| −10 yr | 9% |
| −20 yr | 25% |
| **−30 yr (Gill-scale)** | **56%** |
| −40 yr | 83% |

**56% is a coin flip**, so §5.1 is *not* proof that no rejuvenation exists — it is an underpowered
null. What makes it informative is the **direction**: the point estimate is **+8.2**, not
"negative but short of significance". A true −30 yr effect would have to be masked by ~38 yr of
noise in the wrong direction.

**§5.2 is not subject to this caveat.** The decomposition holds regardless of power.

## 5.4 What this changes

| | before this run | after |
|---|---|---|
| "A1/A3 explain the escalation" | plausible | ❌ **refuted** — fixed, and the null persists |
| "current day-0 labels are fine" | possible | ❌ **no rejuvenation signal at any window** |
| "non-responder control recovers the effect" | proposed | ❌ **the effect *is* the control arm** |

**All three candidate label definitions now fail**, each for its own reason. That is a real result,
and it is the honest end of the transcriptomic-only route on this dataset: Gill needed **methylation**
to establish rejuvenation, and our own §2.4 check found Gill saying existing transcription clocks
*"failed to accurately predict the age of our negative control samples."*

**Recommended next step is no longer a label change.** It is to decide, with the user, between:

1. **An identity-independent anchor** (methylation or a second modality) — the only route that can
   actually settle it, and the one Gill took. Links to Stage 6.
2. **Scope the claim to what is supported** — the fate/safety head remains strong (PR-AUC 0.99) and
   is untouched by any of this; the quantitative rejuvenation claim is not currently supportable on
   this data.

Stage 2's intervention (k≈3 reference cells/donor) remains worth doing under either, because it
attacks the n=1 baseline directly — but its *justification* should now be baseline replication, not
"correcting a known biological offset."

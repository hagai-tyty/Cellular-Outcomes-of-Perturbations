# STAGE 23.2 V5 — ADDENDUM 1: the power CURVE

**Amends:** `STAGE_23_2_ROLE_A_CONFIRMATION_V5.md`, canonical LF sha256
`b576c95aa1036470db1d041e2540f1ca66f4aee1e60e0f088edd34d1b4ba04fb`.
V5 itself is **not edited**. It is frozen at that digest and stays there.

**Why an addendum rather than an edit.** The curve was proposed while the primary arm's 2000-draw
permutation run was still in flight. Editing V5 changed its digest, and the resume cache — which
stamps the protocol it was written under — correctly refused to merge against a protocol that had
moved underneath it. That guard was right and the sequencing was wrong: a frozen protocol must not
be edited while a run against it is executing. Removing exactly this block restored V5's digest to
`b576c95a`, which is itself the proof that this block was the only change. Amendments now go in
dated addenda.

**Status of this addendum:** REPORTING ONLY. It adds no gate, removes no gate, and changes no
gate's threshold, alternative, or decision rule. Gate 18.3 remains exactly as V5 §9.3 defines it.

---

## A1.1 The power curve — reported, never gating

§9.3's gate is evaluated at the pre-registered alternative `oracle AUC = 0.66` and **at no other
value**. That anchor is inherited: 23.2E chose it in advance to match R1's observed R3 ROC-AUC of
`0.6628`, and V5 does not move it.

But a single point answers only "is this design powered for an effect of exactly that size?" The
more informative and equally pre-registrable question is **"for effects of what size IS this design
adequately powered?"** 23.2E already registered `AUC 0.70` as its secondary planning sensitivity,
so extending 0.66 into a curve continues prior registration rather than inventing a new analysis.

```text
  curve points     0.66 (the GATE point), 0.70, 0.75
  null             SHARED across the curve -- generated at beta = 0, independent of the target
  alternative      recomputed per point
  status           REPORTED. Gate 18.3 is evaluated at 0.66 only.
```

**What the curve may not be used for.** §9.4 stands unchanged and is restated because this is
exactly where it would be tempting to cheat:

```text
  evaluating gate 18.3 at any point other than 0.66
  reporting SUPPORTED because power at some larger alternative reaches 0.80
  selecting a curve point after seeing which one clears the threshold
  inferring the true effect size from the observed result and reading power off the curve there
    -- this is the observed-power fallacy and it is forbidden outright
```

The curve's legitimate use is to convert *"the design is underpowered"* into the quantified and
more useful *"the design has >= 0.80 power against effects of AUC >= X; the pre-registered
alternative was 0.66, at which power is 0.64."* That statement is reportable. It does not open
Stage 24.

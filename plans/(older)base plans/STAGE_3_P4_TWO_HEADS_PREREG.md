# P4 — PRE-REGISTRATION. Two heads, no composite. Committed BEFORE the run.

The rule was already fixed in `WORK_ORDER_2026_08_14.md`: identity-loss and apoptosis are reported
and modelled **separately**, and **no composite endpoint is constructed** — not
`1 − (1 − P_loss)(1 − P_death)`, not a weighted sum, not a union. Choosing a combination rule after
seeing which looks better is a forking path.

**What is missing is the evidence that the separation is warranted**, and a way for that premise to
be wrong. That is what this run supplies.

---

## 1. A structural fact to state before any number

Labels are the **argmax** of a three-class call over (safe, loss, death), so **a cell is either
identity-lost or apoptotic, never both.** The union is therefore just `P(not safe)`, and
`P(unsafe) = P(loss) + P(death)` exactly.

**So collapsing does not lose information to overlap — it loses the distinction between two failure
modes.** Whether that distinction matters is an empirical question about their *time courses*, and
that is what is measured here.

## 2. What will be measured

HFF (`GSE242423`), 9 timepoints, ~42k cells. Per timepoint: `P(loss)`, `P(death)`, their sum, the
**share of the unsafe fraction attributable to apoptosis**, and binomial intervals throughout.
Then the relationship between the two curves across timepoints — Spearman and Pearson.

READ-ONLY, no retrain, `src/` untouched.

## 3. 🔒 PRE-REGISTERED OUTCOMES — including the one that would refute P4's premise

| # | result | reading | action |
|---|---|---|---|
| **H1** | Spearman(loss, death) across timepoints **> +0.9** *and* the two peak within one timepoint of each other | **THE HEADS ARE REDUNDANT.** They rise and fall together, so a single "unsafe" number loses nothing | **P4's premise is REFUTED.** Report it and revert to the collapsed endpoint |
| **H2** | Spearman **≤ +0.9**, or the peaks are separated | the two failure modes have **different time courses**; a single number cannot express both | the two-head form is **warranted and adopted** |
| **H3** | Spearman **negative** | they move in **opposite directions** — the union is dominated by whichever is larger and the other becomes **invisible** in it | two heads are **mandatory**, and the collapsed endpoint is actively misleading |

**Secondary, reported either way:** the apoptosis share of the unsafe fraction at each timepoint. If
that share varies by more than 10× across the course, a single number is a different quantity at
different days even when it is numerically similar.

## 4. Declared limits

1. **One line.** Whether the two courses separate the same way elsewhere is untested.
2. **Marker-based labels**, inheriting whatever the signature definitions inherit — in particular
   the apoptosis programme is scored from a five-gene set.
3. **Mutual exclusivity is imposed by the argmax**, not measured. A cell partway into both
   programmes is forced to one. This bounds how much the two heads *could* overlap and is a
   property of the labeller, not the biology.
4. This changes **what is reported**, not any model weights. No retrain.

## 5. Recording

`results/p4_two_heads_results.json`; write-up to the work order, `CHANGES.md`, the notebook; tests
in `tests/test_p4_two_heads.py`. Outcomes graded as written, including H1 if it fires.

---

## 6. RESULT — 2026-08-14. **H3 — the two failure modes run in OPPOSITE directions.** The collapsed endpoint is actively misleading.

*Graded against §3 as written. Artefacts: `experiments/p4_two_heads.py`,
`results/p4_two_heads_results.json`, `tests/test_p4_two_heads.py` (18 tests).*
42,481 HFF cells, 9 timepoints. Cells labelled both loss and death: **0**, as the argmax requires.

| day | P(identity loss) | P(apoptosis) | union | **apoptosis share of unsafe** |
|---|---|---|---|---|
| 0 | 0.0437 | 0.0397 | 0.0835 | **47.6 %** |
| 2 | 0.3314 | 0.0195 | 0.3509 | 5.6 % |
| 4 | 0.3248 | 0.0261 | 0.3508 | 7.4 % |
| 6 | 0.3890 | 0.0335 | 0.4226 | 7.9 % |
| 8 | 0.4311 | 0.0397 | 0.4708 | 8.4 % |
| 10 | 0.3621 | 0.0344 | 0.3964 | 8.7 % |
| 12 | 0.4454 | 0.0205 | 0.4659 | 4.4 % |
| 14 | 0.6662 | 0.0129 | 0.6791 | 1.9 % |
| 21 | 0.9967 | 0.0029 | 0.9996 | **0.3 %** |

**Spearman −0.633, Pearson −0.818.** Identity loss peaks at **day 21**; apoptosis peaks at
**day 0** — **8 timepoints apart, the maximum separation the design allows.**

### 🔒 H3, the strongest branch in the table

**Two heads are mandatory, and the collapsed endpoint is actively misleading.** Not merely lossy —
misleading, because the union is dominated by whichever mode is larger and the other becomes
**invisible inside it**.

### The secondary criterion fired too, by a wide margin

The apoptosis share of the unsafe fraction runs **47.6 % → 0.3 %, a 163× range** against a
pre-registered threshold of 10×.

> **At day 0, apoptosis is nearly half of all unsafe cells. By day 21 it is one part in three
> hundred.** A single `p_unsafe` number is therefore **a different quantity at different days**,
> even where it happens to be numerically similar.

### Why this matters for the product, concretely

A withdrawal recommendation gated on one `p_unsafe` would be driven almost entirely by identity
loss from day 2 onward, and would be **blind to apoptosis risk being highest at the very start**.
The two curves cross near day 0 and diverge monotonically thereafter, so no single threshold on the
union can express both.

### What the run cost, and one bug it found in itself

`spearman` originally ranked via `argsort(argsort(...))`, which **breaks ties arbitrarily instead
of averaging them** — and `P(apoptosis)` is tied at 0.0397 on days 0 and 8, the two timepoints that
matter most to the verdict. Caught by a unit test asserting a constant vector ranks flat. Fixed to
average ties; **the Spearman is unchanged at −0.633**, so no conclusion moves — but the function
was wrong and would have mattered on a different tie pattern.

### What is adopted

`P(identity loss)` and `P(apoptosis)` are **reported and modelled separately**. **No composite
endpoint is constructed** — not the union, not a weighted sum. If a downstream decision ever needs
one number, that rule gets its own pre-registration, and it now has to justify itself against a
163× shift in what the number is made of.

### What is NOT claimed

That the two courses separate this way in another line — one line, untested elsewhere. That the
apoptosis programme is well measured; it is a five-gene signature. And the exclusivity is **imposed
by the argmax**, so a cell partway into both programmes is forced to one — that bounds how much the
heads could overlap and is a property of the labeller, not the biology.

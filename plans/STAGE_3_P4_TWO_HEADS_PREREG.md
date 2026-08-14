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

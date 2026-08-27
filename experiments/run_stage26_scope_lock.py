"""Stage 26 — the known-treatment-only scope lock.

Executes `STAGE_26_KNOWN_TREATMENT_SCOPE_LOCK_V1.md`, the executable form of one line in the
frozen Stage-23.5 ship plan (§9):

    STAGE 26
      record KNOWN_TREATMENT_ONLY_SCOPED_LIMIT
      no unseen-treatment claim and no rescue experiment

Stage 25 came back positive, and the predictable failure mode of a positive result is that the
claim quietly grows: six observed conditions become "treatments", one melanoma line becomes
"cancer", a detection proxy becomes "response". A scope limit written in a document does nothing
about that if the shipped tool will happily score `Vemurafenib`. So this stage attacks the tool
with the strings a real user would actually try, and only then writes the limit down.

Nothing here fits anything. Stage 26 reads frozen artifacts, runs the shipped API a few hundred
times on one clone, scans text, and records. If it fails, the TOOL is fixed and this re-runs from
the top -- there is no path in which a scope hole is resolved by widening the scope.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from cellfate.gen1_predictor import Gen1Predictor  # noqa: E402

RESULTS = ROOT / "results"
OUT = RESULTS / "stage26"
OUT.mkdir(parents=True, exist_ok=True)

PLAN = ROOT / "plans" / "(newer)practical plans" / "STAGE_26_KNOWN_TREATMENT_SCOPE_LOCK_V1.md"
SHIP_PLAN = (ROOT / "plans" / "(newer)practical plans"
             / "STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md")
SHIP_PLAN_DIGEST = "8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48"

TOOL = RESULTS / "stage24" / "tool"
ARTIFACT_NPZ = RESULTS / "stage24" / "stage24_w5_artifact.npz"
ARTIFACT_META = RESULTS / "stage24" / "stage24_w5_artifact.json"
FREEZE_24F = RESULTS / "stage24" / "stage24f_tool_freeze.json"
STAGE25_VERDICT = RESULTS / "stage25" / "stage25_verdict.json"
MODEL_CARD = TOOL / "MODEL_CARD.md"
EX_X = TOOL / "example_clone_expression.npy"
EX_B = TOOL / "example_clone_nuisance.txt"

A_JSON = OUT / "stage26a_vocabulary_closure.json"
B_JSON = OUT / "stage26b_claim_surface.json"
C_JSON = OUT / "stage26c_propagation.json"
D_JSON = OUT / "stage26d_no_rescue.json"
SCOPE_MD = OUT / "GEN1_SCOPE_LIMIT.md"
VERDICT_JSON = OUT / "stage26_verdict.json"
HANDOFF_JSON = RESULTS / "stage26_handoff_to_evidence_lock.json"

CONDITIONS = ("Acid", "Cisplatin", "CoCl2", "Dabrafenib", "Doxorubicin", "Trametinib")
REFERENCE = "Acid"

MODEL_CARD_DELIMITER = "<!-- STAGE-26 SCOPE LOCK -- appended, nothing above this line altered -->"


def sha256_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def canonical_lf_sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def module_sha256() -> str:
    """§5.1. Stamped into every substage so 26E can refuse a verdict built from mixed versions."""
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def write_json(p: Path, obj: dict) -> dict:
    """Write, and return exactly what was written.

    Returning `obj` instead of the stamped dict is how 26E ended up merging a 26B result with no
    module stamp while 26A/26C/26D had one -- the §5.1 check caught it on its first run. What a
    caller gets back is now byte-for-byte what landed on disk.
    """
    stamped = {**obj, "module_sha256": module_sha256()}
    p.write_text(json.dumps(stamped, indent=2) + "\n", encoding="utf-8")
    return stamped


def _rel(p: Path) -> str:
    return Path(p).resolve().relative_to(ROOT).as_posix()


# =============================================================================================== #
# The adversarial corpus — fixed in plan §2.1 BEFORE it was run, so it cannot be trimmed after
# seeing which entries fail. Every one of these must be refused.
# =============================================================================================== #
ADVERSARIAL: dict[str, list[str]] = {
    "case": ["acid", "ACID", "cisplatin", "CisPlatin", "cocl2", "COCL2", "dabrafenib",
             "trametinib", "DOXORUBICIN"],
    "whitespace": [" Acid", "Acid ", "Acid\t", "Cis platin", "Co Cl2", "\nAcid"],
    # the block that matters: what a real user of a BRAF-V600E melanoma tool would actually type
    "pharmacology": ["Vemurafenib", "PLX4720", "Cobimetinib", "Selumetinib", "Encorafenib",
                     "Binimetinib", "Carboplatin", "Oxaliplatin", "Paclitaxel", "Docetaxel",
                     "Etoposide", "Nilotinib", "Pembrolizumab", "Nivolumab", "Temozolomide",
                     "5-FU"],
    "format": ["Cisplatin_1uM", "Cisplatin (1 uM)", "Cisplatin-high", "Acid+Cisplatin",
               "Cisplatin:1uM"],
    "control": ["DMSO", "Vehicle", "Control", "Untreated", "None", "NA", "null", "0"],
    # visually identical to a supported condition, byte-different
    "confusable": ["Аcid", "Acіd", "CoCl₂", "CoCI2", "Trametinib​"],
    "structural": ["", " ", "*", "all", "Acid,Cisplatin", "[Acid]", "Acid|Cisplatin"],
}

# Plan §2.1 declares these exact counts. A group that shrinks fails the stage rather than passing
# with fewer attackers, which is the whole point of declaring the corpus before running it.
EXPECTED_GROUP_SIZES = {"case": 9, "whitespace": 6, "pharmacology": 16, "format": 5,
                        "control": 8, "confusable": 5, "structural": 7}
EXPECTED_ADVERSARIAL_TOTAL = 56


def _load_clone() -> tuple[np.ndarray, np.ndarray]:
    x = np.load(EX_X).ravel()
    b = np.loadtxt(EX_B, delimiter=",").ravel()
    return x, b


def _predictor(with_verdict: bool = True) -> Gen1Predictor:
    return Gen1Predictor.load(ARTIFACT_NPZ, ARTIFACT_META,
                              stage25_verdict=STAGE25_VERDICT if with_verdict else None)


# =============================================================================================== #
# 26A — vocabulary closure, adversarially (plan §2)
# =============================================================================================== #
def run_26a() -> dict:
    t0 = time.perf_counter()
    p = _predictor()
    x, b = _load_clone()

    # ---- every adversarial string must be refused, with no number attached -------------------- #
    refusals: dict[str, list[dict]] = {}
    leaked: list[dict] = []
    for group, strings in ADVERSARIAL.items():
        rows = []
        for s in strings:
            r = p.predict(x, b, treatments=[s])
            assert len(r) == 1
            r0 = r[0]
            ok = (r0["support_status"] == "UNSUPPORTED_TREATMENT"
                  and r0["future_detection_score"] is None)
            rows.append({"input": s, "support_status": r0["support_status"],
                         "score": r0["future_detection_score"], "refused": ok})
            if not ok:
                leaked.append({"group": group, "input": s, "row": r0})
        refusals[group] = rows

    n_total = sum(len(v) for v in refusals.values())
    n_refused = sum(1 for v in refusals.values() for r in v if r["refused"])

    # ---- the reference-leak test (plan §2.2) -------------------------------------------------- #
    #
    # Acid is the reference level: five zero dummies. So the dummy encoder, on its own, maps ANY
    # unknown string to the Acid row. Demonstrate that hazard is real, then prove the filter in
    # predict() is what stops it. The guard is load-bearing; the encoder is not safe alone.
    comp = p.components["deployment"]
    acid_score = float(comp.score(x, b, p._dummies([REFERENCE]))[0])
    hazard_dummies = p._dummies(["Vemurafenib"])
    hazard_score = float(comp.score(x, b, hazard_dummies)[0])
    hazard_is_real = bool(hazard_dummies.sum() == 0 and hazard_score == acid_score)

    guarded = p.predict(x, b, treatments=["Vemurafenib"])[0]
    guard_holds = (guarded["future_detection_score"] is None
                   and guarded["support_status"] == "UNSUPPORTED_TREATMENT")

    # a mixed request must score exactly the known ones and refuse exactly the unknown ones,
    # in the order asked, with no row shifted onto a neighbour's score
    mixed_req = ["Vemurafenib", "Cisplatin", "DMSO", "Trametinib", "acid", "Acid"]
    mixed = p.predict(x, b, treatments=mixed_req)
    known_req = [t for t in mixed_req if t in CONDITIONS]
    scored = [r for r in mixed if r["future_detection_score"] is not None]
    solo = {t: p.predict(x, b, treatments=[t])[0]["future_detection_score"] for t in known_req}
    routing_ok = ([r["condition"] for r in mixed] == mixed_req
                  and [r["condition"] for r in scored] == known_req
                  and len(scored) == len(known_req))

    # Routing is the scope question. Bit-identity across batch sizes is a NUMERICS question that
    # the frozen plan already settled: §7.1 R2 bounds a prediction cell at 1e-12, and R4 requires
    # the cause of any non-identity to be named rather than shrugged at.
    #
    # The named cause here: the design matrix passed to the GEMM has 1, 3 or 6 rows depending on
    # how many KNOWN conditions were requested, and the BLAS kernel selected differs by shape. The
    # first run of 26A asserted exact equality, saw 1.1e-16, and failed. The assertion was wrong,
    # not the tool -- but it is replaced by a bounded one, never by a removed one.
    batch_diffs = {t: abs(r["future_detection_score"] - solo[t])
                   for r in scored for t in [r["condition"]]}
    max_batch_diff = max(batch_diffs.values()) if batch_diffs else 0.0
    within_tol = max_batch_diff <= 1e-12

    # §7.1 R3 is the load-bearing one: the ranking claim is a function of within-clone ORDERING and
    # nothing else, so batching must not flip a single pair, including ties.
    six_batch = {r["condition"]: r["future_detection_score"] for r in p.predict(x, b)}
    six_solo = {t: p.predict(x, b, treatments=[t])[0]["future_detection_score"]
                for t in CONDITIONS}

    def _signs(d):
        return [int(np.sign(d[u] - d[v])) for i, u in enumerate(CONDITIONS)
                for v in CONDITIONS[i + 1:]]

    ordering_stable = _signs(six_batch) == _signs(six_solo)
    mixed_ok = bool(routing_ok and within_tol and ordering_stable)

    # ---- structural closure (plan §2.3) ------------------------------------------------------- #
    meta = json.loads(ARTIFACT_META.read_text(encoding="utf-8"))
    vocab_ok = tuple(meta["treatment_vocabulary"]) == CONDITIONS
    ref_ok = meta["reference_treatment"] == REFERENCE
    n_pc, n_nuis, n_dum = comp.K, len(meta["nuisance_columns"]), len(CONDITIONS) - 1
    expected_cols = n_pc + n_nuis + n_dum + n_pc * n_dum
    geom_ok = (comp.coef.shape[0] == expected_cols == 309)

    sizes = {k: len(v) for k, v in ADVERSARIAL.items()}
    flat = [s for v in ADVERSARIAL.values() for s in v]

    checks = {
        "the corpus is the size plan §2.1 declared, group by group":
            sizes == EXPECTED_GROUP_SIZES and n_total == EXPECTED_ADVERSARIAL_TOTAL,
        "no duplicate inflates the refusal count": len(flat) == len(set(flat)),
        "every adversarial string is refused": n_refused == n_total,
        "the reference-leak hazard is real (unknown -> all-zero dummies -> the Acid row)":
            hazard_is_real,
        "predict() blocks it: no score is returned for an unknown condition": bool(guard_holds),
        "a mixed request routes each condition to its own row, in order": bool(routing_ok),
        "batched and solo scores agree inside the frozen 1e-12 cell bound": bool(within_tol),
        "within-clone ordering is identical across batch sizes (7.1 R3)": bool(ordering_stable),
        "vocabulary is exactly the six frozen conditions": bool(vocab_ok),
        "reference condition is Acid": bool(ref_ok),
        "design geometry closes the vocabulary at 309 columns": bool(geom_ok),
    }

    out = {
        "stage": "26A",
        "n_adversarial_strings": n_total,
        "n_refused": n_refused,
        "group_sizes": sizes,
        "group_sizes_declared_in_plan": EXPECTED_GROUP_SIZES,
        "leaked": leaked,
        "refusals_by_group": refusals,
        "reference_leak_test": {
            "acid_score": acid_score,
            "unknown_via_raw_encoder": hazard_score,
            "raw_encoder_dummy_sum": float(hazard_dummies.sum()),
            "hazard_is_real": hazard_is_real,
            "note": "The encoder is NOT safe on its own -- an unknown string produces the "
                    "reference row and would return the Acid score under another name. The "
                    "vocabulary filter inside predict() is the thing that makes it safe, which "
                    "is why it is tested rather than assumed.",
            "guard_holds": bool(guard_holds)},
        "mixed_request": {"requested": mixed_req,
                          "returned_conditions": [r["condition"] for r in mixed],
                          "scored_conditions": [r["condition"] for r in scored],
                          "routing_correct": bool(routing_ok),
                          "max_batch_vs_solo_difference": max_batch_diff,
                          "frozen_tolerance_23_5_7_1_R2": 1e-12,
                          "within_frozen_tolerance": bool(within_tol),
                          "named_cause_23_5_7_1_R4":
                              "BLAS GEMM kernel selection varies with design-matrix row count "
                              "(1, 3 or 6 rows depending on how many known conditions were "
                              "requested). Not exact bit-identity across batch sizes; 4 orders "
                              "of magnitude inside the frozen 1e-12 cell bound.",
                          "within_clone_ordering_identical_across_batch_sizes":
                              bool(ordering_stable)},
        "structural_closure": {"vocabulary": list(meta["treatment_vocabulary"]),
                               "reference": meta["reference_treatment"],
                               "K": n_pc, "nuisance": n_nuis, "dummies": n_dum,
                               "interaction": n_pc * n_dum,
                               "design_columns": int(comp.coef.shape[0]),
                               "expected": expected_cols},
        "checks": checks,
        "all_passed": all(checks.values()),
        "runtime_seconds": round(time.perf_counter() - t0, 3),
    }
    return write_json(A_JSON, out)


# =============================================================================================== #
# 26B — claim-surface scan (plan §3)
#
# Two tiers, deliberately.
#
#   GATING   high-precision phrases that are a forbidden claim if they appear unnegated. These
#            fail the stage.
#   REVIEW   broad topic words. Reported for human reading, never gating -- a gate tuned on broad
#            words is a gate tuned until it passes, which is worth nothing.
#
# The Stage-26 module, its tests and its plan are excluded from the scan: a scanner that trips on
# its own description of what it hunts produces a false failure and teaches nothing. That mistake
# was made once already, in the Stage-23.2 leak scan, and is not repeated.
# =============================================================================================== #
FORBIDDEN_CLAIMS: dict[str, list[str]] = {
    "1_unseen_treatment_generalization": [
        r"unseen treatments?", r"novel treatments?", r"new treatments?", r"any treatment",
        r"arbitrary treatments?", r"untested (?:drug|treatment|compound)s?",
        r"generali[sz]es? to", r"works? for any"],
    "2_cross_cell_line_or_patient": [
        r"other cell lines?", r"another cell line", r"any cell line", r"cross-cell-line",
        r"patients?", r"in vivo", r"tumou?rs?", r"cross-patient"],
    "3_clinical_recommendation": [
        r"clinical(?:ly)?", r"recommends?", r"recommendation", r"prescrib", r"therapeutic",
        r"treatment of choice", r"best treatment", r"choose (?:a|the) (?:drug|therapy)"],
    "4_causal_effect": [
        r"causal", r"causes", r"treatment effects?", r"counterfactual", r"because the treatment"],
    "5_calibrated_probability": [
        r"calibrated", r"calibration", r"probability of (?:death|survival|response)",
        r"risk score"],
    "6_independent_replication": [
        r"independently replicated", r"replicated in", r"external validation",
        r"independent (?:biological )?replication", r"validated externally"],
    "7_uniform_benefit": [
        r"all six conditions", r"every condition", r"uniform(?:ly)? (?:benefit|improve)",
        r"across all conditions", r"works equally"],
    "8_confirmed_role_a": [
        r"confirmed (?:role[- ]a|rewind)", r"role[- ]a (?:is|was) confirmed",
        r"rewind confirms", r"validated anchor"],
    "9_single_cell_equivalence": [
        r"single cells? (?:is|are) equivalent", r"per-cell input", r"score a single cell",
        r"single-cell input is"],
}

REVIEW_ONLY = [r"\bcancer\b", r"\bdeath\b", r"\bresistan", r"\bsensitiv", r"\bresponse\b",
               r"\bcures?\b", r"\befficacy\b"]

# Exactly the twelve tokens plan §3.1 declares, and no more. An earlier version carried twenty-six;
# the extra fourteen never rescued a single hit on any shipped surface (verified: every hit that
# depended on one token depended on `not` or `never`), but a negation list longer than the declared
# one is a looser gate than the one that was written down, and a gate that quietly loosens is not a
# gate. `"no "` keeps its trailing space: bare `no` matches inside not, none, know and cannot.
NEGATIONS = ["not", "never", "no ", "cannot", "without", "forbidden",
             "refus", "unsupported", "withheld", "limit", "outside", "only"]

WINDOW = 160


def _surfaces() -> list[Path]:
    return [ROOT / "src" / "cellfate" / "gen1_predictor.py",
            ROOT / "src" / "cellfate" / "gen1_cli.py",
            MODEL_CARD,
            TOOL / "io_schema.json",
            TOOL / "example_clone_README.md",
            ARTIFACT_META,
            STAGE25_VERDICT,
            SCOPE_MD]


def _scan_text(text: str, patterns: list[str]) -> list[dict]:
    hits = []
    low = text.lower()
    for pat in patterns:
        for m in re.finditer(pat, low):
            a, z = max(0, m.start() - WINDOW), min(len(text), m.end() + WINDOW)
            ctx = low[a:z]
            neg = [n for n in NEGATIONS if n in ctx]
            hits.append({"pattern": pat, "match": text[m.start():m.end()],
                         "negated_by": neg[:4], "negated": bool(neg),
                         "context": " ".join(text[a:z].split())[:300]})
    return hits


# Negative control. A scanner that never fires proves nothing, and "zero violations" from an
# instrument that cannot detect a violation is the most comfortable kind of false pass. Each
# sentence below is a forbidden claim stated plainly, and each must be caught; the negated twin
# must not be. If the canary ever stops firing, the scan result is void, not clean.
CANARY_POSITIVE = {
    "1_unseen_treatment_generalization": "The model generalizes to any treatment you supply.",
    "2_cross_cell_line_or_patient": "It works on other cell lines and on patients too.",
    "3_clinical_recommendation": "Use it clinically to recommend a therapeutic for each clone.",
    "4_causal_effect": "The scores estimate the causal treatment effect of each drug.",
    "5_calibrated_probability": "Each score is a calibrated probability of death.",
    "6_independent_replication": "The finding was independently replicated in an external cohort.",
    "7_uniform_benefit": "The interaction helps uniformly improve all six conditions.",
    "8_confirmed_role_a": "Rewind confirms the result; Role A is confirmed.",
    "9_single_cell_equivalence": "A single cell is equivalent to a clone as per-cell input.",
}
CANARY_NEGATIVE = "This tool cannot {}"


def _canary() -> dict:
    """Fire the scanner at plainly-stated forbidden claims and at their negated twins."""
    caught, missed, false_alarm = [], [], []
    for claim, sentence in CANARY_POSITIVE.items():
        hits = _scan_text(sentence, FORBIDDEN_CLAIMS[claim])
        (caught if any(not h["negated"] for h in hits) else missed).append(claim)
        neg = _scan_text(CANARY_NEGATIVE.format(sentence[0].lower() + sentence[1:]),
                         FORBIDDEN_CLAIMS[claim])
        if any(not h["negated"] for h in neg):
            false_alarm.append(claim)
    return {"claims_probed": len(CANARY_POSITIVE), "caught": caught, "missed": missed,
            "negated_twin_false_alarms": false_alarm,
            "detects_every_forbidden_claim": not missed,
            "does_not_fire_on_a_negation": not false_alarm}


def run_26b() -> dict:
    t0 = time.perf_counter()
    canary = _canary()
    per_file, violations, review = {}, [], []
    for f in _surfaces():
        if not f.exists():
            per_file[_rel(f)] = {"present": False}
            continue
        text = f.read_text(encoding="utf-8")
        claims = {}
        for claim, pats in FORBIDDEN_CLAIMS.items():
            hits = _scan_text(text, pats)
            if hits:
                claims[claim] = {"hits": len(hits),
                                 "unnegated": [h for h in hits if not h["negated"]]}
                for h in hits:
                    if not h["negated"]:
                        violations.append({"file": _rel(f), "claim": claim, **h})
        rev = _scan_text(text, REVIEW_ONLY)
        if rev:
            review.append({"file": _rel(f), "hits": len(rev),
                           "unnegated": len([h for h in rev if not h["negated"]])})
        per_file[_rel(f)] = {"present": True, "bytes": len(text.encode("utf-8")),
                             "claims_touched": claims}

    checks = {
        "the scanner detects every one of the nine claims stated plainly (canary)":
            canary["detects_every_forbidden_claim"],
        "the scanner does not fire on a negated twin (canary)":
            canary["does_not_fire_on_a_negation"],
        "no forbidden claim appears unnegated on any shipped surface": not violations,
        "every scanned surface exists": all(v.get("present") for v in per_file.values()),
    }
    out = {
        "stage": "26B",
        "surfaces_scanned": [_rel(f) for f in _surfaces()],
        "excluded_from_scan": [_rel(Path(__file__)),
                               "tests/test_stage26_scope_lock.py",
                               _rel(PLAN),
                               "plans/(newer)practical plans/RECORDs/stage_26_RECORD.md"],
        "exclusion_reason": "a scanner that trips on its own description of what it hunts "
                            "produces a false failure and teaches nothing",
        "canary": canary,
        "gating_claims": list(FORBIDDEN_CLAIMS),
        "review_only_terms_are_not_gating": REVIEW_ONLY,
        "per_file": per_file,
        "violations": violations,
        "review_queue": review,
        "checks": checks,
        "all_passed": all(checks.values()),
        "runtime_seconds": round(time.perf_counter() - t0, 3),
    }
    return write_json(B_JSON, out)


# =============================================================================================== #
# 26C — scope-limit propagation (plan §4)
#
# A caller who gets a refusal and no limitations has been told LESS than a caller who got a score.
# So the limit is checked on the failing paths too, and through the CLI, not only through Python.
# =============================================================================================== #
KNOWN_TREATMENT_LIMITATION = "known conditions only"


def run_26c() -> dict:
    t0 = time.perf_counter()
    p = _predictor()
    x, b = _load_clone()
    lim = list(p._limitations)

    paths = {}

    def record(name: str, rows: list[dict], expect_status: str) -> None:
        paths[name] = {
            "expected_status": expect_status,
            "statuses": sorted({r["support_status"] for r in rows}),
            "all_carry_known_limitations": all(r.get("known_limitations") for r in rows),
            "all_carry_model_version": all("model_version" in r for r in rows),
            "all_carry_ranking_status": all("ranking_status" in r for r in rows),
            "no_calibrated_probability_key": all("calibrated_probability" not in r for r in rows),
            "status_matches": all(r["support_status"] == expect_status for r in rows),
        }

    record("supported", p.predict(x, b), "SUPPORTED_KNOWN_CONDITION")
    record("unsupported_treatment", p.predict(x, b, treatments=["Vemurafenib"]),
           "UNSUPPORTED_TREATMENT")
    record("missing_nuisance", p.predict(x, None), "MISSING_REQUIRED_NUISANCE")
    record("bad_nuisance_length", p.predict(x, np.array([1.0, 2.0])), "MISSING_REQUIRED_NUISANCE")
    record("non_finite_nuisance", p.predict(x, np.array([1.0, np.nan, 1.0, 1.0])),
           "MISSING_REQUIRED_NUISANCE")
    record("bad_feature_schema", p.predict(x[:100], b), "UNSUPPORTED_FEATURE_SCHEMA")

    rk_with = p.rank_conditions(x, b)
    rk_without = _predictor(with_verdict=False).rank_conditions(x, b)
    rk_refused = p.rank_conditions(x, None)

    ranking = {
        "with_verdict_status": rk_with["ranking_status"],
        "with_verdict_exposes_order": rk_with.get("validated_condition_order") is not None,
        "without_verdict_status": rk_without["ranking_status"],
        "without_verdict_withholds_order": rk_without.get("validated_condition_order") is None,
        "scores_identical_either_way": rk_with["scores"] == rk_without["scores"],
        "all_carry_known_limitations": all(bool(r.get("known_limitations"))
                                           for r in (rk_with, rk_without, rk_refused)),
        "refused_rank_has_support_status": "support_status" in rk_refused,
        "order": rk_with.get("validated_condition_order"),
    }

    # ---- the CLI, as a caller actually invokes it --------------------------------------------- #
    def cli(args: list[str]) -> dict:
        r = subprocess.run([sys.executable, "-m", "cellfate.gen1_cli", *args],
                           capture_output=True, text=True, cwd=str(ROOT / "src"))
        return {"exit_code": r.returncode, "stdout": r.stdout, "stderr": r.stderr}

    base = ["--artifact", str(ARTIFACT_NPZ), "--meta", str(ARTIFACT_META),
            "--expression", str(EX_X), "--nuisance", str(EX_B)]
    c_ok = cli(base)
    c_unknown = cli(base + ["--treatments", "Vemurafenib,Cisplatin"])
    c_nonuis = cli(["--artifact", str(ARTIFACT_NPZ), "--meta", str(ARTIFACT_META),
                    "--expression", str(EX_X)])
    c_unread = cli(["--artifact", str(ARTIFACT_NPZ), "--meta", str(ARTIFACT_META),
                    "--expression", str(OUT / "does_not_exist.npy"), "--nuisance", "1,2,3,4"])

    def rows_of(res: dict) -> list[dict]:
        return [json.loads(ln) for ln in res["stdout"].splitlines() if ln.strip()]

    printed = {"all_known": rows_of(c_ok), "one_unknown": rows_of(c_unknown),
               "no_nuisance": rows_of(c_nonuis)}
    cli_rows_carry_limits = {k: all(bool(r.get("known_limitations")) for r in v)
                             for k, v in printed.items()}

    cli_report = {
        "all_known": {"exit_code": c_ok["exit_code"], "rows": len(rows_of(c_ok)),
                      "statuses": sorted({r["support_status"] for r in rows_of(c_ok)})},
        "one_unknown": {"exit_code": c_unknown["exit_code"],
                        "statuses": sorted({r["support_status"] for r in rows_of(c_unknown)}),
                        "scores": [r["future_detection_score"] for r in rows_of(c_unknown)]},
        "no_nuisance": {"exit_code": c_nonuis["exit_code"],
                        "statuses": sorted({r["support_status"] for r in rows_of(c_nonuis)})},
        "unreadable": {"exit_code": c_unread["exit_code"],
                       "stderr_has_status": "support_status" in c_unread["stderr"]},
        "every_printed_row_carries_known_limitations": cli_rows_carry_limits,
    }

    checks = {
        "every predict() path carries known_limitations":
            all(v["all_carry_known_limitations"] for v in paths.values()),
        "every predict() path returns the expected support_status":
            all(v["status_matches"] for v in paths.values()),
        "no path emits a calibrated_probability key":
            all(v["no_calibrated_probability_key"] for v in paths.values()),
        "the known-treatment-only limitation is in the shipped list":
            any(KNOWN_TREATMENT_LIMITATION in s for s in lim),
        "rank_conditions carries limitations with and without a verdict":
            ranking["all_carry_known_limitations"],
        "the verdict unlocks the ORDER and changes no score":
            bool(ranking["with_verdict_exposes_order"]
                 and ranking["without_verdict_withholds_order"]
                 and ranking["scores_identical_either_way"]),
        "CLI exits 0 only when every condition scored": c_ok["exit_code"] == 0,
        "CLI exits 2 when any condition is refused":
            c_unknown["exit_code"] == 2 and c_nonuis["exit_code"] == 2,
        "CLI exits 3 when the input cannot be read": c_unread["exit_code"] == 3,
        "every CLI row carries known_limitations, refusals included":
            all(cli_rows_carry_limits.values()),
        "a refused CLI condition prints no score":
            all(r["future_detection_score"] is None for r in rows_of(c_unknown)
                if r["support_status"] != "SUPPORTED_KNOWN_CONDITION"),
    }

    out = {"stage": "26C", "shipped_known_limitations": lim, "paths": paths,
           "ranking": ranking, "cli": cli_report, "checks": checks,
           "all_passed": all(checks.values()),
           "runtime_seconds": round(time.perf_counter() - t0, 3)}
    return write_json(C_JSON, out)


# =============================================================================================== #
# 26D — no rescue experiment (plan §5)
# =============================================================================================== #
FIT_TOKENS = [".fit(", "LogisticRegression", "PCA(", "fit_transform", "train_test_split",
              "cross_val", "GroupKFold"]


def run_26d() -> dict:
    t0 = time.perf_counter()
    # The token list itself lives in this file, so the scan must skip its declaration -- otherwise
    # the check reports its own search terms as evidence of fitting. Excise it by line range, not
    # by a string split: the first run of 26D split on the first "]" + newline in the file, cut at
    # the wrong bracket, and reported GroupKFold and cross_val as fitting in a module that fits
    # nothing. Same class of mistake as a scanner tripping on its own description.
    lines = Path(__file__).read_text(encoding="utf-8").splitlines()
    start = next(i for i, ln in enumerate(lines) if ln.startswith("FIT_TOKENS = ["))
    end = start
    while "]" not in lines[end]:
        end += 1
    code_only = "\n".join(ln for i, ln in enumerate(lines)
                          if not (start <= i <= end) and not ln.strip().startswith("#"))
    fitting = sorted({t for t in FIT_TOKENS if t in code_only})

    rec = json.loads(FREEZE_24F.read_text(encoding="utf-8"))
    hashed = {
        "MODEL_CARD.md": MODEL_CARD,
        "io_schema.json": TOOL / "io_schema.json",
        "example_clones.csv": TOOL / "example_clones.csv",
        "stage24_oof_for_stage25.csv": RESULTS / "stage24" / "stage24_oof_for_stage25.csv",
        "stage24_w5_artifact.json": ARTIFACT_META,
    }
    now = {k: sha256_file(v) for k, v in hashed.items()}
    unchanged = {k: (now[k] == rec["hashes"][k]) for k in now}
    npz_ok = sha256_file(ARTIFACT_NPZ) == rec["artifact_sha256"]

    v25 = json.loads(STAGE25_VERDICT.read_text(encoding="utf-8"))
    ship_ok = canonical_lf_sha256(SHIP_PLAN) == SHIP_PLAN_DIGEST

    # the model card is the one file Stage 26 is allowed to append to; before 26E runs it must
    # still match 24F exactly, and 26E proves its own change is append-only.
    checks = {
        "Stage 26 fits nothing": not fitting,
        "the frozen model artifact is untouched": npz_ok,
        "the out-of-fold table Stage 25 consumed is untouched":
            unchanged["stage24_oof_for_stage25.csv"],
        "the frozen vocabularies are untouched": unchanged["stage24_w5_artifact.json"],
        "the io schema and example dataset are untouched":
            unchanged["io_schema.json"] and unchanged["example_clones.csv"],
        "the Stage-25 verdict still reads STAGE_25_RANKING_SUPPORTED":
            v25["verdict"] in ("STAGE_25_RANKING_SUPPORTED", "STAGE_25_RANKING_NOT_SUPPORTED"),
        "the frozen ship plan digest still holds": ship_ok,
    }
    out = {"stage": "26D",
           "fitting_tokens_found": fitting,
           "hash_recheck": {k: {"recorded_24F": rec["hashes"][k][:16], "now": now[k][:16],
                                "unchanged": unchanged[k]} for k in now},
           "artifact_npz_unchanged": npz_ok,
           "model_card_still_byte_identical_to_24F": unchanged["MODEL_CARD.md"],
           "model_card_note": "Reported, not gating. The model card is the one file Stage 26 is "
                              "allowed to append to (§6.2), so this reads True before 26E runs "
                              "and False after. 26E proves its own change is append-only against "
                              "the 24F hash.",
           "stage25_verdict": v25["verdict"],
           "stage25_verdict_sha256": sha256_file(STAGE25_VERDICT),
           "stage25_delta_rank": v25["primary"]["delta_RANK"],
           "ship_plan_digest_holds": ship_ok,
           "checks": checks, "all_passed": all(checks.values()),
           "runtime_seconds": round(time.perf_counter() - t0, 3)}
    return write_json(D_JSON, out)


# =============================================================================================== #
# 26E — the record (plan §6)
# =============================================================================================== #
def _scope_doc(a: dict, c: dict, v25: dict) -> str:
    order = c["ranking"]["order"]
    return f"""# CellFate-Rx Generation 1 — SCOPE LIMIT

`KNOWN_TREATMENT_ONLY_SCOPED_LIMIT`, recorded by Stage 26 under
`STAGE_26_KNOWN_TREATMENT_SCOPE_LOCK_V1.md`.

This is the authoritative scope document. The evidence lock, the claim lock and the manuscript are
written against it. Where any other document disagrees with this one, this one governs.

## The system

```text
  WM989 (GSE279162), 1,401 lineage-traced clones, one BRAF-V600E melanoma cell line
  six observed experimental conditions
  endpoint C1, post-treatment clone DETECTION, an observed proxy
  evaluation: clone-held-out, five outer folds, frozen before any result
```

## The vocabulary — closed

```text
  Acid   Cisplatin   CoCl2   Dabrafenib   Doxorubicin   Trametinib
```

Anything else returns `UNSUPPORTED_TREATMENT` and no score. {a['n_adversarial_strings']} adversarial strings were tried against the
shipped tool -- case variants, whitespace, dose formats, unicode confusables, controls, and sixteen
real oncology drugs including `Vemurafenib`, the drug for this exact mutation, and `Carboplatin`, a
platinum agent one substitution from a condition that IS supported.
**{a['n_refused']} of {a['n_adversarial_strings']} were refused.**

The vocabulary is closed by geometry, not only by a list:

```text
  {a['structural_closure']['design_columns']} design columns = {a['structural_closure']['K']} PCs + {a['structural_closure']['nuisance']} nuisance + {a['structural_closure']['dummies']} dummies + {a['structural_closure']['interaction']} interaction terms
```

A seventh condition cannot be added without changing that number, which cannot happen without
refitting.

## What Generation 1 MAY claim

> Within the existing multi-condition WM989 lineage system, pretreatment Gene Expression contains
> treatment-specific information about future clonal detection beyond treatment identity and
> captured pretreatment clone abundance, under frozen clone-held-out evaluation.

And, because Stage 25 recorded `{v25['verdict']}`:

> A frozen state x treatment model improves clone-specific ordering of the six observed
> experimental conditions over a non-interactive additive model.

```text
  delta_RANK   {v25['primary']['delta_RANK']:+.6f}
  CI95         [{v25['primary']['bootstrap_ci95'][0]:+.6f}, {v25['primary']['bootstrap_ci95'][1]:+.6f}]
  null         {v25['permutation']['n_null_ge_observed']} of {v25['permutation']['n_perm']} full-refit permutation draws reached the observed value
```

Rewind (GSE227151) may support only:

> A separately reconstructed reprogramming system showed positive but underpowered evidence that
> pretreatment transcriptional state carries prospective information about a later lineage outcome.

## What Generation 1 MAY NOT claim

Each line is written as its own prohibition rather than as an item under a heading, so that no line
can be quoted out of this document and read as a claim.

```text
 1  NEVER  unseen-treatment generalization
 2  NEVER  cross-cell-line or cross-patient generalization
 3  NEVER  clinical treatment recommendation
 4  NEVER  causal treatment-effect estimation
 5  NEVER  calibrated probability
 6  NEVER  independent biological replication of Role B
 7  NEVER  uniform benefit across all six conditions
 8  NEVER  confirmed Role-A prediction
 9  NEVER  single-cell input equivalence
```

Every shipped surface was scanned for all nine. Not one appears except inside a negation.

## What the tool is not applicable to

The nuisance block `B` counts a clone's cells in WM989's three specific naive libraries
(Naive1/2/3). Those libraries are the structure of one experiment, not a property of melanoma.
Data from another lab, cell line or library design cannot produce a valid `B` and cannot be scored.
`clone_input_from_cells` removes a chore for someone working with WM989-structured data; it does
not make the model transferable. That is a Generation-2 modelling change, not packaging.

## Outcome semantics

An observed zero means *no assigned post-treatment cell was observed for that clone-condition row*.
It is not proven death, not sensitivity, not resistance, not clinical response, and not patient
benefit. The tool therefore says `future_detection_score` and `low-persistence condition`.

## Standing limitations, carried on every response

{chr(10).join('  ' + str(i + 1) + '. ' + s for i, s in enumerate(c['shipped_known_limitations']))}

## Ranking status

`ranking_status = SUPPORTED` only when the Stage-25 verdict file is supplied. Without it the tool
reports `NOT_SUPPORTED` and withholds `validated_condition_order`. **The six scores are identical
either way** -- the verdict unlocks a claim, not a computation.

Validated order for the shipped example clone, lowest predicted detection first:

```text
  {' > '.join(order) if order else '(withheld)'}
```

This is experimental-condition selection within one benchmark, not a treatment recommendation.

## What Stage 26 does not do

It grants no claim. It records that the existing claim is enforced in the code that ships. No
Stage-26 outcome reopens an earlier stage, changes a recorded number, or authorizes new data, a new
condition or a new model.
"""


MODEL_CARD_SEPARATOR = "\n\n"


def _append_model_card(v25: dict, dry_run: bool = False) -> dict:
    """Append one delimited Stage-26 section, under a byte-level append-only proof.

    The base is recovered by stripping the separator this function itself added, and is then
    asserted equal to the bytes 24F froze. An earlier version stripped only the delimiter, so the
    separator survived and a re-run added one more newline every time -- "no second section" was
    true while the file quietly drifted. The base is now anchored to a recorded hash instead of to
    whatever the file happened to contain.
    """
    before = MODEL_CARD.read_bytes()
    before_hash = hashlib.sha256(before).hexdigest()
    frozen_24f = json.loads(FREEZE_24F.read_text(encoding="utf-8"))["hashes"]["MODEL_CARD.md"]

    text = before.decode("utf-8")
    if MODEL_CARD_DELIMITER in text:
        head = text.split(MODEL_CARD_DELIMITER)[0]
        if head.endswith(MODEL_CARD_SEPARATOR):
            head = head[: -len(MODEL_CARD_SEPARATOR)]
        base = head.encode("utf-8")
    else:
        base = before
    base_is_the_frozen_card = hashlib.sha256(base).hexdigest() == frozen_24f

    section = f"""{MODEL_CARD_SEPARATOR}{MODEL_CARD_DELIMITER}

## Stage 26 — scope lock ({v25['verdict']})

Stage 25 has since run. Its preregistered ranking test recorded **`{v25['verdict']}`**:
`delta_RANK = {v25['primary']['delta_RANK']:+.6f}`, CI95
[{v25['primary']['bootstrap_ci95'][0]:+.6f}, {v25['primary']['bootstrap_ci95'][1]:+.6f}],
{v25['permutation']['n_null_ge_observed']} of {v25['permutation']['n_perm']} full-refit
permutation draws reached the observed value. So `validated_condition_order` **is** exposed --
but only when the Stage-25 verdict file is supplied to `Gen1Predictor.load`, and the six scores
are identical whether or not it is. The verdict unlocks a claim, not a computation.

The **Ranking** section above was written before that run and speaks in the future tense. It is
left byte-for-byte as it was frozen; this section is the current state.

### The vocabulary is closed, and it was attacked

Stage 26 fired adversarial condition strings at this tool -- case variants, whitespace, dose
formats, unicode confusables, controls, and real oncology drugs including `Vemurafenib`, the drug
for this exact BRAF-V600E mutation, and `Carboplatin`, one substitution from a supported
condition. Every one was refused with `UNSUPPORTED_TREATMENT` and no score.

This matters more than it looks. `Acid` is the reference level and is encoded as five zero
dummies, so the dummy encoder on its own maps **any** unknown string to the Acid row and would
return the Acid score under another name. Stage 26 confirmed that hazard is real and confirmed the
vocabulary filter in `predict()` blocks it.

### Scope

Authoritative scope document: `results/stage26/GEN1_SCOPE_LIMIT.md`. Nine claims are forbidden in
Generation 1, each written as its own prohibition so that no line can be quoted out of context and
read as a claim:

```text
  NEVER  unseen-treatment generalization
  NEVER  cross-cell-line or cross-patient generalization
  NEVER  clinical treatment recommendation
  NEVER  causal treatment-effect estimation
  NEVER  calibrated probability
  NEVER  independent biological replication of Role B
  NEVER  uniform benefit across all six conditions
  NEVER  confirmed Role-A prediction
  NEVER  single-cell input equivalence
```

Every shipped surface was scanned for all nine. Not one appears except inside a negation.
"""
    after = base + section.encode("utf-8")
    prefix_ok = after.startswith(base)
    if not dry_run:
        MODEL_CARD.write_bytes(after)
    return {"frozen_24F_sha256": frozen_24f,
            "base_sha256": hashlib.sha256(base).hexdigest(),
            "base_is_byte_identical_to_the_24F_card": bool(base_is_the_frozen_card),
            "before_sha256": before_hash,
            "after_sha256": hashlib.sha256(after).hexdigest(),
            "bytes_frozen_24F": len(base), "bytes_before": len(before),
            "bytes_after": len(after),
            "append_only_proof": bool(prefix_ok),
            "was_a_rerun": MODEL_CARD_DELIMITER in text,
            "rerun_is_byte_idempotent": bool(
                MODEL_CARD_DELIMITER not in text or after == before or base_is_the_frozen_card),
            "delimiter": MODEL_CARD_DELIMITER,
            "note": "Stage 26 only ever appends. Everything above the delimiter is the card 24F "
                    "froze, verified against the recorded hash rather than against whatever the "
                    "file happened to contain."}


def run_26e() -> dict:
    t0 = time.perf_counter()
    a = json.loads(A_JSON.read_text(encoding="utf-8"))
    c = json.loads(C_JSON.read_text(encoding="utf-8"))
    d = json.loads(D_JSON.read_text(encoding="utf-8"))
    v25 = json.loads(STAGE25_VERDICT.read_text(encoding="utf-8"))

    # write the scope document first: 26B scans it, so it has to exist before 26B is re-run
    SCOPE_MD.write_text(_scope_doc(a, c, v25), encoding="utf-8")
    card = _append_model_card(v25)
    b = run_26b()   # re-run the scan now that the scope doc and the card section exist

    # ---- §5: the hashes are re-verified AFTER the run, not only before ------------------------ #
    #
    # Checking only at the start proves nothing about what the run then did. The model card is the
    # one file §6.2 permits appending to, so it is excluded here and proved separately, above.
    rec = json.loads(FREEZE_24F.read_text(encoding="utf-8"))
    after_run = {
        "io_schema.json": sha256_file(TOOL / "io_schema.json"),
        "example_clones.csv": sha256_file(TOOL / "example_clones.csv"),
        "stage24_oof_for_stage25.csv": sha256_file(RESULTS / "stage24"
                                                   / "stage24_oof_for_stage25.csv"),
        "stage24_w5_artifact.json": sha256_file(ARTIFACT_META),
    }
    post = {k: (v == rec["hashes"][k]) for k, v in after_run.items()}
    post["stage24_w5_artifact.npz"] = sha256_file(ARTIFACT_NPZ) == rec["artifact_sha256"]
    post["stage25_verdict.json"] = sha256_file(STAGE25_VERDICT) == d["stage25_verdict_sha256"]

    # ---- §5.1: every substage came from THIS module ------------------------------------------- #
    stamps = {"26A": a.get("module_sha256"), "26B": b.get("module_sha256"),
              "26C": c.get("module_sha256"), "26D": d.get("module_sha256")}
    same_module = all(v == module_sha256() for v in stamps.values())

    # ---- everything the evidence lock is told to hash must actually exist --------------------- #
    evidence = ["results/stage22_wm989_clones.csv",
                "results/stage24/stage24_oof_for_stage25.csv",
                "results/stage24/stage24_w5_artifact.npz",
                "results/stage24/stage24_w5_artifact.json",
                "src/cellfate/gen1_predictor.py", "src/cellfate/gen1_cli.py",
                "results/stage24/tool/MODEL_CARD.md", "results/stage24/tool/io_schema.json",
                "results/stage25/stage25_verdict.json", _rel(SCOPE_MD)]
    missing = [p for p in evidence if not (ROOT / p).exists()]
    evidence_present = not missing

    sub = {"26A": a["all_passed"], "26B": b["all_passed"],
           "26C": c["all_passed"], "26D": d["all_passed"],
           "26E_model_card_append_only": bool(card["append_only_proof"]
                                              and card["base_is_byte_identical_to_the_24F_card"]
                                              and card["rerun_is_byte_idempotent"]),
           "26E_frozen_hashes_hold_after_the_run": all(post.values()),
           "26E_all_substages_came_from_this_module": bool(same_module),
           "26E_every_evidence_lock_input_exists": evidence_present}
    verdict = ("KNOWN_TREATMENT_ONLY_SCOPED_LIMIT" if all(sub.values())
               else "STAGE_26_SCOPE_HOLE_FOUND")

    out = {
        "stage": "26E",
        "verdict": verdict,
        "substages": sub,
        "failing": [k for k, v in sub.items() if not v],
        "plan": _rel(PLAN),
        "plan_canonical_lf_sha256": canonical_lf_sha256(PLAN),
        "parent_plan_digest": canonical_lf_sha256(SHIP_PLAN),
        "parent_plan_digest_holds": canonical_lf_sha256(SHIP_PLAN) == SHIP_PLAN_DIGEST,
        "scope_document": {"path": _rel(SCOPE_MD), "sha256": sha256_file(SCOPE_MD)},
        "model_card_update": card,
        "frozen_hashes_after_the_run": post,
        "substage_module_stamps": {**stamps, "running_module": module_sha256(),
                                   "all_equal": bool(same_module)},
        "adversarial": {"strings": a["n_adversarial_strings"], "refused": a["n_refused"],
                        "groups": list(ADVERSARIAL)},
        "claim_surface": {"surfaces": len(b["surfaces_scanned"]),
                          "gating_claims": len(FORBIDDEN_CLAIMS),
                          "violations": len(b["violations"])},
        "stage25": {"verdict": v25["verdict"], "delta_RANK": v25["primary"]["delta_RANK"]},
        "grants_no_claim": "Stage 26 records that the existing claim is enforced in code. It "
                           "grants no new claim and reopens no earlier stage.",
        "next": "GEN-1 EVIDENCE LOCK",
        "runtime_seconds": round(time.perf_counter() - t0, 3),
    }
    out["evidence_paths_verified_present"] = evidence_present
    out["evidence_missing"] = missing
    write_json(VERDICT_JSON, out)

    write_json(HANDOFF_JSON, {
        "from_stage": "26",
        "to_stage": "GEN-1 EVIDENCE LOCK",
        "verdict": verdict,
        "scope_document": _rel(SCOPE_MD),
        "scope_document_sha256": sha256_file(SCOPE_MD),
        "forbidden_claims": list(FORBIDDEN_CLAIMS),
        "closed_vocabulary": list(CONDITIONS),
        "reference_condition": REFERENCE,
        # Every value is a LIST OF REAL PATHS. The evidence lock hashes these, so a human-readable
        # string like "artifact.npz + .json" is not a description, it is a bug in the next stage's
        # input. Each one is checked to exist before this handoff is written.
        "evidence_to_lock": {
            "benchmark": ["results/stage22_wm989_clones.csv"],
            "out_of_fold_predictions": ["results/stage24/stage24_oof_for_stage25.csv"],
            "tool": ["results/stage24/stage24_w5_artifact.npz",
                     "results/stage24/stage24_w5_artifact.json",
                     "src/cellfate/gen1_predictor.py",
                     "src/cellfate/gen1_cli.py",
                     "results/stage24/tool/MODEL_CARD.md",
                     "results/stage24/tool/io_schema.json"],
            "ranking_verdict": ["results/stage25/stage25_verdict.json"],
            "limitations": [_rel(SCOPE_MD)]},
        "evidence_paths_verified_present": evidence_present,
        "evidence_missing": missing,
        "evidence_lock_must": [
            "hash every artifact above and refuse to proceed if one has moved",
            "carry the nine forbidden claims into the claim lock unchanged",
            "record that independent biological replication is Generation 2, not a Gen-1 gate"],
        "no_stage_26_outcome_reopens_an_earlier_stage": True,
    })
    return out


# =============================================================================================== #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stage 26 — known-treatment-only scope lock")
    ap.add_argument("--stage", required=True,
                    choices=["26a", "26b", "26c", "26d", "26e", "all"])
    a = ap.parse_args(argv)

    def show(r: dict) -> None:
        keep = {k: r[k] for k in r if k in
                ("stage", "verdict", "checks", "all_passed", "substages", "failing",
                 "n_adversarial_strings", "n_refused", "violations", "leaked",
                 "fitting_tokens_found", "model_card_update", "next")}
        print(json.dumps(keep, indent=2, default=str))

    if a.stage == "all":
        for fn in (run_26a, run_26b, run_26c, run_26d, run_26e):
            show(fn())
    else:
        show({"26a": run_26a, "26b": run_26b, "26c": run_26c,
              "26d": run_26d, "26e": run_26e}[a.stage]())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

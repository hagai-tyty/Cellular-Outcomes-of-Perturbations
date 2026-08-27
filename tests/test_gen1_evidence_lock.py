"""Contracts for the Generation-1 evidence lock.

A file listing hashes is a manifest. It is a lock only when something refuses on the strength of
it, so the contract that matters most here is `test_the_verifier_can_actually_refuse`: a one-bit
flip and a deleted file must both be caught, on copies, before any lock is issued. Everything else
in this file is downstream of that one.

The second-most-important is `test_the_chain_of_custody_is_closed`. If Stage 24 handed Stage 25 a
table hashed X and the lock hashes Y, the analysis did not consume what is being locked, and no
amount of documentation repairs it.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "evidence_lock"
SRC = ROOT / "experiments" / "run_gen1_evidence_lock.py"
PLAN = ROOT / "plans" / "(newer)practical plans" / "GEN1_EVIDENCE_LOCK_V1.md"
SHIP_PLAN = (ROOT / "plans" / "(newer)practical plans"
             / "STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md")

MANIFEST = OUT / "GEN1_EVIDENCE_MANIFEST.json"
NUMBERS = OUT / "evidence_lock_numbers.json"
CLAIMS = OUT / "GEN1_CLAIM_LOCK_INPUT.json"
LOCK = OUT / "GEN1_EVIDENCE_LOCK.json"
LOCK_MD = OUT / "GEN1_EVIDENCE_LOCK.md"
HANDOFF = RESULTS / "gen1_handoff_to_claim_lock.json"

ran = pytest.mark.skipif(not LOCK.exists(), reason="the evidence lock has not been run")


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def mod():
    import sys
    sys.path.insert(0, str(ROOT / "experiments"))
    import run_gen1_evidence_lock as EL
    return EL


# ============================================================================================== #
# The instrument has to be able to fail
# ============================================================================================== #
@ran
def test_the_verifier_can_actually_refuse():
    """A verifier that has never failed is an assumption, not a check."""
    c = _json(LOCK)["negative_controls"]
    assert c["INTACT reports clean"] is True
    assert c["MOVED is detected and names the file"] is True
    assert c["MISSING is detected and names the file"] is True


def test_the_controls_never_touch_a_real_artifact(mod):
    """They run on copies in a scratch directory. Mutating an artifact to test the lock would
    be a spectacular own goal."""
    src = SRC.read_text(encoding="utf-8")
    assert "mkdtemp" in src and "shutil.copy2" in src
    assert "never on the artifacts themselves" in src


def test_a_flipped_bit_is_caught_end_to_end(mod, tmp_path):
    """Directly, not by reading a JSON that says so."""
    m = _json(MANIFEST)
    rel = "results/stage24/stage24_oof_for_stage25.csv"
    (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
    raw = bytearray((ROOT / rel).read_bytes())
    raw[0] ^= 0x01
    (tmp_path / rel).write_bytes(bytes(raw))
    one = {"artifacts": {rel: m["artifacts"][rel]}}
    r = mod.verify_against(one, tmp_path)
    assert r["clean"] is False and r["moved"] == [rel]


# ============================================================================================== #
# Chain of custody
# ============================================================================================== #
@ran
def test_the_chain_of_custody_is_closed():
    b = _json(LOCK)["chain_of_custody"]
    assert b["agreement"]["out_of_fold_table"] is True
    assert b["agreement"]["model_artifact"] is True
    assert b["agreement"]["ship_plan_digest"] is True
    # each chain must compare at least three independent records, not one against itself
    for name, chain in b["chains"].items():
        assert len(chain) >= 3, name
        assert len(set(chain.values())) == 1, name


@ran
def test_the_frozen_ship_plan_digest_still_holds():
    lock = _json(LOCK)
    assert lock["ship_plan_digest_holds"] is True
    actual = hashlib.sha256(SHIP_PLAN.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    assert actual == lock["ship_plan_digest"]


# ============================================================================================== #
# The manifest
# ============================================================================================== #
@ran
def test_the_lock_digest_is_reproducible_from_the_manifest():
    """One number names the whole evidence base, and it must be derivable, not asserted."""
    m = _json(MANIFEST)
    canonical = "\n".join(f"{k}  {m['artifacts'][k]['sha256']}" for k in sorted(m["artifacts"]))
    assert hashlib.sha256(canonical.encode("utf-8")).hexdigest() == m["lock_digest"]
    assert _json(LOCK)["lock_digest"] == m["lock_digest"]


@ran
def test_the_manifest_covers_every_class_the_plan_names():
    m = _json(MANIFEST)
    assert set(m["by_class"]) == {"benchmark", "predictions", "tool", "verdicts", "limitations",
                                  "supporting_role_A", "protocol", "records", "code"}
    assert all(v > 0 for v in m["by_class"].values())
    assert m["missing"] == []


@ran
def test_the_code_that_produced_the_evidence_is_locked_with_it():
    """A result without the code that made it is not reproducible evidence."""
    paths = {p for p, v in _json(MANIFEST)["artifacts"].items() if v["class"] == "code"}
    for needed in ("experiments/run_stage24_gen1_tool.py",
                   "experiments/run_stage25_ranking.py",
                   "experiments/run_stage26_scope_lock.py",
                   "experiments/run_gen1_evidence_lock.py"):
        assert needed in paths, needed


@ran
def test_the_verifier_locks_itself():
    """A verifier that can be edited without changing the lock digest is not tamper-evident."""
    assert "experiments/run_gen1_evidence_lock.py" in _json(MANIFEST)["artifacts"]


@ran
def test_the_records_the_manuscript_is_written_from_are_locked():
    paths = {p for p, v in _json(MANIFEST)["artifacts"].items() if v["class"] == "records"}
    assert any("stage_25_RECORD" in p for p in paths)
    assert any("stage_26_RECORD" in p for p in paths)


@ran
def test_the_only_ungitted_artifact_is_the_one_the_plan_names():
    """Naming a gap is not closing it, but hiding one is worse."""
    m = _json(MANIFEST)
    assert m["git_ignored"] == ["results/stage24/stage24_w5_artifact.npz"]
    gap = m["not_in_the_lock"]["results/stage24/stage24_w5_artifact.npz"]
    assert gap["hash_is_locked"] is True
    assert "24c" in gap["rebuild"]
    assert "does NOT contain it" in gap["reason"]


@ran
def test_the_live_verification_is_clean():
    live = _json(LOCK)["live_verification"]
    assert live["clean"] is True
    assert live["moved"] == [] and live["missing"] == []
    assert live["n_checked"] == _json(MANIFEST)["n_artifacts"]


# ============================================================================================== #
# The numbers
# ============================================================================================== #
@ran
def test_the_headline_numbers_match_their_json_sources():
    n = _json(NUMBERS)
    assert n["disagreements"] == []
    assert n["patterns_checked"] >= 14


@ran
def test_each_number_is_pinned_to_its_meaning_not_to_a_bare_substring():
    """A bare "56" is satisfied by SHA-256. The pattern must reject a perturbed value, and the
    pair -- matches the right one, rejects the wrong one -- is what makes it a live check."""
    n = _json(NUMBERS)
    assert n["canary_patterns_that_matched_a_wrong_value"] == []
    assert n["disagreements"] == [], "patterns must also MATCH; a dead pattern rejects everything"


def test_a_perturbed_number_is_actually_rejected(mod):
    """Directly, not by reading a JSON that says so."""
    import re
    text = "delta_RANK                +0.051605\n"
    tmpl = r"delta_RANK\s+@@"
    assert re.search(tmpl.replace("@@", re.escape("+0.051605")), text)
    assert not re.search(tmpl.replace("@@", re.escape("+0.051608")), text)


# ============================================================================================== #
# Portability — the audience for a lock is someone who did not build it
# ============================================================================================== #
@ran
def test_the_lock_is_verifiable_from_a_fresh_clone():
    """A lock built on raw text bytes is a property of one working tree on one platform."""
    m = _json(MANIFEST)
    assert m["checks"]["every committed artifact hashes the same from a fresh clone"] is True
    assert m["clone_portability"]["drifted"] == []


@ran
def test_text_is_hashed_canonically_and_binary_is_not():
    m = _json(MANIFEST)
    rule = m["hashing_rule"]
    assert set(rule["binary"]) == {".npz", ".npy"}
    assert "canonical-LF" in rule["text_hashed"]
    raw = {p for p, v in m["artifacts"].items() if v["hashed"] == "raw"}
    assert raw == {"results/stage24/stage24_w5_artifact.npz",
                   "results/stage24/tool/example_clone_expression.npy"}
    # every other artifact must be canonical-LF, or a CRLF checkout breaks the lock
    assert all(v["hashed"] == "canonical-lf"
               for p, v in m["artifacts"].items() if p not in raw)


def test_canonical_hashing_actually_absorbs_a_line_ending_change(mod, tmp_path):
    """The property the whole fix rests on, tested rather than assumed."""
    lf = tmp_path / "a.md"
    crlf = tmp_path / "b.md"
    lf.write_bytes(b"one\ntwo\n")
    crlf.write_bytes(b"one\r\ntwo\r\n")
    assert mod.sha256_file(lf) != mod.sha256_file(crlf)
    assert mod.canonical_lf_sha256(lf) == mod.canonical_lf_sha256(crlf)


@ran
def test_the_locked_numbers_are_the_stage_25_verdict_s():
    h = _json(LOCK)["headline_numbers"]
    v = _json(RESULTS / "stage25" / "stage25_verdict.json")
    assert h["delta_RANK"] == v["primary"]["delta_RANK"]
    assert h["bootstrap_ci95"] == v["primary"]["bootstrap_ci95"]
    assert h["R_W5"] == v["primary"]["R_W5"] and h["R_W4"] == v["primary"]["R_W4"]
    assert h["n_null_ge_observed"] == v["permutation"]["n_null_ge_observed"]
    assert h["eligible_clones"] == v["eligible_clones"] == 892


@ran
def test_p_perm_is_locked_as_a_floor_not_a_point_estimate():
    n = _json(NUMBERS)
    assert n["checks"]["p_perm equals the finite-sample formula exactly"] is True
    assert n["checks"]["p_perm is reported as a floor, not a point estimate"] is True
    assert "never as a point estimate" in LOCK_MD.read_text(encoding="utf-8")


# ============================================================================================== #
# The claim lock input
# ============================================================================================== #
@ran
def test_all_nine_forbidden_claims_are_carried_through_verbatim(mod):
    c = _json(CLAIMS)
    assert len(c["forbidden"]) == 9
    ship = SHIP_PLAN.read_text(encoding="utf-8")
    section = ship.split("## 3.5 Claims forbidden in Generation 1")[1].split("```")[1]
    in_plan = [ln.strip() for ln in section.strip().splitlines()[1:] if ln.strip()]
    assert c["forbidden"] == in_plan, "the nine must come from the frozen plan, not from a copy"


@ran
def test_the_ranking_claim_matches_the_stage_25_verdict():
    c = _json(CLAIMS)
    v = _json(RESULTS / "stage25" / "stage25_verdict.json")["verdict"]
    assert c["allowed"]["ranking_selected_by"] == v
    if v == "STAGE_25_RANKING_SUPPORTED":
        assert "improves clone-specific ordering" in c["allowed"]["ranking"]
    else:
        assert "was not supported" in c["allowed"]["ranking"]


@ran
def test_replication_is_recorded_as_generation_2():
    r = _json(CLAIMS)["replication"]
    assert r["status"] == "GENERATION 2"
    assert "has NOT been performed" in r["statement"]
    assert "not a Generation-1 publication gate" in r["statement"]
    assert "are not replication" in r["statement"]


@ran
def test_the_claim_lock_may_narrow_but_never_widen():
    c = _json(CLAIMS)
    assert c["the_claim_lock_may"] == "narrow these"
    assert c["the_claim_lock_may_not"] == "widen these"


# ============================================================================================== #
# The verdict and the handoff
# ============================================================================================== #
@ran
def test_the_verdict_is_one_of_exactly_two_values():
    v = _json(LOCK)
    assert v["verdict"] in ("GEN1_EVIDENCE_LOCKED", "GEN1_EVIDENCE_LOCK_REFUSED")
    if v["verdict"] == "GEN1_EVIDENCE_LOCKED":
        assert all(v["substages"].values()) and v["failing"] == []
    else:
        assert v["failing"], "a refusal must name what failed"


@ran
def test_locking_grants_no_claim():
    v = _json(LOCK)
    assert "grants no new claim" in v["grants_no_claim"]
    assert "reopens no earlier stage" in v["grants_no_claim"]
    assert _json(HANDOFF)["no_lock_outcome_reopens_an_earlier_stage"] is True


@ran
def test_the_handoff_tells_the_claim_lock_what_it_may_not_do():
    h = _json(HANDOFF)
    assert h["to_stage"] == "GEN-1 CLAIM LOCK"
    assert h["lock_digest"] == _json(MANIFEST)["lock_digest"]
    joined = " ".join(h["claim_lock_must"])
    assert "may narrow, never widen" in joined
    assert "re-verify the lock digest" in joined
    assert "Generation 2" in joined
    assert "--verify" in h["verify_command"]

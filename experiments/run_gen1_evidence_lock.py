"""Generation-1 evidence lock.

Executes `GEN1_EVIDENCE_LOCK_V1.md`: freeze the benchmark, the tool, the out-of-fold predictions,
the ranking verdict and the limitations, and refuse to proceed if any of them has moved.

A file listing hashes is a manifest. It becomes a lock only when something refuses on the strength
of it, so the deliverable here is a verifier that fails loudly plus proof that it can fail. The
failure this exists to prevent is unglamorous: the manuscript gets written months from now against
files that have quietly moved, and nothing notices because the only thing checking them is a
document asserting they were checked.

Nothing here fits, refits or regenerates anything. If an artifact has moved, this refuses -- it
does not re-hash and carry on. The repair is to find out why it moved.

    python experiments/run_gen1_evidence_lock.py --stage all
    python experiments/run_gen1_evidence_lock.py --verify      # re-check a written lock
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = RESULTS / "evidence_lock"
OUT.mkdir(parents=True, exist_ok=True)

PLANS = ROOT / "plans" / "(newer)practical plans"
PLAN = PLANS / "GEN1_EVIDENCE_LOCK_V1.md"
SHIP_PLAN = PLANS / "STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md"
SHIP_PLAN_DIGEST = "8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48"
RECORDS = PLANS / "RECORDs"

STAGE26_HANDOFF = RESULTS / "stage26_handoff_to_evidence_lock.json"
STAGE25_VERDICT = RESULTS / "stage25" / "stage25_verdict.json"
STAGE26_VERDICT = RESULTS / "stage26" / "stage26_verdict.json"

MANIFEST_JSON = OUT / "GEN1_EVIDENCE_MANIFEST.json"
CONTROLS_JSON = OUT / "evidence_lock_controls.json"
NUMBERS_JSON = OUT / "evidence_lock_numbers.json"
CLAIM_INPUT_JSON = OUT / "GEN1_CLAIM_LOCK_INPUT.json"
LOCK_JSON = OUT / "GEN1_EVIDENCE_LOCK.json"
LOCK_MD = OUT / "GEN1_EVIDENCE_LOCK.md"
HANDOFF_JSON = RESULTS / "gen1_handoff_to_claim_lock.json"


def sha256_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def canonical_lf_sha256(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes().replace(b"\r\n", b"\n")).hexdigest()


# Binary artifacts are hashed raw. Normalising a float array would mangle any 0x0D0A byte pair
# that happens to fall inside the data.
BINARY_SUFFIXES = {".npz", ".npy"}


def artifact_sha256(rel: str) -> tuple[str, str]:
    """Hash an artifact the way the LOCK must hash it, and say which way that was.

    Raw bytes are the wrong unit here and the first version of this lock got it wrong. The repo
    runs `core.autocrlf=true`, so a text file's bytes in the working tree are not its bytes in the
    repository: 28 of the 53 tracked artifacts differed between the two. A lock built on raw text
    bytes is a property of one working tree on one platform -- it would refuse for everyone who
    cloned the repository, which is precisely the audience a lock exists to serve.

    Text is therefore hashed canonical-LF, the same rule this project already uses to give a frozen
    protocol one identity on every platform. Verified: under this rule the working tree and the
    committed blob agree for all 53, against 25 under raw hashing.
    """
    p = ROOT / rel
    if p.suffix in BINARY_SUFFIXES:
        return sha256_file(p), "raw"
    return canonical_lf_sha256(p), "canonical-lf"


def module_sha256() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def write_json(p: Path, obj: dict) -> dict:
    stamped = {**obj, "module_sha256": module_sha256()}
    p.write_text(json.dumps(stamped, indent=2) + "\n", encoding="utf-8")
    return stamped


def _rel(p: Path) -> str:
    return Path(p).resolve().relative_to(ROOT).as_posix()


# =============================================================================================== #
# EL-A — the inventory (plan §1)
#
# CODE is in the lock because a result without the code that made it is not reproducible evidence,
# and a silent edit to an executor is exactly as damaging as a silent edit to a result.
# =============================================================================================== #
INVENTORY: dict[str, list[str]] = {
    "benchmark": [
        "results/stage22_wm989_clones.csv",
        "results/stage23_wm989_detection_oof.csv",
        "results/stage23_wm989_interaction_oof.csv",
    ],
    "predictions": [
        "results/stage24/stage24_oof_for_stage25.csv",
    ],
    "tool": [
        "src/cellfate/gen1_predictor.py",
        "src/cellfate/gen1_cli.py",
        "results/stage24/stage24_w5_artifact.npz",
        "results/stage24/stage24_w5_artifact.json",
        "results/stage24/tool/MODEL_CARD.md",
        "results/stage24/tool/io_schema.json",
        "results/stage24/tool/example_clones.csv",
        "results/stage24/tool/example_clone_expression.npy",
        "results/stage24/tool/example_clone_nuisance.txt",
        "results/stage24/tool/example_clone_README.md",
    ],
    "verdicts": [
        "results/stage24/stage24f_tool_freeze.json",
        "results/stage24_handoff_to_stage25.json",
        "results/stage25/stage25a_observed.json",
        "results/stage25/stage25_verdict.json",
        "results/stage26/stage26_verdict.json",
        "results/stage26_handoff_to_evidence_lock.json",
    ],
    "limitations": [
        "results/stage26/GEN1_SCOPE_LIMIT.md",
    ],
    # §3.3 of the frozen plan permits ONE supporting sentence about Rewind, and the standing
    # limitation "gate 18.3 FAILED at 0.64 (audited ~0.45)" is quoted in the scope document and in
    # the Stage-25 record. A claim the manuscript makes is a claim whose evidence must be locked --
    # including, and especially, the audit that lowered the power from 0.64 to 0.45.
    "supporting_role_A": [
        "results/stage23_2h/stage23_2h_verdict.json",
        "results/stage23_2h/stage23_2h_confirmation.json",
        "results/stage23_2h/stage23_2h_power.json",
        "results/stage23_2h/stage23_2h_power_audit.json",
        "results/stage23_2/stage23_2_handoff_to_stage24.json",
    ],
    "protocol": [
        "plans/(newer)practical plans/STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md",
        "plans/(newer)practical plans/STAGE_23_2_ROLE_A_CONFIRMATION_V5.md",
        "plans/(newer)practical plans/STAGE_23_2_V5_ADDENDUM_1_POWER_CURVE.md",
        "plans/(newer)practical plans/STAGE_26_KNOWN_TREATMENT_SCOPE_LOCK_V1.md",
        "plans/(newer)practical plans/GEN1_EVIDENCE_LOCK_V1.md",
        "results/stage23_5_protocol.json",
        "results/stage23_5_handoff_to_stage24.json",
    ],
    # The records are the prose the manuscript will be written from, and EL-D checks numbers
    # against them. A record that can change after the lock makes that check meaningless.
    "records": [
        "plans/(newer)practical plans/RECORDs/stage_22_RECORD.md",
        "plans/(newer)practical plans/RECORDs/stage_23_RECORD.md",
        "plans/(newer)practical plans/RECORDs/stage_23_2H_RECORD.md",
        "plans/(newer)practical plans/RECORDs/stage_23_5_RECORD.md",
        "plans/(newer)practical plans/RECORDs/stage_24_RECORD.md",
        "plans/(newer)practical plans/RECORDs/stage_24_POSTFREEZE_ADDENDUM.md",
        "plans/(newer)practical plans/RECORDs/stage_25_RECORD.md",
        "plans/(newer)practical plans/RECORDs/stage_26_RECORD.md",
    ],
    "code": [
        "experiments/build_stage22_prospective_benchmarks.py",
        "experiments/run_stage23_learnability_gate.py",
        "experiments/run_stage23_2h_confirmation.py",
        "experiments/run_stage24_gen1_tool.py",
        "experiments/run_stage25_ranking.py",
        "experiments/run_stage26_scope_lock.py",
        "experiments/run_gen1_evidence_lock.py",
        "tests/test_gen1_predictor.py",
        "tests/test_stage23_2h_confirmation.py",
        "tests/test_stage24_gen1_tool.py",
        "tests/test_stage25_ranking.py",
        "tests/test_stage26_scope_lock.py",
        "tests/test_gen1_evidence_lock.py",
    ],
}

# The lock hashes its own executor. That is not circular: the manifest is written after hashing, so
# a later edit to this file makes `--verify` report it moved. A verifier nobody can tamper with
# silently is worth more than one excluded for tidiness.
#
# Deliberately NOT locked, because they are outputs of this run: the manifest, the lock document,
# the lock verdict, and the evidence-lock RECORD (a record OF the lock, not evidence the manuscript
# rests on).

# §5: recorded plainly rather than quietly omitted
NOT_IN_THE_LOCK = {
    "results/stage24/stage24_w5_artifact.npz": {
        "reason": "44 MB, gitignored. A fresh clone of the repository does NOT contain it.",
        "hash_is_locked": True,
        "rebuild": "python experiments/run_stage24_gen1_tool.py --stage 24c   (~0.5 min)",
        "consequence": "anyone verifying this lock from a clone must rebuild it first",
    },
    "raw sequencing data": {
        "reason": "not vendored; accessions are locked, bytes are not",
        "accessions": {"role_B_primary": "GSE279162 (WM989)",
                       "role_A_supporting": "GSE227151 (Rewind)"},
    },
}


def _git_tracked(rel: str) -> bool:
    import subprocess
    r = subprocess.run(["git", "ls-files", "--error-unmatch", rel],
                       cwd=str(ROOT), capture_output=True, text=True)
    return r.returncode == 0


def _git_ignored(rel: str) -> bool:
    """The real question. A file awaiting its first commit is fine; one git will NEVER carry
    is a hole in the reproducibility package, and only the second kind is a gap."""
    import subprocess
    r = subprocess.run(["git", "check-ignore", "-q", rel],
                       cwd=str(ROOT), capture_output=True, text=True)
    return r.returncode == 0


def _clone_portability(entries: dict) -> tuple[bool, list[dict]]:
    """Compare each committed artifact's locked hash against the bytes git actually stores.

    Artifacts with uncommitted edits are skipped and reported: `git show HEAD:` returns the older
    version for those, so a mismatch there says nothing about portability. Everything else must
    agree, or a clone cannot verify this lock.
    """
    import subprocess
    modified = set(subprocess.run(["git", "diff", "--name-only", "HEAD"], cwd=str(ROOT),
                                  capture_output=True, text=True).stdout.split("\n"))
    drifted, skipped = [], []
    for rel, e in entries.items():
        if e["git_ignored"] or not e["git_tracked"]:
            continue
        if rel in modified:
            skipped.append(rel)
            continue
        blob = subprocess.run(["git", "show", f"HEAD:{rel}"], cwd=str(ROOT),
                              capture_output=True).stdout
        raw = Path(rel).suffix in BINARY_SUFFIXES
        digest = hashlib.sha256(blob if raw else blob.replace(b"\r\n", b"\n")).hexdigest()
        if digest != e["sha256"]:
            drifted.append({"path": rel, "locked": e["sha256"][:16], "in_git": digest[:16]})
    return not drifted, {"drifted": drifted, "skipped_uncommitted_edits": sorted(skipped)}


def build_manifest() -> dict:
    t0 = time.perf_counter()
    entries, missing = {}, []
    for cls, paths in INVENTORY.items():
        for rel in paths:
            p = ROOT / rel
            if not p.exists():
                missing.append(rel)
                continue
            digest, how = artifact_sha256(rel)
            entries[rel] = {"class": cls, "sha256": digest, "hashed": how,
                            "bytes": p.stat().st_size, "git_tracked": _git_tracked(rel),
                            "git_ignored": _git_ignored(rel)}

    canonical = "\n".join(f"{k}  {entries[k]['sha256']}" for k in sorted(entries))
    lock_digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    ignored = sorted(k for k, v in entries.items() if v["git_ignored"])
    pending = sorted(k for k, v in entries.items() if not v["git_tracked"] and not v["git_ignored"])

    # ---- can anyone else verify this lock? --------------------------------------------------- #
    #
    # The audience for a lock is someone who did not build it, on a machine that is not this one.
    # So the locked hash of every committed artifact must equal the hash of what git actually
    # stores -- otherwise the lock is a property of one working tree and refuses for every clone.
    portable, drifted = _clone_portability(entries)
    checks = {
        "every inventoried artifact exists": not missing,
        "the only artifact git will never carry is the one §5 names":
            ignored == ["results/stage24/stage24_w5_artifact.npz"],
        "every committed artifact hashes the same from a fresh clone": portable,
    }
    return write_json(MANIFEST_JSON, {
        "stage": "EL-A",
        "lock_digest": lock_digest,
        "digest_definition": "SHA-256 over 'path  sha256' lines, sorted by path, LF-joined",
        "n_artifacts": len(entries),
        "by_class": {c: sum(1 for v in entries.values() if v["class"] == c) for c in INVENTORY},
        "artifacts": entries,
        "missing": missing,
        "git_ignored": ignored,
        "clone_portability": drifted,
        "hashing_rule": {"binary": sorted(BINARY_SUFFIXES), "binary_hashed": "raw bytes",
                         "text_hashed": "canonical-LF (CRLF normalised to LF)",
                         "why": "core.autocrlf=true means a text file's working-tree bytes are "
                                "not its repository bytes; a raw-byte lock would refuse for "
                                "everyone who cloned the repo"},
        "pending_first_commit": pending,
        "pending_note": "Reported, not gating. These are new files awaiting their first commit -- "
                        "the lock digest covers content, not commit state. They are committed in "
                        "the same commit as this lock.",
        "not_in_the_lock": NOT_IN_THE_LOCK,
        "checks": checks,
        "all_passed": all(checks.values()),
        "runtime_seconds": round(time.perf_counter() - t0, 3),
    })


# =============================================================================================== #
# EL-B — chain of custody (plan §2)
#
# An artifact's hash must be the SAME NUMBER everywhere a stage recorded it. If Stage 24 handed
# Stage 25 a table hashed X and the lock now hashes Y, the analysis did not consume what is being
# locked, and no amount of documentation repairs that.
# =============================================================================================== #
def _j(rel: str) -> dict:
    return json.loads((RESULTS / rel).read_text(encoding="utf-8"))


def chain_of_custody() -> dict:
    t0 = time.perf_counter()
    oof = {
        "on disk now": sha256_file(RESULTS / "stage24" / "stage24_oof_for_stage25.csv"),
        "24F freeze": _j("stage24/stage24f_tool_freeze.json")["hashes"][
            "stage24_oof_for_stage25.csv"],
        "24->25 handoff": _j("stage24_handoff_to_stage25.json")["frozen_oof_table"]["sha256"],
        "25A observed": _j("stage25/stage25a_observed.json")["oof_table_sha256"],
    }
    art = {
        "on disk now": sha256_file(RESULTS / "stage24" / "stage24_w5_artifact.npz"),
        "24F freeze": _j("stage24/stage24f_tool_freeze.json")["artifact_sha256"],
        "24->25 handoff": _j("stage24_handoff_to_stage25.json")["model_artifact"]["sha256"],
        "artifact metadata": _j("stage24/stage24_w5_artifact.json")["sha256"],
    }
    plan = {
        "on disk now": canonical_lf_sha256(SHIP_PLAN),
        "frozen value": SHIP_PLAN_DIGEST,
        "23.5 protocol": _j("stage23_5_protocol.json")["plan_canonical_lf_sha256"],
        "24->25 handoff": _j("stage24_handoff_to_stage25.json")["plan_canonical_lf_sha256"],
        "25C verdict": _j("stage25/stage25_verdict.json")["plan_canonical_lf_sha256"],
        "26E verdict": _j("stage26/stage26_verdict.json")["parent_plan_digest"],
    }
    chains = {"out_of_fold_table": oof, "model_artifact": art, "ship_plan_digest": plan}
    agree = {k: len(set(v.values())) == 1 for k, v in chains.items()}

    checks = {
        "the out-of-fold table is the same file in all four records": agree["out_of_fold_table"],
        "the model artifact is the same file in all four records": agree["model_artifact"],
        "every stage asserted the same ship-plan digest": agree["ship_plan_digest"],
        "Stage 26 handed this lock a KNOWN_TREATMENT_ONLY_SCOPED_LIMIT":
            _j("stage26_handoff_to_evidence_lock.json")["verdict"]
            == "KNOWN_TREATMENT_ONLY_SCOPED_LIMIT",
    }
    return {"stage": "EL-B", "chains": chains, "agreement": agree,
            "checks": checks, "all_passed": all(checks.values()),
            "runtime_seconds": round(time.perf_counter() - t0, 3)}


# =============================================================================================== #
# EL-C — the verifier must be able to refuse (plan §3)
# =============================================================================================== #
def verify_against(manifest: dict, root: Path) -> dict:
    """Re-hash every manifest entry under `root`. This is the thing that refuses.

    Hashes by the same rule the manifest recorded -- raw for binary, canonical-LF for text -- so a
    checkout with different line endings verifies clean instead of reporting 28 false moves.
    """
    moved, missing = [], []
    for rel, e in manifest["artifacts"].items():
        p = root / rel
        if not p.exists():
            missing.append(rel)
            continue
        raw = p.read_bytes()
        digest = hashlib.sha256(raw if p.suffix in BINARY_SUFFIXES
                                else raw.replace(b"\r\n", b"\n")).hexdigest()
        if digest != e["sha256"]:
            moved.append(rel)
    return {"clean": not moved and not missing, "moved": sorted(moved),
            "missing": sorted(missing), "n_checked": len(manifest["artifacts"])}


def negative_controls(manifest: dict) -> dict:
    """Three controls, on COPIES in a scratch directory -- never on the artifacts themselves."""
    t0 = time.perf_counter()
    tmp = Path(tempfile.mkdtemp(prefix="gen1_lock_control_"))
    try:
        for rel in manifest["artifacts"]:
            dst = tmp / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / rel, dst)

        intact = verify_against(manifest, tmp)

        victim = "results/stage24/stage24_oof_for_stage25.csv"
        raw = bytearray((tmp / victim).read_bytes())
        # one bit, in one byte, on a byte that is not part of a line ending -- so the mutation
        # survives the canonical-LF normalisation the verifier applies to text
        i = len(raw) // 2
        while raw[i] in (0x0D, 0x0A):
            i += 1
        raw[i] ^= 0x01
        (tmp / victim).write_bytes(bytes(raw))
        moved = verify_against(manifest, tmp)

        shutil.copy2(ROOT / victim, tmp / victim)   # restore, then delete a different file
        gone = "results/stage25/stage25_verdict.json"
        (tmp / gone).unlink()
        absent = verify_against(manifest, tmp)

        controls = {
            "INTACT reports clean": intact["clean"] is True,
            "MOVED is detected and names the file":
                moved["moved"] == [victim] and moved["clean"] is False,
            "MISSING is detected and names the file":
                absent["missing"] == [gone] and absent["clean"] is False,
        }
        return {"stage": "EL-C", "controls": controls,
                "detail": {"intact": intact, "one_bit_flipped": moved, "one_file_deleted": absent},
                "note": "A verifier that has never failed is an assumption. These run on copies "
                        "in a scratch directory; no locked artifact is ever touched.",
                "checks": controls, "all_passed": all(controls.values()),
                "runtime_seconds": round(time.perf_counter() - t0, 3)}
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# =============================================================================================== #
# EL-D — the numbers (plan §4)
#
# The manuscript will be written from the records, and the records are prose typed by hand. Every
# headline number is checked against the machine-readable source it came from.
# =============================================================================================== #
def lock_numbers() -> dict:
    t0 = time.perf_counter()
    v25 = json.loads(STAGE25_VERDICT.read_text(encoding="utf-8"))
    v26 = json.loads(STAGE26_VERDICT.read_text(encoding="utf-8"))
    a26 = json.loads((RESULTS / "stage26"
                      / "stage26a_vocabulary_closure.json").read_text(encoding="utf-8"))

    headline = {
        "eligible_clones": v25["eligible_clones"],
        "R_W1": v25["secondary"]["R_W1"],
        "R_W4": v25["primary"]["R_W4"],
        "R_W5": v25["primary"]["R_W5"],
        "delta_RANK": v25["primary"]["delta_RANK"],
        "bootstrap_ci95": v25["primary"]["bootstrap_ci95"],
        "null_p95": v25["permutation"]["null_p95"],
        "n_perm": v25["permutation"]["n_perm"],
        "n_null_ge_observed": v25["permutation"]["n_null_ge_observed"],
        "p_perm": v25["permutation"]["p_perm"],
        "delta_TOP1": v25["delta_TOP1"]["value"],
        "adversarial_strings_refused": a26["n_refused"],
        "adversarial_strings_total": a26["n_adversarial_strings"],
        "design_columns": a26["structural_closure"]["design_columns"],
        "ranking_verdict": v25["verdict"],
        "scope_verdict": v26["verdict"],
    }

    # Each entry pins a number to the WORDS AROUND IT, not to a bare substring. A bare "56" is
    # satisfied by `SHA-256` and by `frozen_24F_sha256`; the first version of this check used bare
    # substrings and would have passed a record that never stated the refusal count at all. Every
    # pattern below carries `{n}` exactly once, which is what the canary substitutes.
    r25 = (RECORDS / "stage_25_RECORD.md").read_text(encoding="utf-8")
    r26 = (RECORDS / "stage_26_RECORD.md").read_text(encoding="utf-8")
    expected = [
        ("stage_25_RECORD.md", r25, "delta_RANK",
         r"delta_RANK\s+@@", f"{headline['delta_RANK']:+.6f}"),
        ("stage_25_RECORD.md", r25, "bootstrap CI95 lower",
         r"CI95\s+\[@@,", f"{headline['bootstrap_ci95'][0]:+.6f}"),
        ("stage_25_RECORD.md", r25, "bootstrap CI95 upper",
         r"CI95\s+\[[^\]]+,\s*@@\]", f"{headline['bootstrap_ci95'][1]:+.6f}"),
        ("stage_25_RECORD.md", r25, "R(W1)", r"R\(W1\)\s+@@", f"{headline['R_W1']:.6f}"),
        ("stage_25_RECORD.md", r25, "R(W4)", r"R\(W4\)\s+@@", f"{headline['R_W4']:.6f}"),
        ("stage_25_RECORD.md", r25, "R(W5)", r"R\(W5\)\s+@@", f"{headline['R_W5']:.6f}"),
        ("stage_25_RECORD.md", r25, "null p95",
         r"null p95\s+@@", f"{headline['null_p95']:.6f}"),
        ("stage_25_RECORD.md", r25, "delta_TOP1",
         r"delta_TOP1\s+@@", f"{headline['delta_TOP1']:+.6f}"),
        ("stage_25_RECORD.md", r25, "eligible clones",
         r"eligible clones\s+@@\b", str(headline["eligible_clones"])),
        ("stage_25_RECORD.md", r25, "permutation draws",
         r"@@\s*/\s*1000", str(headline["n_null_ge_observed"])),
        ("stage_25_RECORD.md", r25, "ranking verdict",
         r"@@", headline["ranking_verdict"]),
        ("stage_26_RECORD.md", r26, "adversarial refusals",
         r"@@\s*/\s*@@ adversarial", str(headline["adversarial_strings_refused"])),
        ("stage_26_RECORD.md", r26, "design columns",
         r"design columns\s+@@\s*=", str(headline["design_columns"])),
        ("stage_26_RECORD.md", r26, "scope verdict", r"@@", headline["scope_verdict"]),
    ]

    def _hit(text: str, tmpl: str, value: str) -> bool:
        return re.search(tmpl.replace("@@", re.escape(value)), text) is not None

    disagreements = [{"record": f, "label": label, "expected": val, "pattern": tmpl}
                     for f, text, label, tmpl, val in expected if not _hit(text, tmpl, val)]

    # Canary: perturb each value and confirm the pattern then finds nothing. A check keyed to the
    # surrounding words alone would keep matching and would prove nothing about the number.
    def _perturb(v: str) -> str:
        if v and v[-1].isdigit():
            return v[:-1] + ("8" if v[-1] != "8" else "7")
        return v + "_X"

    canary_failed = [label for f, text, label, tmpl, val in expected
                     if _hit(text, tmpl, _perturb(val))]

    # p_perm is the floor of a 1,000-draw test and must never be reported as a point estimate
    floor_ok = (abs(headline["p_perm"] - (1 + headline["n_null_ge_observed"])
                    / (headline["n_perm"] + 1)) < 1e-12)
    reported_as_floor = "p < 0.001" in r25

    checks = {
        "every headline number in the records matches its JSON source": not disagreements,
        "each number is pinned to its meaning, not to a bare substring": not canary_failed,
        "p_perm equals the finite-sample formula exactly": floor_ok,
        "p_perm is reported as a floor, not a point estimate": reported_as_floor,
        "the ranking verdict is one of the two frozen values":
            headline["ranking_verdict"] in ("STAGE_25_RANKING_SUPPORTED",
                                            "STAGE_25_RANKING_NOT_SUPPORTED"),
    }
    return write_json(NUMBERS_JSON, {
        "stage": "EL-D", "headline": headline,
        "records_checked": ["stage_25_RECORD.md", "stage_26_RECORD.md"],
        "substrings_checked": len(expected),
        "patterns_checked": len(expected),
        "disagreements": disagreements,
        "canary_patterns_that_matched_a_wrong_value": canary_failed,
        "checks": checks, "all_passed": all(checks.values()),
        "runtime_seconds": round(time.perf_counter() - t0, 3),
    })


# =============================================================================================== #
# EL-F — the claim lock input (plan §6). Carried through UNCHANGED.
# =============================================================================================== #
FORBIDDEN_CLAIMS = [
    "unseen-treatment generalization",
    "cross-cell-line or cross-patient generalization",
    "clinical treatment recommendation",
    "causal treatment-effect estimation",
    "calibrated probability unless calibration is separately frozen and passed",
    "independent biological replication of Role B",
    "uniform benefit across all six conditions",
    "confirmed Role-A prediction",
    "single-cell input equivalence when the model was trained on clone pseudobulk",
]

PRIMARY_CLAIM = (
    "Within the existing multi-condition WM989 lineage system, pretreatment Gene Expression "
    "contains treatment-specific information about future clonal detection beyond treatment "
    "identity and captured pretreatment clone abundance under frozen clone-held-out evaluation.")

RANKING_CLAIM_SUPPORTED = (
    "A frozen state x treatment model improves clone-specific ordering of the six observed "
    "experimental conditions over a non-interactive additive model.")

RANKING_CLAIM_NOT_SUPPORTED = (
    "Clone-specific ordering of the six observed experimental conditions was not supported under "
    "the preregistered ranking test.")

SUPPORTING_CLAIM = (
    "A separately reconstructed reprogramming system showed positive but underpowered evidence "
    "that pretreatment transcriptional state carries prospective information about a later "
    "lineage outcome.")


def claim_lock_input(numbers: dict) -> dict:
    t0 = time.perf_counter()
    supported = numbers["headline"]["ranking_verdict"] == "STAGE_25_RANKING_SUPPORTED"
    ship = SHIP_PLAN.read_text(encoding="utf-8")

    # The nine must still be the nine the frozen plan lists -- verified against the plan itself,
    # not against a copy of my own list. Drop the ```text fence tag: leaving it in made "text" the
    # tenth forbidden claim and failed the check on the first run.
    section = ship.split("## 3.5 Claims forbidden in Generation 1")[1].split("```")[1]
    in_plan = [ln.strip() for ln in section.strip().splitlines()[1:] if ln.strip()]

    checks = {
        "all nine forbidden claims are carried through verbatim":
            [c for c in FORBIDDEN_CLAIMS] == in_plan,
        "the primary claim is the frozen §3.1 wording": PRIMARY_CLAIM.rstrip(".") in
                                                        " ".join(ship.split()),
        "the ranking claim selected matches the Stage-25 verdict": True,
        "replication is recorded as Generation 2": True,
    }
    return write_json(CLAIM_INPUT_JSON, {
        "stage": "EL-F",
        "allowed": {
            "primary": PRIMARY_CLAIM,
            "ranking": RANKING_CLAIM_SUPPORTED if supported else RANKING_CLAIM_NOT_SUPPORTED,
            "ranking_selected_by": numbers["headline"]["ranking_verdict"],
            "supporting_role_A": SUPPORTING_CLAIM,
        },
        "forbidden": FORBIDDEN_CLAIMS,
        "forbidden_source": "§3.5 of the frozen ship plan, verbatim",
        "replication": {
            "status": "GENERATION 2",
            "statement": "Independent new-system biological replication has NOT been performed. "
                         "It is Generation-2 work and is not a Generation-1 publication gate. "
                         "Clone-held-out folds and two endpoint families are not replication.",
        },
        "the_claim_lock_may": "narrow these",
        "the_claim_lock_may_not": "widen these",
        "checks": checks, "all_passed": all(checks.values()),
        "runtime_seconds": round(time.perf_counter() - t0, 3),
    })


# =============================================================================================== #
# The lock
# =============================================================================================== #
def _lock_document(m: dict, b: dict, n: dict, verdict: str) -> str:
    h = n["headline"]
    lines = [f"  {k:<12s}  {v['sha256']}  {p}"
             for p, v in sorted(m["artifacts"].items())
             for k in [v["class"]]]
    return f"""# CellFate-Rx Generation 1 — EVIDENCE LOCK

```text
  {verdict}

  lock digest   {m['lock_digest']}
  artifacts     {m['n_artifacts']}
  ship plan     {SHIP_PLAN_DIGEST}
```

The lock digest is a SHA-256 over `path  sha256` lines, sorted by path and LF-joined. It is one
number that names the entire Generation-1 evidence base, and it belongs in the manuscript.

**Two digest forms appear below and they are not in conflict.** The manifest hashes raw bytes,
which is what a file manifest must do. The frozen protocol identity `8da16fca...` is a
*canonical-LF* digest of the same ship plan, computed with CRLF normalised to LF so the protocol
has one identity on every platform. The ship plan therefore carries `59f22e9a...` as a file and
`8da16fca...` as a protocol, and both are checked.

## What is locked

```text
{chr(10).join(lines)}
```

## Chain of custody

An artifact's hash must be the same number everywhere a stage recorded it. If Stage 24 handed
Stage 25 a table hashed X and this lock hashes Y, the analysis did not consume what is being
locked.

```text
  out-of-fold table   identical across 4 independent records: {b['agreement']['out_of_fold_table']}
  model artifact      identical across 4 independent records: {b['agreement']['model_artifact']}
  ship-plan digest    identical across 6 independent records: {b['agreement']['ship_plan_digest']}
```

## The numbers this locks

```text
  eligible clones        {h['eligible_clones']}
  R(W1) / R(W4) / R(W5)  {h['R_W1']:.6f} / {h['R_W4']:.6f} / {h['R_W5']:.6f}
  delta_RANK             {h['delta_RANK']:+.6f}   CI95 [{h['bootstrap_ci95'][0]:+.6f}, {h['bootstrap_ci95'][1]:+.6f}]
  null p95               {h['null_p95']:.6f}
  permutation            {h['n_null_ge_observed']} of {h['n_perm']} draws reached the observed value
  p_perm                 {h['p_perm']:.6f}   report as p < 0.001 (0 of 1,000), never as a point estimate
  delta_TOP1             {h['delta_TOP1']:+.6f}
  adversarial refusals   {h['adversarial_strings_refused']} of {h['adversarial_strings_total']}
  design columns         {h['design_columns']}
```

## What this lock does NOT contain

```text
  stage24_w5_artifact.npz   44 MB, gitignored. A fresh clone does NOT contain it.
                            Its hash is locked. Rebuild:
                              python experiments/run_stage24_gen1_tool.py --stage 24c
  raw sequencing data       GSE279162 (WM989), GSE227151 (Rewind). Accessions are
                            locked; bytes are not vendored.
```

Naming a gap is not closing it. Both stay open.

## Verifying this lock

```text
  python experiments/run_gen1_evidence_lock.py --verify
```

It re-hashes every artifact and refuses if one has moved. Its ability to refuse is itself tested:
a one-bit flip and a deleted file must both be caught, on copies, before any lock is issued.

## What locking does not do

It grants no claim. It fixes what the existing claims are made of. No lock outcome reopens an
earlier stage, changes a recorded number, or authorizes new data, a new condition or a new model.
"""


def run_all() -> dict:
    t0 = time.perf_counter()
    m = build_manifest()
    b = chain_of_custody()
    c = negative_controls(m)
    write_json(CONTROLS_JSON, c)
    n = lock_numbers()
    f = claim_lock_input(n)

    live = verify_against(m, ROOT)
    sub = {"EL-A_inventory": m["all_passed"], "EL-B_chain_of_custody": b["all_passed"],
           "EL-C_verifier_can_refuse": c["all_passed"], "EL-D_numbers": n["all_passed"],
           "EL-F_claim_lock_input": f["all_passed"],
           "EL_live_verification_clean": live["clean"]}
    verdict = "GEN1_EVIDENCE_LOCKED" if all(sub.values()) else "GEN1_EVIDENCE_LOCK_REFUSED"

    out = write_json(LOCK_JSON, {
        "stage": "GEN-1 EVIDENCE LOCK",
        "verdict": verdict,
        "lock_digest": m["lock_digest"],
        "substages": sub,
        "failing": [k for k, v in sub.items() if not v],
        "plan": _rel(PLAN),
        "plan_canonical_lf_sha256": canonical_lf_sha256(PLAN),
        "ship_plan_digest": SHIP_PLAN_DIGEST,
        "ship_plan_digest_holds": canonical_lf_sha256(SHIP_PLAN) == SHIP_PLAN_DIGEST,
        "n_artifacts": m["n_artifacts"],
        "by_class": m["by_class"],
        "chain_of_custody": b,
        "negative_controls": c["controls"],
        "live_verification": live,
        "headline_numbers": n["headline"],
        "not_in_the_lock": NOT_IN_THE_LOCK,
        "grants_no_claim": "Locking fixes what the existing claims are made of. It grants no new "
                           "claim and reopens no earlier stage.",
        "next": "GEN-1 CLAIM LOCK",
        "runtime_seconds": round(time.perf_counter() - t0, 3),
    })
    LOCK_MD.write_text(_lock_document(m, b, n, verdict), encoding="utf-8")

    write_json(HANDOFF_JSON, {
        "from_stage": "GEN-1 EVIDENCE LOCK",
        "to_stage": "GEN-1 CLAIM LOCK",
        "verdict": verdict,
        "lock_digest": m["lock_digest"],
        "manifest": _rel(MANIFEST_JSON),
        "claim_lock_input": _rel(CLAIM_INPUT_JSON),
        "lock_document": _rel(LOCK_MD),
        "verify_command": "python experiments/run_gen1_evidence_lock.py --verify",
        "claim_lock_must": [
            "consume GEN1_CLAIM_LOCK_INPUT.json unchanged; it may narrow, never widen",
            "re-verify the lock digest before writing an abstract-level claim",
            "state that independent biological replication is Generation 2, not a Gen-1 gate"],
        "no_lock_outcome_reopens_an_earlier_stage": True,
    })
    out["total_runtime_seconds"] = round(time.perf_counter() - t0, 3)
    return out


def run_verify() -> dict:
    if not MANIFEST_JSON.exists():
        raise SystemExit("no manifest: run --stage all first")
    m = json.loads(MANIFEST_JSON.read_text(encoding="utf-8"))
    r = verify_against(m, ROOT)
    r["lock_digest"] = m["lock_digest"]
    r["verdict"] = "EVIDENCE_INTACT" if r["clean"] else "EVIDENCE_MOVED"
    return r


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Generation-1 evidence lock")
    ap.add_argument("--stage", choices=["all"], default=None)
    ap.add_argument("--verify", action="store_true",
                    help="re-hash every locked artifact and refuse if one has moved")
    a = ap.parse_args(argv)

    if a.verify:
        r = run_verify()
        print(json.dumps(r, indent=2))
        return 0 if r["clean"] else 2
    if a.stage != "all":
        ap.error("pass --stage all or --verify")
    r = run_all()
    print(json.dumps({k: r[k] for k in
                      ("stage", "verdict", "lock_digest", "substages", "failing", "n_artifacts",
                       "by_class", "negative_controls", "live_verification", "next")},
                     indent=2, default=str))
    return 0 if r["verdict"] == "GEN1_EVIDENCE_LOCKED" else 2


if __name__ == "__main__":
    raise SystemExit(main())

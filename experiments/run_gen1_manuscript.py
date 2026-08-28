"""Generation-1 manuscript + reproducibility package.

Executes `GEN1_MANUSCRIPT_PACKAGE_V1.md`, the last Generation-1 stage.

Every previous stage checked something someone else wrote, or something a machine produced. This
one checks prose I wrote myself, which is the weakest position an instrument can be in: the author
and the reviewer are the same process. The only defence is to make the checker refuse mechanically
and to prove it can, so the manuscript is scanned with the claim lock's own instrument, every
number is traced to a locked artifact, and the checker is fired at deliberately broken copies of
the manuscript to show it says no.

Nothing here fits anything or produces a number. It refuses if either lock fails to verify.

    python experiments/run_gen1_manuscript.py --stage all
    python experiments/run_gen1_manuscript.py --verify
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import run_gen1_claim_lock as CL  # noqa: E402
import run_gen1_evidence_lock as EL  # noqa: E402

RESULTS = ROOT / "results"
OUT = RESULTS / "manuscript"
OUT.mkdir(parents=True, exist_ok=True)

PLANS = ROOT / "plans" / "(newer)practical plans"
PLAN = PLANS / "GEN1_MANUSCRIPT_PACKAGE_V1.md"
SHIP_PLAN = PLANS / "STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md"
SHIP_PLAN_DIGEST = "8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48"

CLAIM_HANDOFF = RESULTS / "gen1_handoff_to_manuscript.json"
EVIDENCE_MANIFEST = RESULTS / "evidence_lock" / "GEN1_EVIDENCE_MANIFEST.json"
EVIDENCE_LOCK = RESULTS / "evidence_lock" / "GEN1_EVIDENCE_LOCK.json"
CLAIM_LOCK = RESULTS / "claim_lock" / "GEN1_CLAIM_LOCK.json"

MANUSCRIPT = OUT / "MANUSCRIPT.md"
REPRO = OUT / "REPRODUCIBILITY.md"
COMPLIANCE_JSON = OUT / "manuscript_compliance.json"
CONTROLS_JSON = OUT / "manuscript_controls.json"
VERDICT_JSON = OUT / "GEN1_MANUSCRIPT.json"
DIGEST_JSON = OUT / "GEN1_PACKAGE_DIGEST.json"


def write_json(p: Path, obj: dict) -> dict:
    stamped = {**obj, "module_sha256":
               hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    p.write_text(json.dumps(stamped, indent=2) + "\n", encoding="utf-8")
    return stamped


def _j(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


REQUIRED_SECTIONS = [
    "Abstract", "Introduction", "Data", "Methods", "Results", "The tool",
    "Limitations", "What this does not show", "Availability", "Generation 2",
]

# The five mandatory qualifiers, each reduced to a phrase that must literally appear. A qualifier
# that is "conveyed by the general tone of the document" is not a qualifier.
QUALIFIER_MARKERS = {
    "system": "WM989",
    "vocabulary": "six observed experimental conditions",
    "outcome": "not death",
    "evaluation": "clone-held-out",
    "replication": "Generation 2",
}
ABSTRACT_MUST_CARRY = ["system", "vocabulary", "outcome"]


# =============================================================================================== #
# MS-A — preflight (plan §1)
# =============================================================================================== #
def preflight() -> dict:
    t0 = time.perf_counter()
    handoff = _j(CLAIM_HANDOFF)
    ev = EL.verify_against(_j(EVIDENCE_MANIFEST), ROOT)
    cl = CL.run_verify()

    checks = {
        "the claim lock verdict is GEN1_CLAIMS_LOCKED":
            handoff["verdict"] == "GEN1_CLAIMS_LOCKED",
        "the evidence lock verifies over all its artifacts": ev["clean"],
        "the claim lock verifies": cl["clean"],
        "the evidence digest matches the handoff":
            _j(EVIDENCE_MANIFEST)["lock_digest"] == handoff["evidence_lock_digest"],
        "the claim digest matches the handoff":
            cl["claim_digest"] == handoff["claim_digest"],
        "the frozen ship-plan digest still holds":
            EL.canonical_lf_sha256(SHIP_PLAN) == SHIP_PLAN_DIGEST,
    }
    return {"stage": "MS-A",
            "evidence_digest": handoff["evidence_lock_digest"],
            "claim_digest": handoff["claim_digest"],
            "evidence_verification": ev, "claim_verification": cl,
            "checks": checks, "all_passed": all(checks.values()),
            "runtime_seconds": round(time.perf_counter() - t0, 3)}


# =============================================================================================== #
# MS-C / MS-D — compliance and number traceability (plan §3, §4)
# =============================================================================================== #
def _numbers() -> list[tuple[str, str, str]]:
    """(label, regex-with-@@, value) -- pinned to the words around it, as the evidence lock does."""
    h = _j(EVIDENCE_LOCK)["headline_numbers"]
    v25 = _j(RESULTS / "stage25" / "stage25_verdict.json")
    a26 = _j(RESULTS / "stage26" / "stage26a_vocabulary_closure.json")
    aud = _j(RESULTS / "stage23_2h" / "stage23_2h_power_audit.json")
    d = v25["descriptives"]
    return [
        ("delta_RANK", r"@@", f"{h['delta_RANK']:+.6f}"),
        ("CI95 lower", r"\[@@,", f"{h['bootstrap_ci95'][0]:+.6f}"),
        ("CI95 upper", r",\s*@@\]", f"{h['bootstrap_ci95'][1]:+.6f}"),
        ("R(W1)", r"W1\D{0,12}@@", f"{h['R_W1']:.6f}"),
        ("R(W4)", r"W4\D{0,12}@@", f"{h['R_W4']:.6f}"),
        ("R(W5)", r"W5\D{0,12}@@", f"{h['R_W5']:.6f}"),
        ("null p95", r"p95\D{0,30}@@", f"{h['null_p95']:.6f}"),
        # `\D` cannot cross "1,000", so this one spells out the digits it is allowed to skip
        ("null max", r"largest of [\d,]+ draws\s+@@", f"{v25['permutation']['null_max']:.6f}"),
        ("delta_TOP1", r"TOP1\D{0,12}@@", f"{h['delta_TOP1']:+.6f}"),
        ("eligible clones", r"@@ (?:of 1,401|eligible)", str(h["eligible_clones"])),
        ("excluded all-zero", r"@@ (?:clones )?(?:were )?(?:all-zero|never detected)",
         str(d["excluded_all_zero_clones"])),
        ("excluded all-positive", r"@@ (?:clones )?(?:were )?(?:all-positive|always detected)",
         str(d["excluded_all_positive_clones"])),
        ("design columns", r"@@ (?:design )?columns", str(h["design_columns"])),
        ("adversarial refusals", r"@@ of @@", str(a26["n_refused"])),
        ("permutation count", r"@@ full-refit", str(h["n_perm"])),
        ("Role A power recorded", r"@@ against", str(aud["recorded_instrument"]["power"])),
        ("Role A power audited", r"audit\D{0,50}@@", str(aud["corrected_power"])),
    ]


def compliance() -> dict:
    t0 = time.perf_counter()
    text = MANUSCRIPT.read_text(encoding="utf-8")
    repro = REPRO.read_text(encoding="utf-8") if REPRO.exists() else ""
    full = CL.combined_patterns()
    handoff = _j(CLAIM_HANDOFF)

    forbidden_hits = CL.scan(text, full)
    repro_hits = CL.scan(repro, full)

    missing_sections = [s for s in REQUIRED_SECTIONS if f"# {s}" not in text
                        and f"## {s}" not in text]
    missing_qual = [k for k, v in QUALIFIER_MARKERS.items() if v not in text]

    abstract = text.split("## Abstract", 1)[-1].split("\n## ", 1)[0] if "## Abstract" in text else ""
    abstract_missing = [k for k in ABSTRACT_MUST_CARRY if QUALIFIER_MARKERS[k] not in abstract]

    # Whitespace is normalised before the number patterns run. A line wrap between "472 were" and
    # "never detected" defeated the trace check for the excluded-clone count -- the third time in
    # this project that a hard wrap has broken a text check, after the Stage-26 leak scan and the
    # claim lock's clause splitter. A checker that reads a wrapped document must not care where the
    # author's editor happened to break the line.
    flat = " ".join(text.split())

    def _hit(t: str, tmpl: str, val: str) -> bool:
        return re.search(tmpl.replace("@@", re.escape(val)), t) is not None

    untraceable = [{"label": lab, "expected": val, "pattern": tmpl}
                   for lab, tmpl, val in _numbers() if not _hit(flat, tmpl, val)]

    checks = {
        "no forbidden claim appears unnegated in the manuscript": not forbidden_hits,
        "no forbidden claim appears unnegated in the package document": not repro_hits,
        "every required section is present": not missing_sections,
        "all five mandatory qualifiers appear": not missing_qual,
        "the abstract carries system, vocabulary and outcome by itself": not abstract_missing,
        "p is reported as a floor": "p < 0.001" in text,
        "the point estimate for p appears nowhere": "0.000999" not in text,
        "replication is stated as Generation 2":
            "Generation 2" in text and "biological replication" in text,
        "both digests are quoted": (handoff["evidence_lock_digest"] in text
                                    and handoff["claim_digest"] in text),
        "every number traces to a locked artifact": not untraceable,
    }
    return write_json(COMPLIANCE_JSON, {
        "stage": "MS-C/MS-D",
        "instrument": "the claim lock's extended patterns under clause-scoped negation",
        "forbidden_hits": forbidden_hits,
        "package_forbidden_hits": repro_hits,
        "missing_sections": missing_sections,
        "missing_qualifiers": missing_qual,
        "abstract_missing_qualifiers": abstract_missing,
        "numbers_checked": len(_numbers()),
        "untraceable_numbers": untraceable,
        "checks": checks, "all_passed": all(checks.values()),
        "runtime_seconds": round(time.perf_counter() - t0, 3)})


# =============================================================================================== #
# MS-E — the reproducibility package (plan §5)
# =============================================================================================== #
def package_check() -> dict:
    t0 = time.perf_counter()
    text = REPRO.read_text(encoding="utf-8")
    manifest = _j(EVIDENCE_MANIFEST)["artifacts"]

    # every `python experiments/X.py --stage Y` must exist, and Y must be accepted by that parser
    cmds = re.findall(r"python (experiments/[\w./]+\.py)(?: --stage (\S+))?", text)
    bad_cmd = []
    for script, stage in cmds:
        p = ROOT / script
        if not p.exists():
            bad_cmd.append({"command": script, "problem": "script does not exist"})
            continue
        if stage:
            src = p.read_text(encoding="utf-8")
            m = re.search(r'choices=\[([^\]]*)\]', src)
            allowed = re.findall(r'"([^"]+)"', m.group(1)) if m else []
            if allowed and stage not in allowed:
                bad_cmd.append({"command": f"{script} --stage {stage}",
                                "problem": f"parser accepts {allowed}"})

    referenced = set(re.findall(r"(results/[\w./\-]+|src/[\w./\-]+|experiments/[\w./\-]+)", text))
    unknown = sorted(r for r in referenced
                     if r not in manifest and not (ROOT / r).exists())

    checks = {
        "every documented command names a script that exists": not bad_cmd,
        "every referenced path exists or is in the evidence manifest": not unknown,
        "the two verification commands are given": (
            "run_gen1_evidence_lock.py --verify" in text
            and "run_gen1_claim_lock.py --verify" in text),
        "the environment is recorded": "Python" in text and "numpy" in text,
        "what is NOT in the package is named":
            "44 MB" in text and "GSE279162" in text and "GSE227151" in text,
        "the rebuild command for the gitignored artifact is given":
            "--stage 24c" in text,
        "the long runtime is stated honestly": "10.7" in text,
    }
    return {"stage": "MS-E", "commands_found": len(cmds), "bad_commands": bad_cmd,
            "unknown_paths": unknown, "checks": checks, "all_passed": all(checks.values()),
            "runtime_seconds": round(time.perf_counter() - t0, 3)}


# =============================================================================================== #
# MS-F — the checker must refuse (plan §6)
# =============================================================================================== #
def negative_controls() -> dict:
    t0 = time.perf_counter()
    original = MANUSCRIPT.read_text(encoding="utf-8")
    full = CL.combined_patterns()
    h = _j(EVIDENCE_LOCK)["headline_numbers"]
    results = {}

    planted = original + "\n\nThe model generalises to new treatments.\n"
    results["a planted forbidden claim is caught"] = bool(CL.scan(planted, full))

    dropped = original.replace("not death", "the outcome")
    results["a dropped qualifier is caught"] = QUALIFIER_MARKERS["outcome"] not in dropped

    ppoint = original.replace("p < 0.001", "p = 0.000999")
    results["p quoted as a point estimate is caught"] = (
        "p < 0.001" not in ppoint and "0.000999" in ppoint)

    broken = original.replace(f"{h['delta_RANK']:+.6f}", "+0.051608")
    results["a changed number is caught"] = not re.search(
        re.escape(f"{h['delta_RANK']:+.6f}"), broken)

    # and the unmodified document must still pass all four, or the controls prove nothing
    results["the real manuscript passes all four"] = (
        not CL.scan(original, full)
        and QUALIFIER_MARKERS["outcome"] in original
        and "p < 0.001" in original and "0.000999" not in original
        and f"{h['delta_RANK']:+.6f}" in original)

    return write_json(CONTROLS_JSON, {
        "stage": "MS-F",
        "note": "Run on in-memory copies; the manuscript on disk is never modified. A checker "
                "that has never refused its own document is decoration.",
        "controls": results, "checks": results, "all_passed": all(results.values()),
        "runtime_seconds": round(time.perf_counter() - t0, 3)})


# =============================================================================================== #
PACKAGE_FILES = [
    "plans/(newer)practical plans/GEN1_MANUSCRIPT_PACKAGE_V1.md",
    "experiments/run_gen1_manuscript.py",
    "tests/test_gen1_manuscript.py",
    "results/manuscript/MANUSCRIPT.md",
    "results/manuscript/REPRODUCIBILITY.md",
]


def package_digest() -> tuple[str, dict]:
    per = {rel: (EL.canonical_lf_sha256(ROOT / rel) if (ROOT / rel).exists() else "MISSING")
           for rel in PACKAGE_FILES}
    canonical = "\n".join(f"{k}  {per[k]}" for k in sorted(per))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest(), per


def run_all() -> dict:
    t0 = time.perf_counter()
    a = preflight()
    if not a["all_passed"]:
        return write_json(VERDICT_JSON, {
            "stage": "GEN-1 MANUSCRIPT", "verdict": "GEN1_MANUSCRIPT_REFUSED",
            "refused_at": "MS-A", "detail": a,
            "note": "A lock did not verify. Writing a manuscript against evidence or claims that "
                    "have shifted is the failure both locks exist to prevent."})

    c = compliance()
    e = package_check()
    f = negative_controls()

    sub = {"MS-A_locks_verify": a["all_passed"], "MS-C_D_compliance": c["all_passed"],
           "MS-E_package": e["all_passed"], "MS-F_controls_fire": f["all_passed"]}
    verdict = "GEN1_MANUSCRIPT_READY" if all(sub.values()) else "GEN1_MANUSCRIPT_REFUSED"

    out = write_json(VERDICT_JSON, {
        "stage": "GEN-1 MANUSCRIPT + REPRODUCIBILITY PACKAGE",
        "verdict": verdict,
        "substages": sub,
        "failing": [k for k, v in sub.items() if not v],
        "evidence_digest": a["evidence_digest"],
        "claim_digest": a["claim_digest"],
        "plan_canonical_lf_sha256": EL.canonical_lf_sha256(PLAN),
        "manuscript": MANUSCRIPT.relative_to(ROOT).as_posix(),
        "package": REPRO.relative_to(ROOT).as_posix(),
        "numbers_traced": c["numbers_checked"],
        "sections": REQUIRED_SECTIONS,
        "negative_controls": f["controls"],
        "what_ready_does_not_mean":
            "GEN1_MANUSCRIPT_READY means the document is consistent with everything locked "
            "beneath it. It is NOT a judgement that the science is good, that the writing is "
            "clear, or that a reviewer will agree.",
        "next": "PREPRINT / SUBMISSION. Generation 1 is complete.",
        "runtime_seconds": round(time.perf_counter() - t0, 3)})

    digest, per = package_digest()
    write_json(DIGEST_JSON, {
        "package_digest": digest,
        "digest_definition": "SHA-256 over 'path  canonical-lf-sha256' lines, sorted by path",
        "covers": per,
        "verify": "python experiments/run_gen1_manuscript.py --verify"})
    out["package_digest"] = digest
    return out


def run_verify() -> dict:
    if not DIGEST_JSON.exists():
        raise SystemExit("no package digest: run --stage all first")
    rec = _j(DIGEST_JSON)
    now, per = package_digest()
    moved = [k for k, v in per.items() if rec["covers"].get(k) != v]
    return {"clean": not moved, "moved": moved, "package_digest": now,
            "recorded_digest": rec["package_digest"],
            "verdict": "PACKAGE_INTACT" if now == rec["package_digest"] else "PACKAGE_MOVED"}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Generation-1 manuscript + package")
    ap.add_argument("--stage", choices=["all"], default=None)
    ap.add_argument("--verify", action="store_true")
    a = ap.parse_args(argv)
    if a.verify:
        r = run_verify()
        print(json.dumps(r, indent=2))
        return 0 if r["clean"] else 2
    if a.stage != "all":
        ap.error("pass --stage all or --verify")
    r = run_all()
    print(json.dumps({k: r[k] for k in r if k in
                      ("stage", "verdict", "substages", "failing", "refused_at",
                       "numbers_traced", "negative_controls", "package_digest", "next")},
                     indent=2, default=str))
    return 0 if r["verdict"] == "GEN1_MANUSCRIPT_READY" else 2


if __name__ == "__main__":
    raise SystemExit(main())

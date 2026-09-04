"""Re-run the three Generation-1 locks in the only order that can succeed, and re-pin the
documents that quote their digests.

The order is not a preference, it is forced by what each layer hashes:

    export_gen1_source_data.py   writes the per-draw source data and the environment lock
    evidence lock                hashes those, and hashes its own module
    GEN1_CLAIM_LOCK_V1.md        pins the evidence digest, and the claim lock hashes that plan
    claim lock                   produces the claim digest
    MANUSCRIPT / README / ...    quote both digests, and the manuscript hashes those documents
    manuscript                   produces the package digest

Run them out of order and a lock refuses against a digest that has already moved. That is the whole
reason this script exists: the sequence was being carried out by hand, and re-pinning four
documents by hand is exactly the "maintained beside the primary rather than checked against it"
failure that the release-verification pass spent its time removing. It went wrong once already -- a
cascade refused at MS-C/MS-D because the README had not been re-pinned.

The exporter is deliberately NOT run here. It is slow, it is the one step that touches the derived
source data, and regenerating that should be an explicit decision rather than a side effect of
re-pinning some digests. Run it first, yourself, when the numbers have actually changed:

    python experiments/export_gen1_source_data.py

The bundle is likewise not built here. It must be built from a COMMITTED tree, because it records
the commit it was cut from and refuses a dirty one. So the release sequence is:

    python experiments/cascade_gen1.py     # re-run the locks, re-pin the documents
    git add -A && git commit               # the bundle will record this commit
    python experiments/make_release_bundle.py

`--check` re-pins nothing and runs no stage. It answers one question: does every document that
quotes a digest quote the current one? Use it to find out whether a cascade is needed at all.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLANS = ROOT / "plans" / "(newer)practical plans"

EVIDENCE_MANIFEST = ROOT / "results" / "evidence_lock" / "GEN1_EVIDENCE_MANIFEST.json"
CLAIM_DIGEST = ROOT / "results" / "claim_lock" / "GEN1_CLAIM_DIGEST.json"
PACKAGE_DIGEST = ROOT / "results" / "manuscript" / "GEN1_PACKAGE_DIGEST.json"

# The claim lock's plan pins the evidence digest ONLY. It is written before a claim digest exists,
# and CL-A refuses if the plan does not name the evidence it was written against.
EVIDENCE_ONLY = [PLANS / "GEN1_CLAIM_LOCK_V1.md"]

# Everything that quotes both. The manuscript plan is included because it pins both in its Entry
# block and nothing else checks it -- it was found stale during the verification pass.
BOTH = [
    ROOT / "results" / "manuscript" / "MANUSCRIPT.md",
    ROOT / "results" / "manuscript" / "REPRODUCIBILITY.md",
    ROOT / "README.md",
    PLANS / "GEN1_MANUSCRIPT_PACKAGE_V1.md",
]

EVIDENCE_PIN = re.compile(r"(evidence(?: lock)?(?: digest)?\s+`?)[a-f0-9]{64}")
CLAIM_PIN = re.compile(r"(claim(?: lock)?(?: digest)?\s+`?)[a-f0-9]{64}")
PLAN_PIN = re.compile(r"lock digest `([a-f0-9]{64})`")


def _run(script: str) -> int:
    """Run one stage with this interpreter, reporting its verdict line and nothing else."""
    proc = subprocess.run([sys.executable, str(ROOT / "experiments" / script), "--stage", "all"],
                          cwd=ROOT, capture_output=True, text=True)
    verdict = ""
    for line in proc.stdout.splitlines():
        if '"verdict"' in line:
            verdict = line.strip().rstrip(",")
            break
    print(f"    {script:32s} exit={proc.returncode}  {verdict}")
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout[-2000:] + "\n" + proc.stderr[-2000:] + "\n")
    return proc.returncode


def _digest(path: Path, key: str) -> str:
    return json.loads(path.read_text(encoding="utf-8"))[key]


def digests() -> tuple[str, str, str]:
    return (_digest(EVIDENCE_MANIFEST, "lock_digest"),
            _digest(CLAIM_DIGEST, "claim_digest"),
            _digest(PACKAGE_DIGEST, "package_digest"))


def _repin_plan(evidence: str) -> bool:
    """The claim-lock plan names the evidence digest in prose as well as in its Entry block, so
    every occurrence of the currently-pinned value moves together."""
    p = EVIDENCE_ONLY[0]
    text = p.read_text(encoding="utf-8")
    m = PLAN_PIN.search(text)
    if not m:
        raise SystemExit(f"{p.name} no longer pins a digest in its Entry block")
    current = m.group(1)
    if current == evidence:
        return False
    p.write_text(text.replace(current, evidence), encoding="utf-8")
    return True


def _repin(path: Path, evidence: str, claim: str) -> bool:
    text = path.read_text(encoding="utf-8")
    new = CLAIM_PIN.sub(lambda m: m.group(1) + claim,
                        EVIDENCE_PIN.sub(lambda m: m.group(1) + evidence, text))
    if new == text:
        return False
    path.write_text(new, encoding="utf-8")
    return True


def check() -> int:
    """Does every document that quotes a digest quote the current one?"""
    evidence, claim, package = digests()
    stale: list[str] = []
    for p in EVIDENCE_ONLY:
        if evidence not in p.read_text(encoding="utf-8"):
            stale.append(f"{p.name}: does not name the current evidence digest")
    for p in BOTH:
        text = p.read_text(encoding="utf-8")
        if evidence not in text:
            stale.append(f"{p.name}: does not name the current evidence digest")
        if claim not in text:
            stale.append(f"{p.name}: does not name the current claim digest")
    print(json.dumps({
        "evidence_digest": evidence,
        "claim_digest": claim,
        "package_digest": package,
        "documents_checked": [p.name for p in EVIDENCE_ONLY + BOTH],
        "stale": stale,
        "verdict": "PINS_CURRENT" if not stale else "PINS_STALE",
    }, indent=2))
    return 0 if not stale else 2


def cascade() -> int:
    print("  1. evidence lock")
    if _run("run_gen1_evidence_lock.py"):
        print("  REFUSED at the evidence lock; nothing downstream was run.")
        return 2
    evidence = _digest(EVIDENCE_MANIFEST, "lock_digest")

    print("  2. re-pin the claim-lock plan")
    moved = "re-pinned" if _repin_plan(evidence) else "already current"
    print(f"    {'GEN1_CLAIM_LOCK_V1.md':32s} {moved}  {evidence[:12]}")

    print("  3. claim lock")
    if _run("run_gen1_claim_lock.py"):
        print("  REFUSED at the claim lock. A CL-A failure means the plan pin and the evidence "
              "digest disagree, which this script has just set -- so read the refusal rather "
              "than re-running.")
        return 2
    claim = _digest(CLAIM_DIGEST, "claim_digest")

    print("  4. re-pin the documents that quote both digests")
    for p in BOTH:
        print(f"    {p.name:32s} {'re-pinned' if _repin(p, evidence, claim) else 'already current'}")

    print("  5. manuscript + package")
    if _run("run_gen1_manuscript.py"):
        print("  REFUSED at the manuscript. Its `failing` list names the substage; MS-C/MS-D "
              "covers the digest quotes, the submitted abstract and the secondary documents.")
        return 2

    print()
    if check():
        print("  A document is still stale AFTER cascading, which means it quotes a digest in a "
              "form these patterns do not recognise. Fix that pin by hand and re-run --check.")
        return 2
    print()
    print("  Next, in this order. The bundle records the commit it is cut from and refuses a "
          "dirty tree:")
    print("    git add -A && git commit")
    print("    python experiments/make_release_bundle.py")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Cascade the Generation-1 locks and re-pin documents")
    ap.add_argument("--check", action="store_true",
                    help="run no stage; report whether every quoted digest is current")
    a = ap.parse_args(argv)
    return check() if a.check else cascade()


if __name__ == "__main__":
    raise SystemExit(main())

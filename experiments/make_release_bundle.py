"""Assemble the Zenodo / GitHub-release bundle for Generation 1.

The repository deliberately gitignores the 44 MB model artifact, so a clone cannot verify the
evidence lock without rebuilding it. An archive has no such excuse: this bundle **includes** the
artifact, which is the one thing the repository cannot carry and the archive must.

Output goes to `dist/`, which is gitignored — a 44 MB zip does not belong in git history. The bundle
is deterministic apart from the zip's own timestamps, and every member is listed with its SHA-256 so
a downloader can check the contents rather than trusting the container.

    python experiments/make_release_bundle.py            # build + checksums
    python experiments/make_release_bundle.py --check     # re-verify an existing bundle
"""

from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIST = ROOT / "dist"
RESULTS = ROOT / "results"

BUNDLE = DIST / "cellfate-rx-gen1-bundle.zip"
SUMS = DIST / "SHA256SUMS.txt"
CONTENTS = DIST / "BUNDLE_CONTENTS.json"

# The artifact the repository cannot carry. Present here, absent from git, and the reason this
# bundle exists at all rather than pointing people at the GitHub tarball.
GITIGNORED_BUT_REQUIRED = ["results/stage24/stage24_w5_artifact.npz"]

# Directories taken whole. Everything under them is evidence, a record, or the tooling that checks
# both; nothing here is scratch.
TREES = [
    "src/cellfate",
    "experiments",
    "tests",
    "plans/(newer)practical plans",
    "results/evidence_lock",
    "results/claim_lock",
    "results/manuscript",
    "results/stage24/tool",
    "results/stage25",
    "results/stage26",
]

FILES = [
    "README.md",
    "ARCHITECTURE.md",
    "CITATION.cff",
    "pyproject.toml",
    ".gitattributes",
    ".github/workflows/ci.yml",
    "results/stage22_wm989_clones.csv",
    "results/stage23_wm989_detection_oof.csv",
    "results/stage23_wm989_interaction_oof.csv",
    "results/stage24/stage24_oof_for_stage25.csv",
    "results/stage24/stage24_w5_artifact.json",
    "results/stage24/stage24f_tool_freeze.json",
    "results/stage23_5_protocol.json",
    "results/stage23_5_handoff_to_stage24.json",
    "results/stage24_handoff_to_stage25.json",
    "results/stage26_handoff_to_evidence_lock.json",
    "results/gen1_handoff_to_claim_lock.json",
    "results/gen1_handoff_to_manuscript.json",
    "results/stage23_2h/stage23_2h_verdict.json",
    "results/stage23_2h/stage23_2h_confirmation.json",
    "results/stage23_2h/stage23_2h_power.json",
    "results/stage23_2h/stage23_2h_power_audit.json",
    "results/stage23_2/stage23_2_handoff_to_stage24.json",
]

SKIP_SUFFIX = {".pyc"}
SKIP_DIR = {"__pycache__", ".pytest_cache", ".cache"}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _members() -> list[str]:
    seen: dict[str, None] = {}
    for rel in FILES + GITIGNORED_BUT_REQUIRED:
        if (ROOT / rel).is_file():
            seen[rel] = None
    for tree in TREES:
        base = ROOT / tree
        if not base.is_dir():
            continue
        for p in sorted(base.rglob("*")):
            if not p.is_file() or p.suffix in SKIP_SUFFIX:
                continue
            if any(part in SKIP_DIR for part in p.relative_to(ROOT).parts):
                continue
            seen[p.relative_to(ROOT).as_posix()] = None
    return sorted(seen)


def build() -> int:
    DIST.mkdir(exist_ok=True)
    members = _members()

    missing = [m for m in GITIGNORED_BUT_REQUIRED if not (ROOT / m).is_file()]
    if missing:
        print("REFUSED — the archive must carry what the repository cannot:")
        for m in missing:
            print(f"  {m}")
        print("\nRebuild it first:\n  python experiments/run_stage24_gen1_tool.py --stage 24c")
        return 2

    digests = {rel: sha256(ROOT / rel) for rel in members}
    total = sum((ROOT / m).stat().st_size for m in members)

    with zipfile.ZipFile(BUNDLE, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for rel in members:
            z.write(ROOT / rel, arcname=f"cellfate-rx-gen1/{rel}")

    SUMS.write_text("".join(f"{digests[r]}  {r}\n" for r in members), encoding="utf-8")

    locks = {}
    for name, path, key in [
        ("evidence", "results/evidence_lock/GEN1_EVIDENCE_MANIFEST.json", "lock_digest"),
        ("claim", "results/claim_lock/GEN1_CLAIM_DIGEST.json", "claim_digest"),
        ("package", "results/manuscript/GEN1_PACKAGE_DIGEST.json", "package_digest"),
    ]:
        p = ROOT / path
        if p.is_file():
            locks[name] = json.loads(p.read_text(encoding="utf-8"))[key]

    CONTENTS.write_text(json.dumps({
        "bundle": BUNDLE.name,
        "bundle_sha256": sha256(BUNDLE),
        "n_files": len(members),
        "uncompressed_bytes": total,
        "compressed_bytes": BUNDLE.stat().st_size,
        "lock_digests": locks,
        "includes_gitignored": GITIGNORED_BUT_REQUIRED,
        "verify": [
            "python experiments/run_gen1_evidence_lock.py --verify",
            "python experiments/run_gen1_claim_lock.py --verify",
            "python experiments/run_gen1_manuscript.py --verify",
        ],
        "note": "Unpack, then run the three verify commands from the unpacked root. The model "
                "artifact is included here precisely because git does not carry it, so the "
                "archive verifies without a rebuild.",
        "files": digests,
    }, indent=2) + "\n", encoding="utf-8")

    print(f"  bundle       {BUNDLE.relative_to(ROOT).as_posix()}")
    print(f"  files        {len(members)}")
    print(f"  uncompressed {total / 1e6:.1f} MB")
    print(f"  compressed   {BUNDLE.stat().st_size / 1e6:.1f} MB")
    print(f"  sha256       {sha256(BUNDLE)}")
    for k, v in locks.items():
        print(f"  {k + ' digest':<13s}{v}")
    print(f"\n  checksums    {SUMS.relative_to(ROOT).as_posix()}")
    return 0


def check() -> int:
    if not SUMS.is_file() or not BUNDLE.is_file():
        print("no bundle to check — run without --check first")
        return 2
    recorded = {}
    for line in SUMS.read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, rel = line.split("  ", 1)
            recorded[rel] = digest

    bad, gone = [], []
    for rel, digest in recorded.items():
        p = ROOT / rel
        if not p.is_file():
            gone.append(rel)
        elif sha256(p) != digest:
            bad.append(rel)

    with zipfile.ZipFile(BUNDLE) as z:
        inside = {n[len("cellfate-rx-gen1/"):] for n in z.namelist()
                  if n.startswith("cellfate-rx-gen1/")}
    not_archived = sorted(set(recorded) - inside)

    ok = not bad and not gone and not not_archived
    print(json.dumps({
        "checked": len(recorded),
        "changed_since_build": sorted(bad),
        "missing_from_working_tree": sorted(gone),
        "listed_but_not_in_archive": not_archived,
        "bundle_sha256": sha256(BUNDLE),
        "verdict": "BUNDLE_INTACT" if ok else "BUNDLE_MISMATCH",
    }, indent=2))
    return 0 if ok else 2


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build the Generation-1 release bundle")
    ap.add_argument("--check", action="store_true",
                    help="re-verify an existing bundle against its checksums")
    a = ap.parse_args(argv)
    return check() if a.check else build()


if __name__ == "__main__":
    raise SystemExit(main())

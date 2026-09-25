"""Every reference and every table must be cited in the body of the manuscript, in order.

BMC Bioinformatics' technical check returned the first submission (2026-09-25) for exactly this:
Tables 2-5 had captions but no in-text citation, and references 7-13 were cited only under
"Availability of data and materials". That section is replaced in the published article by the
statement entered in the submission form, which carries no reference numbers, so those references
would have been cited nowhere. The body is therefore everything before "## Declarations".

This file is deliberately outside every digest, so tightening it never moves a lock.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "results" / "manuscript" / "MANUSCRIPT.md"

# [3], [3, 5], [3-5] or [3–5]. At most three digits, and not followed by "(", so an ORCID link such
# as [0009-0004-4503-6629](https://...) is never read as a citation.
CITATION = re.compile(r"\[(\d{1,3}(?:\s*[,–-]\s*\d{1,3})*)\](?!\()")
CAPTION = re.compile(r"^\*\*Table (\d+)\*\*")


def _text() -> str:
    return MANUSCRIPT.read_text(encoding="utf-8")


def _body() -> str:
    return _text().split("\n## Declarations", 1)[0]


def _reference_numbers() -> list[int]:
    refs = _text().split("\n## References", 1)[1]
    return [int(m.group(1)) for m in re.finditer(r"^(\d+)\.\s", refs, flags=re.M)]


def _cited_in_order(text: str) -> list[int]:
    """Reference numbers in order of first citation, with ranges expanded."""
    seen: list[int] = []
    for m in CITATION.finditer(text):
        parts = re.split(r"\s*([,–-])\s*", m.group(1))
        nums = [int(parts[0])]
        for sep, n in zip(parts[1::2], parts[2::2]):
            nums += list(range(nums[-1] + 1, int(n) + 1)) if sep != "," else [int(n)]
        seen += [n for n in nums if n not in seen]
    return seen


def test_every_reference_is_cited_in_the_body():
    listed = _reference_numbers()
    assert listed == list(range(1, len(listed) + 1)), "reference list is not numbered 1..N"
    cited = set(_cited_in_order(_body()))
    missing = [n for n in listed if n not in cited]
    assert not missing, f"references cited only outside the body, or not at all: {missing}"


def test_references_are_numbered_in_order_of_first_citation():
    order = _cited_in_order(_body())
    assert order == sorted(order), f"first citations out of order (Vancouver): {order}"


def test_every_table_is_cited_in_the_prose_in_order():
    lines = _body().splitlines()
    captions = [int(m.group(1)) for l in lines if (m := CAPTION.match(l))]
    assert captions == list(range(1, len(captions) + 1)), f"tables not numbered 1..N: {captions}"
    first: dict[int, int] = {}
    for i, line in enumerate(lines):
        if CAPTION.match(line):
            continue                       # a caption, or a mention inside one, is not a citation
        for m in re.finditer(r"\bTable (\d+)\b", line):
            first.setdefault(int(m.group(1)), i)
    uncited = [n for n in captions if n not in first]
    assert not uncited, f"tables with a caption but no in-text citation: {uncited}"
    order = sorted(captions, key=first.get)
    assert order == captions, f"tables first cited out of order: {order}"
    late = [n for n in captions
            if first[n] > next(i for i, l in enumerate(lines) if CAPTION.match(l)
                               and int(CAPTION.match(l).group(1)) == n)]
    assert not late, f"tables whose first citation comes after the table itself: {late}"

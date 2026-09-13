"""Render the Generation-1 submission files from the manuscript, without ever editing it.

    python experiments/render_gen1_submission.py           # refuses while any FILL marker remains
    python experiments/render_gen1_submission.py --draft   # placeholders allowed; writes to draft/

Outputs, in results/manuscript/submission/ (or its draft/ subfolder):

    figure_1_design.pdf, figure_2_primary.pdf, figure_3_robustness.pdf   separate figure files
    MANUSCRIPT_BMC.docx       the manuscript as written: legends in the text, figures separate
    MANUSCRIPT_bioRxiv.docx   the same text with the three figures embedded above their legends
    MANUSCRIPT_bioRxiv.pdf    that .docx exported through Microsoft Word
    COVER_LETTER.docx         section 6 of SUBMISSION.md
    RENDER_MANIFEST.json      SHA-256 of every output, and the tools that made them

Why two manuscript files. BMC wants figures as separate files with the legends in the text; bioRxiv
wants one PDF with the figures in it. The Markdown manuscript is the single source for both.

Why Word makes the PDF. Pandoc writes PDF through LaTeX by default, and no LaTeX is installed here.
Word is -- Microsoft 365, which renders SVG -- so the figures stay vector all the way to the PDF.

None of these outputs is byte-deterministic: .docx and .pdf carry timestamps. They are therefore not
locked. The manifest's hashes identify exactly which files were submitted.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_DIR = ROOT / "results" / "manuscript"
MANUSCRIPT = MANUSCRIPT_DIR / "MANUSCRIPT.md"
SUBMIT = MANUSCRIPT_DIR / "SUBMISSION.md"
FIGURE_DIR = MANUSCRIPT_DIR / "figures"
OUT = MANUSCRIPT_DIR / "submission"

FIGURES = {1: "figure_1_design", 2: "figure_2_primary", 3: "figure_3_robustness"}
FILL_MARKER = re.compile(r"<<FILL\b[^>]*>>")


def find_pandoc() -> str | None:
    found = shutil.which("pandoc")
    if found:
        return found
    for candidate in (Path(os.environ.get("LOCALAPPDATA", "")) / "Pandoc" / "pandoc.exe",
                      Path("C:/Program Files/Pandoc/pandoc.exe")):
        if candidate.is_file():
            return str(candidate)
    return None


def have_svg_converter() -> bool:
    return (importlib.util.find_spec("svglib") is not None
            and importlib.util.find_spec("reportlab") is not None)


def embed_figures(text: str) -> str:
    """The manuscript with each figure embedded directly above its legend.

    Refuses rather than guessing: a missing legend, or figures that are already embedded, is an error.
    """
    head, sep, legends = text.partition("\n## Figure legends\n")
    if not sep:
        raise ValueError("the manuscript has no Figure legends section")
    if "![" in legends:
        raise ValueError("figures are already embedded")
    for n, stem in FIGURES.items():
        m = re.search(rf"^\*\*Figure {n}\.", legends, re.M)
        if not m:
            raise ValueError(f"no legend for Figure {n}")
        legends = legends[:m.start()] + f"![](figures/{stem}.svg){{width=100%}}\n\n" + legends[m.start():]
    return head + sep + legends


def cover_letter(submission: str) -> str:
    """Section 6 of the submission pack, unquoted, with inline-code marks removed."""
    if "## 6. Cover letter" not in submission:
        raise ValueError("the submission pack has no section 6 cover letter")
    section = submission.split("## 6. Cover letter", 1)[1].split("\n## 7.", 1)[0]
    kept = []
    for line in section.splitlines():
        if line.startswith("> "):
            kept.append(line[2:])
        elif line.strip() == ">":
            kept.append("")
    body = "\n".join(kept).strip()
    return re.sub(r"`([^`]*)`", r"\1", body) + "\n"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _pandoc_docx(pandoc: str, markdown: str, out: Path) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "source.md"
        src.write_text(markdown, encoding="utf-8")
        subprocess.run([pandoc, str(src), "--from", "markdown", "--to", "docx", "--output", str(out),
                        f"--resource-path={MANUSCRIPT_DIR}"], check=True, capture_output=True, text=True)


def _svg_to_pdf(svg: Path, pdf: Path) -> None:
    from reportlab.graphics import renderPDF
    from svglib.svglib import svg2rlg
    drawing = svg2rlg(str(svg))
    if drawing is None:
        raise RuntimeError(f"svglib could not read {svg.name}")
    renderPDF.drawToFile(drawing, str(pdf))


def _word_pdf(docx: Path, pdf: Path) -> tuple[bool, str]:
    """Export through Word by COM automation. Returns (ok, detail); never raises."""
    def quote(p: Path) -> str:
        return str(p.resolve()).replace("'", "''")
    script = ("$ErrorActionPreference = 'Stop'; "
              "$word = New-Object -ComObject Word.Application; $word.Visible = $false; "
              f"try {{ $doc = $word.Documents.Open('{quote(docx)}', $false, $true); "
              f"$doc.ExportAsFixedFormat('{quote(pdf)}', 17); $doc.Close($false) }} "
              "finally { $word.Quit() }")
    try:
        r = subprocess.run(["powershell", "-NoProfile", "-NonInteractive", "-Command", script],
                           capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as e:
        return False, str(e)
    return (r.returncode == 0 and pdf.is_file()), (r.stderr or r.stdout).strip()[:400]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Render the Generation-1 submission files")
    ap.add_argument("--draft", action="store_true",
                    help="allow FILL placeholders and write to submission/draft/ instead")
    a = ap.parse_args(argv)

    sys.path.insert(0, str(ROOT / "experiments"))
    import run_gen1_manuscript as M

    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    submission = SUBMIT.read_text(encoding="utf-8")

    # ---- every precondition before a single file is written ----------------------------------------
    problems: list[str] = []
    pandoc = find_pandoc()
    if not pandoc:
        problems.append("pandoc is not installed:  winget install --id JohnMacFarlane.Pandoc -e")
    if not have_svg_converter():
        problems.append("svglib and reportlab are not installed:  pip install svglib reportlab")
    shape = {k: v for k, v in M.structure_problems(manuscript).items() if v}
    if shape:
        problems.append(f"the manuscript does not have the required structure: {shape}")
    missing_figures = [f"{stem}.svg" for stem in FIGURES.values() if not (FIGURE_DIR / f"{stem}.svg").is_file()]
    if missing_figures:
        problems.append(f"figures missing: {missing_figures}")
    try:
        letter = cover_letter(submission)
        biorxiv_text = embed_figures(manuscript)
    except ValueError as e:
        problems.append(str(e))
        letter, biorxiv_text = "", ""
    fills = FILL_MARKER.findall(manuscript) + FILL_MARKER.findall(letter)
    if fills and not a.draft:
        problems.append(f"{len(fills)} FILL marker(s) remain; fill them, or render a --draft")

    if problems:
        print("REFUSED -- nothing was written:")
        for p in problems:
            print(f"  {p}")
        return 2

    out_dir = OUT / "draft" if a.draft else OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}

    for stem in FIGURES.values():
        pdf = out_dir / f"{stem}.pdf"
        _svg_to_pdf(FIGURE_DIR / f"{stem}.svg", pdf)
        outputs[pdf.name] = pdf

    for name, text in (("MANUSCRIPT_BMC.docx", manuscript),
                       ("MANUSCRIPT_bioRxiv.docx", biorxiv_text),
                       ("COVER_LETTER.docx", letter)):
        _pandoc_docx(pandoc, text, out_dir / name)
        outputs[name] = out_dir / name

    pdf = out_dir / "MANUSCRIPT_bioRxiv.pdf"
    word_ok, word_detail = _word_pdf(out_dir / "MANUSCRIPT_bioRxiv.docx", pdf)
    if word_ok:
        outputs[pdf.name] = pdf

    version = subprocess.run([pandoc, "--version"], capture_output=True, text=True).stdout.splitlines()
    manifest = {
        "draft": a.draft,
        "fill_markers_remaining": len(fills),
        "manuscript_canonical_lf_sha256": hashlib.sha256(
            manuscript.replace("\r\n", "\n").encode("utf-8")).hexdigest(),
        "pandoc": version[0] if version else "unknown",
        "pdf_exported_by_word": word_ok,
        "word_detail": "" if word_ok else word_detail,
        "outputs": {name: sha256(path) for name, path in sorted(outputs.items())},
        "note": "Not byte-deterministic (.docx and .pdf carry timestamps), so not locked. These hashes "
                "identify the exact files that were submitted.",
    }
    (out_dir / "RENDER_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    for name in sorted(outputs):
        print(f"  wrote {out_dir.relative_to(ROOT).as_posix()}/{name}")
    if not word_ok:
        print()
        print("  Word could not export the PDF automatically. Open MANUSCRIPT_bioRxiv.docx in Word and")
        print("  use File > Save As > PDF, into the same folder, then re-run this script.")
        print(f"  detail: {word_detail}")
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

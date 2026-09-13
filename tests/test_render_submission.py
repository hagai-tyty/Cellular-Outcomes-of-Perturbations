"""Contracts for the submission renderer.

None of these needs pandoc, svglib or Word: they test the parts that decide WHAT is rendered -- that
the bioRxiv copy embeds all three figures and changes nothing else, that the cover letter comes out
clean, and that the renderer refuses before writing anything when a precondition fails. The rendering
itself is checked by inspecting the produced pages at release time, not here.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import render_gen1_submission as R  # noqa: E402

MANUSCRIPT = ROOT / "results" / "manuscript" / "MANUSCRIPT.md"
SUBMIT = ROOT / "results" / "manuscript" / "SUBMISSION.md"


def test_the_biorxiv_copy_embeds_each_figure_above_its_own_legend():
    text = MANUSCRIPT.read_text(encoding="utf-8")
    embedded = R.embed_figures(text)
    for n, stem in R.FIGURES.items():
        image = f"![](figures/{stem}.svg){{width=100%}}"
        assert embedded.count(image) == 1, stem
        assert embedded.index(image) < embedded.index(f"**Figure {n}.")
    positions = [embedded.index(f"figures/{stem}.svg") for stem in R.FIGURES.values()]
    assert positions == sorted(positions), "figures must appear in order"


def test_embedding_changes_nothing_but_the_three_image_lines():
    text = MANUSCRIPT.read_text(encoding="utf-8")
    embedded = R.embed_figures(text)
    stripped = embedded
    for stem in R.FIGURES.values():
        stripped = stripped.replace(f"![](figures/{stem}.svg){{width=100%}}\n\n", "", 1)
    assert stripped == text


def test_embedding_refuses_rather_than_guessing():
    text = MANUSCRIPT.read_text(encoding="utf-8")
    with pytest.raises(ValueError, match="already embedded"):
        R.embed_figures(R.embed_figures(text))
    with pytest.raises(ValueError, match="no legend for Figure 2"):
        R.embed_figures(text.replace("**Figure 2.", "**Fig 2.", 1))
    with pytest.raises(ValueError, match="no Figure legends section"):
        R.embed_figures(text.replace("\n## Figure legends\n", "\n## Legends\n", 1))


def test_the_cover_letter_comes_out_unquoted_and_without_code_marks():
    letter = R.cover_letter(SUBMIT.read_text(encoding="utf-8"))
    assert letter.startswith("Dear Editor,")
    assert not any(line.startswith(">") for line in letter.splitlines())
    assert "`" not in letter
    assert "**Preprint.**" in letter and "**Software licence.**" in letter


def test_the_renderer_refuses_and_writes_nothing_without_its_tools(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(R, "OUT", tmp_path / "out")
    monkeypatch.setattr(R, "find_pandoc", lambda: None)
    monkeypatch.setattr(R, "have_svg_converter", lambda: False)
    assert R.main(["--draft"]) == 2
    assert not (tmp_path / "out").exists()
    out = capsys.readouterr().out
    assert "pandoc is not installed" in out and "svglib" in out


def test_the_renderer_refuses_a_fill_marker_unless_drafting(tmp_path, monkeypatch, capsys):
    manuscript = tmp_path / "MANUSCRIPT.md"
    manuscript.write_text(MANUSCRIPT.read_text(encoding="utf-8") + "\n<<FILL: something>>\n",
                          encoding="utf-8")
    monkeypatch.setattr(R, "MANUSCRIPT", manuscript)
    monkeypatch.setattr(R, "OUT", tmp_path / "out")
    monkeypatch.setattr(R, "find_pandoc", lambda: "pandoc")
    monkeypatch.setattr(R, "have_svg_converter", lambda: True)
    assert R.main([]) == 2
    assert not (tmp_path / "out").exists()
    assert "FILL marker(s) remain" in capsys.readouterr().out


def test_the_renderer_never_writes_the_manuscript():
    src = (ROOT / "experiments" / "render_gen1_submission.py").read_text(encoding="utf-8")
    assert "MANUSCRIPT.write_text" not in src
    assert "SUBMIT.write_text" not in src

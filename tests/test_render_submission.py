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


# ---- code blocks must not wrap on the page -----------------------------------------------------------
STYLES = ('<w:styles><w:style w:type="paragraph" w:styleId="SourceCode"><w:name w:val="Source Code" />'
          '</w:style><w:style w:customStyle="1" w:styleId="VerbatimChar" w:type="character">'
          '<w:rPr><w:rFonts w:ascii="Consolas" w:hAnsi="Consolas" /><w:sz w:val="22" /></w:rPr></w:style>'
          '<w:style w:styleId="Normal"><w:rPr><w:sz w:val="24" /></w:rPr></w:style></w:styles>')


def test_the_code_capacity_reproduces_what_word_did_at_11_pt():
    # the first render, at pandoc's 11 pt, broke both digest lines after their 71st character
    assert R.code_line_capacity(22) == 71
    assert R.code_line_capacity(R.CODE_HALF_POINTS) == 98


def test_every_code_line_in_the_manuscript_fits_at_the_rendered_size():
    text = MANUSCRIPT.read_text(encoding="utf-8")
    assert 0 < R.longest_code_line(text) <= R.code_line_capacity(R.CODE_HALF_POINTS)


def test_the_code_size_patch_changes_only_the_code_style():
    patched = R.set_code_size(STYLES, 16)
    assert patched.count('<w:sz w:val="16" />') == 1 and '<w:sz w:val="22" />' not in patched
    assert patched.replace('<w:sz w:val="16" />', '<w:sz w:val="22" />') == STYLES


def test_the_code_size_patch_refuses_rather_than_guessing():
    with pytest.raises(ValueError, match="one VerbatimChar style"):
        R.set_code_size(STYLES.replace("VerbatimChar", "Verbatim"), 16)
    with pytest.raises(ValueError, match="one size"):
        R.set_code_size(STYLES.replace('<w:sz w:val="22" />', ""), 16)


def test_the_renderer_refuses_a_code_line_that_would_wrap(tmp_path, monkeypatch, capsys):
    manuscript = tmp_path / "MANUSCRIPT.md"
    too_long = "x" * (R.code_line_capacity(R.CODE_HALF_POINTS) + 1)
    manuscript.write_text(MANUSCRIPT.read_text(encoding="utf-8") + f"\n```text\n{too_long}\n```\n",
                          encoding="utf-8")
    monkeypatch.setattr(R, "MANUSCRIPT", manuscript)
    monkeypatch.setattr(R, "OUT", tmp_path / "out")
    monkeypatch.setattr(R, "find_pandoc", lambda: "pandoc")
    monkeypatch.setattr(R, "have_svg_converter", lambda: True)
    assert R.main(["--draft"]) == 2
    assert not (tmp_path / "out").exists()
    assert "would wrap" in capsys.readouterr().out


def test_an_unmeasured_word_step_is_not_read_as_no_wrapping():
    assert R.wrapped_count("WRAPPED=0\n") == 0
    assert R.wrapped_count("noise\nWRAPPED=3\r\n") == 3
    assert R.wrapped_count("") is None
    assert R.wrapped_count("the export failed") is None

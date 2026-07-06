import pytest

from cellfate.common.panel import GenePanel


def test_hash_stable_and_order_sensitive():
    p1 = GenePanel(["A", "B", "C"])
    p2 = GenePanel(["A", "B", "C"])
    p3 = GenePanel(["C", "B", "A"])
    assert p1.hash() == p2.hash()
    assert p1.hash() != p3.hash()      # feature ORDER matters
    assert len(p1) == 3


def test_rejects_empty_and_duplicates():
    with pytest.raises(ValueError):
        GenePanel([])
    with pytest.raises(ValueError):
        GenePanel(["A", "A"])


def test_save_load_roundtrip(tmp_path):
    p = GenePanel([f"GENE{i}" for i in range(20)])
    f = tmp_path / "panel.txt"
    p.save(f)
    loaded = GenePanel.load(f)
    assert loaded == p and loaded.hash() == p.hash()


def test_load_ignores_comments(tmp_path):
    f = tmp_path / "p.txt"
    f.write_text("# header\nA\n\nB\n# mid\nC\n")
    p = GenePanel.load(f)
    assert p.genes == ("A", "B", "C")

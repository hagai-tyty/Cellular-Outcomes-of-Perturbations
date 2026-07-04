from cellfate.common import constants as C


def test_class_order_is_fixed():
    assert C.CLASSES == ("safe", "loss", "death")
    assert C.N_CLASSES == 3
    assert C.SAFE_IDX == 0 and C.LOSS_IDX == 1 and C.DEATH_IDX == 2


def test_idx_maps_consistent():
    for i, c in enumerate(C.CLASSES):
        assert C.CLASS_TO_IDX[c] == i
        assert C.IDX_TO_CLASS[i] == c


def test_cancer_sources_mask_policy():
    assert "tahoe" in C.CANCER_SOURCES


def test_cell_cycle_lists_nonempty_and_unique():
    assert len(C.S_GENES) > 10 and len(set(C.S_GENES)) == len(C.S_GENES)
    assert len(C.G2M_GENES) > 10 and len(set(C.G2M_GENES)) == len(C.G2M_GENES)


import pytest
from pydantic import ValidationError

from cellfate.common.schemas import (
    BundleMeta,
    ConformalParams,
    ManifestRow,
    Modality,
    ResParams,
    ScalerParams,
)
from conftest import G_TEST, make_sample


def test_valid_sample_roundtrips():
    s = make_sample()
    assert s.u_modality is Modality.CHEM
    s2 = type(s).model_validate(s.model_dump())
    assert s2 == s


def test_y_cls_must_sum_to_one():
    with pytest.raises(ValidationError):
        make_sample(y_cls=[0.5, 0.1, 0.1])


def test_y_cls_wrong_length():
    with pytest.raises(ValidationError):
        make_sample(y_cls=[0.5, 0.5])


def test_chem_requires_fingerprint():
    with pytest.raises(ValidationError):
        make_sample(u_chem_fp=None)


def test_fingerprint_must_be_bits():
    with pytest.raises(ValidationError):
        make_sample(u_chem_fp=[2] * 2048)


def test_age_mask_requires_finite_age():
    with pytest.raises(ValidationError):
        make_sample(age_mask=True, y_age=None)


def test_masked_age_allows_none():
    s = make_sample(age_mask=False, y_age=None)
    assert s.y_age is None


def test_masked_age_allows_nan():
    # age_mask=False must tolerate NaN (how masked ages round-trip through Parquet)
    s = make_sample(age_mask=False, y_age=float("nan"))
    assert s.age_mask is False


def test_masked_age_rejects_finite_value():
    # reverse invariant: a masked-out row must not carry a spurious finite label
    with pytest.raises(ValidationError):
        make_sample(age_mask=False, y_age=5.0)


def test_genetic_sample_rejects_stray_fingerprint():
    # the chemical fingerprint must be absent for non-chemical modalities
    with pytest.raises(ValidationError):
        make_sample(u_modality=Modality.GENETIC, u_gene_emb=[0.1, 0.2],
                    u_chem_fp=[0] * 2048, scaffold_id=None)


def test_chem_requires_scaffold():
    with pytest.raises(ValidationError):
        make_sample(scaffold_id=None)


def test_genetic_sample_ok_without_fingerprint():
    s = make_sample(u_modality=Modality.GENETIC, u_chem_fp=None,
                    u_gene_emb=[0.1, 0.2], scaffold_id=None)
    assert s.u_gene_emb == [0.1, 0.2]


def test_manifest_row_from_sample():
    s = make_sample()
    r = ManifestRow.from_sample(s, shard_id="tahoe_AAA", row_idx=3)
    assert r.cell_id == s.cell_id and r.shard_id == "tahoe_AAA" and r.row_idx == 3


def test_scaler_params_validation():
    sp = ScalerParams(x_mean=[0.0] * G_TEST, x_std=[1.0] * G_TEST,
                      dt_mean=[0.0, 0.0], dt_std=[1.0, 1.0],
                      proliferation_coef=[0.0, 0.0], gene_panel_hash="abc")
    assert len(sp.x_mean) == G_TEST
    with pytest.raises(ValidationError):
        ScalerParams(x_mean=[0.0] * G_TEST, x_std=[1.0] * (G_TEST - 1),
                     dt_mean=[0.0, 0.0], dt_std=[1.0, 1.0],
                     proliferation_coef=[0.0, 0.0], gene_panel_hash="abc")


def test_res_params_bounds():
    ResParams()  # defaults valid
    with pytest.raises(ValidationError):
        ResParams(k=0.5)            # k must be >= 1
    with pytest.raises(ValidationError):
        ResParams(tau_safe=1.5)     # must be in (0,1)


def test_conformal_requires_quantile_per_level():
    ConformalParams(levels=[0.9], q={"0.9": 1.5})
    with pytest.raises(ValidationError):
        ConformalParams(levels=[0.9], q={})


def test_bundle_meta_class_order_enforced():
    BundleMeta(n_members=3, gene_panel_hash="h")
    with pytest.raises(ValidationError):
        BundleMeta(n_members=3, gene_panel_hash="h", classes=["loss", "safe", "death"])

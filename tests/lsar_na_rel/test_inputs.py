"""Failing-first contract tests for canonical S1 input loading (#27)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lsar_na_rel.config import load_talk_v0_config, profile_feature_columns
from lsar_na_rel.inputs import (
    EXPECTED_PRIMARY_HEADER,
    EXPECTED_SOURCE_IDS,
    FORBIDDEN_FEATURE_COLUMNS,
    RISK_ALLOWED_VARIANTS,
    load_canonical_cohort,
    parse_variant,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_ROOT = REPO_ROOT / "workflows" / "lsar-na-rel"
PUBLIC_DIR = REPO_ROOT / "data" / "public" / "lsar-na-rel"
PRIMARY_CSV = PUBLIC_DIR / "na-rel-primary.csv"
LABELS_CSV = PUBLIC_DIR / "labels.csv"
CONFIG_PATH = WORKFLOW_ROOT / "config" / "lsar-na-rel.json"


@pytest.fixture(scope="module")
def talk_v0_config():
    return load_talk_v0_config(CONFIG_PATH)


@pytest.fixture(scope="module")
def registered_feature_columns(talk_v0_config):
    return profile_feature_columns(talk_v0_config, "talk-v0-no-foldx")


@pytest.fixture(scope="module")
def cohort(registered_feature_columns):
    return load_canonical_cohort(
        primary_path=PRIMARY_CSV,
        labels_path=LABELS_CSV,
        feature_column_names=registered_feature_columns,
    )


def test_public_inputs_exist_without_duplication():
    assert PRIMARY_CSV.is_file()
    assert LABELS_CSV.is_file()
    assert not (WORKFLOW_ROOT / "data" / "public").exists()


def test_exact_primary_header():
    header = PRIMARY_CSV.read_text(encoding="utf-8").splitlines()[0]
    assert header == EXPECTED_PRIMARY_HEADER
    assert header == "variant,LSAR_mut,LSAR_WT,Na_rel,ln(Na_rel)"


def test_exactly_sixteen_unique_variants(cohort):
    variants = cohort.primary["variant"].tolist()
    assert len(variants) == 16
    assert len(set(variants)) == 16


def test_variant_parser_accepts_canonical_labels(cohort):
    for label in cohort.primary["variant"]:
        parsed = parse_variant(label)
        assert parsed.label == label
        assert parsed.wt.isalpha() and parsed.wt.isupper() and len(parsed.wt) == 1
        assert parsed.mut.isalpha() and parsed.mut.isupper() and len(parsed.mut) == 1
        assert parsed.position >= 1


def test_positive_numeric_na_rel(cohort):
    na_rel = cohort.primary["Na_rel"].astype(float)
    assert np.isfinite(na_rel).all()
    assert (na_rel > 0).all()


def test_recomputed_log_matches_public_ln_within_1e_10(cohort):
    na_rel = cohort.primary["Na_rel"].astype(float)
    public_ln = cohort.primary["ln(Na_rel)"].astype(float)
    recomputed = np.log(na_rel.to_numpy())
    assert np.max(np.abs(recomputed - public_ln.to_numpy())) < 1e-10
    assert np.max(np.abs(cohort.target.to_numpy() - public_ln.to_numpy())) < 1e-10


def test_one_to_one_join_to_labels(cohort):
    assert len(cohort.primary) == 16
    assert len(cohort.groups) == 16
    assert cohort.groups.index.equals(cohort.primary.index)
    assert cohort.groups.notna().all()
    assert cohort.eligibility.index.equals(cohort.primary.index)


def test_exactly_six_source_id_groups(cohort):
    groups = set(cohort.groups.tolist())
    assert groups == EXPECTED_SOURCE_IDS
    assert len(groups) == 6


def test_normalization_tier_counts(cohort):
    tiers = cohort.eligibility["normalization_tier"]
    assert (tiers == "Primary").sum() == 13
    assert (tiers == "Risky allowed").sum() == 3


def test_risky_set_exactly_three_named_variants(cohort):
    risky = set(
        cohort.primary.loc[
            cohort.eligibility["normalization_tier"] == "Risky allowed", "variant"
        ]
    )
    assert risky == RISK_ALLOWED_VARIANTS


def test_g256e_absent(cohort):
    assert "G256E" not in set(cohort.primary["variant"])


def test_leakage_fields_cannot_enter_estimator_feature_columns(cohort, talk_v0_config):
    for profile_name in ("talk-v0-foldx", "talk-v0-no-foldx"):
        columns = profile_feature_columns(talk_v0_config, profile_name)
        overlap = set(columns) & FORBIDDEN_FEATURE_COLUMNS
        assert not overlap, f"{profile_name} leaked {overlap}"

    assert set(cohort.feature_column_names).isdisjoint(FORBIDDEN_FEATURE_COLUMNS)
    assert "source_id" not in cohort.feature_column_names
    assert cohort.groups.name == "source_id"
    # Structural separation: target/groups/eligibility are not feature values.
    assert cohort.feature_values is None or cohort.feature_values.empty
    assert list(cohort.target.index) == list(cohort.primary.index)
    assert set(cohort.eligibility.columns).isdisjoint(set(cohort.feature_column_names))


def test_approved_profile_schemas_are_encoded(talk_v0_config):
    assert profile_feature_columns(talk_v0_config, "talk-v0-foldx") == (
        "foldx_ddg_bh_fh",
        "wt_site_rsa_8act_mean_ab",
        "ihm_interface_flag_8act",
    )
    assert profile_feature_columns(talk_v0_config, "talk-v0-no-foldx") == (
        "delta_charge",
        "wt_site_rsa_8act_mean_ab",
        "ihm_interface_flag_8act",
    )


def test_forbidden_metadata_rejected_as_feature_columns():
    with pytest.raises(ValueError, match="Forbidden columns"):
        load_canonical_cohort(
            primary_path=PRIMARY_CSV,
            labels_path=LABELS_CSV,
            feature_column_names=("source_id", "delta_charge"),
        )


def test_loader_keeps_target_groups_eligibility_and_features_separate(cohort):
    assert isinstance(cohort.target, pd.Series)
    assert isinstance(cohort.groups, pd.Series)
    assert isinstance(cohort.eligibility, pd.DataFrame)
    assert isinstance(cohort.feature_column_names, tuple)
    assert cohort.target.name == "y"
    assert "y" not in cohort.feature_column_names
    assert "source_id" not in cohort.primary.columns
    assert "normalization_tier" not in cohort.primary.columns

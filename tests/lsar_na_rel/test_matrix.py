"""Matrix assembly from committed features.csv (no Na_rel-based score asserts)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lsar_na_rel.config import load_talk_v0_config, ladder_feature_columns
from lsar_na_rel.inputs import FORBIDDEN_FEATURE_COLUMNS, load_canonical_cohort
from lsar_na_rel.matrix import (
    assemble_model_matrix,
    build_talk_v0_matrix,
    load_committed_features,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_ROOT = REPO_ROOT / "workflows" / "lsar-na-rel"
PUBLIC_DIR = REPO_ROOT / "data" / "public" / "lsar-na-rel"
PRIMARY_CSV = PUBLIC_DIR / "na-rel-primary.csv"
LABELS_CSV = PUBLIC_DIR / "labels.csv"
FEATURES_CSV = PUBLIC_DIR / "features.csv"
CONFIG_PATH = WORKFLOW_ROOT / "config" / "lsar-na-rel.json"

IHM_POSITIVES = frozenset({"D382Y", "R403Q"})


@pytest.fixture(scope="module")
def config():
    return load_talk_v0_config(CONFIG_PATH)


def test_load_committed_features_does_not_require_labels(config):
    columns = ladder_feature_columns(config, "ladder-3-nofoldx")
    frame = load_committed_features(FEATURES_CSV, columns)
    assert list(frame.columns) == ["variant", *columns]
    assert "Na_rel" not in frame.columns
    assert "y" not in frame.columns
    assert "source_id" not in frame.columns


def test_ihm_flag_positives_are_exactly_d382y_and_r403q(config):
    columns = ladder_feature_columns(config, "ladder-2-rsa-ihm")
    frame = load_committed_features(FEATURES_CSV, columns)
    positives = set(
        frame.loc[frame["ihm_interface_flag_8act"].astype(int) == 1, "variant"]
    )
    assert positives == IHM_POSITIVES


def test_assemble_restores_primary_order_from_shuffled_features(config):
    columns = ladder_feature_columns(config, "ladder-1-rsa")
    cohort = load_canonical_cohort(
        primary_path=PRIMARY_CSV,
        labels_path=LABELS_CSV,
        feature_column_names=columns,
    )
    features = load_committed_features(FEATURES_CSV, columns)
    shuffled = features.sample(frac=1.0, random_state=0).reset_index(drop=True)
    matrix = assemble_model_matrix(
        run_id="ladder-1-rsa",
        feature_columns=columns,
        features=shuffled,
        cohort=cohort,
    )
    expected = pd.read_csv(PRIMARY_CSV)["variant"].astype(str).tolist()
    assert matrix.variants.tolist() == expected
    assert list(matrix.X.columns) == list(columns)
    assert len(matrix.X) == 16


def test_x_excludes_target_groups_and_forbidden_columns(config):
    matrix = build_talk_v0_matrix(
        run_id="ladder-3-nofoldx",
        config=config,
        primary_path=PRIMARY_CSV,
        labels_path=LABELS_CSV,
        features_path=FEATURES_CSV,
    )
    assert set(matrix.X.columns).isdisjoint(FORBIDDEN_FEATURE_COLUMNS)
    assert "source_id" not in matrix.X.columns
    assert "y" not in matrix.X.columns
    assert matrix.groups.name == "source_id"
    assert matrix.y.name == "y"
    assert np.isfinite(matrix.X.to_numpy(dtype=float)).all()
    assert len(matrix.groups.unique()) == 6


def test_missing_feature_column_is_a_hard_error(config, tmp_path):
    columns = ladder_feature_columns(config, "ladder-3-nofoldx")
    broken = pd.read_csv(FEATURES_CSV).drop(columns=["delta_charge"])
    path = tmp_path / "features.csv"
    broken.to_csv(path, index=False)
    with pytest.raises(ValueError, match="delta_charge"):
        load_committed_features(path, columns)


def test_non_finite_feature_is_a_hard_error(config, tmp_path):
    columns = ladder_feature_columns(config, "ladder-1-rsa")
    broken = pd.read_csv(FEATURES_CSV)
    broken.loc[0, "wt_site_rsa_8act_mean_ab"] = np.nan
    path = tmp_path / "features.csv"
    broken.to_csv(path, index=False)
    features = load_committed_features(path, columns)
    cohort = load_canonical_cohort(
        primary_path=PRIMARY_CSV,
        labels_path=LABELS_CSV,
        feature_column_names=columns,
    )
    with pytest.raises(ValueError, match="non-finite"):
        assemble_model_matrix(
            run_id="ladder-1-rsa",
            feature_columns=columns,
            features=features,
            cohort=cohort,
        )

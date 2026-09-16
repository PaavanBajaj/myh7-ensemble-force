"""Official Talk-v0 labeled feature table (variant, features, y)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lsar_na_rel.config import (
    APPROVED_LADDER_RUN_IDS,
    load_talk_v0_config,
    run_feature_columns,
)
from lsar_na_rel.model_table import (
    IHM_MIN_DIST_COLUMN,
    MODEL_TABLE_COLUMNS,
    build_talk_v0_model_table,
    write_talk_v0_model_table,
)
from lsar_na_rel.structure import CONTACT_CUTOFF_ANGSTROM

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_ROOT = REPO_ROOT / "workflows" / "lsar-na-rel"
PUBLIC_DIR = REPO_ROOT / "data" / "public" / "lsar-na-rel"
PRIMARY_CSV = PUBLIC_DIR / "na-rel-primary.csv"
LABELS_CSV = PUBLIC_DIR / "labels.csv"
FEATURES_CSV = PUBLIC_DIR / "features.csv"
DISTANCE_CSV = PUBLIC_DIR / "ihm_sparse_flag_variant_summary.csv"
CONFIG_PATH = WORKFLOW_ROOT / "config" / "lsar-na-rel.json"


@pytest.fixture(scope="module")
def table() -> pd.DataFrame:
    return build_talk_v0_model_table(
        primary_path=PRIMARY_CSV,
        labels_path=LABELS_CSV,
        features_path=FEATURES_CSV,
        distance_path=DISTANCE_CSV,
    )


def test_model_table_has_sixteen_rows_in_primary_order(table: pd.DataFrame):
    expected = pd.read_csv(PRIMARY_CSV)["variant"].astype(str).tolist()
    assert table["variant"].tolist() == expected
    assert list(table.columns) == list(MODEL_TABLE_COLUMNS)


def test_label_y_matches_log_na_rel_and_public_ln(table: pd.DataFrame):
    primary = pd.read_csv(PRIMARY_CSV)
    recomputed = np.log(primary["Na_rel"].to_numpy(dtype=float))
    public_ln = primary["ln(Na_rel)"].to_numpy(dtype=float)
    assert np.max(np.abs(table["y"].to_numpy(dtype=float) - recomputed)) < 1e-10
    assert np.max(np.abs(table["y"].to_numpy(dtype=float) - public_ln)) < 1e-10
    assert np.max(np.abs(table["ln_na_rel_public"].to_numpy(dtype=float) - public_ln)) < 1e-10


def test_ihm_min_dist_is_the_angstrom_quantity_thresholded_at_4p5(table: pd.DataFrame):
    dist = table[IHM_MIN_DIST_COLUMN].to_numpy(dtype=float)
    flag = table["ihm_interface_flag_8act"].to_numpy(dtype=int)
    expected_flag = (dist <= CONTACT_CUTOFF_ANGSTROM).astype(int)
    assert CONTACT_CUTOFF_ANGSTROM == 4.5
    assert np.array_equal(flag, expected_flag)
    by_variant = table.set_index("variant")
    assert int(by_variant.loc["D382Y", "ihm_interface_flag_8act"]) == 1
    assert int(by_variant.loc["R403Q", "ihm_interface_flag_8act"]) == 1
    assert float(by_variant.loc["D382Y", IHM_MIN_DIST_COLUMN]) == pytest.approx(
        3.132480, abs=1e-6
    )
    assert float(by_variant.loc["R403Q", IHM_MIN_DIST_COLUMN]) == pytest.approx(
        2.678151, abs=1e-6
    )
    assert int(by_variant.loc["H251N", "ihm_interface_flag_8act"]) == 0
    assert float(by_variant.loc["H251N", IHM_MIN_DIST_COLUMN]) > 4.5


def test_write_model_table_roundtrip(table: pd.DataFrame, tmp_path: Path):
    path = tmp_path / "talk-v0-model-table.csv"
    write_talk_v0_model_table(path, table)
    loaded = pd.read_csv(path)
    assert loaded["variant"].tolist() == table["variant"].tolist()
    assert IHM_MIN_DIST_COLUMN in loaded.columns
    assert "y" in loaded.columns


def test_exploratory_probe_is_not_a_runnable_model():
    config = load_talk_v0_config(CONFIG_PATH)
    assert APPROVED_LADDER_RUN_IDS == (
        "ladder-1-rsa",
        "ladder-2-rsa-ihm",
        "ladder-3-nofoldx",
    )
    with pytest.raises(KeyError, match="probe-ihm-min-dist"):
        run_feature_columns(config, "probe-ihm-min-dist")

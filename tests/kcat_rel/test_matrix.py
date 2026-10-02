from __future__ import annotations

import ast
import csv
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kcat_rel.features import build_feature_table
from kcat_rel.matrix import build_model_table
from kcat_rel.registry import ACTIVE_FEATURE_COLUMNS, FROZEN_PRIMARY_VARIANTS, load_registry


WORKFLOW_ROOT = Path(__file__).resolve().parents[2] / "results/kcat-rel"
DATA_ROOT = Path(__file__).resolve().parents[2] / "data/public/kcat-rel"
CONFIG = WORKFLOW_ROOT / "config" / "kcat-rel-v0.json"


@pytest.fixture
def reference_sequence() -> str:
    source = ast.parse((Path(__file__).with_name("test_chemistry.py")).read_text(encoding="utf-8"))
    return next(
        ast.literal_eval(node.value)
        for node in source.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "REFERENCE_SEQUENCE" for target in node.targets)
    )


@pytest.fixture
def valid_features(reference_sequence: str) -> pd.DataFrame:
    """A real prior, label-blind feature artifact used for the post-audit join."""
    return build_feature_table(
        FROZEN_PRIMARY_VARIANTS,
        reference_sequence,
        Path(__file__).with_name("fixtures") / "8act-mini.cif",
        load_registry(CONFIG),
    ).table


@pytest.fixture
def canonical_path() -> Path:
    return DATA_ROOT / "canonical-labels.csv"


@pytest.fixture
def measurements_path() -> Path:
    return DATA_ROOT / "measurement-evidence.csv"


def test_model_table_joins_only_after_feature_audit(
    valid_features: pd.DataFrame, canonical_path: Path, measurements_path: Path
) -> None:
    """Changing row order or joining labels without validated evidence must fail."""
    table = build_model_table(valid_features, canonical_path, measurements_path)

    assert tuple(table.variant) == FROZEN_PRIMARY_VARIANTS
    assert np.allclose(table.y, pd.read_csv(canonical_path).ln_kcat_rel)
    assert tuple(table.columns) == ("variant", "y", "y_error", "source_id", *ACTIVE_FEATURE_COLUMNS)
    assert table.loc[:, ACTIVE_FEATURE_COLUMNS].equals(valid_features.loc[:, ACTIVE_FEATURE_COLUMNS])


def test_model_table_rejects_duplicate_feature_rows(
    valid_features: pd.DataFrame, canonical_path: Path, measurements_path: Path
) -> None:
    """Permitting duplicate feature variants would make the label join many-to-one."""
    invalid = pd.concat([valid_features, valid_features.iloc[[0]]], ignore_index=True)

    with pytest.raises(ValueError, match="one-to-one|frozen cohort"):
        build_model_table(invalid, canonical_path, measurements_path)


def test_model_table_rejects_tampered_measurement_evidence_before_label_use(
    valid_features: pd.DataFrame, canonical_path: Path, measurements_path: Path, tmp_path: Path
) -> None:
    """Skipping the evidence audit would wrongly construct a table from broken evidence."""
    tampered = tmp_path / "measurement-evidence.csv"
    with measurements_path.open(newline="", encoding="utf-8") as source:
        rows = list(csv.DictReader(source))
        fieldnames = list(rows[0])
    rows[0]["wt_context_match"] = "false"
    with tampered.open("w", newline="", encoding="utf-8") as destination:
        writer = csv.DictWriter(destination, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with pytest.raises(Exception, match="matched WT"):
        build_model_table(valid_features, canonical_path, tampered)

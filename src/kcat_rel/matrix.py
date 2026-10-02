"""Post-audit assembly of ``kcat_rel`` labels, metadata, and fixed features."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from kcat_rel.evidence import audit_evidence
from kcat_rel.registry import ACTIVE_FEATURE_COLUMNS, FROZEN_PRIMARY_VARIANTS


_FEATURE_COLUMNS = ("variant", *ACTIVE_FEATURE_COLUMNS)
_MODEL_COLUMNS = ("variant", "y", "y_error", "source_id", *ACTIVE_FEATURE_COLUMNS)


def _audit_features(features: pd.DataFrame) -> pd.DataFrame:
    """Accept only the exact finite, ordered feature artifact before evidence opens."""
    if not isinstance(features, pd.DataFrame) or tuple(features.columns) != _FEATURE_COLUMNS:
        raise ValueError("feature table schema does not match the active profile")
    if len(features) != len(FROZEN_PRIMARY_VARIANTS) or tuple(features["variant"]) != FROZEN_PRIMARY_VARIANTS:
        raise ValueError("feature table does not match the frozen cohort and order")
    if features["variant"].duplicated().any():
        raise ValueError("feature table must have a one-to-one frozen cohort")
    values = features.loc[:, ACTIVE_FEATURE_COLUMNS].to_numpy(dtype=float)
    if values.shape != (len(FROZEN_PRIMARY_VARIANTS), len(ACTIVE_FEATURE_COLUMNS)) or not np.isfinite(values).all():
        raise ValueError("feature table contains non-finite active-profile values")
    return features.copy()


def _read_audited_labels(canonical_path: str | Path) -> pd.DataFrame:
    """Load only the fields permitted after the independent evidence audit."""
    try:
        canonical = pd.read_csv(canonical_path)
    except (OSError, pd.errors.ParserError) as exc:
        raise ValueError(f"could not read canonical labels: {canonical_path}") from exc
    required = ("variant", "ln_kcat_rel", "ln_kcat_rel_error", "source_id")
    if any(column not in canonical.columns for column in required):
        raise ValueError("canonical labels lack required model-table fields")
    labels = canonical.loc[:, required].rename(
        columns={"ln_kcat_rel": "y", "ln_kcat_rel_error": "y_error"}
    )
    if len(labels) != len(FROZEN_PRIMARY_VARIANTS) or tuple(labels["variant"]) != FROZEN_PRIMARY_VARIANTS:
        raise ValueError("canonical labels do not match the frozen cohort and order")
    if labels["variant"].duplicated().any():
        raise ValueError("canonical labels must have a one-to-one frozen cohort")
    values = labels.loc[:, ["y", "y_error"]].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (labels["y_error"] < 0).any():
        raise ValueError("canonical labels contain invalid outcome uncertainty")
    if labels["source_id"].isna().any() or (labels["source_id"].astype(str).str.strip() == "").any():
        raise ValueError("canonical labels contain invalid source grouping")
    return labels


def build_model_table(
    features: pd.DataFrame, canonical_path: str | Path, measurements_path: str | Path
) -> pd.DataFrame:
    """Audit evidence, then join its labels one-to-one to the fixed feature table.

    The returned ``y_error`` and ``source_id`` are retained metadata, not
    predictors; consumers must select only ``ACTIVE_FEATURE_COLUMNS`` for X.
    """
    audited_features = _audit_features(features)
    # This must precede any direct read of outcome-adjacent canonical fields.
    audit_evidence(canonical_path, measurements_path)
    labels = _read_audited_labels(canonical_path)
    table = labels.merge(audited_features, on="variant", how="inner", sort=False, validate="one_to_one")
    if len(table) != len(FROZEN_PRIMARY_VARIANTS) or tuple(table["variant"]) != FROZEN_PRIMARY_VARIANTS:
        raise ValueError("model table one-to-one join does not match the frozen cohort")
    return table.loc[:, _MODEL_COLUMNS]

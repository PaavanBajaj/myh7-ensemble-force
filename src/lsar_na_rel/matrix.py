"""Deterministic Talk-v0 feature matrix from the committed features snapshot."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from lsar_na_rel.config import TalkV0Config, ladder_feature_columns
from lsar_na_rel.inputs import FORBIDDEN_FEATURE_COLUMNS, CohortBundle, load_canonical_cohort


@dataclass(frozen=True)
class ModelMatrix:
    run_id: str
    feature_columns: tuple[str, ...]
    X: pd.DataFrame
    y: pd.Series
    groups: pd.Series
    variants: pd.Series
    features_path: Path


def load_committed_features(
    features_path: str | Path, feature_columns: Sequence[str]
) -> pd.DataFrame:
    path = Path(features_path)
    columns = tuple(str(name) for name in feature_columns)
    if not columns:
        raise ValueError("feature_columns must be non-empty")
    leaked = set(columns) & FORBIDDEN_FEATURE_COLUMNS
    if leaked:
        raise ValueError(f"Forbidden columns selected as features: {sorted(leaked)}")
    frame = pd.read_csv(path)
    if "variant" not in frame.columns:
        raise ValueError(f"{path} missing variant column")
    missing = [name for name in columns if name not in frame.columns]
    if missing:
        raise ValueError(f"{path} missing feature columns: {missing}")
    return frame.loc[:, ["variant", *columns]].copy()


def assemble_model_matrix(
    *,
    run_id: str,
    feature_columns: Sequence[str],
    features: pd.DataFrame,
    cohort: CohortBundle,
    features_path: str | Path | None = None,
) -> ModelMatrix:
    columns = tuple(str(name) for name in feature_columns)
    if tuple(cohort.feature_column_names) != columns:
        raise ValueError("cohort.feature_column_names must match frozen run columns")
    variants = cohort.primary["variant"].astype(str).rename("variant")
    feat = features.copy()
    feat["variant"] = feat["variant"].astype(str)
    if feat["variant"].duplicated().any():
        raise ValueError("features table has duplicate variant keys")
    present = set(feat["variant"])
    missing = [v for v in variants.tolist() if v not in present]
    if missing:
        raise ValueError(f"Missing features for variants: {missing}")
    merged = pd.DataFrame({"variant": variants}).merge(
        feat, on="variant", how="left", validate="one_to_one"
    )
    X = merged.loc[:, list(columns)].apply(pd.to_numeric, errors="coerce")
    values = X.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("non-finite values in selected feature columns")
    leaked = set(X.columns) & FORBIDDEN_FEATURE_COLUMNS
    if leaked:
        raise ValueError(f"Forbidden columns present in X: {sorted(leaked)}")
    return ModelMatrix(
        run_id=str(run_id),
        feature_columns=columns,
        X=X.reset_index(drop=True),
        y=cohort.target.reset_index(drop=True),
        groups=cohort.groups.reset_index(drop=True),
        variants=variants.reset_index(drop=True),
        features_path=Path(features_path) if features_path is not None else Path("."),
    )


def build_talk_v0_matrix(
    *,
    run_id: str,
    config: TalkV0Config,
    primary_path: str | Path,
    labels_path: str | Path,
    features_path: str | Path,
) -> ModelMatrix:
    columns = ladder_feature_columns(config, run_id)
    features = load_committed_features(features_path, columns)
    cohort = load_canonical_cohort(
        primary_path=primary_path,
        labels_path=labels_path,
        feature_column_names=columns,
    )
    return assemble_model_matrix(
        run_id=run_id,
        feature_columns=columns,
        features=features,
        cohort=cohort,
        features_path=features_path,
    )

"""Read-only validation and descriptive checks for the frozen MYH7 diagnostic.

The run manifest pins artifact bytes before this module accepts their contents;
semantic checks then make malformed hash-bypassed fixtures diagnosable. This
module intentionally imports no fitting, validation, or sensitivity runtime.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from kcat_rel.features import _audit_feature_table, _audit_manifest_contract
from kcat_rel.registry import ACTIVE_FEATURE_COLUMNS, FROZEN_PRIMARY_VARIANTS, FROZEN_SOURCE_STUDIES, load_registry


_MODEL_TABLE_COLUMNS = ("variant", "y", "y_error", "source_id", *ACTIVE_FEATURE_COLUMNS)
_OOF_COLUMNS = ("model_kind", "variant", "source_id", "truth", "model_prediction", "baseline_prediction", "outer_fold", "selected_inner_hyperparameters", "scaler_mean", "scaler_scale")
_MODEL_KINDS = ("ridge", "gpr")
_OOF_ROW_MODEL_KINDS = ("gpr", "ridge")
_MANIFEST_ARTIFACTS = {"AUDIT.md", "inner-fold-audit.json", "metrics.json", "model-table.csv", "run.log", "sensitivity.json", "study-oof.csv", "tuning-records.csv", "variant-oof.csv"}
FROZEN_RUN_MANIFEST_SHA256 = "db10cbbbbb61cb7e12f1a2953074f342c582132fe44720ae9b58e553ef669901"


@dataclass(frozen=True)
class SummaryPaths:
    """Explicit paths to the immutable source artifacts for a frozen-run summary."""

    config: Path | str
    features: Path | str
    feature_manifest: Path | str
    model_table: Path | str
    study_oof: Path | str
    variant_oof: Path | str
    metrics: Path | str
    sensitivity: Path | str
    tuning_records: Path | str
    run_manifest: Path | str


@dataclass(frozen=True)
class SummaryInputs:
    """Validated frozen artifacts, retaining data only (never estimator objects)."""

    paths: SummaryPaths
    registry: Any
    features: pd.DataFrame
    model_table: pd.DataFrame
    study_oof: pd.DataFrame
    variant_oof: pd.DataFrame
    tuning_records: pd.DataFrame | None
    metrics: Mapping[str, Any]
    sensitivity: Mapping[str, Any]
    run_manifest: Mapping[str, Any]
    artifact_sha256: Mapping[str, str] = field(default_factory=dict)
    run_manifest_root_sha256: str = ""


@dataclass(frozen=True)
class SummaryChecks:
    """Descriptive values reconstructed from a byte-pinned frozen run."""

    feature_correlations: Mapping[str, float]
    study_oof: pd.DataFrame
    variant_oof: pd.DataFrame
    model_metrics: Mapping[str, Mapping[str, Mapping[str, float]]]
    baseline_metrics: Mapping[str, Mapping[str, Mapping[str, float]]]
    mae_improvement: Mapping[str, Mapping[str, float]]
    rmse_degradation: Mapping[str, Mapping[str, float]]
    failed_rules: Mapping[str, tuple[str, ...]]
    pass_gate: Mapping[str, bool]
    weighted_sensitivity: Mapping[str, Mapping[str, Mapping[str, Mapping[str, float]]]]
    perturbation_sensitivity: Mapping[str, Mapping[str, int | float]]
    label_uncertainties: Mapping[str, float]
    cohort_size: int
    source_study_count: int
    artifact_sha256: Mapping[str, str]
    run_manifest_root_sha256: str


def _sha256_file(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _summary_artifact_hashes(paths: SummaryPaths) -> dict[str, str]:
    """Hash every exact source byte stream retained in a rendered summary."""
    sources = {
        "config": paths.config,
        "features_csv": paths.features,
        "feature_manifest": paths.feature_manifest,
        "model_table": paths.model_table,
        "study_oof": paths.study_oof,
        "variant_oof": paths.variant_oof,
        "metrics": paths.metrics,
        "sensitivity": paths.sensitivity,
        "tuning_records": paths.tuning_records,
        "run_manifest": paths.run_manifest,
    }
    try:
        return {name: _sha256_file(path) for name, path in sources.items()}
    except OSError as exc:
        raise ValueError("could not hash a summary source artifact") from exc


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(character in "0123456789abcdef" for character in value)


def _read_csv(path: Path | str, name: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError) as exc:
        raise ValueError(f"could not read {name}") from exc


def _read_mapping(path: Path | str, name: str) -> Mapping[str, Any]:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"could not read {name}") from exc
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a JSON mapping")
    return value


def _finite_vector(value: Any, *, expected_length: int, name: str, positive: bool = False) -> np.ndarray:
    try:
        values = np.asarray(value, dtype=float).reshape(-1)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite numeric vector") from exc
    if len(values) != expected_length or not np.isfinite(values).all() or (positive and np.any(values <= 0.0)):
        raise ValueError(f"{name} must be a finite numeric vector")
    return values


def _finite_number(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
        raise ValueError(f"{name} must be a finite numeric value")
    return float(value)


def _unique_mapping(pairs: list[tuple[str, Any]], name: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"{name} contains duplicate key {key!r}")
        result[key] = value
    return result


def _parse_mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, str):
        raise ValueError(f"{name} must be a JSON object")
    try:
        parsed = json.loads(value, object_pairs_hook=lambda pairs: _unique_mapping(pairs, name))
    except json.JSONDecodeError as exc:
        raise ValueError(f"{name} must be valid JSON") from exc
    if not isinstance(parsed, Mapping):
        raise ValueError(f"{name} must be a JSON object")
    return parsed


def _parse_vector(value: Any, width: int, name: str, *, positive: bool = False) -> np.ndarray:
    if not isinstance(value, str):
        raise ValueError(f"{name} must be a JSON numeric vector")
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{name} must be valid JSON") from exc
    if not isinstance(parsed, list):
        raise ValueError(f"{name} must be a JSON numeric vector")
    return _finite_vector(parsed, expected_length=width, name=name, positive=positive)


def _validate_run_manifest(
    manifest: Mapping[str, Any], paths: SummaryPaths, registry: Any, observed_hashes: Mapping[str, str]
) -> None:
    required = {"active_profile", "artifact_sha256", "command", "fallback_reason", "git_revision", "input_sha256", "post_execution_provenance_repair", "schema_version", "seed", "timestamp_utc", "tool_versions", "variants"}
    if set(manifest) != required or manifest["schema_version"] != 1:
        raise ValueError("run manifest has an unexpected frozen schema")
    if manifest["active_profile"] != registry.active_profile or manifest["fallback_reason"] != registry.fallback_reason:
        raise ValueError("run manifest does not identify the frozen active profile")
    if manifest["variants"] != list(FROZEN_PRIMARY_VARIANTS) or manifest["seed"] != registry.seed:
        raise ValueError("run manifest does not preserve frozen cohort or seed")
    if not isinstance(manifest["git_revision"], str) or len(manifest["git_revision"]) != 40 or not isinstance(manifest["command"], list) or "diagnostic" not in manifest["command"]:
        raise ValueError("run manifest lacks historical diagnostic execution provenance")
    hashes = manifest["input_sha256"]
    expected_input = {"canonical_labels", "feature_manifest", "features_csv", "measurement_evidence", "registry"}
    if not isinstance(hashes, Mapping) or set(hashes) != expected_input or not all(_is_sha256(value) for value in hashes.values()):
        raise ValueError("run manifest has incomplete input SHA-256 provenance")
    artifacts = manifest["artifact_sha256"]
    if not isinstance(artifacts, Mapping) or set(artifacts) != _MANIFEST_ARTIFACTS or not all(_is_sha256(value) for value in artifacts.values()):
        raise ValueError("run manifest has incomplete artifact SHA-256 provenance")
    repair = manifest["post_execution_provenance_repair"]
    repair_keys = {"repair_type", "execution_git_head", "execution_worktree_state", "source_snapshot_method", "source_sha256", "inner_fold_audit", "derived_model_table"}
    if not isinstance(repair, Mapping) or set(repair) != repair_keys or repair.get("repair_type") != "artifact_only_no_diagnostic_rerun" or repair.get("execution_git_head") != manifest["git_revision"] or not isinstance(repair.get("execution_worktree_state"), str) or not repair["execution_worktree_state"] or not isinstance(repair.get("source_snapshot_method"), str) or not repair["source_snapshot_method"]:
        raise ValueError("run manifest post-execution provenance is inconsistent")
    source_hashes = repair["source_sha256"]
    if not isinstance(source_hashes, Mapping) or not all(isinstance(key, str) and _is_sha256(value) for key, value in source_hashes.items()) or source_hashes.get("config/kcat-rel-v0.json") != hashes["registry"]:
        raise ValueError("run manifest source snapshot does not bind the frozen registry")
    inner, derived = repair["inner_fold_audit"], repair["derived_model_table"]
    if not isinstance(inner, Mapping) or inner.get("method") != "post_execution_reconstruction" or inner.get("original_execution_revision") != manifest["git_revision"] or inner.get("diagnostic_rerun") is not False or inner.get("scientific_result_files_changed") is not False:
        raise ValueError("run manifest asserts invalid diagnostic rerun provenance")
    if not isinstance(derived, Mapping) or derived.get("path") != "data/derived/kcat-rel-v0-no-msa-model-table.csv" or derived.get("sha256") != artifacts["model-table.csv"] or derived.get("byte_identical_to") != "results/diagnostic-v0/model-table.csv":
        raise ValueError("run manifest does not bind the derived model table")
    expected_hashes = {
        "config": hashes["registry"],
        "features_csv": hashes["features_csv"],
        "feature_manifest": hashes["feature_manifest"],
        "model_table": artifacts["model-table.csv"],
        "study_oof": artifacts["study-oof.csv"],
        "variant_oof": artifacts["variant-oof.csv"],
        "metrics": artifacts["metrics.json"],
        "sensitivity": artifacts["sensitivity.json"],
        "tuning_records": artifacts["tuning-records.csv"],
        "run_manifest": FROZEN_RUN_MANIFEST_SHA256,
    }
    display_names = {
        "config": "registry", "features_csv": "features", "feature_manifest": "feature manifest",
        "model_table": "model table", "study_oof": "study OOF", "variant_oof": "variant OOF",
        "metrics": "metrics", "sensitivity": "sensitivity", "tuning_records": "tuning records",
        "run_manifest": "run manifest",
    }
    for name, expected in expected_hashes.items():
        if observed_hashes.get(name) != expected:
            raise ValueError(f"{display_names[name]} SHA-256 mismatch against frozen run manifest")


def _load_audited_features(paths: SummaryPaths, registry: Any) -> pd.DataFrame:
    features = _read_csv(paths.features, "features")
    feature_manifest = _read_mapping(paths.feature_manifest, "feature manifest")
    try:
        _audit_feature_table(features, registry)
        _audit_manifest_contract(feature_manifest, features, Path(paths.features).read_bytes(), Path(paths.features))
    except (OSError, TypeError, ValueError) as exc:
        raise ValueError("feature artifact pair fails the frozen audit") from exc
    structure = feature_manifest.get("structure")
    if not isinstance(structure, Mapping) or structure.get("compressed_checksum_verified") is not True or structure.get("input_sha256") != registry.structure["compressed_sha256"] or structure.get("expected_compressed_sha256") != registry.structure["compressed_sha256"]:
        raise ValueError("feature artifact pair lacks pinned 8ACT provenance")
    return features.copy()


def _validate_model_table(table: pd.DataFrame, features: pd.DataFrame) -> None:
    if tuple(table.columns) != _MODEL_TABLE_COLUMNS or tuple(table["variant"]) != FROZEN_PRIMARY_VARIANTS:
        raise ValueError("model table does not preserve the exact frozen schema and cohort order")
    if table["variant"].duplicated().any() or any(not isinstance(value, str) or not value for value in table["source_id"]):
        raise ValueError("model table has invalid variant or source identities")
    if set(table["source_id"]) != set(FROZEN_SOURCE_STUDIES):
        raise ValueError("model table sources do not match the frozen source studies")
    try:
        numeric = table.loc[:, ["y", "y_error", *ACTIVE_FEATURE_COLUMNS]].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("model table labels and predictors must be numeric") from exc
    if not np.isfinite(numeric).all():
        raise ValueError("model table contains non-finite labels, errors, or predictors")
    if not table.loc[:, ["variant", *ACTIVE_FEATURE_COLUMNS]].equals(features.loc[:, ["variant", *ACTIVE_FEATURE_COLUMNS]]):
        raise ValueError("model table predictors differ from audited frozen features")


def _validate_parameters(value: Any, model_kind: str, registry: Any, name: str) -> None:
    params = _parse_mapping(value, name) if isinstance(value, str) else value
    if not isinstance(params, Mapping):
        raise ValueError(f"{name} must be a parameter mapping")
    if model_kind == "ridge":
        candidates = tuple({"alpha": float(alpha)} for alpha in registry.model_grids["ridge_alphas"])
    else:
        gpr = registry.model_grids["gpr"]
        candidates = tuple({"length_scale": float(length), "noise_variance": float(noise)} for length in gpr["length_scales"] for noise in gpr["noise_variances"])
    if any(type(item) not in (int, float) or not np.isfinite(item) for item in params.values()) or dict(params) not in candidates:
        raise ValueError(f"{name} is not an approved frozen hyperparameter candidate")


def _expected_oof_rows(model_table: pd.DataFrame, view: str) -> list[tuple[str, str, str, int]]:
    expected: list[tuple[str, str, str, int]] = []
    for model_kind in _OOF_ROW_MODEL_KINDS:
        for fold in range(9 if view == "study" else 17):
            indices = [fold] if view == "variant" else [index for index, source in enumerate(model_table["source_id"]) if source == FROZEN_SOURCE_STUDIES[fold]]
            expected.extend((model_kind, model_table.iloc[index]["variant"], model_table.iloc[index]["source_id"], fold) for index in indices)
    return expected


def _validate_oof(table: pd.DataFrame, model_table: pd.DataFrame, view: str, registry: Any) -> None:
    if tuple(table.columns) != _OOF_COLUMNS:
        raise ValueError(f"{view} OOF has an unexpected schema")
    actual = list(table.loc[:, ["model_kind", "variant", "source_id", "outer_fold"]].itertuples(index=False, name=None))
    if actual != _expected_oof_rows(model_table, view):
        raise ValueError(f"{view} OOF row order, identities, or outer-fold semantics differ from the frozen protocol")
    expected = model_table.set_index("variant")
    for model_kind in _MODEL_KINDS:
        rows = table.loc[table["model_kind"] == model_kind].set_index("variant").loc[list(FROZEN_PRIMARY_VARIANTS)]
        if not np.array_equal(rows["truth"].to_numpy(dtype=float), expected["y"].to_numpy(dtype=float)):
            raise ValueError(f"{view} OOF truth does not match the frozen model table")
        try:
            values = rows.loc[:, ["truth", "model_prediction", "baseline_prediction"]].to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{view} OOF predictions must be numeric") from exc
        if not np.isfinite(values).all():
            raise ValueError(f"{view} OOF contains non-finite predictions")
    for row in table.itertuples(index=False):
        _validate_parameters(row.selected_inner_hyperparameters, row.model_kind, registry, f"{view} OOF hyperparameters")
        _parse_vector(row.scaler_mean, len(ACTIVE_FEATURE_COLUMNS), f"{view} OOF scaler mean")
        _parse_vector(row.scaler_scale, len(ACTIVE_FEATURE_COLUMNS), f"{view} OOF scaler scale", positive=True)


def _parse_indices(value: Any, name: str) -> tuple[int, ...]:
    if not isinstance(value, str):
        raise ValueError(f"{name} must be a JSON integer vector")
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{name} must be valid JSON") from exc
    if not isinstance(parsed, list) or not parsed or any(type(item) is not int for item in parsed) or len(set(parsed)) != len(parsed):
        raise ValueError(f"{name} must be a non-empty unique integer vector")
    return tuple(parsed)


def _candidate_name(params: Mapping[str, Any], model_kind: str) -> str:
    if model_kind == "ridge":
        return f"alpha={float(params['alpha']):g}"
    return f"length_scale={float(params['length_scale']):g},noise_variance={float(params['noise_variance']):g}"


def _expected_fold_indices(model_table: pd.DataFrame, view: str, fold: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if view == "study":
        test = tuple(index for index, source in enumerate(model_table["source_id"]) if source == FROZEN_SOURCE_STUDIES[fold])
    else:
        test = (fold,)
    return tuple(index for index in range(len(model_table)) if index not in test), test


def _validate_tuning_records(
    tuning: pd.DataFrame, model_table: pd.DataFrame, study_oof: pd.DataFrame, variant_oof: pd.DataFrame, registry: Any
) -> None:
    columns = ("model_kind", "view", "outer_fold", "train_indices", "test_indices", "selected_params", "inner_scores", "scaler_mean", "scaler_scale", "baseline_mean")
    if tuple(tuning.columns) != columns or len(tuning) != 52:
        raise ValueError("tuning records do not have the exact frozen schema or row count")
    for model_kind in _MODEL_KINDS:
        for view, folds, oof in (("study", 9, study_oof), ("variant", 17, variant_oof)):
            rows = tuning.loc[(tuning["model_kind"] == model_kind) & (tuning["view"] == view)]
            if len(rows) != folds or tuple(rows["outer_fold"]) != tuple(range(folds)):
                raise ValueError("tuning records lack complete ordered model/view/fold coverage")
            for record in rows.itertuples(index=False):
                train, test = _parse_indices(record.train_indices, "tuning train indices"), _parse_indices(record.test_indices, "tuning test indices")
                expected_train, expected_test = _expected_fold_indices(model_table, view, record.outer_fold)
                if train != expected_train or test != expected_test:
                    raise ValueError("tuning train/test indices do not match the frozen outer-fold protocol")
                params = _parse_mapping(record.selected_params, "tuning selected hyperparameters")
                _validate_parameters(params, model_kind, registry, "tuning selected hyperparameters")
                scores = _parse_mapping(record.inner_scores, "tuning inner scores")
                candidates = (
                    ({"alpha": float(alpha)} for alpha in registry.model_grids["ridge_alphas"])
                    if model_kind == "ridge"
                    else ({"length_scale": float(length), "noise_variance": float(noise)} for length in registry.model_grids["gpr"]["length_scales"] for noise in registry.model_grids["gpr"]["noise_variances"])
                )
                candidate_map = {_candidate_name(candidate, model_kind): candidate for candidate in candidates}
                if set(scores) != set(candidate_map) or any(_finite_number(value, "tuning inner score") < 0.0 for value in scores.values()):
                    raise ValueError("tuning inner scores do not cover the exact frozen parameter grid")
                best_score = min(float(value) for value in scores.values())
                tied = [candidate_map[name] for name, value in scores.items() if float(value) == best_score]
                selected = max(tied, key=lambda item: tuple(item.values()))
                if dict(params) != selected:
                    raise ValueError("tuning selected hyperparameters are not the deterministic frozen minimizer")
                mean = _parse_vector(record.scaler_mean, len(ACTIVE_FEATURE_COLUMNS), "tuning scaler mean")
                scale = _parse_vector(record.scaler_scale, len(ACTIVE_FEATURE_COLUMNS), "tuning scaler scale", positive=True)
                train_table = model_table.iloc[list(train)]
                expected_mean = train_table.loc[:, ACTIVE_FEATURE_COLUMNS].to_numpy(float).mean(axis=0)
                expected_scale = train_table.loc[:, ACTIVE_FEATURE_COLUMNS].to_numpy(float).std(axis=0, ddof=0)
                if not np.allclose(mean, expected_mean, rtol=0, atol=1e-12) or not np.allclose(scale, expected_scale, rtol=0, atol=1e-12):
                    raise ValueError("tuning scaler state is not reconstructed from frozen training rows")
                baseline = _finite_number(record.baseline_mean, "tuning baseline mean")
                if not np.isclose(baseline, train_table["y"].mean(), rtol=0, atol=1e-12):
                    raise ValueError("tuning baseline mean is not the frozen outer-training mean")
                oof_rows = oof.loc[(oof["model_kind"] == model_kind) & (oof["outer_fold"] == record.outer_fold)]
                if any(_parse_mapping(value, "OOF selected hyperparameters") != params for value in oof_rows["selected_inner_hyperparameters"]) or any(not np.allclose(_parse_vector(value, len(ACTIVE_FEATURE_COLUMNS), "OOF scaler mean"), mean, rtol=0, atol=1e-12) for value in oof_rows["scaler_mean"]) or any(not np.allclose(_parse_vector(value, len(ACTIVE_FEATURE_COLUMNS), "OOF scaler scale", positive=True), scale, rtol=0, atol=1e-12) for value in oof_rows["scaler_scale"]) or not np.allclose(oof_rows["baseline_prediction"].to_numpy(float), baseline, rtol=0, atol=1e-12):
                    raise ValueError("OOF fold state does not match the reconstructed tuning protocol")


def load_summary_inputs(paths: SummaryPaths) -> SummaryInputs:
    """Load byte-pinned artifacts, then validate their frozen semantic contract."""
    artifact_hashes = _summary_artifact_hashes(paths)
    if artifact_hashes["run_manifest"] != FROZEN_RUN_MANIFEST_SHA256:
        raise ValueError("run manifest immutable SHA-256 mismatch")
    registry = load_registry(paths.config)
    run_manifest = _read_mapping(paths.run_manifest, "run manifest")
    _validate_run_manifest(run_manifest, paths, registry, artifact_hashes)
    features = _load_audited_features(paths, registry)
    model_table = _read_csv(paths.model_table, "model table")
    _validate_model_table(model_table, features)
    study_oof, variant_oof = _read_csv(paths.study_oof, "study OOF"), _read_csv(paths.variant_oof, "variant OOF")
    tuning_records = _read_csv(paths.tuning_records, "tuning records")
    _validate_oof(study_oof, model_table, "study", registry)
    _validate_oof(variant_oof, model_table, "variant", registry)
    _validate_tuning_records(tuning_records, model_table, study_oof, variant_oof, registry)
    return SummaryInputs(
        paths, registry, features, model_table, study_oof, variant_oof, tuning_records,
        _read_mapping(paths.metrics, "metrics"), _read_mapping(paths.sensitivity, "sensitivity"), run_manifest,
        artifact_hashes, artifact_hashes["run_manifest"],
    )


def compute_feature_correlations(features: pd.DataFrame, y: np.ndarray) -> dict[str, float]:
    """Calculate descriptive Pearson correlations in the frozen feature order."""
    if tuple(features.columns) != ACTIVE_FEATURE_COLUMNS:
        raise ValueError("summary features do not match the frozen fallback registry")
    target = _finite_vector(y, expected_length=len(features), name="y")
    if np.ptp(target) == 0.0:
        raise ValueError("y: Pearson r is undefined for a constant target")
    result = {}
    for name in ACTIVE_FEATURE_COLUMNS:
        values = _finite_vector(features[name], expected_length=len(features), name=name)
        if np.ptp(values) == 0.0:
            raise ValueError(f"{name}: Pearson r is undefined for a constant feature")
        result[name] = float(np.corrcoef(values, target)[0, 1])
    return result


def _metrics_value(value: Any, name: str) -> dict[str, float]:
    if not isinstance(value, Mapping) or set(value) != {"mae", "rmse"}:
        raise ValueError(f"{name} must have MAE and RMSE")
    return {key: _finite_number(value[key], f"{name} {key}") for key in ("mae", "rmse")}


def _variant_metrics(truth: Any, prediction: Any) -> dict[str, float]:
    y, estimate = np.asarray(truth, dtype=float).reshape(-1), np.asarray(prediction, dtype=float).reshape(-1)
    if not len(y) or len(y) != len(estimate) or not np.isfinite(y).all() or not np.isfinite(estimate).all():
        raise ValueError("variant metrics require equal non-empty finite arrays")
    residual = y - estimate
    return {"mae": float(np.mean(np.abs(residual))), "rmse": float(np.sqrt(np.mean(np.square(residual))))}


def _study_macro_metrics(truth: Any, prediction: Any, groups: Any) -> dict[str, float]:
    y, estimate, source = np.asarray(truth, dtype=float).reshape(-1), np.asarray(prediction, dtype=float).reshape(-1), np.asarray(groups, dtype=object).reshape(-1)
    if not len(y) or len(y) != len(estimate) or len(y) != len(source) or not np.isfinite(y).all() or not np.isfinite(estimate).all() or any(not isinstance(value, str) or not value for value in source):
        raise ValueError("study metrics require equal finite arrays and non-empty source identities")
    per_study = [_variant_metrics(y[source == study], estimate[source == study]) for study in sorted(set(source), key=str)]
    return {key: float(np.mean([metric[key] for metric in per_study])) for key in ("mae", "rmse")}


def _view_metrics(rows: pd.DataFrame, view: str) -> tuple[dict[str, float], dict[str, float]]:
    if view == "study":
        return _study_macro_metrics(rows["truth"], rows["model_prediction"], rows["source_id"]), _study_macro_metrics(rows["truth"], rows["baseline_prediction"], rows["source_id"])
    return _variant_metrics(rows["truth"], rows["model_prediction"]), _variant_metrics(rows["truth"], rows["baseline_prediction"])


def _fraction(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else (0.0 if numerator == 0.0 else float(np.copysign(np.inf, numerator)))


def _gate_values(model: Mapping[str, Mapping[str, float]], baseline: Mapping[str, Mapping[str, float]], study_rows: pd.DataFrame, registry: Any) -> tuple[dict[str, float], dict[str, float], tuple[str, ...]]:
    improvement = {view: _fraction(baseline[view]["mae"] - model[view]["mae"], baseline[view]["mae"]) for view in ("study", "variant")}
    degradation = {view: _fraction(model[view]["rmse"] - baseline[view]["rmse"], baseline[view]["rmse"]) for view in ("study", "variant")}
    thresholds = getattr(registry, "gate_thresholds", {}) if registry is not None else {}
    failures = []
    if improvement["study"] < float(thresholds.get("minimum_mae_improvement_fraction", 0.05)): failures.append("study_macro_mae")
    if improvement["variant"] < float(thresholds.get("minimum_mae_improvement_fraction", 0.05)): failures.append("variant_mae")
    if degradation["study"] > float(thresholds.get("maximum_rmse_degradation_fraction", 0.05)): failures.append("study_macro_rmse")
    if degradation["variant"] > float(thresholds.get("maximum_rmse_degradation_fraction", 0.05)): failures.append("variant_rmse")
    sources = np.asarray(study_rows["source_id"], dtype=object)
    influence = {str(source): _study_macro_metrics(study_rows.loc[sources != source, "truth"], study_rows.loc[sources != source, "baseline_prediction"], sources[sources != source])["mae"] - _study_macro_metrics(study_rows.loc[sources != source, "truth"], study_rows.loc[sources != source, "model_prediction"], sources[sources != source])["mae"] for source in sorted(set(sources), key=str)}
    expected_sources = set(getattr(registry, "source_studies", FROZEN_SOURCE_STUDIES))
    if set(influence) != expected_sources or not all(np.isfinite(value) and value > 0.0 for value in influence.values()): failures.append("single_study_influence")
    return improvement, degradation, tuple(failures)


def _strict_numeric_list(value: Any, length: int, name: str, *, positive: bool = False) -> np.ndarray:
    if not isinstance(value, list) or len(value) != length or any(isinstance(item, bool) or not isinstance(item, (int, float)) for item in value):
        raise ValueError(f"{name} must be an exact finite numeric sequence")
    return _finite_vector(value, expected_length=length, name=name, positive=positive)


def _summary_sensitivity(
    sensitivity: Mapping[str, Any], model_table: pd.DataFrame, tuning: pd.DataFrame, registry: Any
) -> tuple[dict[str, Any], dict[str, Any]]:
    if not isinstance(sensitivity, Mapping) or set(sensitivity) != set(_MODEL_KINDS):
        raise ValueError("sensitivity must contain exact Ridge and GPR model keys")
    truth, sources = model_table["y"].to_numpy(dtype=float), model_table["source_id"].to_numpy()
    weighted_summary: dict[str, Any] = {}; perturbation_summary: dict[str, Any] = {}
    for model_kind in _MODEL_KINDS:
        model = sensitivity[model_kind]
        if not isinstance(model, Mapping) or set(model) != {"weighted", "perturbation"} or not isinstance(model["weighted"], Mapping) or set(model["weighted"]) != {"study", "variant"}:
            raise ValueError("sensitivity must retain exact model and view provenance")
        weighted_summary[model_kind] = {}
        for view in ("study", "variant"):
            item = model["weighted"][view]
            required = {"model_kind", "view", "prediction", "baseline_prediction", "metrics", "baseline_metrics", "fold_hyperparameters", "fold_sample_weights", "fold_reported_variances"}
            if not isinstance(item, Mapping) or set(item) != required or item["model_kind"] != model_kind or item["view"] != view:
                raise ValueError("weighted sensitivity model/view provenance is invalid")
            prediction = _finite_vector(item["prediction"], expected_length=len(truth), name="weighted sensitivity prediction")
            baseline = _finite_vector(item["baseline_prediction"], expected_length=len(truth), name="weighted sensitivity baseline")
            expected_model = _study_macro_metrics(truth, prediction, sources) if view == "study" else _variant_metrics(truth, prediction)
            expected_baseline = _study_macro_metrics(truth, baseline, sources) if view == "study" else _variant_metrics(truth, baseline)
            model_metrics, baseline_metrics = _metrics_value(item["metrics"], "weighted sensitivity model metrics"), _metrics_value(item["baseline_metrics"], "weighted sensitivity baseline metrics")
            if any(not np.isclose(model_metrics[key], expected_model[key], rtol=0, atol=1e-12) or not np.isclose(baseline_metrics[key], expected_baseline[key], rtol=0, atol=1e-12) for key in ("mae", "rmse")):
                raise ValueError("weighted sensitivity metrics do not match stored predictions")
            fold_count = 9 if view == "study" else 17
            if any(not isinstance(item[key], list) or len(item[key]) != fold_count for key in ("fold_hyperparameters", "fold_sample_weights", "fold_reported_variances")):
                raise ValueError("weighted sensitivity lacks frozen fold provenance")
            primary_folds = tuning.loc[(tuning["model_kind"] == model_kind) & (tuning["view"] == view)].sort_values("outer_fold")
            if tuple(primary_folds["outer_fold"]) != tuple(range(fold_count)):
                raise ValueError("weighted sensitivity primary fold coverage is invalid")
            expected_baseline = np.empty(len(truth), dtype=float)
            for position, primary in enumerate(primary_folds.itertuples(index=False)):
                params = item["fold_hyperparameters"][position]
                _validate_parameters(params, model_kind, registry, "weighted sensitivity hyperparameters")
                if dict(params) != dict(_parse_mapping(primary.selected_params, "tuning selected hyperparameters")):
                    raise ValueError("weighted sensitivity hyperparameters differ from the primary fold")
                train, test = _parse_indices(primary.train_indices, "tuning train indices"), _parse_indices(primary.test_indices, "tuning test indices")
                expected_baseline[list(test)] = _finite_number(primary.baseline_mean, "tuning baseline mean")
                errors = model_table.iloc[list(train)]["y_error"].to_numpy(dtype=float)
                weights, variances = item["fold_sample_weights"][position], item["fold_reported_variances"][position]
                if model_kind == "ridge":
                    if np.any(errors <= 0.0):
                        raise ValueError("weighted sensitivity Ridge training errors must be positive")
                    expected_weights = 1.0 / np.square(errors)
                    expected_weights /= expected_weights.mean()
                    observed_weights = _strict_numeric_list(weights, len(train), "weighted sensitivity Ridge weights", positive=True)
                    if not np.allclose(observed_weights, expected_weights, rtol=0, atol=1e-12):
                        raise ValueError("weighted sensitivity Ridge weights differ from frozen training uncertainty")
                    if not isinstance(variances, list) or variances:
                        raise ValueError("weighted sensitivity Ridge reported variances must be empty lists")
                else:
                    if not isinstance(weights, list) or weights:
                        raise ValueError("weighted sensitivity GPR weights must be empty lists")
                    observed_variances = _strict_numeric_list(variances, len(train), "weighted sensitivity GPR variances")
                    if np.any(observed_variances < 0.0) or not np.allclose(observed_variances, np.square(errors), rtol=0, atol=1e-12):
                        raise ValueError("weighted sensitivity GPR variances differ from frozen training uncertainty")
            if not np.allclose(baseline, expected_baseline, rtol=0, atol=1e-12):
                raise ValueError("weighted sensitivity baseline differs from primary outer-training means")
            weighted_summary[model_kind][view] = {"model_metrics": model_metrics, "baseline_metrics": baseline_metrics}
        perturbation = model["perturbation"]
        required = {"draws", "seed", "metric_quantiles", "gate_pass_count", "gate_pass_rate", "fold_hyperparameters"}
        if not isinstance(perturbation, Mapping) or set(perturbation) != required or type(perturbation["draws"]) is not int or type(perturbation["seed"]) is not int or type(perturbation["gate_pass_count"]) is not int or type(perturbation["gate_pass_rate"]) not in (int, float):
            raise ValueError("perturbation sensitivity headline types are invalid")
        draws, count, rate = perturbation["draws"], perturbation["gate_pass_count"], float(perturbation["gate_pass_rate"])
        if draws != 1000 or perturbation["seed"] != registry.seed or count < 0 or count > draws or not np.isfinite(rate) or not np.isclose(rate, count / draws, rtol=0, atol=1e-12):
            raise ValueError("perturbation sensitivity headline is inconsistent")
        if not isinstance(perturbation["fold_hyperparameters"], Mapping) or set(perturbation["fold_hyperparameters"]) != {"study", "variant"}:
            raise ValueError("perturbation sensitivity lacks model/view provenance")
        for view, folds in (("study", 9), ("variant", 17)):
            params = perturbation["fold_hyperparameters"][view]
            if not isinstance(params, list) or len(params) != folds: raise ValueError("perturbation sensitivity lacks frozen fold provenance")
            primary = tuning.loc[(tuning["model_kind"] == model_kind) & (tuning["view"] == view)].sort_values("outer_fold")
            expected = [_parse_mapping(value, "tuning selected hyperparameters") for value in primary["selected_params"]]
            for item in params: _validate_parameters(item, model_kind, registry, "perturbation sensitivity hyperparameters")
            if params != expected:
                raise ValueError("perturbation sensitivity hyperparameters differ from primary tuning")
        perturbation_summary[model_kind] = {"draws": draws, "seed": perturbation["seed"], "gate_pass_count": count, "gate_pass_rate": rate}
    return weighted_summary, perturbation_summary


def compute_summary_checks(inputs: SummaryInputs) -> SummaryChecks:
    """Recompute frozen summary scores with local pure arithmetic only."""
    metrics = inputs.metrics
    if not isinstance(metrics, Mapping) or set(metrics) != {"winner", "models"} or metrics.get("winner") is not None or not isinstance(metrics["models"], Mapping) or set(metrics["models"]) != set(_MODEL_KINDS):
        raise ValueError("frozen metrics must retain exact models and a null winner")
    model_metrics: dict[str, Any] = {}; baseline_metrics: dict[str, Any] = {}; improvements: dict[str, Any] = {}; degradations: dict[str, Any] = {}; rules: dict[str, tuple[str, ...]] = {}; gates: dict[str, bool] = {}
    for model_kind in _MODEL_KINDS:
        decision = metrics["models"][model_kind]
        required = {"pass_gate", "failed_rules", "mae_improvement", "rmse_degradation", "model_views", "baseline_views"}
        if not isinstance(decision, Mapping) or set(decision) != required or decision["pass_gate"] is not False or not isinstance(decision["failed_rules"], list) or set(decision["model_views"]) != {"study", "variant"} or set(decision["baseline_views"]) != {"study", "variant"}:
            raise ValueError("frozen primary model gate schema is invalid")
        model_metrics[model_kind] = {}; baseline_metrics[model_kind] = {}
        for view, table in (("study", inputs.study_oof), ("variant", inputs.variant_oof)):
            actual_model, actual_baseline = _view_metrics(table.loc[table["model_kind"] == model_kind], view)
            stored_model, stored_baseline = _metrics_value(decision["model_views"][view], f"{model_kind} {view} model metrics"), _metrics_value(decision["baseline_views"][view], f"{model_kind} {view} baseline metrics")
            if any(not np.isclose(stored_model[key], actual_model[key], rtol=0, atol=1e-12) or not np.isclose(stored_baseline[key], actual_baseline[key], rtol=0, atol=1e-12) for key in ("mae", "rmse")):
                raise ValueError(f"{model_kind} {view} metrics do not match frozen OOF predictions")
            model_metrics[model_kind][view], baseline_metrics[model_kind][view] = stored_model, stored_baseline
        improvement, degradation, failed = _gate_values(model_metrics[model_kind], baseline_metrics[model_kind], inputs.study_oof.loc[inputs.study_oof["model_kind"] == model_kind], inputs.registry)
        for name, observed, expected in (("mae_improvement", decision["mae_improvement"], improvement), ("rmse_degradation", decision["rmse_degradation"], degradation)):
            if not isinstance(observed, Mapping) or set(observed) != {"study", "variant"} or any(not np.isclose(_finite_number(observed[view], name), expected[view], rtol=0, atol=1e-12) for view in expected):
                raise ValueError("metrics relative-score summary does not match frozen OOF predictions")
        if decision["failed_rules"] != list(failed): raise ValueError("metrics failed rules do not match frozen OOF predictions")
        improvements[model_kind], degradations[model_kind], rules[model_kind], gates[model_kind] = improvement, degradation, failed, False
    if inputs.tuning_records is None:
        weighted, perturbation = {}, {}
    else:
        weighted, perturbation = _summary_sensitivity(inputs.sensitivity, inputs.model_table, inputs.tuning_records, inputs.registry)
    errors = {
        str(variant): float(error)
        for variant, error in zip(inputs.model_table["variant"], inputs.model_table["y_error"], strict=True)
    }
    return SummaryChecks(
        compute_feature_correlations(inputs.features.loc[:, ACTIVE_FEATURE_COLUMNS], inputs.model_table["y"].to_numpy(dtype=float)),
        inputs.study_oof.copy(deep=True), inputs.variant_oof.copy(deep=True),
        model_metrics, baseline_metrics, improvements, degradations, rules, gates, weighted, perturbation,
        errors, len(inputs.model_table), len(set(inputs.model_table["source_id"])),
        dict(inputs.artifact_sha256), inputs.run_manifest_root_sha256,
    )


def summary_payload(checks: SummaryChecks) -> dict[str, Any]:
    """Return only JSON-safe, descriptive values retained by the summary checks."""
    return {
        "status": "posthoc_descriptive_only",
        "cohort_size": checks.cohort_size,
        "source_study_count": checks.source_study_count,
        "run_manifest_root_sha256": checks.run_manifest_root_sha256,
        "artifact_sha256": dict(checks.artifact_sha256),
        "feature_names": list(checks.feature_correlations),
        "feature_correlations": dict(checks.feature_correlations),
        "label_uncertainties": dict(checks.label_uncertainties),
        "study_oof": checks.study_oof.to_dict(orient="records"),
        "variant_oof": checks.variant_oof.to_dict(orient="records"),
        "model_metrics": checks.model_metrics,
        "baseline_metrics": checks.baseline_metrics,
        "mae_improvement": checks.mae_improvement,
        "rmse_degradation": checks.rmse_degradation,
        "failed_rules": {model: list(value) for model, value in checks.failed_rules.items()},
        "pass_gate": checks.pass_gate,
        "weighted_sensitivity": checks.weighted_sensitivity,
        "perturbation_sensitivity": checks.perturbation_sensitivity,
    }


def preflight_summary_output(output_dir: Path | str) -> Path:
    """Require an output destination that is missing or already empty."""
    destination = Path(output_dir)
    if destination.exists():
        if not destination.is_dir():
            raise ValueError("summary output path must be a directory")
        if any(destination.iterdir()):
            raise ValueError("summary output directory must be missing or empty")
    if not destination.parent.is_dir():
        raise ValueError("summary output parent directory must already exist")
    return destination


def _fsync_file(path: Path) -> None:
    with path.open("rb") as handle:
        os.fsync(handle.fileno())


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write_summary_json(path: Path, checks: SummaryChecks) -> None:
    text = json.dumps(summary_payload(checks), indent=2, sort_keys=True, allow_nan=False) + "\n"
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())


def write_summary_outputs(checks: SummaryChecks, output_dir: Path | str) -> tuple[Path, tuple[Path, Path, Path]]:
    """Stage and atomically publish JSON plus deterministic figures from checks only."""
    from kcat_rel.summary_figures import (
        render_error_comparison,
        render_feature_correlations,
        render_study_oof,
    )

    destination = preflight_summary_output(output_dir)
    staging = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging-", dir=destination.parent))
    filenames = (
        "feature-correlations.png",
        "study-oof-observed-vs-predicted.png",
        "held-out-error-vs-baseline.png",
    )
    try:
        json_path = staging / "modeling-summary.json"
        _write_summary_json(json_path, checks)
        figure_paths = tuple(staging / name for name in filenames)
        render_feature_correlations(checks, figure_paths[0])
        render_study_oof(checks, figure_paths[1])
        render_error_comparison(checks, figure_paths[2])
        for path in (json_path, *figure_paths):
            _fsync_file(path)
        _fsync_directory(staging)
        os.replace(staging, destination)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    return destination / "modeling-summary.json", tuple(destination / name for name in filenames)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Render a read-only summary of frozen kcat_rel artifacts.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--feature-manifest", required=True)
    parser.add_argument("--model-table", required=True)
    parser.add_argument("--study-oof", required=True)
    parser.add_argument("--variant-oof", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--sensitivity", required=True)
    parser.add_argument("--tuning-records")
    parser.add_argument("--run-manifest")
    parser.add_argument("--output-dir", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Validate explicit frozen artifacts and render their descriptive summary."""
    args = _parser().parse_args(argv)
    destination = preflight_summary_output(args.output_dir)
    diagnostic_dir = Path(args.metrics).parent
    paths = SummaryPaths(
        config=args.config,
        features=args.features,
        feature_manifest=args.feature_manifest,
        model_table=args.model_table,
        study_oof=args.study_oof,
        variant_oof=args.variant_oof,
        metrics=args.metrics,
        sensitivity=args.sensitivity,
        tuning_records=args.tuning_records or diagnostic_dir / "tuning-records.csv",
        run_manifest=args.run_manifest or diagnostic_dir / "run-manifest.json",
    )
    write_summary_outputs(compute_summary_checks(load_summary_inputs(paths)), destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

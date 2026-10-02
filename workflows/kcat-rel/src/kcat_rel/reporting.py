"""Deterministic, failure-clean reporting for frozen ``kcat_rel`` runs."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import tempfile
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from kcat_rel.registry import ACTIVE_FEATURE_COLUMNS, FROZEN_PRIMARY_VARIANTS, FROZEN_SOURCE_STUDIES
from kcat_rel.validation import Metrics, evaluate_gate, inner_group_splits, select_winner, single_study_influence, study_macro_metrics, variant_metrics


_ARTIFACT_FILENAMES = {
    "run_manifest": "run-manifest.json",
    "model_table": "model-table.csv",
    "tuning_records": "tuning-records.csv",
    "study_oof": "study-oof.csv",
    "variant_oof": "variant-oof.csv",
    "inner_fold_audit": "inner-fold-audit.json",
    "metrics": "metrics.json",
    "sensitivity": "sensitivity.json",
    "audit_markdown": "AUDIT.md",
    "run_log": "run.log",
}
_INPUT_HASH_KEYS = {"features_csv", "feature_manifest", "canonical_labels", "measurement_evidence", "registry"}
_HASHED_OUTPUT_FILENAMES = set(_ARTIFACT_FILENAMES.values()) - {"run-manifest.json"}
_FROZEN_EXECUTION_REVISION = "91affadb0f5bdb3fbf8f09b22209836f38ddf9c5"


def required_artifacts(output: str | Path) -> set[Path]:
    """Return every required completed-run artifact, including the final marker."""
    directory = Path(output)
    return {directory / name for name in (*_ARTIFACT_FILENAMES.values(), "SUCCESS")}


def preflight_diagnostic_output(output: str | Path) -> None:
    """Reject an existing non-empty destination without changing it."""
    target = Path(output)
    if target.exists() and (not target.is_dir() or any(target.iterdir())):
        raise ValueError("diagnostic output must be a new or empty directory")


def _json_value(value: Any) -> Any:
    if is_dataclass(value):
        return _json_value(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, Path):
        return str(value)
    return value


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(_json_value(value), sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def _csv_bytes(table: pd.DataFrame) -> bytes:
    if not isinstance(table, pd.DataFrame):
        raise ValueError("diagnostic table artifacts must be pandas DataFrames")
    stable = table.copy()
    for column in stable.columns:
        if stable[column].dtype == object:
            stable[column] = stable[column].map(
                lambda value: json.dumps(_json_value(value), sort_keys=True, separators=(",", ":"))
                if isinstance(value, (dict, list, tuple, np.ndarray))
                else value
            )
    return stable.to_csv(index=False, lineterminator="\n").encode("utf-8")


def _sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _normalize_sensitivity_containers(value: Any, path: tuple[str, ...] = ()) -> Any:
    """Canonicalize only the native containers retained by CLI ``asdict``.

    Scalar values and mapping keys are deliberately untouched so validation can
    reject bool/numeric aliases instead of coercing them into valid JSON values.
    NumPy vectors are supported only for the two weighted prediction fields.
    """
    if isinstance(value, Mapping):
        return {key: _normalize_sensitivity_containers(item, (*path, key)) for key, item in value.items()}
    if type(value) in (list, tuple):
        return [_normalize_sensitivity_containers(item, path) for item in value]
    if type(value) is np.ndarray and len(path) == 4 and path[1] == "weighted" and path[3] in {"prediction", "baseline_prediction"}:
        if value.ndim != 1 or value.dtype.kind not in "iuf":
            raise ValueError("weighted prediction arrays must be one-dimensional numeric vectors")
        return value.tolist()
    return value


def _validate_result(result: Mapping[str, Any]) -> None:
    if not isinstance(result, Mapping) or set(result) != set(_ARTIFACT_FILENAMES):
        raise ValueError("diagnostic result must contain every required report artifact")
    for name in ("run_manifest", "inner_fold_audit", "metrics", "sensitivity"):
        if not isinstance(result[name], Mapping):
            raise ValueError(f"diagnostic {name} must be a mapping")
    for name in ("model_table", "tuning_records", "study_oof", "variant_oof"):
        if not isinstance(result[name], pd.DataFrame):
            raise ValueError(f"diagnostic {name} must be a DataFrame")
    for name in ("audit_markdown", "run_log"):
        if not isinstance(result[name], str) or not result[name]:
            raise ValueError(f"diagnostic {name} must be non-empty text")
    _validate_run_manifest(result["run_manifest"])
    _validate_model_table(result["model_table"])
    _validate_tuning(result["tuning_records"])
    _validate_oof(result["study_oof"], "study")
    _validate_oof(result["variant_oof"], "variant")
    _validate_metrics(result["metrics"])
    _validate_sensitivity(result["sensitivity"])
    _validate_inner_fold_audit(result["inner_fold_audit"], result["model_table"], result["tuning_records"], result["run_manifest"])
    _validate_cross_artifacts(result)


def _require_columns(table: pd.DataFrame, columns: tuple[str, ...], name: str) -> None:
    if not set(columns).issubset(table.columns):
        raise ValueError(f"{name} lacks required columns")


def _validate_run_manifest(manifest: Mapping[str, Any]) -> None:
    required = {"schema_version", "active_profile", "fallback_reason", "variants", "seed", "input_sha256", "tool_versions", "command", "git_revision", "timestamp_utc"}
    if not required.issubset(manifest) or manifest["schema_version"] != 1 or manifest["variants"] != list(FROZEN_PRIMARY_VARIANTS):
        raise ValueError("run manifest lacks the frozen diagnostic contract")
    hashes = manifest["input_sha256"]
    if not isinstance(hashes, Mapping) or set(hashes) != _INPUT_HASH_KEYS or any(not isinstance(value, str) or len(value) != 64 or any(char not in "0123456789abcdef" for char in value) for value in hashes.values()):
        raise ValueError("run manifest lacks complete input hashes")
    if not isinstance(manifest["tool_versions"], Mapping) or not manifest["tool_versions"] or not isinstance(manifest["command"], list):
        raise ValueError("run manifest lacks command or version provenance")
    output_hashes = manifest.get("artifact_sha256")
    if not isinstance(output_hashes, Mapping) or set(output_hashes) != _HASHED_OUTPUT_FILENAMES or any(not isinstance(value, str) or len(value) != 64 or any(char not in "0123456789abcdef" for char in value) for value in output_hashes.values()):
        raise ValueError("run manifest lacks complete output hashes")


def _validate_model_table(table: pd.DataFrame) -> None:
    columns = ("variant", "y", "y_error", "source_id", *ACTIVE_FEATURE_COLUMNS)
    if tuple(table.columns) != columns or tuple(table["variant"]) != FROZEN_PRIMARY_VARIANTS:
        raise ValueError("model table does not preserve the exact frozen schema and cohort")
    if table["variant"].duplicated().any() or not np.isfinite(table.loc[:, ["y", "y_error", *ACTIVE_FEATURE_COLUMNS]].to_numpy(dtype=float)).all():
        raise ValueError("model table contains duplicate or non-finite values")


def _validate_tuning(table: pd.DataFrame) -> None:
    _require_columns(table, ("model_kind", "view", "outer_fold", "selected_params", "inner_scores", "scaler_mean", "scaler_scale", "baseline_mean"), "tuning records")
    if len(table) != 52 or set(table["model_kind"]) != {"ridge", "gpr"} or any(isinstance(value, bool) or not isinstance(value, (int, np.integer)) for value in table["outer_fold"]):
        raise ValueError("tuning records do not cover both models and all outer folds")
    for model in ("ridge", "gpr"):
        for view, count, folds in (("study", 9, set(range(9))), ("variant", 17, set(range(17)))):
            subset = table[(table.model_kind == model) & (table.view == view)]
            if len(subset) != count or set(subset.outer_fold) != folds:
                raise ValueError("tuning records have incomplete outer-fold coverage")
    numeric = table["baseline_mean"].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError("tuning records contain non-finite values")


def _validate_oof(table: pd.DataFrame, view: str) -> None:
    _require_columns(table, ("model_kind", "variant", "source_id", "truth", "model_prediction", "baseline_prediction", "outer_fold", "selected_inner_hyperparameters", "scaler_mean", "scaler_scale"), f"{view} OOF")
    if len(table) != 34 or set(table.model_kind) != {"ridge", "gpr"} or any(isinstance(value, bool) or not isinstance(value, (int, np.integer)) for value in table["outer_fold"]):
        raise ValueError(f"{view} OOF does not cover both candidate models")
    expected_folds = set(range(9 if view == "study" else 17))
    for model in ("ridge", "gpr"):
        subset = table[table.model_kind == model]
        if len(subset) != 17 or tuple(sorted(subset.variant)) != tuple(sorted(FROZEN_PRIMARY_VARIANTS)) or set(subset.outer_fold) != expected_folds:
            raise ValueError(f"{view} OOF has incomplete cohort or fold coverage")
    if not np.isfinite(table.loc[:, ["truth", "model_prediction", "baseline_prediction"]].to_numpy(dtype=float)).all():
        raise ValueError(f"{view} OOF contains non-finite predictions")
    for row in table.itertuples(index=False):
        _validate_parameters(row.model_kind, _object(row.selected_inner_hyperparameters, "OOF hyperparameters"))


def _validate_metrics(metrics: Mapping[str, Any]) -> None:
    if set(metrics) != {"winner", "models"} or metrics["winner"] not in {None, "ridge", "gpr"} or set(metrics["models"]) != {"ridge", "gpr"}:
        raise ValueError("metrics lacks complete candidate decisions")
    for decision in metrics["models"].values():
        required = {"pass_gate", "failed_rules", "mae_improvement", "rmse_degradation", "model_views", "baseline_views"}
        if not isinstance(decision, Mapping) or not required.issubset(decision) or set(decision["model_views"]) != {"study", "variant"} or set(decision["baseline_views"]) != {"study", "variant"}:
            raise ValueError("metrics has incomplete gate decision")
        values = [metric[key] for views in (decision["model_views"], decision["baseline_views"]) for metric in views.values() for key in ("mae", "rmse")]
        if not np.isfinite(np.asarray(values, dtype=float)).all():
            raise ValueError("metrics contains non-finite scores")


def _validate_sensitivity(sensitivity: Mapping[str, Any]) -> None:
    if set(sensitivity) != {"ridge", "gpr"}:
        raise ValueError("sensitivity lacks both candidate models")
    for model_kind, model in sensitivity.items():
        if not isinstance(model, Mapping) or set(model) != {"weighted", "perturbation"} or not isinstance(model["weighted"], Mapping) or not isinstance(model["perturbation"], Mapping):
            raise ValueError("sensitivity lacks required frozen secondary views")
        if set(model["weighted"]) != {"study", "variant"}:
            raise ValueError("sensitivity lacks complete weighted outer views")
        for view, weighted in model["weighted"].items():
            required = {"prediction", "baseline_prediction", "metrics", "baseline_metrics", "fold_hyperparameters", "fold_sample_weights", "fold_reported_variances"}
            if not isinstance(weighted, Mapping) or set(weighted) != required | {"model_kind", "view"} or weighted.get("view") != view or weighted.get("model_kind") not in {"ridge", "gpr"}:
                raise ValueError("sensitivity lacks complete weighted OOF values")
            fold_count = 9 if view == "study" else 17
            if any(type(weighted[key]) is not list or len(weighted[key]) != fold_count for key in ("fold_hyperparameters", "fold_sample_weights", "fold_reported_variances")):
                raise ValueError("sensitivity lacks frozen fold detail")
            for params in weighted["fold_hyperparameters"]:
                _validate_parameters(model_kind, params)
            _numeric_vector(weighted["prediction"], 17, "weighted prediction")
            _numeric_vector(weighted["baseline_prediction"], 17, "weighted baseline prediction")
            if not isinstance(weighted["fold_sample_weights"], (list, tuple)) or not isinstance(weighted["fold_reported_variances"], (list, tuple)):
                raise ValueError("sensitivity contains non-finite predictions")
        perturbation = model["perturbation"]
        if set(perturbation) != {"draws", "seed", "metric_quantiles", "gate_pass_count", "gate_pass_rate", "fold_hyperparameters"} or type(perturbation.get("draws")) is not int or perturbation["draws"] != 1000 or type(perturbation.get("seed")) is not int or perturbation.get("seed") != 20260916 or not isinstance(perturbation.get("metric_quantiles"), Mapping) or set(perturbation["metric_quantiles"]) != {"study", "variant"} or type(perturbation.get("gate_pass_count")) is not int or type(perturbation.get("gate_pass_rate")) not in (int, float) or not np.isfinite(perturbation["gate_pass_rate"]) or not 0 <= perturbation["gate_pass_rate"] <= 1:
            raise ValueError("sensitivity lacks the required 1,000-draw perturbation summary")
        for view_values in perturbation["metric_quantiles"].values():
            if not isinstance(view_values, Mapping) or set(view_values) != {"mae", "rmse", "baseline_mae", "baseline_rmse"} or any(not isinstance(summary, Mapping) or set(summary) != {"q025", "q500", "q975"} or any(type(value) not in (int, float) or not np.isfinite(value) or value < 0 for value in summary.values()) or not (summary["q025"] <= summary["q500"] <= summary["q975"]) for summary in view_values.values()):
                raise ValueError("sensitivity lacks complete perturbation quantiles")
        if perturbation["gate_pass_count"] < 0 or perturbation["gate_pass_count"] > perturbation["draws"] or not np.isclose(perturbation["gate_pass_rate"], perturbation["gate_pass_count"] / perturbation["draws"], rtol=0, atol=1e-12):
            raise ValueError("sensitivity perturbation gate summary is inconsistent")
        fold_params = perturbation["fold_hyperparameters"]
        if not isinstance(fold_params, Mapping) or set(fold_params) != {"study", "variant"}:
            raise ValueError("perturbation lacks both fold parameter sequences")
        for view, count in (("study", 9), ("variant", 17)):
            if type(fold_params[view]) is not list or len(fold_params[view]) != count:
                raise ValueError("perturbation lacks complete fold parameter sequences")
            for params in fold_params[view]:
                _validate_parameters(model_kind, params)


def _numeric_vector(value: Any, length: int, name: str, *, positive: bool = False) -> np.ndarray:
    if not isinstance(value, (list, tuple)) or len(value) != length or any(type(item) not in (int, float) for item in value):
        raise ValueError(f"{name} must be a strict one-dimensional numeric sequence")
    array = np.asarray(value, dtype=float)
    if array.ndim != 1 or not np.isfinite(array).all() or (positive and np.any(array <= 0)):
        raise ValueError(f"{name} must be a strict one-dimensional numeric sequence")
    return array


def _expected_outer_indices(model_table: pd.DataFrame, view: str, fold: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if view == "variant":
        test = (fold,)
    elif view == "study":
        test = tuple(index for index, source in enumerate(model_table["source_id"]) if source == FROZEN_SOURCE_STUDIES[fold])
    else:
        raise ValueError("unknown outer view")
    train = tuple(index for index in range(len(model_table)) if index not in test)
    return train, test


def reconstruct_inner_fold_audit(model_table: pd.DataFrame, tuning_records: pd.DataFrame) -> dict[str, Any]:
    """Reconstruct deterministic nested splits/scalers without fitting or reading predictions."""
    _validate_model_table(model_table)
    _validate_tuning(tuning_records)
    records: list[dict[str, Any]] = []
    features = model_table.loc[:, ACTIVE_FEATURE_COLUMNS].to_numpy(float)
    sources = model_table["source_id"].to_numpy(object)
    for row in tuning_records.itertuples(index=False):
        view, model, outer_fold = str(row.view), str(row.model_kind), int(row.outer_fold)
        outer_train, outer_test = _expected_outer_indices(model_table, view, outer_fold)
        if _indices(row.train_indices, "train_indices") != outer_train or _indices(row.test_indices, "test_indices") != outer_test:
            raise ValueError("cannot reconstruct from noncanonical outer folds")
        outer = np.asarray(outer_train, dtype=int)
        for inner_fold, split in enumerate(inner_group_splits(sources[outer])):
            inner_test = tuple(int(index) for index in outer[split.test])
            inner_train = tuple(int(index) for index in outer[split.train])
            held_source = sources[inner_test[0]]
            values = features[list(inner_train)]
            records.append({
                "model_kind": model, "view": view, "outer_fold": outer_fold, "inner_fold": inner_fold,
                "held_out_source_id": str(held_source), "train_indices": list(inner_train), "test_indices": list(inner_test),
                "scaler_mean": [float(value) for value in values.mean(axis=0)],
                "scaler_scale": [float(value) for value in values.std(axis=0, ddof=0)],
            })
    return {
        "schema_version": 1,
        "provenance": {
            "method": "post_execution_reconstruction",
            "original_execution_revision": _FROZEN_EXECUTION_REVISION,
            "diagnostic_rerun": False,
        },
        "records": records,
    }


def _validate_inner_fold_audit(audit: Mapping[str, Any], model_table: pd.DataFrame, tuning: pd.DataFrame, manifest: Mapping[str, Any]) -> None:
    if not isinstance(audit, Mapping) or set(audit) != {"schema_version", "provenance", "records"} or type(audit["schema_version"]) is not int or audit["schema_version"] != 1:
        raise ValueError("inner-fold audit has an invalid schema")
    provenance = audit["provenance"]
    provenance_fields = {"method", "original_execution_revision", "diagnostic_rerun"}
    if not isinstance(provenance, Mapping) or set(provenance) != provenance_fields or provenance["diagnostic_rerun"] is not False:
        raise ValueError("inner-fold audit lacks exact provenance")
    if type(provenance["method"]) is not str or type(provenance["original_execution_revision"]) is not str:
        raise ValueError("inner-fold provenance method and revision must be strings")
    expected = reconstruct_inner_fold_audit(model_table, tuning)
    repair = manifest.get("post_execution_provenance_repair")
    if repair is not None or manifest["git_revision"] == _FROZEN_EXECUTION_REVISION:
        if not isinstance(repair, Mapping) or not isinstance(repair.get("inner_fold_audit"), Mapping):
            raise ValueError("frozen inner-fold reconstruction lacks manifest provenance")
        declared = repair["inner_fold_audit"]
        if set(declared) != provenance_fields | {"scientific_result_files_changed"} or declared["scientific_result_files_changed"] is not False:
            raise ValueError("frozen inner-fold reconstruction changed scientific files")
        if provenance != expected["provenance"] or manifest["git_revision"] != _FROZEN_EXECUTION_REVISION or repair.get("execution_git_head") != _FROZEN_EXECUTION_REVISION:
            raise ValueError("inner-fold reconstruction does not retain the original execution")
    else:
        declared = manifest.get("inner_fold_audit")
        if not isinstance(declared, Mapping) or set(declared) != provenance_fields:
            raise ValueError("inner-fold execution lacks exact manifest provenance")
        if provenance["method"] != "recorded_during_execution" or provenance["original_execution_revision"] != manifest["git_revision"]:
            raise ValueError("inner-fold execution provenance does not match run")
    if declared["diagnostic_rerun"] is not False or any(declared[key] != provenance[key] for key in provenance_fields):
        raise ValueError("inner-fold audit and manifest provenance disagree")
    if type(audit["records"]) is not list:
        raise ValueError("inner-fold records must be an actual record sequence")
    fields = {"model_kind", "view", "outer_fold", "inner_fold", "held_out_source_id", "train_indices", "test_indices", "scaler_mean", "scaler_scale"}
    for record in audit["records"]:
        if not isinstance(record, Mapping) or set(record) != fields:
            raise ValueError("inner-fold record has an invalid field schema")
        if any(type(record[key]) is not str for key in ("model_kind", "view", "held_out_source_id")) or any(type(record[key]) is not int for key in ("outer_fold", "inner_fold")):
            raise ValueError("inner-fold identities must have exact string/integer types")
        for key in ("train_indices", "test_indices"):
            if type(record[key]) is not list:
                raise ValueError("inner-fold indices must be actual integer sequences")
            _indices(record[key], key)
        _numeric_vector(record["scaler_mean"], 3, "inner scaler mean")
        _numeric_vector(record["scaler_scale"], 3, "inner scaler scale", positive=True)
    if audit["records"] != expected["records"]:
        raise ValueError("inner-fold audit differs from exact training-only reconstruction")


def _object(value: Any, name: str) -> Mapping[str, Any]:
    if isinstance(value, Mapping):
        return value
    if isinstance(value, str):
        try:
            def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
                result: dict[str, Any] = {}
                for key, item in pairs:
                    if key in result:
                        raise ValueError(f"{name} contains duplicate key {key!r}")
                    result[key] = item
                return result
            parsed = json.loads(value, object_pairs_hook=unique_object)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{name} is not valid JSON") from exc
        if isinstance(parsed, Mapping):
            return parsed
    raise ValueError(f"{name} must be a non-empty mapping")


def _indices(value: Any, name: str) -> tuple[int, ...]:
    if isinstance(value, str):
        value = json.loads(value)
    if not isinstance(value, (list, tuple)) or not value or any(type(item) is not int for item in value) or len(set(value)) != len(value):
        raise ValueError(f"{name} must be non-empty integer indices")
    return tuple(int(item) for item in value)


def _candidate_grid(model: str) -> dict[str, dict[str, float]]:
    if model == "ridge":
        return {f"alpha={alpha:g}": {"alpha": alpha} for alpha in (0.01, 0.1, 1.0, 10.0, 100.0)}
    return {
        f"length_scale={length:g},noise_variance={noise:g}": {"length_scale": length, "noise_variance": noise}
        for length in (0.5, 1.0, 2.0) for noise in (0.0025, 0.01, 0.04)
    }


def _validate_parameters(model: str, params: Any) -> None:
    grid = _candidate_grid(model)
    if not isinstance(params, Mapping) or set(params) != set(next(iter(grid.values()))):
        raise ValueError("hyperparameters do not have the exact frozen field schema")
    if any(type(value) not in (int, float) or not np.isfinite(value) for value in params.values()):
        raise ValueError("hyperparameters must have finite builtin numeric values")
    if dict(params) not in grid.values():
        raise ValueError("hyperparameters are not exactly one frozen candidate")


def _validate_tuning_choice(model: str, params: Mapping[str, Any], scores: Mapping[str, Any]) -> None:
    grid = _candidate_grid(model)
    if set(scores) != set(grid) or any(type(value) not in (int, float) or not np.isfinite(value) or value < 0 for value in scores.values()):
        raise ValueError("inner scores are not a one-to-one complete frozen grid")
    _validate_parameters(model, params)
    best_score = min(scores.values())
    tied = [grid[name] for name in grid if scores[name] == best_score]
    selected = max(tied, key=lambda candidate: (candidate.get("alpha", -np.inf), candidate.get("length_scale", -np.inf), candidate.get("noise_variance", -np.inf)))
    if dict(params) != selected:
        raise ValueError("selected parameters are not the frozen deterministic minimizer")


def _validate_cross_artifacts(result: Mapping[str, Any]) -> None:
    model_table = result["model_table"].set_index("variant")
    if not model_table.index.is_unique:
        raise ValueError("model table variant identities must be unique")
    metrics = result["metrics"]
    recomputed_gates = {}
    for view in ("study", "variant"):
        oof = result[f"{view}_oof"]
        for model_kind in ("ridge", "gpr"):
            model_rows = oof[oof.model_kind == model_kind]
            for fold in range(9 if view == "study" else 17):
                _, expected_test = _expected_outer_indices(result["model_table"], view, fold)
                fold_rows = model_rows[model_rows.outer_fold == fold]
                expected_variants = tuple(result["model_table"].iloc[list(expected_test)]["variant"])
                expected_sources = tuple(result["model_table"].iloc[list(expected_test)]["source_id"])
                expected_truth = result["model_table"].iloc[list(expected_test)]["y"].to_numpy(float)
                if tuple(fold_rows["variant"]) != expected_variants or tuple(fold_rows["source_id"]) != expected_sources or not np.array_equal(fold_rows["truth"].to_numpy(float), expected_truth):
                    raise ValueError("OOF variant, source, or truth does not match the exact canonical outer fold")
            rows = model_rows.set_index("variant")
            if not rows.index.is_unique:
                raise ValueError("OOF variant identities must be unique per model")
            rows = rows.loc[list(FROZEN_PRIMARY_VARIANTS)]
            if not np.array_equal(rows.index.to_numpy(), model_table.index.to_numpy()) or not np.allclose(rows["truth"].to_numpy(float), model_table["y"].to_numpy(float), rtol=0, atol=1e-12) or tuple(rows["source_id"]) != tuple(model_table["source_id"]):
                raise ValueError("OOF truth, source, or variant identities drift from the model table")
            scorer = study_macro_metrics if view == "study" else variant_metrics
            model_score = scorer(rows["truth"], rows["model_prediction"], rows["source_id"]) if view == "study" else scorer(rows["truth"], rows["model_prediction"])
            baseline_score = scorer(rows["truth"], rows["baseline_prediction"], rows["source_id"]) if view == "study" else scorer(rows["truth"], rows["baseline_prediction"])
            decision = metrics["models"][model_kind]
            for key, observed, expected in (("model_views", decision["model_views"][view], model_score), ("baseline_views", decision["baseline_views"][view], baseline_score)):
                if not np.isclose(float(observed["mae"]), expected.mae, rtol=0, atol=1e-12) or not np.isclose(float(observed["rmse"]), expected.rmse, rtol=0, atol=1e-12):
                    raise ValueError("metrics do not match OOF predictions")

    for model_kind in ("ridge", "gpr"):
        study_rows = result["study_oof"][result["study_oof"].model_kind == model_kind].set_index("variant").loc[list(FROZEN_PRIMARY_VARIANTS)]
        decision = metrics["models"][model_kind]
        recomputed = evaluate_gate(
            {view: Metrics(**decision["model_views"][view]) for view in ("study", "variant")},
            {view: Metrics(**decision["baseline_views"][view]) for view in ("study", "variant")},
            single_study_influence(study_rows["truth"], study_rows["model_prediction"], study_rows["baseline_prediction"], study_rows["source_id"]),
        )
        if recomputed.pass_gate != decision["pass_gate"] or tuple(decision["failed_rules"]) != recomputed.failed_rules or any(not np.isclose(float(decision[key][view]), expected[view], rtol=0, atol=1e-12) for key, expected in (("mae_improvement", recomputed.mae_improvement), ("rmse_degradation", recomputed.rmse_degradation)) for view in ("study", "variant")):
            raise ValueError("metrics gate decision does not match OOF predictions")
        recomputed_gates[model_kind] = recomputed
    if metrics["winner"] != select_winner(recomputed_gates):
        raise ValueError("metrics winner does not match gate decisions")
    for _, row in result["tuning_records"].iterrows():
        params = _object(row["selected_params"], "selected_params")
        _validate_parameters(row["model_kind"], params)
        scores = _object(row["inner_scores"], "inner_scores")
        train, test = _indices(row["train_indices"], "train_indices"), _indices(row["test_indices"], "test_indices")
        mean, scale = row["scaler_mean"], row["scaler_scale"]
        if isinstance(mean, str): mean = json.loads(mean)
        if isinstance(scale, str): scale = json.loads(scale)
        if not params or not scores or len(mean) != 3 or len(scale) != 3 or not np.isfinite(np.asarray([*mean, *scale], dtype=float)).all() or np.any(np.asarray(scale, dtype=float) <= 0) or set(train).intersection(test):
            raise ValueError("tuning records lack complete valid fold audit fields")
        view, model, fold = str(row["view"]), str(row["model_kind"]), int(row["outer_fold"])
        records = result[f"{view}_oof"]
        records = records[(records.model_kind == model) & (records.outer_fold == fold)]
        expected_train, expected_test = _expected_outer_indices(result["model_table"], view, fold)
        if tuple(test) != expected_test or tuple(train) != expected_train:
            raise ValueError("tuning records do not match OOF fold assignments")
        expected_params = _object(records.iloc[0]["selected_inner_hyperparameters"], "OOF selected_inner_hyperparameters")
        if any(_object(value, "OOF selected_inner_hyperparameters") != expected_params for value in records["selected_inner_hyperparameters"]) or params != expected_params:
            raise ValueError("tuning selected parameters do not match OOF records")
        if any(tuple(np.asarray(json.loads(value) if isinstance(value, str) else value, dtype=float)) != tuple(np.asarray(mean, dtype=float)) for value in records["scaler_mean"]) or any(tuple(np.asarray(json.loads(value) if isinstance(value, str) else value, dtype=float)) != tuple(np.asarray(scale, dtype=float)) for value in records["scaler_scale"]):
            raise ValueError("tuning scaler audit does not match OOF records")
        expected_mean = model_table.iloc[list(train)][list(ACTIVE_FEATURE_COLUMNS)].to_numpy(float).mean(axis=0)
        expected_scale = model_table.iloc[list(train)][list(ACTIVE_FEATURE_COLUMNS)].to_numpy(float).std(axis=0, ddof=0)
        if not np.allclose(mean, expected_mean, rtol=0, atol=1e-12) or not np.allclose(scale, expected_scale, rtol=0, atol=1e-12) or not np.isclose(float(row["baseline_mean"]), model_table.iloc[list(train)]["y"].mean(), rtol=0, atol=1e-12):
            raise ValueError("tuning scaler or baseline is not the training-only model-table value")
        if not np.allclose(records["baseline_prediction"].to_numpy(float), float(row["baseline_mean"]), rtol=0, atol=1e-12):
            raise ValueError("primary OOF baseline does not match outer-train mean")
        _validate_tuning_choice(model, params, scores)
    for model, values in result["sensitivity"].items():
        for view, weighted in values["weighted"].items():
            if weighted["model_kind"] != model or weighted["view"] != view:
                raise ValueError("weighted sensitivity provenance does not match enclosing model/view")
            source = model_table["source_id"]
            truth = model_table["y"]
            scorer = study_macro_metrics if view == "study" else variant_metrics
            predicted = scorer(truth, weighted["prediction"], source) if view == "study" else scorer(truth, weighted["prediction"])
            baseline = scorer(truth, weighted["baseline_prediction"], source) if view == "study" else scorer(truth, weighted["baseline_prediction"])
            if not np.isclose(float(weighted["metrics"]["mae"]), predicted.mae, rtol=0, atol=1e-12) or not np.isclose(float(weighted["metrics"]["rmse"]), predicted.rmse, rtol=0, atol=1e-12) or not np.isclose(float(weighted["baseline_metrics"]["mae"]), baseline.mae, rtol=0, atol=1e-12) or not np.isclose(float(weighted["baseline_metrics"]["rmse"]), baseline.rmse, rtol=0, atol=1e-12):
                raise ValueError("weighted sensitivity metrics do not match predictions")
            folds = result["tuning_records"][(result["tuning_records"].model_kind == model) & (result["tuning_records"].view == view)].sort_values("outer_fold")
            if list(weighted["fold_hyperparameters"]) != [_object(row.selected_params, "selected_params") for row in folds.itertuples()] or len(folds) != len(weighted["fold_sample_weights"]):
                raise ValueError("weighted sensitivity does not reuse primary fold parameters")
            expected_baseline = np.empty(17, dtype=float)
            for fold in folds.itertuples():
                expected_baseline[list(_indices(fold.test_indices, "test_indices"))] = float(fold.baseline_mean)
            if not np.allclose(_numeric_vector(weighted["baseline_prediction"], 17, "weighted baseline"), expected_baseline, rtol=0, atol=1e-12):
                raise ValueError("weighted baseline does not match exact outer-training means")
            for record, fold in zip(weighted["fold_sample_weights"], folds.itertuples()):
                train = _indices(fold.train_indices, "train_indices")
                errors = model_table.iloc[list(train)]["y_error"].to_numpy(float)
                if model == "ridge":
                    expected = 1 / np.square(errors); expected /= expected.mean()
                    observed = _numeric_vector(record, len(train), "Ridge fold weights", positive=True)
                    if not np.allclose(observed, expected, rtol=0, atol=1e-12): raise ValueError("Ridge sensitivity weights differ from frozen uncertainty")
                elif not isinstance(record, (list, tuple)) or len(record) != 0: raise ValueError("GPR sensitivity must carry an actual empty weight sequence")
            for record, fold in zip(weighted["fold_reported_variances"], folds.itertuples()):
                errors = model_table.iloc[list(_indices(fold.train_indices, "train_indices"))]["y_error"].to_numpy(float)
                if model == "gpr":
                    observed = _numeric_vector(record, len(errors), "GPR fold variances")
                    if np.any(observed < 0) or not np.allclose(observed, np.square(errors), rtol=0, atol=1e-12): raise ValueError("GPR sensitivity variances differ from frozen uncertainty")
                elif not isinstance(record, (list, tuple)) or len(record) != 0: raise ValueError("Ridge sensitivity must carry an actual empty variance sequence")
        perturbation = values["perturbation"]
        expected_hyperparameters = {
            view: [_object(row.selected_params, "selected_params") for row in result["tuning_records"][(result["tuning_records"].model_kind == model) & (result["tuning_records"].view == view)].sort_values("outer_fold").itertuples()]
            for view in ("study", "variant")
        }
        if perturbation["fold_hyperparameters"] != expected_hyperparameters:
            raise ValueError("perturbation sensitivity does not reuse both primary fold parameter sets")


def _write_bytes(path: Path, content: bytes) -> None:
    with path.open("wb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())


def write_diagnostic_artifacts(output: str | Path, result: Mapping[str, Any]) -> None:
    """Write a complete deterministic report and create ``SUCCESS`` strictly last.

    Callers must provide a new or empty output directory.  Files are assembled
    in a sibling staging directory; a failed write cleans both staging files and
    every newly created output file, leaving no partial diagnostic to mistake
    for a real run.
    """
    if isinstance(result, Mapping) and "sensitivity" in result:
        result = {**result, "sensitivity": _normalize_sensitivity_containers(result["sensitivity"])}
    _validate_result(result)
    target = Path(output)
    # Recheck immediately before writing: the CLI's early preflight cannot
    # protect against another process populating the destination while fitting.
    preflight_diagnostic_output(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{target.name}.", dir=target.parent))
    created_target = False
    copied: list[Path] = []
    try:
        content = {
            "model_table": _csv_bytes(result["model_table"]),
            "tuning_records": _csv_bytes(result["tuning_records"]),
            "study_oof": _csv_bytes(result["study_oof"]),
            "variant_oof": _csv_bytes(result["variant_oof"]),
            "inner_fold_audit": _json_bytes(result["inner_fold_audit"]),
            "metrics": _json_bytes(result["metrics"]),
            "sensitivity": _json_bytes(result["sensitivity"]),
            "audit_markdown": result["audit_markdown"].encode("utf-8"),
            "run_log": result["run_log"].encode("utf-8"),
        }
        manifest = dict(result["run_manifest"])
        if "artifact_sha256" in manifest:
            manifest["artifact_sha256"] = {
                _ARTIFACT_FILENAMES[name]: _sha256(content[name]) for name in sorted(content)
            }
        content["run_manifest"] = _json_bytes(manifest)
        for name, filename in _ARTIFACT_FILENAMES.items():
            _write_bytes(staging / filename, content[name])
        if not target.exists():
            target.mkdir()
            created_target = True
        for filename in _ARTIFACT_FILENAMES.values():
            destination = target / filename
            os.replace(staging / filename, destination)
            copied.append(destination)
        # This marker is deliberately the only post-completion write.
        _write_bytes(staging / "SUCCESS", b"completed\n")
        os.replace(staging / "SUCCESS", target / "SUCCESS")
    except BaseException:
        for path in copied:
            path.unlink(missing_ok=True)
        (target / "SUCCESS").unlink(missing_ok=True)
        if created_target:
            target.rmdir()
        raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)

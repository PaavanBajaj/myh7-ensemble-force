"""Secondary reported-uncertainty checks for a frozen primary diagnostic.

These routines deliberately refit only the primary model kind and its recorded
outer-fold hyperparameters.  They are descriptive: no candidate search, model
selection, profile selection, or primary gate decision is performed here.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.linear_model import Ridge

from kcat_rel.registry import ACTIVE_FEATURE_COLUMNS, FROZEN_PRIMARY_VARIANTS
from kcat_rel.validation import (
    FoldResult,
    Metrics,
    OofResult,
    _input_digest,
    evaluate_gate,
    single_study_influence,
    study_macro_metrics,
    variant_metrics,
)


@dataclass(frozen=True)
class SensitivityResult:
    """OOF result of a secondary, reported-uncertainty-weighted refit."""

    view: str
    model_kind: str
    prediction: np.ndarray
    baseline_prediction: np.ndarray
    metrics: Metrics
    baseline_metrics: Metrics
    fold_hyperparameters: tuple[dict[str, float], ...]
    fold_sample_weights: tuple[tuple[float, ...], ...]
    fold_reported_variances: tuple[tuple[float, ...], ...]


@dataclass(frozen=True)
class PerturbationSummary:
    """Fixed-seed quantile and primary-gate pass summaries across label draws."""

    draws: int
    seed: int
    metric_quantiles: Mapping[str, Mapping[str, Mapping[str, float]]]
    gate_pass_count: int | None
    gate_pass_rate: float | None
    fold_hyperparameters: Any


def _as_features(X: Any, primary: OofResult) -> np.ndarray:
    if not isinstance(X, pd.DataFrame) or tuple(X.columns) != ("variant", *ACTIVE_FEATURE_COLUMNS):
        raise ValueError("X feature schema must be exactly variant plus the active frozen features")
    try:
        features = X.loc[:, ACTIVE_FEATURE_COLUMNS].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("X predictor columns must be numeric") from exc
    if features.ndim != 2 or not np.isfinite(features).all():
        raise ValueError("X predictor columns must be finite")
    variants = tuple(X["variant"])
    if variants != FROZEN_PRIMARY_VARIANTS or variants != _variant_order(primary):
        raise ValueError("X variants do not match the frozen primary cohort snapshot")
    return features


def _as_target(y: Any, primary: OofResult, n_rows: int) -> np.ndarray:
    values = np.asarray(y, dtype=float).reshape(-1)
    if len(values) != n_rows or not np.isfinite(values).all():
        raise ValueError("y must be finite and match X rows")
    if not np.array_equal(values, _truth_order(primary)):
        raise ValueError("y must exactly match frozen primary OOF truth in variant order")
    return values


def _as_error(y_error: Any, n_rows: int) -> np.ndarray:
    values = np.asarray(y_error, dtype=float).reshape(-1)
    if len(values) != n_rows or not np.isfinite(values).all() or np.any(values < 0.0):
        raise ValueError("y_error must be finite, non-negative, and match X rows")
    return values


def _as_groups(groups: Any, primary: OofResult, n_rows: int) -> np.ndarray:
    values = np.asarray(groups, dtype=object).reshape(-1)
    if len(values) != n_rows or any(not isinstance(value, str) or not value.strip() for value in values):
        raise ValueError("groups must be non-empty strings and match X rows")
    expected = _source_order(primary)
    if tuple(values) != expected:
        raise ValueError("groups do not match the frozen primary folds")
    return values


def _require_primary_input_digest(
    primary: OofResult, features: np.ndarray, target: np.ndarray, source: np.ndarray
) -> str:
    digest = _input_digest(
        FROZEN_PRIMARY_VARIANTS,
        source,
        target,
        ACTIVE_FEATURE_COLUMNS,
        features,
    )
    if digest != primary.input_digest:
        raise ValueError("supplied inputs do not match the frozen primary input digest")
    return digest


def _records_for_fold(primary: OofResult, fold: FoldResult) -> pd.DataFrame:
    required = {"outer_fold", "variant", "source_id", "truth"}
    if not isinstance(primary.records, pd.DataFrame) or not required.issubset(primary.records.columns):
        raise ValueError("primary_result lacks complete OOF audit records")
    records = primary.records.loc[primary.records["outer_fold"] == fold.outer_fold]
    if len(records) != len(fold.test_indices):
        raise ValueError("primary_result OOF records do not match frozen folds")
    indices = np.asarray(fold.test_indices, dtype=int)
    try:
        record_truth = records["truth"].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("primary_result has invalid OOF truth") from exc
    if (
        tuple(records["variant"]) != tuple(primary.variants[index] for index in indices)
        or tuple(records["source_id"]) != tuple(primary.source_ids[index] for index in indices)
        or not np.array_equal(record_truth, np.asarray(primary.truth, dtype=float)[indices])
    ):
        raise ValueError("primary_result OOF records differ from the frozen primary snapshot")
    return records


def _validate_folds(primary: OofResult, n_rows: int) -> None:
    if primary.view not in {"study", "variant"} or primary.model_kind not in {"ridge", "gpr"}:
        raise ValueError("primary_result must be a primary Ridge or GPR outer view")
    if not primary.folds:
        raise ValueError("primary_result must retain frozen outer folds")
    if (
        tuple(primary.variants) != FROZEN_PRIMARY_VARIANTS
        or len(primary.source_ids) != n_rows
        or len(primary.truth) != n_rows
        or not isinstance(primary.input_digest, str)
        or len(primary.input_digest) != 64
    ):
        raise ValueError("primary_result has an invalid frozen primary input snapshot")
    if any(not isinstance(value, str) or not value for value in primary.source_ids):
        raise ValueError("primary_result has invalid frozen primary source identities")
    if not np.isfinite(np.asarray(primary.truth, dtype=float)).all():
        raise ValueError("primary_result has invalid frozen primary truth")
    seen: list[int] = []
    for fold in primary.folds:
        train = np.asarray(fold.train_indices, dtype=int)
        test = np.asarray(fold.test_indices, dtype=int)
        if not len(train) or not len(test) or np.intersect1d(train, test).size:
            raise ValueError("primary_result has invalid frozen fold assignments")
        if np.any(train < 0) or np.any(test < 0) or np.any(train >= n_rows) or np.any(test >= n_rows):
            raise ValueError("primary_result has out-of-range frozen fold assignments")
        if len(fold.scaler_mean) == 0 or len(fold.scaler_mean) != len(fold.scaler_scale):
            raise ValueError("primary_result has invalid frozen scaler state")
        scale = np.asarray(fold.scaler_scale, dtype=float)
        mean = np.asarray(fold.scaler_mean, dtype=float)
        if not np.isfinite(mean).all() or not np.isfinite(scale).all() or np.any(scale <= 0.0):
            raise ValueError("primary_result has invalid frozen scaler state")
        _records_for_fold(primary, fold)
        seen.extend(test.tolist())
    if sorted(seen) != list(range(n_rows)):
        raise ValueError("primary_result folds must provide one OOF row per input row")


def _variant_order(primary: OofResult) -> tuple[str, ...]:
    n_rows = len(primary.prediction)
    _validate_folds(primary, n_rows)
    return primary.variants


def _source_order(primary: OofResult) -> tuple[str, ...]:
    n_rows = len(primary.prediction)
    _validate_folds(primary, n_rows)
    return primary.source_ids


def _truth_order(primary: OofResult) -> np.ndarray:
    n_rows = len(primary.prediction)
    _validate_folds(primary, n_rows)
    return np.asarray(primary.truth, dtype=float)


def _kernel(params: Mapping[str, float]):
    try:
        length_scale = float(params["length_scale"])
        noise_variance = float(params["noise_variance"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("primary GPR fold has invalid frozen hyperparameters") from exc
    if not np.isfinite(length_scale) or length_scale <= 0.0 or not np.isfinite(noise_variance) or noise_variance < 0.0:
        raise ValueError("primary GPR fold has invalid frozen hyperparameters")
    return (
        ConstantKernel(1.0, constant_value_bounds="fixed")
        * Matern(length_scale=length_scale, length_scale_bounds="fixed", nu=1.5)
        + WhiteKernel(noise_level=noise_variance, noise_level_bounds="fixed")
    )


def _fit_frozen_fold(
    primary: OofResult,
    fold: FoldResult,
    features: np.ndarray,
    target: np.ndarray,
    *,
    weights: np.ndarray | None = None,
    reported_variance: np.ndarray | None = None,
) -> np.ndarray:
    train = np.asarray(fold.train_indices, dtype=int)
    test = np.asarray(fold.test_indices, dtype=int)
    mean = np.asarray(fold.scaler_mean, dtype=float)
    scale = np.asarray(fold.scaler_scale, dtype=float)
    if features.shape[1] != len(mean):
        raise ValueError("X feature width does not match the frozen primary scaler")
    train_X = (features[train] - mean) / scale
    test_X = (features[test] - mean) / scale
    if primary.model_kind == "ridge":
        try:
            alpha = float(fold.selected_params["alpha"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("primary Ridge fold has invalid frozen hyperparameters") from exc
        if not np.isfinite(alpha) or alpha < 0.0:
            raise ValueError("primary Ridge fold has invalid frozen hyperparameters")
        estimator = Ridge(alpha=alpha, fit_intercept=True)
        estimator.fit(train_X, target[train], sample_weight=weights)
    else:
        estimator = GaussianProcessRegressor(
            kernel=_kernel(fold.selected_params),
            alpha=0.0 if reported_variance is None else reported_variance,
            optimizer=None,
            normalize_y=False,
            random_state=20260916,
        )
        estimator.fit(train_X, target[train])
    prediction = np.asarray(estimator.predict(test_X), dtype=float)
    if not np.isfinite(prediction).all():
        raise RuntimeError("frozen sensitivity fold did not produce finite predictions")
    return prediction


def _metrics_for_view(view: str, target: np.ndarray, prediction: np.ndarray, groups: np.ndarray) -> Metrics:
    if view == "study":
        return study_macro_metrics(target, prediction, groups)
    return variant_metrics(target, prediction)


def _refit(
    primary: OofResult,
    features: np.ndarray,
    target: np.ndarray,
    errors: np.ndarray | None,
    *,
    weighted: bool,
) -> tuple[np.ndarray, np.ndarray, tuple[tuple[float, ...], ...], tuple[tuple[float, ...], ...]]:
    prediction = np.empty(len(target), dtype=float)
    baseline = np.empty(len(target), dtype=float)
    all_weights: list[tuple[float, ...]] = []
    all_variances: list[tuple[float, ...]] = []
    for fold in primary.folds:
        train = np.asarray(fold.train_indices, dtype=int)
        test = np.asarray(fold.test_indices, dtype=int)
        weights: np.ndarray | None = None
        variance: np.ndarray | None = None
        if weighted and primary.model_kind == "ridge":
            assert errors is not None
            train_errors = errors[train]
            if np.any(train_errors <= 0.0):
                raise ValueError("positive y_error is required for Ridge inverse-variance weighting")
            weights = 1.0 / np.square(train_errors)
            weights /= float(np.mean(weights))
            all_weights.append(tuple(float(value) for value in weights))
            all_variances.append(())
        elif weighted and primary.model_kind == "gpr":
            assert errors is not None
            variance = np.square(errors[train])
            all_weights.append(())
            all_variances.append(tuple(float(value) for value in variance))
        else:
            all_weights.append(())
            all_variances.append(())
        prediction[test] = _fit_frozen_fold(
            primary, fold, features, target, weights=weights, reported_variance=variance
        )
        baseline[test] = float(np.mean(target[train]))
    if not np.isfinite(prediction).all() or not np.isfinite(baseline).all():
        raise RuntimeError("frozen sensitivity did not produce complete finite OOF predictions")
    return prediction, baseline, tuple(all_weights), tuple(all_variances)


def run_weighted_sensitivity(
    primary_result: OofResult, X: Any, y: Any, y_error: Any, groups: Any
) -> SensitivityResult:
    """Refit one primary outer view with reported uncertainty as a secondary check.

    Ridge receives per-training-fold inverse-variance weights normalized to mean
    one.  GPR keeps the selected WhiteKernel noise and adds each training row's
    reported variance to its covariance diagonal.
    """
    if not isinstance(primary_result, OofResult):
        raise ValueError("weighted sensitivity requires one frozen OofResult")
    features = _as_features(X, primary_result)
    target = _as_target(y, primary_result, len(features))
    errors = _as_error(y_error, len(features))
    source = _as_groups(groups, primary_result, len(features))
    _require_primary_input_digest(primary_result, features, target, source)
    prediction, baseline, weights, variances = _refit(
        primary_result, features, target, errors, weighted=True
    )
    return SensitivityResult(
        view=primary_result.view,
        model_kind=primary_result.model_kind,
        prediction=prediction,
        baseline_prediction=baseline,
        metrics=_metrics_for_view(primary_result.view, target, prediction, source),
        baseline_metrics=_metrics_for_view(primary_result.view, target, baseline, source),
        fold_hyperparameters=primary_result.fold_hyperparameters,
        fold_sample_weights=weights,
        fold_reported_variances=variances,
    )


def _primary_views(primary_result: OofResult | Mapping[str, OofResult]) -> dict[str, OofResult]:
    if isinstance(primary_result, OofResult):
        return {primary_result.view: primary_result}
    if not isinstance(primary_result, Mapping) or set(primary_result) != {"study", "variant"}:
        raise ValueError("primary_result must be one OofResult or frozen study and variant OofResults")
    views = dict(primary_result)
    if not all(isinstance(result, OofResult) and result.view == view for view, result in views.items()):
        raise ValueError("primary_result views must match their frozen OofResults")
    if views["study"].model_kind != views["variant"].model_kind:
        raise ValueError("primary_result views must retain one selected model kind")
    return views


def _quantiles(values: list[float]) -> Mapping[str, float]:
    return {
        "q025": float(np.quantile(values, 0.025)),
        "q500": float(np.quantile(values, 0.5)),
        "q975": float(np.quantile(values, 0.975)),
    }


def run_perturbation_sensitivity(
    primary_result: OofResult | Mapping[str, OofResult],
    X: Any,
    y: Any,
    y_error: Any,
    groups: Any,
    draws: int = 1000,
    seed: int = 20260916,
) -> PerturbationSummary:
    """Summarize fixed-seed, unweighted label perturbations with frozen fits.

    A two-view primary mapping additionally yields the fraction of perturbations
    satisfying the unchanged primary gate.  A single view cannot define that
    gate and therefore reports ``None`` for its pass count and rate.
    """
    if isinstance(draws, bool) or not isinstance(draws, (int, np.integer)) or draws <= 0:
        raise ValueError("draws must be a positive integer")
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)):
        raise ValueError("seed must be an integer")
    views = _primary_views(primary_result)
    first = next(iter(views.values()))
    features = _as_features(X, first)
    target = _as_target(y, first, len(features))
    errors = _as_error(y_error, len(features))
    source = _as_groups(groups, first, len(features))
    supplied_digest = _require_primary_input_digest(first, features, target, source)
    for view, result in views.items():
        if (
            _variant_order(result) != _variant_order(first)
            or _source_order(result) != _source_order(first)
            or not np.array_equal(_truth_order(result), _truth_order(first))
            or not np.array_equal(target, _truth_order(result))
            or result.input_digest != first.input_digest
            or result.input_digest != supplied_digest
        ):
            raise ValueError(f"{view} primary snapshot does not match the frozen canonical input")

    values = {
        view: {"mae": [], "rmse": [], "baseline_mae": [], "baseline_rmse": []}
        for view in views
    }
    passed = 0
    rng = np.random.default_rng(int(seed))
    for _ in range(int(draws)):
        perturbed = target + rng.normal(loc=0.0, scale=errors, size=len(target))
        model_metrics: dict[str, Metrics] = {}
        baseline_metrics: dict[str, Metrics] = {}
        predictions: dict[str, np.ndarray] = {}
        baselines: dict[str, np.ndarray] = {}
        for view, result in views.items():
            prediction, baseline, _, _ = _refit(result, features, perturbed, None, weighted=False)
            predictions[view] = prediction
            baselines[view] = baseline
            model_metrics[view] = _metrics_for_view(view, perturbed, prediction, source)
            baseline_metrics[view] = _metrics_for_view(view, perturbed, baseline, source)
            values[view]["mae"].append(model_metrics[view].mae)
            values[view]["rmse"].append(model_metrics[view].rmse)
            values[view]["baseline_mae"].append(baseline_metrics[view].mae)
            values[view]["baseline_rmse"].append(baseline_metrics[view].rmse)
        if set(views) == {"study", "variant"}:
            influence = single_study_influence(
                perturbed, predictions["study"], baselines["study"], source
            )
            if evaluate_gate(model_metrics, baseline_metrics, influence).pass_gate:
                passed += 1

    quantiles = {
        view: {name: _quantiles(metric_values) for name, metric_values in metrics.items()}
        for view, metrics in values.items()
    }
    fold_hyperparameters: Any
    if set(views) == {"study", "variant"}:
        fold_hyperparameters = {view: result.fold_hyperparameters for view, result in views.items()}
        gate_pass_count: int | None = passed
        gate_pass_rate: float | None = passed / int(draws)
    else:
        fold_hyperparameters = first.fold_hyperparameters
        gate_pass_count = None
        gate_pass_rate = None
    return PerturbationSummary(
        draws=int(draws),
        seed=int(seed),
        metric_quantiles=quantiles,
        gate_pass_count=gate_pass_count,
        gate_pass_rate=gate_pass_rate,
        fold_hyperparameters=fold_hyperparameters,
    )

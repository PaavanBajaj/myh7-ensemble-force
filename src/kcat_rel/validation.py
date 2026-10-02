"""Leakage-safe nested validation and model gates for the frozen diagnostic."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Any, Iterator, Mapping

import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from kcat_rel.registry import FROZEN_SOURCE_STUDIES


@dataclass(frozen=True)
class Metrics:
    """Prediction errors on the original, unscaled target scale."""

    mae: float
    rmse: float


@dataclass(frozen=True)
class IndexSplit:
    """One deterministic, disjoint train/test index split."""

    train: np.ndarray
    test: np.ndarray


@dataclass(frozen=True)
class InnerFoldResult:
    """Training-only predictor-scaler state for one nested inner split."""

    train_indices: np.ndarray
    test_indices: np.ndarray
    scaler_mean: tuple[float, ...]
    scaler_scale: tuple[float, ...]


@dataclass(frozen=True)
class FoldResult:
    """The training-only state and selected hyperparameters for one outer fold."""

    outer_fold: int
    train_indices: np.ndarray
    test_indices: np.ndarray
    scaler_mean: tuple[float, ...]
    scaler_scale: tuple[float, ...]
    baseline_mean: float
    selected_params: dict[str, float]
    inner_scores: dict[str, float]
    inner_folds: tuple[InnerFoldResult, ...]


@dataclass(frozen=True)
class OofResult:
    """Out-of-fold predictions plus complete outer-fold audit information."""

    view: str
    model_kind: str
    prediction: np.ndarray
    baseline_prediction: np.ndarray
    folds: tuple[FoldResult, ...]
    records: pd.DataFrame
    variants: tuple[str, ...]
    source_ids: tuple[str, ...]
    truth: tuple[float, ...]
    input_digest: str

    @property
    def oof(self) -> pd.DataFrame:
        """Alias for the OOF audit rows used by reporting consumers."""
        return self.records

    @property
    def fold_hyperparameters(self) -> tuple[dict[str, float], ...]:
        """Frozen outer-fold selections for downstream sensitivity analyses."""
        return tuple(dict(fold.selected_params) for fold in self.folds)


@dataclass(frozen=True)
class GateDecision:
    """Whether one model clears every predeclared forward-prediction gate."""

    pass_gate: bool
    failed_rules: tuple[str, ...]
    mae_improvement: dict[str, float]
    rmse_degradation: dict[str, float]
    model_views: Mapping[str, Metrics]
    baseline_views: Mapping[str, Metrics]


def _as_numeric_matrix(X: Any, registry: Any) -> tuple[np.ndarray, tuple[str, ...]]:
    """Accept only the exact frozen predictor schema and primary variant order."""
    if not isinstance(X, pd.DataFrame):
        raise ValueError("X must be a DataFrame carrying audited frozen variant identities")
    feature_columns = tuple(getattr(registry, "feature_columns", ()))
    columns = tuple(X.columns)
    if columns == ("variant", *feature_columns):
        variants = tuple(X["variant"])
        values = X.loc[:, feature_columns].to_numpy(dtype=float)
    elif columns == feature_columns:
        variants = tuple(X.index)
        values = X.to_numpy(dtype=float)
    else:
        raise ValueError("X feature schema must contain exactly variant plus the active frozen features")
    expected_variants = tuple(getattr(registry, "variants", ()))
    if (
        not expected_variants
        or any(not isinstance(variant, str) or not variant.strip() for variant in variants)
        or variants != expected_variants
    ):
        raise ValueError("X variants do not match the ordered frozen primary cohort")
    if values.ndim != 2 or values.shape[0] == 0 or values.shape[1] == 0 or not np.isfinite(values).all():
        raise ValueError("X must be a non-empty finite two-dimensional predictor matrix")
    return values, variants


def _input_digest(
    variants: tuple[str, ...],
    source_ids: Any,
    truth: Any,
    feature_columns: tuple[str, ...],
    values: np.ndarray,
) -> str:
    """Fingerprint all canonical primary inputs with strict, unambiguous encoding."""
    matrix = np.ascontiguousarray(values, dtype="<f8")
    target = np.ascontiguousarray(np.asarray(truth, dtype="<f8").reshape(-1))
    sources = tuple(source_ids)
    if (
        matrix.ndim != 2
        or matrix.shape[0] != len(variants)
        or len(sources) != len(variants)
        or len(target) != len(variants)
        or not np.isfinite(matrix).all()
        or not np.isfinite(target).all()
        or any(not isinstance(value, str) for value in (*feature_columns, *variants, *sources))
    ):
        raise ValueError("cannot digest non-canonical primary inputs")

    def update_strings(items: tuple[str, ...]) -> None:
        digest.update(len(items).to_bytes(8, byteorder="big"))
        for item in items:
            encoded = item.encode("utf-8")
            digest.update(len(encoded).to_bytes(8, byteorder="big"))
            digest.update(encoded)

    digest = hashlib.sha256()
    digest.update(b"kcat_rel_primary_input_v2\0")
    update_strings(feature_columns)
    update_strings(variants)
    update_strings(sources)
    digest.update(np.asarray(matrix.shape, dtype="<i8").tobytes())
    digest.update(len(target).to_bytes(8, byteorder="big"))
    digest.update(target.tobytes())
    digest.update(matrix.tobytes())
    return digest.hexdigest()


def _as_target(y: Any, n_rows: int) -> np.ndarray:
    values = np.asarray(y, dtype=float).reshape(-1)
    if len(values) != n_rows or not np.isfinite(values).all():
        raise ValueError("y must be finite and match X rows")
    return values


def _as_groups(groups: Any, n_rows: int) -> np.ndarray:
    values = np.asarray(groups, dtype=object).reshape(-1)
    if len(values) != n_rows or any(not isinstance(value, str) or not value.strip() for value in values):
        raise ValueError("groups must be non-empty strings and match X rows")
    return values


def _group_labels(groups: np.ndarray) -> tuple[object, ...]:
    """Provide a deterministic ordering without relying on mixed-type NumPy sorting."""
    return tuple(sorted(set(groups.tolist()), key=lambda value: str(value)))


def inner_group_splits(groups: Any) -> Iterator[IndexSplit]:
    """Yield leave-one-source-study-out splits for nested tuning."""
    values = _as_groups(groups, len(np.asarray(groups).reshape(-1)))
    labels = _group_labels(values)
    if len(labels) < 2:
        raise ValueError("nested source-study CV requires at least two source studies")
    all_indices = np.arange(len(values), dtype=int)
    for label in labels:
        test = all_indices[values == label]
        train = all_indices[values != label]
        yield IndexSplit(train=train, test=test)


def study_macro_metrics(y: Any, prediction: Any, groups: Any) -> Metrics:
    """Average each source study's MAE and RMSE equally, regardless of row count."""
    truth = np.asarray(y, dtype=float).reshape(-1)
    estimate = np.asarray(prediction, dtype=float).reshape(-1)
    source = _as_groups(groups, len(truth))
    if len(estimate) != len(truth) or not np.isfinite(truth).all() or not np.isfinite(estimate).all():
        raise ValueError("y and prediction must be finite arrays of equal length")
    maes: list[float] = []
    rmses: list[float] = []
    for label in _group_labels(source):
        residual = truth[source == label] - estimate[source == label]
        maes.append(float(np.mean(np.abs(residual))))
        rmses.append(float(np.sqrt(np.mean(np.square(residual)))))
    return Metrics(mae=float(np.mean(maes)), rmse=float(np.mean(rmses)))


def variant_metrics(y: Any, prediction: Any) -> Metrics:
    """Score variant LOOCV as one pooled set of held-out variants."""
    truth = np.asarray(y, dtype=float).reshape(-1)
    estimate = np.asarray(prediction, dtype=float).reshape(-1)
    if len(truth) == 0 or len(estimate) != len(truth) or not np.isfinite(truth).all() or not np.isfinite(estimate).all():
        raise ValueError("y and prediction must be non-empty finite arrays of equal length")
    residual = truth - estimate
    return Metrics(mae=float(np.mean(np.abs(residual))), rmse=float(np.sqrt(np.mean(np.square(residual)))))


def _candidate_parameters(model_kind: str, registry: Any) -> tuple[dict[str, float], ...]:
    grids = getattr(registry, "model_grids")
    if model_kind == "ridge":
        return tuple({"alpha": float(alpha)} for alpha in grids["ridge_alphas"])
    if model_kind == "gpr":
        gpr = grids["gpr"]
        return tuple(
            {"length_scale": float(length_scale), "noise_variance": float(noise_variance)}
            for length_scale in gpr["length_scales"]
            for noise_variance in gpr["noise_variances"]
        )
    raise ValueError("model_kind must be 'ridge' or 'gpr'")


def _build_estimator(model_kind: str, params: Mapping[str, float], seed: int):
    if model_kind == "ridge":
        return Ridge(alpha=params["alpha"], fit_intercept=True)
    kernel = (
        ConstantKernel(1.0, constant_value_bounds="fixed")
        * Matern(length_scale=params["length_scale"], length_scale_bounds="fixed", nu=1.5)
        + WhiteKernel(noise_level=params["noise_variance"], noise_level_bounds="fixed")
    )
    return GaussianProcessRegressor(kernel=kernel, optimizer=None, normalize_y=False, random_state=seed)


def _fit_predict(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    model_kind: str,
    params: Mapping[str, float],
    seed: int,
) -> tuple[np.ndarray, tuple[float, ...], tuple[float, ...]]:
    """Fit a fresh predictor scaler and estimator entirely inside one fit."""
    scaler = StandardScaler().fit(X_train)
    estimator = _build_estimator(model_kind, params, seed)
    estimator.fit(scaler.transform(X_train), y_train)
    return (
        np.asarray(estimator.predict(scaler.transform(X_test)), dtype=float),
        tuple(float(value) for value in scaler.mean_),
        tuple(float(value) for value in scaler.scale_),
    )


def _params_key(model_kind: str, params: Mapping[str, float]) -> tuple[float, ...]:
    # Secondary fields encode the approved deterministic preference on exact MAE ties.
    if model_kind == "ridge":
        return (-params["alpha"],)
    return (-params["length_scale"], -params["noise_variance"])


def _choose_inner_params(
    X: np.ndarray,
    y: np.ndarray,
    groups: np.ndarray,
    model_kind: str,
    registry: Any,
    outer_train_indices: np.ndarray,
) -> tuple[dict[str, float], dict[str, float], tuple[InnerFoldResult, ...]]:
    candidates = _candidate_parameters(model_kind, registry)
    scored: list[tuple[float, tuple[float, ...], dict[str, float], str, tuple[InnerFoldResult, ...]]] = []
    for params in candidates:
        prediction = np.empty(len(y), dtype=float)
        inner_folds: list[InnerFoldResult] = []
        for split in inner_group_splits(groups):
            prediction[split.test], scaler_mean, scaler_scale = _fit_predict(
                X[split.train], y[split.train], X[split.test], model_kind, params, int(registry.seed)
            )
            inner_folds.append(
                InnerFoldResult(
                    train_indices=outer_train_indices[split.train],
                    test_indices=outer_train_indices[split.test],
                    scaler_mean=scaler_mean,
                    scaler_scale=scaler_scale,
                )
            )
        score = study_macro_metrics(y, prediction, groups).mae
        name = ",".join(f"{key}={value:g}" for key, value in params.items())
        scored.append((score, _params_key(model_kind, params), params, name, tuple(inner_folds)))
    _, _, params, _, selected_inner_folds = min(scored, key=lambda candidate: (candidate[0], candidate[1]))
    return (
        dict(params),
        {name: float(candidate_score) for candidate_score, _, _, name, _ in scored},
        selected_inner_folds,
    )


def _outer_splits(view: str, groups: np.ndarray) -> Iterator[IndexSplit]:
    indices = np.arange(len(groups), dtype=int)
    if view == "study":
        yield from inner_group_splits(groups)
        return
    if view == "variant":
        for index in indices:
            yield IndexSplit(train=indices[indices != index], test=np.array([index], dtype=int))
        return
    raise ValueError("view must be 'study' or 'variant'")


def run_outer_validation(X: Any, y: Any, groups: Any, view: str, model_kind: str, registry: Any) -> OofResult:
    """Run nested source-study tuning inside each leakage-safe outer fold.

    Feature scaling is always learned from that fit's training rows; targets are
    deliberately never scaled.  Baseline predictions are the same outer-fold
    training mean used for the candidate comparison.
    """
    features, variants = _as_numeric_matrix(X, registry)
    target = _as_target(y, len(features))
    source = _as_groups(groups, len(features))
    if set(source) != set(getattr(registry, "source_studies", FROZEN_SOURCE_STUDIES)):
        raise ValueError("groups do not match the complete frozen source-study set")
    prediction = np.empty(len(target), dtype=float)
    baseline_prediction = np.empty(len(target), dtype=float)
    folds: list[FoldResult] = []
    records: list[dict[str, Any]] = []

    for outer_fold, split in enumerate(_outer_splits(view, source)):
        train_X, train_y = features[split.train], target[split.train]
        params, inner_scores, inner_folds = _choose_inner_params(
            train_X, train_y, source[split.train], model_kind, registry, split.train
        )
        scaler = StandardScaler().fit(train_X)
        estimator = _build_estimator(model_kind, params, int(registry.seed))
        estimator.fit(scaler.transform(train_X), train_y)
        prediction[split.test] = estimator.predict(scaler.transform(features[split.test]))
        baseline_mean = float(np.mean(train_y))
        baseline_prediction[split.test] = baseline_mean
        fold = FoldResult(
            outer_fold=outer_fold,
            train_indices=split.train,
            test_indices=split.test,
            scaler_mean=tuple(float(value) for value in scaler.mean_),
            scaler_scale=tuple(float(value) for value in scaler.scale_),
            baseline_mean=baseline_mean,
            selected_params=params,
            inner_scores=inner_scores,
            inner_folds=inner_folds,
        )
        folds.append(fold)
        for index in split.test:
            records.append(
                {
                    "variant": variants[index],
                    "source_id": str(source[index]),
                    "truth": float(target[index]),
                    "model_prediction": float(prediction[index]),
                    "baseline_prediction": float(baseline_prediction[index]),
                    "outer_fold": outer_fold,
                    "selected_inner_hyperparameters": dict(params),
                    "scaler_mean": fold.scaler_mean,
                    "scaler_scale": fold.scaler_scale,
                }
            )
    if len(records) != len(target) or not np.isfinite(prediction).all() or not np.isfinite(baseline_prediction).all():
        raise RuntimeError("outer validation did not produce complete finite OOF predictions")
    if len({row["variant"] for row in records}) != len(target):
        raise RuntimeError("outer validation did not produce one OOF row per frozen variant")
    return OofResult(
        view=view,
        model_kind=model_kind,
        prediction=prediction,
        baseline_prediction=baseline_prediction,
        folds=tuple(folds),
        records=pd.DataFrame(records),
        variants=variants,
        source_ids=tuple(str(value) for value in source),
        truth=tuple(float(value) for value in target),
        input_digest=_input_digest(variants, source, target, tuple(registry.feature_columns), features),
    )


def single_study_influence(y: Any, model_prediction: Any, baseline_prediction: Any, groups: Any) -> dict[str, float]:
    """Return baseline-minus-model macro-MAE after deletion of each source study."""
    truth = np.asarray(y, dtype=float).reshape(-1)
    model = np.asarray(model_prediction, dtype=float).reshape(-1)
    baseline = np.asarray(baseline_prediction, dtype=float).reshape(-1)
    source = _as_groups(groups, len(truth))
    if len(model) != len(truth) or len(baseline) != len(truth):
        raise ValueError("model and baseline predictions must match y rows")
    if len(_group_labels(source)) < 2:
        raise ValueError("single-study influence requires at least two source studies")
    gains: dict[str, float] = {}
    for label in _group_labels(source):
        keep = source != label
        gains[str(label)] = float(
            study_macro_metrics(truth[keep], baseline[keep], source[keep]).mae
            - study_macro_metrics(truth[keep], model[keep], source[keep]).mae
        )
    return gains


def _fraction(numerator: float, denominator: float) -> float:
    """Preserve the gate's intended ordering when a degenerate baseline is perfect."""
    if denominator != 0.0:
        return numerator / denominator
    if numerator == 0.0:
        return 0.0
    return float(np.copysign(np.inf, numerator))


def _require_views(views: Mapping[str, Metrics]) -> None:
    if set(views) != {"study", "variant"} or not all(isinstance(metric, Metrics) for metric in views.values()):
        raise ValueError("views must contain exactly 'study' and 'variant' Metrics")
    try:
        values = tuple(
            float(value)
            for metric in views.values()
            for value in (metric.mae, metric.rmse)
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("gate metrics must be finite and non-negative") from exc
    if any(not np.isfinite(value) or value < 0.0 for value in values):
        raise ValueError("gate metrics must be finite and non-negative")


def evaluate_gate(
    model_views: Mapping[str, Metrics], baseline_views: Mapping[str, Metrics], influence_gains: Mapping[str, float], registry: Any | None = None
) -> GateDecision:
    """Evaluate every fixed primary gate for one candidate model."""
    _require_views(model_views)
    _require_views(baseline_views)
    thresholds = getattr(registry, "gate_thresholds", None) if registry is not None else None
    minimum_improvement = float((thresholds or {}).get("minimum_mae_improvement_fraction", 0.05))
    maximum_degradation = float((thresholds or {}).get("maximum_rmse_degradation_fraction", 0.05))
    mae_improvement = {
        view: _fraction(baseline_views[view].mae - model_views[view].mae, baseline_views[view].mae)
        for view in ("study", "variant")
    }
    rmse_degradation = {
        view: _fraction(model_views[view].rmse - baseline_views[view].rmse, baseline_views[view].rmse)
        for view in ("study", "variant")
    }
    failures: list[str] = []
    if mae_improvement["study"] < minimum_improvement:
        failures.append("study_macro_mae")
    if mae_improvement["variant"] < minimum_improvement:
        failures.append("variant_mae")
    if rmse_degradation["study"] > maximum_degradation:
        failures.append("study_macro_rmse")
    if rmse_degradation["variant"] > maximum_degradation:
        failures.append("variant_rmse")
    expected_studies = set(getattr(registry, "source_studies", FROZEN_SOURCE_STUDIES))
    if (
        set(influence_gains) != expected_studies
        or not all(np.isfinite(float(gain)) and float(gain) > 0.0 for gain in influence_gains.values())
    ):
        failures.append("single_study_influence")
    return GateDecision(
        pass_gate=not failures,
        failed_rules=tuple(failures),
        mae_improvement=mae_improvement,
        rmse_degradation=rmse_degradation,
        model_views=dict(model_views),
        baseline_views=dict(baseline_views),
    )


def select_winner(decisions: Mapping[str, GateDecision], registry: Any | None = None) -> str | None:
    """Select a passing candidate, preferring Ridge inside the approved tie band."""
    passing = {name: decision for name, decision in decisions.items() if decision.pass_gate}
    if not passing:
        return None
    if len(passing) == 1:
        return next(iter(passing))
    if "ridge" in passing and "gpr" in passing:
        ridge_mae = passing["ridge"].model_views["study"].mae
        gpr_mae = passing["gpr"].model_views["study"].mae
        better = min(ridge_mae, gpr_mae)
        thresholds = getattr(registry, "gate_thresholds", None) if registry is not None else None
        tie_fraction = float((thresholds or {}).get("winner_effective_tie_fraction", 0.02))
        if _fraction(abs(ridge_mae - gpr_mae), better) <= tie_fraction:
            return "ridge"
    return min(passing, key=lambda name: (passing[name].model_views["study"].mae, name))

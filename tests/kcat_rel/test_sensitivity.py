"""Secondary uncertainty sensitivities preserve the frozen primary diagnostic."""

from __future__ import annotations

import copy
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.linear_model import Ridge

from kcat_rel.registry import ACTIVE_FEATURE_COLUMNS, load_registry
from kcat_rel.sensitivity import run_perturbation_sensitivity, run_weighted_sensitivity
from kcat_rel.validation import (
    evaluate_gate,
    run_outer_validation,
    single_study_influence,
    study_macro_metrics,
    variant_metrics,
)


WORKFLOW_ROOT = Path(__file__).resolve().parents[2] / "results/kcat-rel"
CONFIG = WORKFLOW_ROOT / "config" / "kcat-rel-v0.json"


@pytest.fixture
def registry():
    return load_registry(CONFIG)


@pytest.fixture
def primary_inputs(registry):
    X = pd.DataFrame(
        {
            "variant": registry.variants,
            "delta_charge": np.linspace(-2.0, 2.0, 17),
            "grantham_distance": np.arange(17, dtype=float) * 11.0 + 5.0,
            "catalytic_motif_min_ca_distance_8act_mean_ab": np.arange(17, dtype=float) ** 2 + 3.0,
        }
    )
    y = np.array(
        [-1.0, 1.0, 4.0, 6.0, 20.0, 22.0, 3.0, 5.0, -2.0, 8.0, 9.0, -3.0, 2.0, -4.0, 7.0, 0.0, 11.0]
    )
    groups = np.array(
        [
            "nandwani2025", "adhikari2016", "adhikari2019", "adhikari2016", "nag2015", "sommese2013",
            "adhikari2019", "nandwani2025", "sarkar2020", "kawana2017", "kawana2017", "pathak2026",
            "morck2022", "morck2022", "morck2022", "morck2022", "morck2022",
        ]
    )
    errors = np.linspace(0.1, 0.9, len(y))
    return X, y, errors, groups


@pytest.fixture
def ridge_primary(primary_inputs, registry):
    X, y, _, groups = primary_inputs
    return run_outer_validation(X, y, groups, "study", "ridge", registry)


def _manual_frozen_prediction(primary, X, y, *, errors=None):
    """Independent estimator oracle using only the primary fold audit state."""
    features = X.loc[:, ACTIVE_FEATURE_COLUMNS].to_numpy(dtype=float)
    prediction = np.empty(len(y), dtype=float)
    baseline = np.empty(len(y), dtype=float)
    for fold in primary.folds:
        train, test = fold.train_indices, fold.test_indices
        scaled_train = (features[train] - np.asarray(fold.scaler_mean)) / np.asarray(fold.scaler_scale)
        scaled_test = (features[test] - np.asarray(fold.scaler_mean)) / np.asarray(fold.scaler_scale)
        if primary.model_kind == "ridge":
            weights = None
            if errors is not None:
                weights = 1.0 / np.square(errors[train])
                weights /= weights.mean()
            estimator = Ridge(alpha=fold.selected_params["alpha"], fit_intercept=True)
            estimator.fit(scaled_train, y[train], sample_weight=weights)
        else:
            kernel = (
                ConstantKernel(1.0, constant_value_bounds="fixed")
                * Matern(
                    length_scale=fold.selected_params["length_scale"],
                    length_scale_bounds="fixed",
                    nu=1.5,
                )
                + WhiteKernel(
                    noise_level=fold.selected_params["noise_variance"],
                    noise_level_bounds="fixed",
                )
            )
            estimator = GaussianProcessRegressor(
                kernel=kernel,
                alpha=0.0 if errors is None else np.square(errors[train]),
                optimizer=None,
                normalize_y=False,
                random_state=20260916,
            ).fit(scaled_train, y[train])
        prediction[test] = estimator.predict(scaled_test)
        baseline[test] = y[train].mean()
    return prediction, baseline


def _assert_primary_unchanged(before, after) -> None:
    assert before.view == after.view
    assert before.model_kind == after.model_kind
    assert before.variants == after.variants
    assert before.source_ids == after.source_ids
    assert before.truth == after.truth
    assert before.input_digest == after.input_digest
    assert before.fold_hyperparameters == after.fold_hyperparameters
    assert np.array_equal(before.prediction, after.prediction)
    assert np.array_equal(before.baseline_prediction, after.baseline_prediction)
    pd.testing.assert_frame_equal(before.records, after.records)
    for before_fold, after_fold in zip(before.folds, after.folds, strict=True):
        assert before_fold.outer_fold == after_fold.outer_fold
        assert np.array_equal(before_fold.train_indices, after_fold.train_indices)
        assert np.array_equal(before_fold.test_indices, after_fold.test_indices)
        assert before_fold.scaler_mean == after_fold.scaler_mean
        assert before_fold.scaler_scale == after_fold.scaler_scale
        assert before_fold.baseline_mean == after_fold.baseline_mean
        assert before_fold.selected_params == after_fold.selected_params
        assert before_fold.inner_scores == after_fold.inner_scores
        assert len(before_fold.inner_folds) == len(after_fold.inner_folds)
        for before_inner, after_inner in zip(before_fold.inner_folds, after_fold.inner_folds, strict=True):
            assert np.array_equal(before_inner.train_indices, after_inner.train_indices)
            assert np.array_equal(before_inner.test_indices, after_inner.test_indices)
            assert before_inner.scaler_mean == after_inner.scaler_mean
            assert before_inner.scaler_scale == after_inner.scaler_scale


def test_weighted_sensitivity_reuses_primary_hyperparameters_and_cannot_select_winner(
    primary_inputs, ridge_primary
) -> None:
    """A secondary fit must not re-tune or expose a replacement primary choice."""
    X, y, errors, groups = primary_inputs

    result = run_weighted_sensitivity(ridge_primary, X, y, errors, groups)

    assert result.fold_hyperparameters == ridge_primary.fold_hyperparameters
    assert not hasattr(result, "selected_model")
    assert result.model_kind == "ridge"
    assert np.isfinite(result.prediction).all()


def test_ridge_inverse_variance_weights_are_normalized_inside_each_training_fold(
    primary_inputs, ridge_primary
) -> None:
    """Using cohort-wide normalization would let held-out uncertainty affect a fit."""
    X, y, errors, groups = primary_inputs

    result = run_weighted_sensitivity(ridge_primary, X, y, errors, groups)

    for primary_fold, weights in zip(ridge_primary.folds, result.fold_sample_weights, strict=True):
        expected = 1.0 / np.square(errors[primary_fold.train_indices])
        expected /= expected.mean()
        assert weights == pytest.approx(expected)
        assert np.mean(weights) == pytest.approx(1.0)


def test_weighted_ridge_prediction_uses_the_frozen_scaler_and_hyperparameters(
    primary_inputs, ridge_primary
) -> None:
    """Refitting a scaler or alpha would diverge from this fold-audit oracle."""
    X, y, errors, groups = primary_inputs

    result = run_weighted_sensitivity(ridge_primary, X, y, errors, groups)
    expected, expected_baseline = _manual_frozen_prediction(ridge_primary, X, y, errors=errors)

    assert result.prediction == pytest.approx(expected)
    assert result.baseline_prediction == pytest.approx(expected_baseline)


def test_gpr_uncertainty_sensitivity_adds_reported_variance_to_frozen_noise(
    primary_inputs, registry
) -> None:
    """Replacing, rather than adding to, selected homoscedastic noise changes the model."""
    X, y, errors, groups = primary_inputs
    primary = run_outer_validation(X, y, groups, "study", "gpr", registry)

    result = run_weighted_sensitivity(primary, X, y, errors, groups)
    fold = primary.folds[0]
    features = X.loc[:, registry.feature_columns].to_numpy(dtype=float)
    scaled_train = (features[fold.train_indices] - np.asarray(fold.scaler_mean)) / np.asarray(fold.scaler_scale)
    scaled_test = (features[fold.test_indices] - np.asarray(fold.scaler_mean)) / np.asarray(fold.scaler_scale)
    kernel = (
        ConstantKernel(1.0, constant_value_bounds="fixed")
        * Matern(length_scale=fold.selected_params["length_scale"], length_scale_bounds="fixed", nu=1.5)
        + WhiteKernel(noise_level=fold.selected_params["noise_variance"], noise_level_bounds="fixed")
    )
    expected = GaussianProcessRegressor(
        kernel=kernel,
        alpha=np.square(errors[fold.train_indices]),
        optimizer=None,
        normalize_y=False,
        random_state=registry.seed,
    ).fit(scaled_train, y[fold.train_indices]).predict(scaled_test)

    assert result.prediction[fold.test_indices] == pytest.approx(expected)
    assert result.fold_reported_variances[0] == pytest.approx(np.square(errors[fold.train_indices]))


def test_sensitivity_rejects_any_schema_or_identity_drift(primary_inputs, ridge_primary) -> None:
    """Width-only checking would admit outcome-adjacent or replacement predictors."""
    X, y, errors, groups = primary_inputs
    invalid_tables = (
        X.rename(columns={"delta_charge": "renamed_feature"}),
        X.loc[:, ["variant", *reversed(ACTIVE_FEATURE_COLUMNS)]],
        X.assign(outcome_adjacent=y),
        X.assign(variant=list(reversed(X.variant))),
    )

    for invalid in invalid_tables:
        with pytest.raises(ValueError, match="schema|variant"):
            run_weighted_sensitivity(ridge_primary, invalid, y, errors, groups)


def test_sensitivity_rejects_features_or_cohort_not_identical_to_the_primary_snapshot(
    primary_inputs, ridge_primary
) -> None:
    """A matching schema is insufficient when predictor values or mutable OOF IDs differ."""
    X, y, errors, groups = primary_inputs
    changed_values = X.copy()
    changed_values.loc[0, "grantham_distance"] += 0.01
    forged_primary = copy.deepcopy(ridge_primary)
    forged_primary.records.loc[0, "variant"] = "D382Y"
    forged_X = X.copy()
    for fold in forged_primary.folds:
        rows = forged_primary.records.loc[
            forged_primary.records.outer_fold == fold.outer_fold, "variant"
        ].tolist()
        forged_X.loc[fold.test_indices, "variant"] = rows

    with pytest.raises(ValueError, match="snapshot|primary|cohort"):
        run_weighted_sensitivity(ridge_primary, changed_values, y, errors, groups)
    with pytest.raises(ValueError, match="snapshot|primary|cohort"):
        run_weighted_sensitivity(forged_primary, forged_X, y, errors, groups)


def test_sensitivity_rejects_a_coordinated_source_and_truth_forgery(primary_inputs, ridge_primary) -> None:
    """A digest excluding source/target values would accept coordinated mutable snapshots."""
    X, y, errors, groups = primary_inputs
    forged_groups = groups.astype(object)
    forged_groups[0] = "forged_source"
    forged_truth = y.copy()
    forged_truth[0] += 0.01
    forged = replace(
        ridge_primary,
        source_ids=tuple(forged_groups),
        truth=tuple(forged_truth),
    )
    forged_records = forged.records.copy()
    for fold in forged.folds:
        indices = fold.test_indices
        mask = forged_records.outer_fold == fold.outer_fold
        forged_records.loc[mask, "source_id"] = forged_groups[indices]
        forged_records.loc[mask, "truth"] = forged_truth[indices]
    forged = replace(forged, records=forged_records)

    with pytest.raises(ValueError, match="digest|snapshot|primary"):
        run_weighted_sensitivity(forged, X, forged_truth, errors, forged_groups)


def test_sensitivity_binds_target_to_frozen_oof_truth_and_does_not_mutate_primary(
    primary_inputs, ridge_primary
) -> None:
    """Changing the supplied labels must not create a new sensitivity target."""
    X, y, errors, groups = primary_inputs
    before = copy.deepcopy(ridge_primary)

    run_weighted_sensitivity(ridge_primary, X, y, errors, groups)
    with pytest.raises(ValueError, match="target|truth"):
        run_weighted_sensitivity(ridge_primary, X, y + 0.01, errors, groups)
    with pytest.raises(ValueError, match="target|truth"):
        run_perturbation_sensitivity(ridge_primary, X, y[::-1], errors, groups, draws=2)

    _assert_primary_unchanged(before, ridge_primary)


def test_two_view_perturbation_requires_every_primary_snapshot_to_agree(primary_inputs, registry) -> None:
    """Checking only the first mapping entry would allow a forged second validation view."""
    X, y, errors, groups = primary_inputs
    study = run_outer_validation(X, y, groups, "study", "ridge", registry)
    variant = run_outer_validation(X, y, groups, "variant", "ridge", registry)
    forged_records = variant.records.copy()
    forged_records.loc[0, "truth"] += 0.01
    forged_variant = replace(variant, records=forged_records)

    with pytest.raises(ValueError, match="truth|primary|snapshot"):
        run_perturbation_sensitivity(
            {"study": study, "variant": forged_variant}, X, y, errors, groups, draws=2
        )


@pytest.mark.parametrize("model_kind", ["ridge", "gpr"])
def test_perturbation_matches_independent_frozen_refit_quantiles_and_gate_rate(
    model_kind,
    primary_inputs, registry
) -> None:
    """The RNG draws must drive frozen refits, metrics, quantiles, and the unchanged gate."""
    X, y, errors, groups = primary_inputs
    primary = {
        "study": run_outer_validation(X, y, groups, "study", model_kind, registry),
        "variant": run_outer_validation(X, y, groups, "variant", model_kind, registry),
    }
    before = copy.deepcopy(primary)

    a = run_perturbation_sensitivity(primary, X, y, errors, groups, draws=4, seed=20260916)
    b = run_perturbation_sensitivity(primary, X, y, errors, groups, draws=4, seed=20260916)
    expected_values = {
        view: {"mae": [], "rmse": [], "baseline_mae": [], "baseline_rmse": []}
        for view in primary
    }
    expected_passes = 0
    rng = np.random.default_rng(20260916)
    for _ in range(4):
        perturbed = y + rng.normal(0.0, errors, size=len(y))
        model_metrics, baseline_metrics, prediction, baseline = {}, {}, {}, {}
        for view, result in primary.items():
            prediction[view], baseline[view] = _manual_frozen_prediction(result, X, perturbed)
            if view == "study":
                model_metrics[view] = study_macro_metrics(perturbed, prediction[view], groups)
                baseline_metrics[view] = study_macro_metrics(perturbed, baseline[view], groups)
            else:
                model_metrics[view] = variant_metrics(perturbed, prediction[view])
                baseline_metrics[view] = variant_metrics(perturbed, baseline[view])
            expected_values[view]["mae"].append(model_metrics[view].mae)
            expected_values[view]["rmse"].append(model_metrics[view].rmse)
            expected_values[view]["baseline_mae"].append(baseline_metrics[view].mae)
            expected_values[view]["baseline_rmse"].append(baseline_metrics[view].rmse)
        influence = single_study_influence(perturbed, prediction["study"], baseline["study"], groups)
        expected_passes += evaluate_gate(model_metrics, baseline_metrics, influence).pass_gate

    assert a == b
    assert a.fold_hyperparameters == {view: result.fold_hyperparameters for view, result in primary.items()}
    for view, metrics in expected_values.items():
        for metric_name, draw_values in metrics.items():
            assert a.metric_quantiles[view][metric_name] == pytest.approx(
                {
                    "q025": np.quantile(draw_values, 0.025),
                    "q500": np.quantile(draw_values, 0.5),
                    "q975": np.quantile(draw_values, 0.975),
                }
            )
    assert a.gate_pass_count == expected_passes
    assert a.gate_pass_rate == pytest.approx(expected_passes / 4.0)
    assert not hasattr(a, "selected_model")
    for view in primary:
        _assert_primary_unchanged(before[view], primary[view])


@pytest.mark.parametrize("bad_errors", [
    np.array([0.1] * 16),
    np.array([-0.1] + [0.1] * 16),
    np.array([np.nan] + [0.1] * 16),
])
def test_sensitivities_reject_invalid_uncertainty_shapes_and_values(
    primary_inputs, ridge_primary, bad_errors
) -> None:
    """Malformed uncertainty cannot silently alter sample weights or perturbations."""
    X, y, _, groups = primary_inputs

    with pytest.raises(ValueError, match="y_error"):
        run_weighted_sensitivity(ridge_primary, X, y, bad_errors, groups)
    with pytest.raises(ValueError, match="y_error"):
        run_perturbation_sensitivity(ridge_primary, X, y, bad_errors, groups, draws=2)


def test_ridge_weighted_sensitivity_rejects_zero_uncertainty(primary_inputs, ridge_primary) -> None:
    """A zero reported standard deviation has undefined inverse-variance weight."""
    X, y, errors, groups = primary_inputs
    errors[0] = 0.0

    with pytest.raises(ValueError, match="positive"):
        run_weighted_sensitivity(ridge_primary, X, y, errors, groups)

"""Behavioral tests for leakage-safe nested ``kcat_rel`` validation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from dataclasses import replace
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.preprocessing import StandardScaler

from kcat_rel.registry import load_registry
from kcat_rel.validation import (
    Metrics,
    evaluate_gate,
    inner_group_splits,
    run_outer_validation,
    select_winner,
    single_study_influence,
    study_macro_metrics,
    variant_metrics,
)


WORKFLOW_ROOT = Path(__file__).resolve().parents[2] / "workflows/kcat-rel"
CONFIG = WORKFLOW_ROOT / "config" / "kcat-rel-v0.json"
CANONICAL_SOURCE_IDS = (
    "adhikari2016", "adhikari2019", "kawana2017", "morck2022", "nag2015", "nandwani2025", "pathak2026",
    "sarkar2020", "sommese2013",
)


@pytest.fixture
def registry():
    return load_registry(CONFIG)


@pytest.fixture
def primary_inputs(registry):
    """The complete, ordered primary cohort and its nine canonical sources."""
    X = pd.DataFrame(
        {
            "variant": registry.variants,
            "delta_charge": np.linspace(-2.0, 2.0, 17),
            "grantham_distance": np.arange(17, dtype=float) * 11.0 + 5.0,
            "catalytic_motif_min_ca_distance_8act_mean_ab": np.arange(17, dtype=float) ** 2 + 3.0,
        }
    )
    y = np.array([-1.0, 1.0, 4.0, 6.0, 20.0, 22.0, 3.0, 5.0, -2.0, 8.0, 9.0, -3.0, 2.0, -4.0, 7.0, 0.0, 11.0])
    groups = np.array(
        [
            "nandwani2025", "adhikari2016", "adhikari2019", "adhikari2016", "nag2015", "sommese2013",
            "adhikari2019", "nandwani2025", "sarkar2020", "kawana2017", "kawana2017", "pathak2026",
            "morck2022", "morck2022", "morck2022", "morck2022", "morck2022",
        ]
    )
    return X, y, groups


def test_scaler_and_mean_baseline_use_outer_training_rows_only(primary_inputs, registry) -> None:
    """A fit that includes the held-out study would leak its feature/target mean."""
    X, y, groups = primary_inputs

    result = run_outer_validation(X, y, groups, "study", "ridge", registry)
    first = result.folds[0]

    assert first.scaler_mean == pytest.approx(
        X.loc[first.train_indices, registry.feature_columns].to_numpy().mean(axis=0)
    )
    assert np.all(result.baseline_prediction[first.test_indices] == y[first.train_indices].mean())


def test_inner_splits_never_cross_source_studies(primary_inputs) -> None:
    """An inner split that shares a source across train/test would tune on leakage."""
    _, _, groups = primary_inputs

    for split in inner_group_splits(groups):
        assert set(groups[split.train]).isdisjoint(set(groups[split.test]))


def test_exact_inner_ties_choose_the_preferred_approved_hyperparameters(primary_inputs, registry) -> None:
    """Choosing the first grid value on a score tie would violate frozen tie rules."""
    X, _, groups = primary_inputs
    y = np.zeros(len(X), dtype=float)

    ridge = run_outer_validation(X, y, groups, "study", "ridge", registry)
    gpr = run_outer_validation(X, y, groups, "study", "gpr", registry)

    assert all(fold.selected_params == {"alpha": 100.0} for fold in ridge.folds)
    assert all(
        fold.selected_params == {"length_scale": 2.0, "noise_variance": 0.04} for fold in gpr.folds
    )


def test_oof_result_keeps_fold_hyperparameters_for_frozen_sensitivities(primary_inputs, registry) -> None:
    """A sensitivity that re-tunes instead of reusing primary choices would alter the primary result."""
    X, y, groups = primary_inputs

    result = run_outer_validation(X, y, groups, "study", "ridge", registry)

    assert result.fold_hyperparameters == tuple(fold.selected_params for fold in result.folds)
    assert set(result.oof.columns) >= {
        "variant", "source_id", "truth", "model_prediction", "baseline_prediction", "outer_fold",
        "selected_inner_hyperparameters", "scaler_mean", "scaler_scale",
    }


def test_oof_result_keeps_immutable_primary_input_snapshots(primary_inputs, registry) -> None:
    """Secondary work must bind to values captured by the primary run, not mutable OOF rows."""
    X, y, groups = primary_inputs

    result = run_outer_validation(X, y, groups, "study", "ridge", registry)

    assert result.variants == registry.variants
    assert result.source_ids == tuple(groups)
    assert result.truth == tuple(y)
    assert isinstance(result.input_digest, str)
    assert len(result.input_digest) == 64
    changed_target = run_outer_validation(X, y + 0.01, groups, "study", "ridge", registry)
    swapped_sources = groups.copy()
    swapped_sources[[0, 1]] = swapped_sources[[1, 0]]
    changed_sources = run_outer_validation(X, y, swapped_sources, "study", "ridge", registry)
    assert changed_target.input_digest != result.input_digest
    assert changed_sources.input_digest != result.input_digest


def test_inner_scalers_fit_only_their_own_nested_training_rows(primary_inputs, registry) -> None:
    """Fitting an inner scaler on an outer-training holdout would leak its feature values into tuning."""
    X, y, groups = primary_inputs

    result = run_outer_validation(X, y, groups, "study", "ridge", registry)
    outer = result.folds[0]
    inner = outer.inner_folds[0]

    expected = X.loc[inner.train_indices, registry.feature_columns].to_numpy().mean(axis=0)
    assert inner.scaler_mean == pytest.approx(expected)


def test_gpr_predictions_use_the_original_unscaled_target(primary_inputs, registry) -> None:
    """Enabling GPR target normalization changes fixed-kernel predictions and is prohibited."""
    X, y, groups = primary_inputs

    result = run_outer_validation(X, y, groups, "study", "gpr", registry)
    fold = result.folds[0]
    features = X.loc[:, registry.feature_columns].to_numpy()
    scaler = StandardScaler().fit(features[fold.train_indices])
    kernel = (
        ConstantKernel(1.0, constant_value_bounds="fixed")
        * Matern(length_scale=fold.selected_params["length_scale"], length_scale_bounds="fixed", nu=1.5)
        + WhiteKernel(noise_level=fold.selected_params["noise_variance"], noise_level_bounds="fixed")
    )
    raw_target_fit = GaussianProcessRegressor(kernel=kernel, optimizer=None, normalize_y=False).fit(
        scaler.transform(features[fold.train_indices]), y[fold.train_indices]
    ).predict(scaler.transform(features[fold.test_indices]))
    normalized_target_fit = GaussianProcessRegressor(kernel=kernel, optimizer=None, normalize_y=True).fit(
        scaler.transform(features[fold.train_indices]), y[fold.train_indices]
    ).predict(scaler.transform(features[fold.test_indices]))

    assert result.prediction[fold.test_indices] == pytest.approx(raw_target_fit)
    assert result.prediction[fold.test_indices] != pytest.approx(normalized_target_fit)
    assert fold.baseline_mean == pytest.approx(y[fold.train_indices].mean())


def test_validation_accepts_only_the_frozen_feature_schema_and_variant_identities(primary_inputs, registry) -> None:
    """Accepting arbitrary columns or identities would reintroduce outcome-adjacent predictors or a wrong cohort."""
    X, y, groups = primary_inputs
    indexed_features = X.loc[:, registry.feature_columns].copy()
    indexed_features.index = X.variant

    accepted = run_outer_validation(indexed_features, y, groups, "study", "ridge", registry)
    assert set(accepted.oof.variant) == set(registry.variants)
    with pytest.raises(ValueError, match="DataFrame|variant"):
        run_outer_validation(X.loc[:, registry.feature_columns].to_numpy(), y, groups, "study", "ridge", registry)
    with pytest.raises(ValueError, match="schema"):
        run_outer_validation(X.assign(y=y), y, groups, "study", "ridge", registry)
    with pytest.raises(ValueError, match="feature"):
        run_outer_validation(X.drop(columns=registry.feature_columns[-1]), y, groups, "study", "ridge", registry)
    with pytest.raises(ValueError, match="cohort|variant"):
        run_outer_validation(X.assign(variant=["D382Y", *registry.variants[1:]]), y, groups, "study", "ridge", registry)


@pytest.mark.parametrize("bad_group", [np.nan, 7, ""])
def test_validation_rejects_non_string_or_empty_source_ids(primary_inputs, registry, bad_group) -> None:
    """Coercing malformed source IDs into labels would silently corrupt grouped validation."""
    X, y, groups = primary_inputs
    groups = groups.astype(object)
    groups[0] = bad_group

    with pytest.raises(ValueError, match="groups"):
        run_outer_validation(X, y, groups, "study", "ridge", registry)


def test_validation_produces_one_oof_row_for_every_frozen_variant(primary_inputs, registry) -> None:
    """A skipped outer fold would leave an incomplete OOF result available for scoring."""
    X, y, groups = primary_inputs

    result = run_outer_validation(X, y, groups, "study", "ridge", registry)

    assert len(result.oof) == len(registry.variants)
    assert tuple(sorted(result.oof.variant)) == tuple(sorted(registry.variants))
    assert result.oof.variant.is_unique
    assert np.isfinite(result.prediction).all()
    assert np.isfinite(result.baseline_prediction).all()


def test_study_macro_metrics_weight_studies_equally() -> None:
    """Pooling rows would let a larger source study dominate the diagnostic score."""
    y = np.array([0.0, 10.0, 10.0, 10.0])
    prediction = np.array([2.0, 10.0, 10.0, 0.0])
    groups = np.array(["small", "large", "large", "large"])

    metrics = study_macro_metrics(y, prediction, groups)

    assert metrics.mae == pytest.approx((2.0 + (10.0 / 3.0)) / 2.0)
    assert metrics.rmse == pytest.approx((2.0 + np.sqrt(100.0 / 3.0)) / 2.0)


def test_variant_metrics_pool_all_held_out_variants() -> None:
    """Treating variant LOOCV as grouped data would score the wrong outer view."""
    metrics = variant_metrics(np.array([0.0, 2.0]), np.array([1.0, -1.0]))

    assert metrics == Metrics(mae=2.0, rmse=np.sqrt(5.0))


def test_single_study_influence_recomputes_equal_study_macro_mae() -> None:
    """Deleting rows rather than re-aggregating studies can hide a single-study failure."""
    y = np.array([0.0, 10.0, 10.0, 10.0])
    model = np.array([1.0, 10.0, 10.0, 2.0])
    baseline = np.array([4.0, 10.0, 10.0, 0.0])
    groups = np.array(["small", "large", "large", "large"])

    gains = single_study_influence(y, model, baseline, groups)

    assert gains == pytest.approx({"large": 3.0, "small": 2.0 / 3.0})


def test_gate_requires_every_threshold_and_each_deletion_positive(registry) -> None:
    """Passing aggregate scores despite a harmful deletion would wrongly authorize prediction."""
    baseline = {"study": Metrics(mae=10.0, rmse=10.0), "variant": Metrics(mae=10.0, rmse=10.0)}
    model = {"study": Metrics(mae=9.0, rmse=10.4), "variant": Metrics(mae=9.0, rmse=10.4)}

    gains = {source: 0.2 for source in CANONICAL_SOURCE_IDS}
    gains["adhikari2016"] = 0.0
    decision = evaluate_gate(model, baseline, gains, registry)

    assert decision.pass_gate is False
    assert decision.failed_rules == ("single_study_influence",)


def test_gate_rejects_an_incomplete_source_study_deletion_set(registry) -> None:
    """Omitting a source from deletion analysis would leave a required influence gate untested."""
    baseline = {"study": Metrics(mae=10.0, rmse=10.0), "variant": Metrics(mae=10.0, rmse=10.0)}
    model = {"study": Metrics(mae=9.5, rmse=10.5), "variant": Metrics(mae=9.5, rmse=10.5)}
    incomplete_gains = {source: 0.1 for source in CANONICAL_SOURCE_IDS[:-1]}

    decision = evaluate_gate(model, baseline, incomplete_gains, registry)

    assert decision.pass_gate is False
    assert decision.failed_rules == ("single_study_influence",)


def test_gate_accepts_exact_five_percent_boundaries_and_rejects_wrong_sign(registry) -> None:
    """Strict inequalities at 5% would change the predeclared forward-prediction decision."""
    baseline = {"study": Metrics(mae=10.0, rmse=10.0), "variant": Metrics(mae=10.0, rmse=10.0)}
    gains = {source: 0.1 for source in CANONICAL_SOURCE_IDS}

    passing = evaluate_gate(
        {"study": Metrics(mae=9.5, rmse=10.5), "variant": Metrics(mae=9.5, rmse=10.5)}, baseline, gains, registry
    )
    worse_mae = evaluate_gate(
        {"study": Metrics(mae=10.5, rmse=10.5), "variant": Metrics(mae=10.5, rmse=10.5)}, baseline, gains, registry
    )
    worse_rmse = evaluate_gate(
        {"study": Metrics(mae=9.5, rmse=10.5001), "variant": Metrics(mae=9.5, rmse=10.5001)}, baseline, gains, registry
    )

    assert passing.pass_gate is True
    assert worse_mae.failed_rules == ("study_macro_mae", "variant_mae")
    assert worse_rmse.failed_rules == ("study_macro_rmse", "variant_rmse")


@pytest.mark.parametrize(
    ("owner", "view", "field", "invalid"),
    [
        ("model", "study", "mae", np.nan),
        ("model", "variant", "rmse", np.inf),
        ("baseline", "study", "rmse", -np.inf),
        ("baseline", "variant", "mae", -0.01),
    ],
)
def test_gate_rejects_non_finite_or_negative_metrics(registry, owner, view, field, invalid) -> None:
    """NaN/Inf comparison semantics must not let an invalid diagnostic metric pass a gate."""
    baseline = {"study": Metrics(mae=10.0, rmse=10.0), "variant": Metrics(mae=10.0, rmse=10.0)}
    model = {"study": Metrics(mae=9.0, rmse=10.0), "variant": Metrics(mae=9.0, rmse=10.0)}
    gains = {source: 0.1 for source in CANONICAL_SOURCE_IDS}

    (model if owner == "model" else baseline)[view] = replace(
        (model if owner == "model" else baseline)[view], **{field: invalid}
    )

    with pytest.raises(ValueError, match="finite and non-negative"):
        evaluate_gate(model, baseline, gains, registry)


def test_winner_prefers_ridge_within_effective_study_mae_tie(registry) -> None:
    """Selecting the lower raw score inside the 2% tie band violates the conservative winner rule."""
    baseline = {"study": Metrics(mae=10.0, rmse=10.0), "variant": Metrics(mae=10.0, rmse=10.0)}
    ridge = evaluate_gate(
        {"study": Metrics(mae=9.18, rmse=10.1), "variant": Metrics(mae=9.18, rmse=10.1)},
        baseline,
        {source: 0.1 for source in CANONICAL_SOURCE_IDS}, registry,
    )
    gpr = evaluate_gate(
        {"study": Metrics(mae=9.0, rmse=10.1), "variant": Metrics(mae=9.0, rmse=10.1)},
        baseline,
        {source: 0.1 for source in CANONICAL_SOURCE_IDS}, registry,
    )

    assert select_winner({"ridge": ridge, "gpr": gpr}, registry) == "ridge"


def test_winner_selects_lower_study_mae_outside_the_two_percent_tie(registry) -> None:
    """Preferring Ridge beyond the effective-tie band would ignore a materially better GPR score."""
    baseline = {"study": Metrics(mae=10.0, rmse=10.0), "variant": Metrics(mae=10.0, rmse=10.0)}
    gains = {source: 0.1 for source in CANONICAL_SOURCE_IDS}
    ridge = evaluate_gate(
        {"study": Metrics(mae=9.19, rmse=10.1), "variant": Metrics(mae=9.19, rmse=10.1)}, baseline, gains, registry
    )
    gpr = evaluate_gate(
        {"study": Metrics(mae=9.0, rmse=10.1), "variant": Metrics(mae=9.0, rmse=10.1)}, baseline, gains, registry
    )

    assert select_winner({"ridge": ridge, "gpr": gpr}, registry) == "gpr"

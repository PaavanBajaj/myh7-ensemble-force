"""Synthetic LOSO GPR/baseline tests. Do not import or use S1 Na_rel values."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import LeaveOneGroupOut

from lsar_na_rel.config import KernelSettings, load_talk_v0_config
from lsar_na_rel.validation import (
    ScaledGPR,
    mae,
    make_baseline_estimator,
    rmse,
    run_loso_gpr_and_baseline,
)

CONFIG = load_talk_v0_config(
    Path(__file__).resolve().parents[2] / "workflows" / "lsar-na-rel" / "config" / "lsar-na-rel.json"
)
SEED = CONFIG.seed


def run_validation(X, y, groups, *, random_state=SEED, kernel_settings=None):
    return run_loso_gpr_and_baseline(
        X,
        y,
        groups,
        random_state=random_state,
        kernel_settings=CONFIG.kernel if kernel_settings is None else kernel_settings,
    )

# Hand-calculated: (|1-1| + |2-3| + |4-2|) / 3 = 1.0
SIMPLE_Y_TRUE = np.array([1.0, 2.0, 4.0])
SIMPLE_Y_PRED = np.array([1.0, 3.0, 2.0])
SIMPLE_MAE = 1.0
SIMPLE_MSE = 5.0 / 3.0

# Hand-calculated LOSO mean baseline on y = 0..11, groups of 2.
BASELINE_Y = np.arange(12, dtype=float)
BASELINE_GROUPS = np.array([f"g{i // 2}" for i in range(12)])
BASELINE_PRED = np.array(
    [6.5, 6.5, 6.1, 6.1, 5.7, 5.7, 5.3, 5.3, 4.9, 4.9, 4.5, 4.5]
)
BASELINE_MAE = 3.6
BASELINE_MSE = 17.05


def test_mae_matches_hand_calculated_simple_fixture():
    assert mae(SIMPLE_Y_TRUE, SIMPLE_Y_PRED) == SIMPLE_MAE


def test_rmse_squared_matches_hand_calculated_simple_mse():
    value = rmse(SIMPLE_Y_TRUE, SIMPLE_Y_PRED)
    assert value**2 == pytest.approx(SIMPLE_MSE, abs=1e-12)


def test_loso_baseline_predictions_match_hand_calculated_train_means():
    X = np.zeros((12, 1), dtype=float)
    result = run_validation(X, BASELINE_Y, BASELINE_GROUPS)
    assert result.n_folds == 6
    assert result.splitter_name == "LeaveOneGroupOut"
    assert result.y_pred_baseline == pytest.approx(BASELINE_PRED, abs=1e-12)
    assert mae(result.y_true, result.y_pred_baseline) == pytest.approx(
        BASELINE_MAE, abs=1e-12
    )
    assert rmse(result.y_true, result.y_pred_baseline) ** 2 == pytest.approx(
        BASELINE_MSE, abs=1e-12
    )


def test_oof_arrays_restore_original_row_order():
    X = np.arange(12, dtype=float).reshape(12, 1)
    y = BASELINE_Y + 0.25
    result = run_validation(X, y, BASELINE_GROUPS)
    assert np.array_equal(result.y_true, y)
    assert result.y_pred_gpr.shape == (12,)
    assert result.fold_id.shape == (12,)


def test_scaler_mean_excludes_held_out_group():
    X = np.array([[0.0, 2.0], [10.0, 4.0], [100.0, 50.0]])
    y = np.array([0.0, 1.0, 2.0])
    groups = np.array(["a", "a", "b"])
    result = run_validation(X, y, groups)
    held_b = [fold for fold in result.folds if fold.group == "b"]
    assert len(held_b) == 1
    assert held_b[0].scaler_mean == pytest.approx((5.0, 3.0), abs=1e-12)


def test_gpr_returns_finite_mean_and_std():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(12, 2))
    y = rng.normal(size=12)
    groups = BASELINE_GROUPS
    result = run_validation(X, y, groups)
    assert np.isfinite(result.y_pred_gpr).all()
    assert np.isfinite(result.y_std_gpr).all()
    assert (result.y_std_gpr >= 0).all()


def test_fixed_seed_is_stable():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(12, 2))
    y = np.linspace(-1.0, 1.0, 12)
    a = run_validation(X, y, BASELINE_GROUPS, random_state=SEED)
    b = run_validation(X, y, BASELINE_GROUPS, random_state=SEED)
    assert a.y_pred_gpr == pytest.approx(b.y_pred_gpr, abs=1e-12)
    assert a.y_std_gpr == pytest.approx(b.y_std_gpr, abs=1e-12)


def test_single_group_falls_back_to_loocv_and_names_splitter():
    X = np.array([[0.0], [1.0], [2.0]])
    y = np.array([1.0, 3.0, 5.0])
    groups = np.array(["only", "only", "only"])
    result = run_validation(X, y, groups)
    assert result.splitter_name == "LeaveOneOut"
    assert result.n_folds == 3
    # Leave-one-out mean baseline: for held-out row i, mean of the other two.
    # Hand-calculated: hold 0 -> mean(3,5)=4; hold 1 -> mean(1,5)=3; hold 2 -> mean(1,3)=2
    assert result.y_pred_baseline == pytest.approx(np.array([4.0, 3.0, 2.0]), abs=1e-12)


def test_kernel_and_seed_constants():
    model = ScaledGPR(CONFIG.kernel, random_state=SEED)
    assert model.gpr.random_state == 20260811
    assert model.gpr.normalize_y is True
    assert model.gpr.n_restarts_optimizer == 10
    assert make_baseline_estimator().strategy == "mean"
    logo = LeaveOneGroupOut()
    assert logo.get_n_splits(groups=BASELINE_GROUPS) == 6


def test_scaled_gpr_uses_injected_kernel_settings():
    settings = KernelSettings(
        type=CONFIG.kernel.type,
        constant_value=2.0,
        length_scale=3.0,
        nu=1.5,
        noise_level=0.25,
        normalize_y=False,
        n_restarts_optimizer=0,
        random_state=17,
    )
    model = ScaledGPR(settings, random_state=17)
    params = model.gpr.kernel.get_params()
    assert params["k1__k1__constant_value"] == 2.0
    assert params["k1__k2__length_scale"] == 3.0
    assert params["k2__noise_level"] == 0.25
    assert model.gpr.normalize_y is False
    assert model.gpr.n_restarts_optimizer == 0

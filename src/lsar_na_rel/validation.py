"""Fixed Talk-v0 LOSO GPR and training-fold mean baseline."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.dummy import DummyRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.model_selection import LeaveOneGroupOut, LeaveOneOut
from sklearn.preprocessing import StandardScaler

from lsar_na_rel.config import KernelSettings


def mae(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    return float(np.mean(np.abs(y_true - y_pred)))


def rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def _kernel(settings: KernelSettings):
    return (
        ConstantKernel(settings.constant_value)
        * Matern(length_scale=settings.length_scale, nu=settings.nu)
        + WhiteKernel(noise_level=settings.noise_level)
    )


class ScaledGPR:
    """StandardScaler (train only) then isotropic Matern GPR with WhiteKernel."""

    def __init__(self, settings: KernelSettings, random_state: int | None = None) -> None:
        self.scaler = StandardScaler()
        seed = settings.random_state if random_state is None else int(random_state)
        self.gpr = GaussianProcessRegressor(
            kernel=_kernel(settings),
            normalize_y=settings.normalize_y,
            n_restarts_optimizer=settings.n_restarts_optimizer,
            random_state=seed,
        )

    def fit(self, X, y):
        Xs = self.scaler.fit_transform(np.asarray(X, dtype=float))
        self.gpr.fit(Xs, np.asarray(y, dtype=float))
        return self

    def predict(self, X, return_std: bool = False):
        Xs = self.scaler.transform(np.asarray(X, dtype=float))
        return self.gpr.predict(Xs, return_std=return_std)


def make_baseline_estimator() -> DummyRegressor:
    return DummyRegressor(strategy="mean")


@dataclass(frozen=True)
class FoldRecord:
    fold_id: int
    group: str
    n_train: int
    n_test: int
    scaler_mean: tuple[float, ...]
    kernel_: str


@dataclass(frozen=True)
class OofResult:
    y_true: np.ndarray
    y_pred_gpr: np.ndarray
    y_std_gpr: np.ndarray
    y_pred_baseline: np.ndarray
    fold_id: np.ndarray
    n_folds: int
    splitter_name: str
    folds: tuple[FoldRecord, ...]


def _as_numpy_xyg(X, y, groups):
    if isinstance(X, pd.DataFrame):
        X_np = X.to_numpy(dtype=float)
    else:
        X_np = np.asarray(X, dtype=float)
    y_np = np.asarray(y, dtype=float).reshape(-1)
    g_np = np.asarray(groups).astype(str).reshape(-1)
    if not (len(X_np) == len(y_np) == len(g_np)):
        raise ValueError("X, y, and groups must have the same length")
    return X_np, y_np, g_np


def _split_is_usable(splitter, X, y, groups) -> bool:
    try:
        splits = list(splitter.split(X, y, groups))
    except (ValueError, TypeError):
        return False
    if not splits:
        return False
    for train_idx, test_idx in splits:
        if len(train_idx) == 0 or len(test_idx) == 0:
            return False
    return True


def run_loso_gpr_and_baseline(
    X,
    y,
    groups,
    *,
    kernel_settings: KernelSettings,
    random_state: int | None = None,
) -> OofResult:
    X_np, y_np, g_np = _as_numpy_xyg(X, y, groups)
    n = len(y_np)
    logo = LeaveOneGroupOut()
    if _split_is_usable(logo, X_np, y_np, g_np):
        splitter = logo
        splitter_name = "LeaveOneGroupOut"
        split_kwargs = {"groups": g_np}
    else:
        splitter = LeaveOneOut()
        splitter_name = "LeaveOneOut"
        split_kwargs = {}

    y_pred_gpr = np.empty(n, dtype=float)
    y_std_gpr = np.empty(n, dtype=float)
    y_pred_baseline = np.empty(n, dtype=float)
    fold_id = np.empty(n, dtype=int)
    records: list[FoldRecord] = []

    for i, (train_idx, test_idx) in enumerate(splitter.split(X_np, y_np, **split_kwargs)):
        gpr = ScaledGPR(kernel_settings, random_state=random_state)
        gpr.fit(X_np[train_idx], y_np[train_idx])
        mean, std = gpr.predict(X_np[test_idx], return_std=True)
        y_pred_gpr[test_idx] = np.asarray(mean, dtype=float)
        y_std_gpr[test_idx] = np.asarray(std, dtype=float)

        baseline = make_baseline_estimator()
        baseline.fit(X_np[train_idx], y_np[train_idx])
        y_pred_baseline[test_idx] = baseline.predict(X_np[test_idx])
        fold_id[test_idx] = i

        held = g_np[test_idx]
        group_label = str(held[0]) if splitter_name == "LeaveOneGroupOut" else f"row{test_idx[0]}"
        scaler_mean = tuple(float(v) for v in np.asarray(gpr.scaler.mean_).tolist())
        records.append(
            FoldRecord(
                fold_id=i,
                group=group_label,
                n_train=int(len(train_idx)),
                n_test=int(len(test_idx)),
                scaler_mean=scaler_mean,
                kernel_=str(gpr.gpr.kernel_),
            )
        )

    return OofResult(
        y_true=y_np.copy(),
        y_pred_gpr=y_pred_gpr,
        y_std_gpr=y_std_gpr,
        y_pred_baseline=y_pred_baseline,
        fold_id=fold_id,
        n_folds=len(records),
        splitter_name=splitter_name,
        folds=tuple(records),
    )

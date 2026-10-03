"""Euclidean distance dependence estimators with explicit degeneracy statuses."""
from __future__ import annotations

import math

import numpy as np


def fmt(value: float | None) -> str:
    return "" if value is None else f"{value:.12g}"


def distance_scores(xs: np.ndarray, ys: np.ndarray) -> dict:
    """Conventional dCor and signed U-centered squared dCor (never sqrt U)."""
    xs, ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    if xs.ndim == 1:
        xs = xs[:, None]
    if ys.ndim == 1:
        ys = ys[:, None]
    if xs.ndim != 2 or ys.ndim != 2 or len(xs) != len(ys):
        raise ValueError("distance inputs must be paired two-dimensional matrices")
    if not np.isfinite(xs).all() or not np.isfinite(ys).all():
        raise ValueError("nonfinite distance inputs")
    n = len(xs)
    empty = {"distance_correlation": "", "u_centered_squared_distance_correlation": ""}
    if n < 2:
        return empty | {"status": "n_lt2", "u_status": "n_lt4"}
    if xs.shape[1] == 0 or ys.shape[1] == 0 or np.all(xs == xs[0]) or np.all(ys == ys[0]):
        return empty | {"status": "constant_side", "u_status": "n_lt4" if n < 4 else "constant_side"}
    a = np.linalg.norm(xs[:, None, :] - xs[None, :, :], axis=2)
    b = np.linalg.norm(ys[:, None, :] - ys[None, :, :], axis=2)
    def centered(d: np.ndarray) -> np.ndarray:
        return d - d.mean(axis=0)[None, :] - d.mean(axis=1)[:, None] + d.mean()
    aa, bb = centered(a), centered(b)
    den = math.sqrt(float(np.sum(aa * aa)) * float(np.sum(bb * bb)))
    dcor_squared = float(np.sum(aa * bb)) / den
    # Floating-point noise only: this population/sample V statistic is nonnegative.
    if dcor_squared < -1e-12 or dcor_squared > 1 + 1e-12:
        raise ValueError("conventional distance-correlation numeric failure")
    result = {"distance_correlation": fmt(math.sqrt(min(1.0, max(0.0, dcor_squared)))),
              "u_centered_squared_distance_correlation": "", "status": "estimated", "u_status": "n_lt4"}
    if n < 4:
        return result
    def u_centered(d: np.ndarray) -> np.ndarray:
        u = d - d.sum(axis=0)[None, :] / (n - 2) - d.sum(axis=1)[:, None] / (n - 2)
        u += d.sum() / ((n - 1) * (n - 2))
        np.fill_diagonal(u, 0.0)
        return u
    ua, ub = u_centered(a), u_centered(b)
    # The common n(n-3) inner-product denominator cancels in the normalized ratio.
    norm_ua, norm_ub = float(np.linalg.norm(ua)), float(np.linalg.norm(ub))
    # Singleton binary groups can have exactly zero U variance algebraically but
    # leave rounding residuals. Detect self-norm degeneracy relative to distances,
    # without thresholding the cross-product or clipping small signed U scores.
    zero_a = norm_ua <= 64 * np.finfo(float).eps * float(np.linalg.norm(a))
    zero_b = norm_ub <= 64 * np.finfo(float).eps * float(np.linalg.norm(b))
    unorm = norm_ua * norm_ub
    if zero_a or zero_b:
        result["u_status"] = "zero_u_distance_variance"
    else:
        result["u_centered_squared_distance_correlation"] = fmt(float(np.sum(ua * ub)) / unorm)
        result["u_status"] = "estimated"
    return result

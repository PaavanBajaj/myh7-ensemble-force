"""Deterministic, noninteractive figures for validated frozen summaries."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

if TYPE_CHECKING:
    from kcat_rel.modeling_summary import SummaryChecks


_WARM_BROWN = "#9B4F3B"
_NEUTRAL_GRAY = "#8A8782"
_DPI = 160
_MODEL_NAMES = {"ridge": "Ridge", "gpr": "GPR"}
_FEATURE_NAMES = {
    "delta_charge": "Charge change",
    "grantham_distance": "Grantham distance",
    "catalytic_motif_min_ca_distance_8act_mean_ab": "8ACT catalytic-motif\ndistance",
}


def _save(figure: plt.Figure, path: Path) -> None:
    """Save without date metadata, then close the figure before publication."""
    try:
        figure.savefig(
            path,
            dpi=_DPI,
            facecolor="white",
            metadata={"Software": "kcat_rel modeling summary"},
        )
    finally:
        plt.close(figure)


def _reserve_caption_space(figure: plt.Figure) -> None:
    """Keep the bottom figure caption outside the constrained axes rectangle."""
    figure.get_layout_engine().set(rect=(0.0, 0.15, 1.0, 0.83))


def render_feature_correlations(checks: SummaryChecks, path: Path) -> None:
    """Render registry-ordered full-cohort Pearson associations."""
    names = list(checks.feature_correlations)
    values = [checks.feature_correlations[name] for name in names]
    figure, axis = plt.subplots(figsize=(6.4, 4.5), layout="constrained", facecolor="white")
    _reserve_caption_space(figure)
    positions = np.arange(len(names))
    bars = axis.bar(positions, values, color=_WARM_BROWN, width=0.64)
    axis.axhline(0.0, color=_NEUTRAL_GRAY, linewidth=0.8)
    axis.set_xticks(positions, [_FEATURE_NAMES.get(name, name) for name in names])
    axis.set_ylabel("Pearson r with ln(kcat_rel)")
    axis.set_title("Descriptive full-cohort feature correlations")
    axis.set_ylim(-1.12, 1.12)
    for bar, value in zip(bars, values, strict=True):
        vertical = value + (0.045 if value >= 0 else -0.045)
        axis.text(bar.get_x() + bar.get_width() / 2, vertical, f"{value:.3f}", ha="center", va="bottom" if value >= 0 else "top", fontsize=9)
    figure.text(
        0.5,
        0.035,
        f"Full cohort (n = {checks.cohort_size}); descriptive only—not held-out performance or feature selection.",
        ha="center",
        color=_NEUTRAL_GRAY,
        fontsize=8,
    )
    _save(figure, path)


def render_study_oof(checks: SummaryChecks, path: Path) -> None:
    """Render source-study-held-out predictions with horizontal label uncertainty."""
    rows_by_model = {
        model: checks.study_oof.loc[checks.study_oof["model_kind"] == model].copy()
        for model in ("ridge", "gpr")
    }
    all_values = np.concatenate(
        [
            np.concatenate(
                [
                    rows["truth"].to_numpy(dtype=float)
                    - np.asarray([checks.label_uncertainties[str(variant)] for variant in rows["variant"]]),
                    rows["truth"].to_numpy(dtype=float)
                    + np.asarray([checks.label_uncertainties[str(variant)] for variant in rows["variant"]]),
                    rows["model_prediction"].to_numpy(dtype=float),
                    rows["baseline_prediction"].to_numpy(dtype=float),
                ]
            )
            for rows in rows_by_model.values()
        ]
    )
    minimum, maximum = float(all_values.min()), float(all_values.max())
    padding = max((maximum - minimum) * 0.06, 0.05)
    limits = (minimum - padding, maximum + padding)
    figure, axes = plt.subplots(1, 2, figsize=(9.5, 4.5), sharex=True, sharey=True, layout="constrained", facecolor="white")
    _reserve_caption_space(figure)
    for axis, model in zip(axes, ("ridge", "gpr"), strict=True):
        rows = rows_by_model[model]
        truth = rows["truth"].to_numpy(dtype=float)
        errors = np.asarray([checks.label_uncertainties[str(variant)] for variant in rows["variant"]], dtype=float)
        prediction = rows["model_prediction"].to_numpy(dtype=float)
        baseline = rows["baseline_prediction"].to_numpy(dtype=float)
        axis.errorbar(truth, prediction, xerr=errors, fmt="o", color=_WARM_BROWN, ecolor=_WARM_BROWN, elinewidth=0.8, capsize=2.0, markersize=4.8, label="Model prediction")
        axis.scatter(truth, baseline, facecolors="none", edgecolors=_NEUTRAL_GRAY, linewidths=1.0, s=35, label="Training-mean baseline")
        axis.plot(limits, limits, color=_NEUTRAL_GRAY, linewidth=0.9, linestyle="--", label="Identity")
        axis.set(xlim=limits, ylim=limits, aspect="equal")
        axis.set_title(f"{_MODEL_NAMES[model]} source-study OOF")
        axis.set_xlabel("Observed ln(kcat_rel)")
        axis.grid(alpha=0.18, linewidth=0.6)
    axes[0].set_ylabel("Predicted ln(kcat_rel)")
    axes[1].legend(loc="best", fontsize=7.5, frameon=False)
    figure.text(0.5, 0.035, "Horizontal error bars show reported observed-label uncertainty; n = 17 across 9 source-study folds.", ha="center", color=_NEUTRAL_GRAY, fontsize=8)
    _save(figure, path)


def render_error_comparison(checks: SummaryChecks, path: Path) -> None:
    """Render study-macro and variant MAE relative to the frozen baselines."""
    views = ("study", "variant")
    labels = ("Source-study macro MAE\n(9 folds)", "Leave-one-variant MAE\n(17 folds)")
    series = {
        "Training-mean baseline": [_baseline_value(checks, view) for view in views],
        "Ridge": [checks.model_metrics["ridge"][view]["mae"] for view in views],
        "GPR": [checks.model_metrics["gpr"][view]["mae"] for view in views],
    }
    colors = (_NEUTRAL_GRAY, _WARM_BROWN, "#C98D7E")
    figure, axis = plt.subplots(figsize=(7.5, 4.6), layout="constrained", facecolor="white")
    _reserve_caption_space(figure)
    positions = np.arange(len(views))
    width = 0.23
    maximum = max(value for values in series.values() for value in values)
    label_offset = max(maximum * 0.018, 0.012)
    for index, ((name, values), color) in enumerate(zip(series.items(), colors, strict=True)):
        bars = axis.bar(positions + (index - 1) * width, values, width=width, label=name, color=color)
        for bar, value in zip(bars, values, strict=True):
            axis.text(bar.get_x() + bar.get_width() / 2, value + label_offset, f"{value:.3f}", ha="center", va="bottom", fontsize=8, rotation=90)
    axis.set_xticks(positions, labels)
    axis.set_ylabel("Mean absolute error, ln(kcat_rel)")
    axis.set_title("Held-out error versus training-mean baseline")
    axis.set_ylim(0.0, maximum + 4 * label_offset)
    axis.legend(frameon=False, fontsize=8)
    axis.grid(axis="y", alpha=0.18, linewidth=0.6)
    figure.text(0.5, 0.035, "Smaller is better. Neither model passed the predeclared gate.", ha="center", color=_NEUTRAL_GRAY, fontsize=8)
    _save(figure, path)


def _baseline_value(checks: SummaryChecks, view: str) -> float:
    """Both models share an audited baseline; reject a divergent rendering input."""
    ridge = checks.baseline_metrics["ridge"][view]["mae"]
    gpr = checks.baseline_metrics["gpr"][view]["mae"]
    if not np.isclose(ridge, gpr, rtol=0.0, atol=1e-12):
        raise ValueError(f"{view} baseline MAE differs between frozen models")
    return float(ridge)

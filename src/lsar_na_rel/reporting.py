"""Honesty-packet reporting for Talk-v0 OOF runs (no SUCCESS marker)."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lsar_na_rel.validation import OofResult, mae, rmse

_FLOAT_TOKEN = "0.33333333333333331"

_NEXT_STEPS_BODY = """# Next steps (Talk-v0 diagnostic ladder)

This run is an exploratory small-n private diagnostic. Risky allowed rows
R403Q, R663H, and G768R are included. Do not treat pooled OOF errors as a
clinical or mechanistic claim.

No FoldX-derived values are present in this release. The modeling family is
chemistry plus checksummed 8ACT features from the committed features.csv
snapshot. The IHM interface flag remains 2/16 positives
(D382Y, R403Q) at the frozen 4.5 Å cutoff; do not retune geometry from these
scores.

Re-derive metrics from metrics.json and oof-predictions.csv in this directory.
Ship the other pre-registered ladder run_ids on the same rows and folds.
Do not add a fourth feature combination or drop rows after seeing scores.
"""

def dump_json_full_floats(payload: object, path: str | Path) -> None:
    def convert(node):
        if isinstance(node, float):
            if not math.isfinite(node):
                raise ValueError(f"non-finite float in metrics: {node}")
            return f"__FLOAT__{format(node, '.17g')}__"
        if isinstance(node, dict):
            return {str(k): convert(v) for k, v in node.items()}
        if isinstance(node, (list, tuple)):
            return [convert(v) for v in node]
        return node

    raw = json.dumps(convert(payload), indent=2)
    text = re.sub(r'"__FLOAT__([^"]+)__"', r"\1", raw)
    Path(path).write_text(text + "\n", encoding="utf-8")


def write_oof_csv(
    path: str | Path,
    variants: pd.Series,
    source_ids: pd.Series,
    oof: OofResult,
) -> None:
    frame = pd.DataFrame(
        {
            "variant": variants.to_numpy(),
            "source_id": source_ids.to_numpy(),
            "y_true": np.asarray(oof.y_true, dtype=float),
            "y_pred_gpr": np.asarray(oof.y_pred_gpr, dtype=float),
            "y_std_gpr": np.asarray(oof.y_std_gpr, dtype=float),
            "y_pred_baseline": np.asarray(oof.y_pred_baseline, dtype=float),
            "fold_id": np.asarray(oof.fold_id),
        }
    )
    frame.to_csv(path, index=False)


def build_metrics_payload(
    run_id: str,
    feature_columns: Sequence[str],
    seed: int,
    oof: OofResult,
) -> dict:
    gpr_mae = mae(oof.y_true, oof.y_pred_gpr)
    gpr_rmse = rmse(oof.y_true, oof.y_pred_gpr)
    baseline_mae = mae(oof.y_true, oof.y_pred_baseline)
    baseline_rmse = rmse(oof.y_true, oof.y_pred_baseline)
    return {
        "run_id": run_id,
        "feature_columns": list(feature_columns),
        "seed": seed,
        "n_rows": int(len(oof.y_true)),
        "n_folds": int(oof.n_folds),
        "splitter_name": oof.splitter_name,
        "gpr_mae": gpr_mae,
        "gpr_rmse": gpr_rmse,
        "baseline_mae": baseline_mae,
        "baseline_rmse": baseline_rmse,
        "delta_mae_gpr_minus_baseline": gpr_mae - baseline_mae,
        "delta_rmse_gpr_minus_baseline": gpr_rmse - baseline_rmse,
    }


def write_metrics_json(path: str | Path, payload: object) -> None:
    dump_json_full_floats(payload, path)


def write_observed_vs_predicted_png(
    path: str | Path,
    oof: OofResult,
    run_id: str,
) -> None:
    y_true = np.asarray(oof.y_true, dtype=float)
    y_gpr = np.asarray(oof.y_pred_gpr, dtype=float)
    y_base = np.asarray(oof.y_pred_baseline, dtype=float)

    fig, ax = plt.subplots()
    ax.scatter(y_true, y_gpr, label="GPR")
    ax.scatter(y_true, y_base, label="training-fold mean baseline")
    lo = float(np.min([y_true.min(), y_gpr.min(), y_base.min()]))
    hi = float(np.max([y_true.max(), y_gpr.max(), y_base.max()]))
    ax.plot([lo, hi], [lo, hi], linestyle="--", color="gray", label="y=x")
    ax.set_xlabel("y_true")
    ax.set_ylabel("y_pred")
    ax.set_title(f"Observed vs predicted (OOF) — {run_id}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, format="png")
    plt.close(fig)


def write_next_steps(path: str | Path, *, run_id: str) -> None:
    Path(path).write_text(_NEXT_STEPS_BODY, encoding="utf-8")


def write_reports(
    run_dir: str | Path,
    *,
    run_id: str,
    feature_columns: Sequence[str],
    seed: int,
    variants: pd.Series,
    source_ids: pd.Series,
    oof: OofResult,
) -> None:
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    write_oof_csv(run_dir / "oof-predictions.csv", variants, source_ids, oof)
    payload = build_metrics_payload(run_id, feature_columns, seed, oof)
    write_metrics_json(run_dir / "metrics.json", payload)
    write_observed_vs_predicted_png(
        run_dir / "observed-vs-predicted.png", oof, run_id=run_id
    )
    write_next_steps(run_dir / "next-steps.md", run_id=run_id)

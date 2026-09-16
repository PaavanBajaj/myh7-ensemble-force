"""Reporting artifacts from synthetic OOF (not S1 Na_rel)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lsar_na_rel.reporting import (
    build_metrics_payload,
    dump_json_full_floats,
    write_next_steps,
    write_observed_vs_predicted_png,
    write_oof_csv,
    write_reports,
)
from lsar_na_rel.validation import OofResult

SYNTHETIC = OofResult(
    y_true=np.array([0.0, 2.0, 4.0]),
    # Abs errors 1,1,1 → MAE/MSE/RMSE = 1.0 (brief listed [1,2,3] but hand-calc is 1.0)
    y_pred_gpr=np.array([1.0, 1.0, 3.0]),
    y_std_gpr=np.array([0.1, 0.1, 0.1]),
    y_pred_baseline=np.array([2.0, 2.0, 2.0]),
    fold_id=np.array([0, 0, 1]),
    n_folds=2,
    splitter_name="LeaveOneGroupOut",
    folds=(),
)


def test_oof_csv_schema_and_primary_order(tmp_path: Path):
    path = tmp_path / "oof-predictions.csv"
    write_oof_csv(
        path,
        variants=pd.Series(["A1C", "D2E", "F3G"]),
        source_ids=pd.Series(["s1", "s1", "s2"]),
        oof=SYNTHETIC,
    )
    frame = pd.read_csv(path)
    assert list(frame.columns) == [
        "variant",
        "source_id",
        "y_true",
        "y_pred_gpr",
        "y_std_gpr",
        "y_pred_baseline",
        "fold_id",
    ]
    assert frame["variant"].tolist() == ["A1C", "D2E", "F3G"]


def test_metrics_match_hand_calculated_synthetic_oof():
    payload = build_metrics_payload(
        run_id="ladder-1-rsa",
        feature_columns=("wt_site_rsa_8act_mean_ab",),
        seed=20260811,
        oof=SYNTHETIC,
    )
    assert payload["n_rows"] == 3
    assert payload["n_folds"] == 2
    assert payload["gpr_mae"] == 1.0
    assert payload["gpr_rmse"] == 1.0
    assert payload["baseline_mae"] == pytest.approx(4.0 / 3.0, abs=1e-12)
    assert payload["baseline_rmse"] ** 2 == pytest.approx(8.0 / 3.0, abs=1e-12)
    assert payload["delta_mae_gpr_minus_baseline"] == pytest.approx(
        1.0 - (4.0 / 3.0), abs=1e-12
    )


def test_dump_json_full_floats_writes_seventeen_sig_digits(tmp_path: Path):
    path = tmp_path / "metrics.json"
    dump_json_full_floats({"gpr_mae": 1.0 / 3.0}, path)
    text = path.read_text(encoding="utf-8")
    assert "0.33333333333333331" in text or "0.3333333333333333" in text
    loaded = json.loads(text)
    assert loaded["gpr_mae"] == pytest.approx(1.0 / 3.0, abs=1e-16)


def test_plot_is_nonempty_png(tmp_path: Path):
    path = tmp_path / "observed-vs-predicted.png"
    write_observed_vs_predicted_png(path, SYNTHETIC, run_id="ladder-1-rsa")
    data = path.read_bytes()
    assert data.startswith(b"\x89PNG")
    assert path.stat().st_size > 0


def test_next_steps_mentions_required_limitations_without_ranking(tmp_path: Path):
    path = tmp_path / "next-steps.md"
    write_next_steps(path, run_id="ladder-3-nofoldx")
    text = path.read_text(encoding="utf-8").lower()
    words = text.split()
    assert len(words) <= 600
    assert "foldx" in text
    assert "no foldx-derived values" in text
    assert "operator" not in text
    assert "d382y" in text and "r403q" in text
    assert "2/16" in text or "2 of 16" in text
    assert "risky" in text
    assert "best" not in text
    assert "winner" not in text
    assert "significant" not in text


def test_write_reports_does_not_create_success(tmp_path: Path):
    write_reports(
        tmp_path,
        run_id="ladder-1-rsa",
        feature_columns=("wt_site_rsa_8act_mean_ab",),
        seed=20260811,
        variants=pd.Series(["A1C", "D2E", "F3G"]),
        source_ids=pd.Series(["s1", "s1", "s2"]),
        oof=SYNTHETIC,
    )
    assert (tmp_path / "oof-predictions.csv").is_file()
    assert (tmp_path / "metrics.json").is_file()
    assert (tmp_path / "observed-vs-predicted.png").is_file()
    assert (tmp_path / "next-steps.md").is_file()
    assert not (tmp_path / "SUCCESS").exists()

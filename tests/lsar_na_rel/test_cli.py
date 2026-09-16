"""CLI wiring tests. Do not assert S1 Na_rel MAE numbers."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lsar_na_rel.cli import main, select_run_ids, set_blas_threads_one, write_success_marker
from lsar_na_rel.config import load_talk_v0_config
from lsar_na_rel.validation import FoldRecord, OofResult

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_ROOT = REPO_ROOT / "workflows" / "lsar-na-rel"
CONFIG = load_talk_v0_config(WORKFLOW_ROOT / "config" / "lsar-na-rel.json")


def test_select_run_ids_one_or_all():
    assert select_run_ids("ladder-1-rsa", False) == ("ladder-1-rsa",)
    assert select_run_ids(None, True) == (
        "ladder-1-rsa",
        "ladder-2-rsa-ihm",
        "ladder-3-nofoldx",
    )
    with pytest.raises(ValueError, match="exactly one"):
        select_run_ids("ladder-1-rsa", True)
    with pytest.raises(ValueError, match="exactly one"):
        select_run_ids(None, False)
    with pytest.raises(ValueError, match="registered run_id"):
        select_run_ids("probe-ihm-min-dist", False)


def test_set_blas_threads_one(monkeypatch):
    for key in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        monkeypatch.delenv(key, raising=False)
    set_blas_threads_one()
    assert os.environ["OMP_NUM_THREADS"] == "1"
    assert os.environ["MKL_NUM_THREADS"] == "1"
    assert os.environ["OPENBLAS_NUM_THREADS"] == "1"
    assert os.environ["VECLIB_MAXIMUM_THREADS"] == "1"
    assert os.environ["NUMEXPR_NUM_THREADS"] == "1"


def test_success_marker_requires_other_artifacts(tmp_path: Path):
    with pytest.raises(FileNotFoundError):
        write_success_marker(tmp_path, CONFIG.required_artifacts)
    for name in (
        "run-manifest.json",
        "metrics.json",
        "oof-predictions.csv",
        "observed-vs-predicted.png",
        "run.log",
        "next-steps.md",
    ):
        (tmp_path / name).write_text("x", encoding="utf-8")
    write_success_marker(tmp_path, CONFIG.required_artifacts)
    assert (tmp_path / "SUCCESS").read_text(encoding="utf-8").strip() == "SUCCESS"


def test_cli_one_run_writes_success_last(tmp_path, monkeypatch):
    oof = OofResult(
        y_true=np.array([0.0, 1.0]),
        y_pred_gpr=np.array([0.2, 0.8]),
        y_std_gpr=np.array([0.1, 0.1]),
        y_pred_baseline=np.array([0.5, 0.5]),
        fold_id=np.array([0, 1]),
        n_folds=2,
        splitter_name="LeaveOneGroupOut",
        folds=(
            FoldRecord(0, "g0", 1, 1, (0.0,), "k"),
            FoldRecord(1, "g1", 1, 1, (1.0,), "k"),
        ),
    )

    class FakeMatrix:
        run_id = "ladder-1-rsa"
        feature_columns = ("wt_site_rsa_8act_mean_ab",)
        X = pd.DataFrame({"wt_site_rsa_8act_mean_ab": [0.1, 0.2]})
        y = pd.Series([0.0, 1.0], name="y")
        groups = pd.Series(["g0", "g1"], name="source_id")
        variants = pd.Series(["A1C", "D2E"], name="variant")
        features_path = tmp_path / "features.csv"

    monkeypatch.setattr("lsar_na_rel.cli.build_talk_v0_matrix", lambda **kwargs: FakeMatrix())
    monkeypatch.setattr("lsar_na_rel.cli.run_loso_gpr_and_baseline", lambda *a, **k: oof)
    monkeypatch.setattr("lsar_na_rel.cli.utc_timestamp", lambda: "20260816T000000Z")

    features = tmp_path / "features.csv"
    primary = tmp_path / "primary.csv"
    labels = tmp_path / "labels.csv"
    features.write_text("variant,wt_site_rsa_8act_mean_ab\nA1C,0.1\nD2E,0.2\n", encoding="utf-8")
    primary.write_text("variant,LSAR_mut,LSAR_WT,Na_rel,ln(Na_rel)\n", encoding="utf-8")
    labels.write_text("variant,source_id\n", encoding="utf-8")
    out_root = tmp_path / "results"
    argv = [
        "talk-v0",
        "--run-id",
        "ladder-1-rsa",
        "--primary",
        str(primary),
        "--labels",
        str(labels),
        "--features",
        str(features),
        "--output-root",
        str(out_root),
        "--config",
        str(WORKFLOW_ROOT / "config" / "lsar-na-rel.json"),
    ]
    assert main(argv) == 0
    run_dir = out_root / "ladder-1-rsa-20260816T000000Z"
    assert (run_dir / "SUCCESS").is_file()
    names = [
        "run-manifest.json",
        "metrics.json",
        "oof-predictions.csv",
        "observed-vs-predicted.png",
        "run.log",
        "next-steps.md",
        "SUCCESS",
    ]
    mtimes = [(run_dir / name).stat().st_mtime_ns for name in names]
    assert mtimes[-1] >= max(mtimes[:-1])
    manifest = (run_dir / "run-manifest.json").read_text(encoding="utf-8")
    assert "unavailable_no_values" in manifest
    assert "/Users/" not in (run_dir / "run.log").read_text(encoding="utf-8")


def test_cli_failure_does_not_write_success(tmp_path, monkeypatch):
    def boom(**kwargs):
        raise RuntimeError("injected failure")

    monkeypatch.setattr("lsar_na_rel.cli.build_talk_v0_matrix", boom)
    monkeypatch.setattr("lsar_na_rel.cli.utc_timestamp", lambda: "20260816T000001Z")
    out_root = tmp_path / "results"
    dummy = tmp_path / "x.csv"
    dummy.write_text("x\n", encoding="utf-8")
    argv = [
        "talk-v0",
        "--run-id",
        "ladder-1-rsa",
        "--primary",
        str(dummy),
        "--labels",
        str(dummy),
        "--features",
        str(dummy),
        "--output-root",
        str(out_root),
        "--config",
        str(WORKFLOW_ROOT / "config" / "lsar-na-rel.json"),
    ]
    assert main(argv) == 1
    run_dir = out_root / "ladder-1-rsa-20260816T000001Z"
    assert not (run_dir / "SUCCESS").exists()

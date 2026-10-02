from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kcat_rel.reporting import required_artifacts, write_diagnostic_artifacts
from kcat_rel.registry import FROZEN_PRIMARY_VARIANTS


def _complete_result() -> dict:
    """Load the frozen completed run as a realistic cross-artifact fixture."""
    root = Path(__file__).resolve().parents[2] / "workflows/kcat-rel/results" / "diagnostic-v0"
    return {
        "run_manifest": json.loads((root / "run-manifest.json").read_text()),
        "model_table": pd.read_csv(root / "model-table.csv"),
        "tuning_records": pd.read_csv(root / "tuning-records.csv"),
        "study_oof": pd.read_csv(root / "study-oof.csv"),
        "variant_oof": pd.read_csv(root / "variant-oof.csv"),
        "inner_fold_audit": json.loads((root / "inner-fold-audit.json").read_text()),
        "metrics": json.loads((root / "metrics.json").read_text()),
        "sensitivity": json.loads((root / "sensitivity.json").read_text()),
        "audit_markdown": (root / "AUDIT.md").read_text(),
        "run_log": (root / "run.log").read_text(),
    }


def _obsolete_complete_result() -> dict:
    """Retained only to make this fixture replacement explicit in test history."""
    variants = list(FROZEN_PRIMARY_VARIANTS)
    model_table = pd.DataFrame({"variant": variants, "y": [0.0] * 17, "y_error": [0.1] * 17,
                                "source_id": ["study"] * 17, "delta_charge": [0.0] * 17,
                                "grantham_distance": [1.0] * 17,
                                "catalytic_motif_min_ca_distance_8act_mean_ab": [2.0] * 17})
    def oof(view: str) -> pd.DataFrame:
        folds = (list(range(9)) + [8] * 8) if view == "study" else list(range(17))
        return pd.DataFrame({"model_kind": [kind for kind in ("ridge", "gpr") for _ in variants],
                             "variant": variants * 2, "source_id": ["study"] * 34, "truth": [0.0] * 34,
                             "model_prediction": [0.0] * 34, "baseline_prediction": [0.0] * 34,
                             "outer_fold": folds * 2, "selected_inner_hyperparameters": ["{}"] * 34,
                             "scaler_mean": ["[]"] * 34, "scaler_scale": ["[]"] * 34})
    metric = {"mae": 0.0, "rmse": 0.0}
    decision = {"pass_gate": False, "failed_rules": ["study_macro_mae"], "mae_improvement": {"study": 0.0, "variant": 0.0}, "rmse_degradation": {"study": 0.0, "variant": 0.0}, "model_views": {"study": metric, "variant": metric}, "baseline_views": {"study": metric, "variant": metric}}
    tuning = pd.DataFrame({"model_kind": ["ridge"] * 26 + ["gpr"] * 26, "view": ["study"] * 9 + ["variant"] * 17 + ["study"] * 9 + ["variant"] * 17, "outer_fold": list(range(9)) + list(range(17)) + list(range(9)) + list(range(17)), "selected_params": ["{}"] * 52, "inner_scores": ["{}"] * 52, "scaler_mean": ["[]"] * 52, "scaler_scale": ["[]"] * 52, "baseline_mean": [0.0] * 52})
    return {"run_manifest": {"schema_version": 1, "active_profile": "kcat-rel-v0-no-msa", "fallback_reason": "x", "variants": variants, "seed": 1, "input_sha256": {"a": "a" * 64}, "tool_versions": {"python": "x"}, "command": ["diagnostic"], "git_revision": None, "timestamp_utc": "2026-01-01T00:00:00+00:00"}, "model_table": model_table, "tuning_records": tuning, "study_oof": oof("study"), "variant_oof": oof("variant"), "metrics": {"winner": None, "models": {"ridge": decision, "gpr": decision}}, "sensitivity": {"ridge": {"weighted": {}, "perturbation": {"draws": 1000, "seed": 1, "gate_pass_rate": 0.0}}, "gpr": {"weighted": {}, "perturbation": {"draws": 1000, "seed": 1, "gate_pass_rate": 0.0}}}, "audit_markdown": "# Audit\n", "run_log": "complete\n"}


def test_success_marker_is_written_last_and_only_after_complete_artifacts(tmp_path: Path) -> None:
    """A completed diagnostic must be auditable before it is marked successful."""
    result = _complete_result()

    write_diagnostic_artifacts(tmp_path, result)

    assert (tmp_path / "SUCCESS").exists()
    assert required_artifacts(tmp_path).issubset(set(tmp_path.iterdir()))
    assert (tmp_path / "SUCCESS").stat().st_mtime_ns >= max(
        path.stat().st_mtime_ns for path in required_artifacts(tmp_path) - {tmp_path / "SUCCESS"}
    )
    manifest = json.loads((tmp_path / "run-manifest.json").read_text())
    assert all(
        hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() == digest
        for name, digest in manifest["artifact_sha256"].items()
    )


def test_reporting_failure_leaves_no_partial_output_or_success_marker(tmp_path: Path) -> None:
    """An invalid incomplete result must not create an apparently completed run."""
    result = {"run_manifest": {"schema_version": 1}}

    try:
        write_diagnostic_artifacts(tmp_path, result)
    except ValueError:
        pass
    else:  # pragma: no cover - guards the intentional failure assertion
        raise AssertionError("incomplete result unexpectedly wrote artifacts")

    assert not (tmp_path / "SUCCESS").exists()
    assert not any(tmp_path.iterdir())


@pytest.mark.parametrize("field", ["model_table", "study_oof", "metrics", "sensitivity"])
def test_reporting_rejects_empty_or_incomplete_required_artifacts(tmp_path: Path, field: str) -> None:
    """Empty tables and missing decision/sensitivity data cannot receive SUCCESS."""
    result = _complete_result()
    result[field] = pd.DataFrame() if field in {"model_table", "study_oof"} else {}
    with pytest.raises(ValueError):
        write_diagnostic_artifacts(tmp_path, result)
    assert not any(tmp_path.iterdir())


def test_reporting_rejects_oof_truth_or_source_drift_from_model_table(tmp_path: Path) -> None:
    """OOF values cannot be marked complete if they are not the frozen model inputs."""
    result = _complete_result()
    result["study_oof"].loc[0, "truth"] += 0.01
    with pytest.raises(ValueError, match="truth|source|variant"):
        write_diagnostic_artifacts(tmp_path, result)
    assert not any(tmp_path.iterdir())


@pytest.mark.parametrize("mutation", ["metrics", "winner", "tuning", "scaler", "fold", "params", "sensitivity", "weights", "perturbation", "hashes"])
def test_reporting_rejects_cross_artifact_inconsistency(tmp_path: Path, mutation: str) -> None:
    """Changing a gate, fold audit, sensitivity detail, or hash must block SUCCESS."""
    result = _complete_result()
    if mutation == "metrics":
        result["metrics"]["models"]["ridge"]["model_views"]["study"]["mae"] += 0.01
    elif mutation == "winner":
        result["metrics"]["winner"] = "ridge"
    elif mutation == "tuning":
        result["tuning_records"].loc[0, "selected_params"] = "{}"
    elif mutation == "scaler":
        result["tuning_records"].loc[0, "scaler_scale"] = "[]"
    elif mutation == "fold":
        result["tuning_records"].loc[0, "train_indices"] = result["tuning_records"].loc[0, "test_indices"]
    elif mutation == "params":
        result["tuning_records"].loc[0, "selected_params"] = '{"length_scale":99,"noise_variance":0.04}'
    elif mutation == "sensitivity":
        result["sensitivity"]["ridge"]["weighted"]["study"]["fold_sample_weights"] = []
    elif mutation == "weights":
        result["sensitivity"]["ridge"]["weighted"]["study"]["fold_sample_weights"][0][0] += 1
    elif mutation == "perturbation":
        result["sensitivity"]["gpr"]["perturbation"]["seed"] = True
    else:
        result["run_manifest"]["input_sha256"].pop("registry")
    with pytest.raises(ValueError):
        write_diagnostic_artifacts(tmp_path, result)
    assert not any(tmp_path.iterdir())

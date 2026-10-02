from __future__ import annotations

import copy
import json
from pathlib import Path

import pandas as pd
import pytest

from kcat_rel.reporting import reconstruct_inner_fold_audit, write_diagnostic_artifacts


ROOT = Path(__file__).resolve().parents[2] / "workflows/kcat-rel/results" / "diagnostic-v0"


def _result() -> dict:
    return {
        "run_manifest": json.loads((ROOT / "run-manifest.json").read_text()),
        "model_table": pd.read_csv(ROOT / "model-table.csv"),
        "tuning_records": pd.read_csv(ROOT / "tuning-records.csv"),
        "study_oof": pd.read_csv(ROOT / "study-oof.csv"),
        "variant_oof": pd.read_csv(ROOT / "variant-oof.csv"),
        "inner_fold_audit": json.loads((ROOT / "inner-fold-audit.json").read_text()),
        "metrics": json.loads((ROOT / "metrics.json").read_text()),
        "sensitivity": json.loads((ROOT / "sensitivity.json").read_text()),
        "audit_markdown": (ROOT / "AUDIT.md").read_text(),
        "run_log": (ROOT / "run.log").read_text(),
    }


def _rejects(tmp_path: Path, result: dict) -> None:
    with pytest.raises(ValueError):
        write_diagnostic_artifacts(tmp_path, result)
    assert not any(tmp_path.iterdir())


def test_reconstructed_inner_audit_is_exact_and_training_only() -> None:
    result = _result()
    assert result["inner_fold_audit"] == reconstruct_inner_fold_audit(
        result["model_table"], result["tuning_records"]
    )


@pytest.mark.parametrize(
    "mutation",
    ["variant_identity", "study_identity", "oof_fold", "oof_fold_bool", "tuning_fold_bool", "outer_train_order", "outer_test_order"],
)
def test_reporting_rejects_outer_fold_identity_mutations(tmp_path: Path, mutation: str) -> None:
    result = _result()
    if mutation == "variant_identity":
        result["variant_oof"].loc[0, ["variant", "source_id", "truth"]] = result["variant_oof"].loc[1, ["variant", "source_id", "truth"]].to_list()
    elif mutation == "study_identity":
        result["study_oof"].loc[0, "outer_fold"] = 1
    elif mutation == "oof_fold":
        result["variant_oof"].loc[0, "outer_fold"] = 1
    elif mutation == "oof_fold_bool":
        result["variant_oof"]["outer_fold"] = result["variant_oof"]["outer_fold"].astype(object)
        result["variant_oof"].loc[0, "outer_fold"] = False
    elif mutation == "tuning_fold_bool":
        result["tuning_records"]["outer_fold"] = result["tuning_records"]["outer_fold"].astype(object)
        result["tuning_records"].loc[0, "outer_fold"] = False
    elif mutation == "outer_train_order":
        record = json.loads(result["tuning_records"].loc[0, "train_indices"])
        record[:2] = reversed(record[:2])
        result["tuning_records"].loc[0, "train_indices"] = json.dumps(record)
    else:
        row = result["tuning_records"].query("view == 'study'").index[0]
        record = json.loads(result["tuning_records"].loc[row, "test_indices"])
        record.reverse()
        result["tuning_records"].loc[row, "test_indices"] = json.dumps(record)
    _rejects(tmp_path, result)


@pytest.mark.parametrize(
    "mutation",
    ["ridge_alias", "ridge_bool", "ridge_extra", "ridge_missing", "gpr_alias", "gpr_off_grid", "gpr_missing", "gpr_extra", "negative_score", "string_score", "missing_score", "extra_score", "duplicate_score", "duplicate_param"],
)
def test_reporting_rejects_noncanonical_grid_encodings(tmp_path: Path, mutation: str) -> None:
    result = _result()
    table = result["tuning_records"]
    row = table.index[table.model_kind == ("ridge" if mutation.startswith("ridge") else "gpr")][0]
    params = json.loads(table.loc[row, "selected_params"])
    scores = json.loads(table.loc[row, "inner_scores"])
    if mutation == "ridge_alias": params = {"alpha": 1.0, "a": 1.0}
    elif mutation == "ridge_bool": params = {"alpha": True}
    elif mutation == "ridge_extra": params["noise"] = 0.1
    elif mutation == "ridge_missing": params = {}
    elif mutation == "gpr_alias": params = {"length": params["length_scale"], "noise_variance": params["noise_variance"]}
    elif mutation == "gpr_off_grid": params["length_scale"] = 0.75
    elif mutation == "gpr_missing": params.pop("noise_variance")
    elif mutation == "gpr_extra": params["amplitude"] = 1.0
    elif mutation == "negative_score": scores[next(iter(scores))] = -0.1
    elif mutation == "string_score": scores[next(iter(scores))] = "0.1"
    elif mutation == "missing_score": scores.pop(next(iter(scores)))
    elif mutation == "extra_score": scores["alias=1"] = 0.0
    elif mutation == "duplicate_score":
        encoded = json.dumps(scores)
        key = next(iter(scores))
        table.loc[row, "inner_scores"] = encoded[:-1] + f',"{key}":{scores[key]}' + "}"
    else:
        encoded = json.dumps(params)
        key = next(iter(params))
        table.loc[row, "selected_params"] = encoded[:-1] + f',"{key}":{params[key]}' + "}"
    if mutation != "duplicate_param": table.loc[row, "selected_params"] = json.dumps(params)
    if mutation != "duplicate_score": table.loc[row, "inner_scores"] = json.dumps(scores)
    _rejects(tmp_path, result)


@pytest.mark.parametrize("mutation", ["provenance", "revision", "missing_record", "inner_test", "inner_train", "inner_scaler", "inner_model"])
def test_reporting_rejects_inner_fold_audit_mutations(tmp_path: Path, mutation: str) -> None:
    result = _result()
    audit = result["inner_fold_audit"]
    if mutation == "provenance": audit["provenance"]["method"] = "recorded_at_execution"
    elif mutation == "revision": audit["provenance"]["original_execution_revision"] = "ba9b69b"
    elif mutation == "missing_record": audit["records"].pop()
    elif mutation == "inner_test": audit["records"][0]["test_indices"] = audit["records"][1]["test_indices"]
    elif mutation == "inner_train": audit["records"][0]["train_indices"].reverse()
    elif mutation == "inner_scaler": audit["records"][0]["scaler_mean"][0] += 1
    else: audit["records"][0]["model_kind"] = "ridge" if audit["records"][0]["model_kind"] == "gpr" else "gpr"
    _rejects(tmp_path, result)


@pytest.mark.parametrize(
    "mutation",
    ["top_model", "view", "weighted_model", "prediction_column", "prediction_bool", "baseline_value", "weights_nested", "weights_bool", "weights_empty_false", "variances_negative", "variances_empty_none"],
)
def test_reporting_rejects_weighted_sensitivity_mutations(tmp_path: Path, mutation: str) -> None:
    result = _result()
    weighted = result["sensitivity"]["ridge"]["weighted"]["variant"]
    if mutation == "top_model": result["sensitivity"]["ridge"], result["sensitivity"]["gpr"] = result["sensitivity"]["gpr"], result["sensitivity"]["ridge"]
    elif mutation == "view": weighted["view"] = "study"
    elif mutation == "weighted_model": weighted["model_kind"] = "gpr"
    elif mutation == "prediction_column": weighted["prediction"] = [[value] for value in weighted["prediction"]]
    elif mutation == "prediction_bool": weighted["prediction"][0] = True
    elif mutation == "baseline_value": weighted["baseline_prediction"][0] += 0.01
    elif mutation == "weights_nested": weighted["fold_sample_weights"][0] = [[value] for value in weighted["fold_sample_weights"][0]]
    elif mutation == "weights_bool": weighted["fold_sample_weights"][0][0] = True
    elif mutation == "weights_empty_false": result["sensitivity"]["gpr"]["weighted"]["variant"]["fold_sample_weights"][0] = False
    elif mutation == "variances_negative": result["sensitivity"]["gpr"]["weighted"]["variant"]["fold_reported_variances"][0][0] = -1.0
    else: weighted["fold_reported_variances"][0] = None
    _rejects(tmp_path, result)


@pytest.mark.parametrize(
    "mutation",
    ["seed_float", "draws_bool", "missing_fold_params", "wrong_fold_params", "quantile_string", "quantile_bool", "quantile_negative", "quantile_order", "count_float", "rate_string", "rate_mismatch"],
)
def test_reporting_rejects_perturbation_mutations(tmp_path: Path, mutation: str) -> None:
    result = _result()
    perturbation = result["sensitivity"]["ridge"]["perturbation"]
    summary = perturbation["metric_quantiles"]["study"]["mae"]
    if mutation == "seed_float": perturbation["seed"] = 20260916.0
    elif mutation == "draws_bool": perturbation["draws"] = True
    elif mutation == "missing_fold_params": perturbation["fold_hyperparameters"].pop("variant")
    elif mutation == "wrong_fold_params": perturbation["fold_hyperparameters"]["variant"][0] = {"alpha": 99.0}
    elif mutation == "quantile_string": summary["q025"] = "0.1"
    elif mutation == "quantile_bool": summary["q025"] = False
    elif mutation == "quantile_negative": summary["q025"] = -0.1
    elif mutation == "quantile_order": summary["q025"] = summary["q975"] + 1
    elif mutation == "count_float": perturbation["gate_pass_count"] = 4.0
    elif mutation == "rate_string": perturbation["gate_pass_rate"] = "0.004"
    else: perturbation["gate_pass_rate"] += 0.001
    _rejects(tmp_path, result)

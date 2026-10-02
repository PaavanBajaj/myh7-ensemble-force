"""Reject Python equality aliases before cross-artifact comparisons."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from kcat_rel.reporting import write_diagnostic_artifacts
from .test_reporting_producers import _producer_result


def _rejects(tmp_path: Path, result: dict) -> None:
    with pytest.raises(ValueError):
        write_diagnostic_artifacts(tmp_path, result)
    assert not any(tmp_path.iterdir())


@pytest.mark.parametrize("copy", ["oof", "weighted", "perturbation"])
@pytest.mark.parametrize("model,view,fold,parameter", [
    ("ridge", "study", 2, "alpha"),
    ("ridge", "variant", 2, "alpha"),
    ("gpr", "variant", 1, "length_scale"),
])
@pytest.mark.parametrize("mutation", ["bool", "string", "numpy", "nan", "infinity", "extra", "missing", "renamed", "nested"])
def test_each_hyperparameter_copy_requires_exact_schema_and_scalar_types(
    tmp_path: Path, copy: str, model: str, view: str, fold: int, parameter: str, mutation: str,
) -> None:
    """True and NumPy scalars must not compare equal to a frozen numeric 1.0."""
    result = _producer_result()
    if copy == "oof":
        table = result[f"{view}_oof"]
        table["selected_inner_hyperparameters"] = table["selected_inner_hyperparameters"].astype(object)
        index = table.index[(table.model_kind == model) & (table.outer_fold == fold)][0]
        params = json.loads(table.at[index, "selected_inner_hyperparameters"])
        table.at[index, "selected_inner_hyperparameters"] = params
    elif copy == "weighted":
        params = result["sensitivity"][model]["weighted"][view]["fold_hyperparameters"][fold]
    else:
        params = result["sensitivity"][model]["perturbation"]["fold_hyperparameters"][view][fold]
    assert params[parameter] == 1.0  # Hand-checked frozen fold, so bool aliases its value.
    if mutation == "bool": params[parameter] = True
    elif mutation == "string": params[parameter] = "1.0"
    elif mutation == "numpy": params[parameter] = np.float64(1.0)
    elif mutation == "nan": params[parameter] = float("nan")
    elif mutation == "infinity": params[parameter] = float("inf")
    elif mutation == "extra": params["amplitude"] = 1.0
    elif mutation == "missing": params.pop(parameter)
    elif mutation == "renamed": params["alias"] = params.pop(parameter)
    else: params[parameter] = [1.0]
    _rejects(tmp_path, result)


@pytest.mark.parametrize("field", ["outer_fold", "inner_fold", "train_indices", "test_indices"])
@pytest.mark.parametrize("alias", [False, 0.0, "0", np.int64(0), True, 1.0, "1", np.int64(1)], ids=["false", "zero_float", "zero_string", "zero_numpy", "true", "one_float", "one_string", "one_numpy"])
def test_inner_audit_rejects_nonbuiltin_integer_aliases(tmp_path: Path, field: str, alias: object) -> None:
    """A correctly valued bool, float, string, or NumPy index cannot earn SUCCESS."""
    result = _producer_result()
    expected = int(alias)
    for record in result["inner_fold_audit"]["records"]:
        if field.endswith("indices"):
            if expected in record[field]:
                record[field][record[field].index(expected)] = alias
                break
        elif record[field] == expected:
            record[field] = alias
            break
    else:
        pytest.fail("fixture has no requested canonical index")
    _rejects(tmp_path, result)


@pytest.mark.parametrize("mutation", ["schema_bool", "schema_float", "extra_field", "missing_field", "indices_string", "indices_nested", "scaler_bool", "scaler_numpy", "scaler_nested", "records_mapping"])
def test_inner_audit_requires_typed_nested_record_schema(tmp_path: Path, mutation: str) -> None:
    result = _producer_result()
    audit = result["inner_fold_audit"]
    record = audit["records"][9]  # Its first scaler mean is exactly zero.
    if mutation == "schema_bool": audit["schema_version"] = True
    elif mutation == "schema_float": audit["schema_version"] = 1.0
    elif mutation == "extra_field": record["extra"] = 0
    elif mutation == "missing_field": record.pop("inner_fold")
    elif mutation == "indices_string": record["train_indices"] = json.dumps(record["train_indices"])
    elif mutation == "indices_nested": record["train_indices"] = [[index] for index in record["train_indices"]]
    elif mutation == "scaler_bool": record["scaler_mean"][0] = False
    elif mutation == "scaler_numpy": record["scaler_mean"][0] = np.float64(0.0)
    elif mutation == "scaler_nested": record["scaler_mean"] = [[value] for value in record["scaler_mean"]]
    else: audit["records"] = dict(enumerate(audit["records"]))
    _rejects(tmp_path, result)


@pytest.mark.parametrize("copy", ["weighted", "perturbation"])
@pytest.mark.parametrize("malformed", [False, None, "invalid", {0: {"alpha": 1.0}}, np.array([1.0])], ids=["bool", "none", "string", "mapping", "array"])
def test_fold_parameter_copies_require_real_sequences(tmp_path: Path, copy: str, malformed: object) -> None:
    result = _producer_result()
    model = result["sensitivity"]["ridge"]
    if copy == "weighted": model[copy]["study"]["fold_hyperparameters"] = malformed
    else: model[copy]["fold_hyperparameters"]["study"] = malformed
    _rejects(tmp_path, result)

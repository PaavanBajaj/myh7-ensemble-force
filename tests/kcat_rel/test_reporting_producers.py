"""Exercise report completion with native sensitivity producers, without fitting."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from kcat_rel import cli
from kcat_rel.reporting import write_diagnostic_artifacts
from kcat_rel.sensitivity import PerturbationSummary, SensitivityResult
from kcat_rel.validation import Metrics
from .test_reporting import _complete_result


def _producer_result() -> dict:
    """Rehydrate recorded values into the real producer dataclasses and CLI asdict path."""
    result = _complete_result()
    for model in result["sensitivity"].values():
        for view, record in model["weighted"].items():
            weighted = SensitivityResult(
                view=view,
                model_kind=record["model_kind"],
                prediction=np.array(record["prediction"], dtype=float),
                baseline_prediction=np.array(record["baseline_prediction"], dtype=float),
                metrics=Metrics(**record["metrics"]),
                baseline_metrics=Metrics(**record["baseline_metrics"]),
                fold_hyperparameters=tuple(record["fold_hyperparameters"]),
                fold_sample_weights=tuple(tuple(row) for row in record["fold_sample_weights"]),
                fold_reported_variances=tuple(tuple(row) for row in record["fold_reported_variances"]),
            )
            model["weighted"][view] = cli.asdict(weighted)
        record = model["perturbation"]
        perturbation = PerturbationSummary(
            draws=record["draws"], seed=record["seed"],
            metric_quantiles=record["metric_quantiles"],
            gate_pass_count=record["gate_pass_count"],
            gate_pass_rate=record["gate_pass_rate"],
            fold_hyperparameters={view: tuple(params) for view, params in record["fold_hyperparameters"].items()},
        )
        model["perturbation"] = cli.asdict(perturbation)
    return result


@pytest.mark.parametrize("container", ["arrays", "tuples", "both"])
def test_real_sensitivity_producers_complete_report(tmp_path: Path, container: str) -> None:
    """Rejecting producer vectors or tuple fold copies must prevent this SUCCESS."""
    result = _producer_result()
    if container == "arrays":
        for model, record in _complete_result()["sensitivity"].items():
            result["sensitivity"][model]["perturbation"] = record["perturbation"]
    elif container == "tuples":
        for model in result["sensitivity"].values():
            for record in model["weighted"].values():
                record["prediction"] = record["prediction"].tolist()
                record["baseline_prediction"] = record["baseline_prediction"].tolist()

    write_diagnostic_artifacts(tmp_path, result)

    assert (tmp_path / "SUCCESS").read_text() == "completed\n"
    written = json.loads((tmp_path / "sensitivity.json").read_text())
    assert written["ridge"]["perturbation"]["fold_hyperparameters"]["variant"][0] == {"alpha": 0.1}
    assert len(written["gpr"]["weighted"]["study"]["prediction"]) == 17
    assert written["ridge"]["perturbation"]["draws"] == 1000
    assert json.loads((tmp_path / "metrics.json").read_text())["winner"] is None


@pytest.mark.parametrize("field", ["prediction", "baseline_prediction"])
@pytest.mark.parametrize("mutation", ["column", "row", "scalar", "bool", "string", "object_nested", "object_numeric"])
def test_native_prediction_container_rejects_malformed_vectors(tmp_path: Path, field: str, mutation: str) -> None:
    """Array support must never flatten, broadcast, or coerce malformed producer values."""
    result = _producer_result()
    record = result["sensitivity"]["gpr"]["weighted"]["study"]
    values = record[field]
    if mutation == "column": record[field] = values.reshape(17, 1)
    elif mutation == "row": record[field] = values.reshape(1, 17)
    elif mutation == "scalar": record[field] = np.array(values[0])
    elif mutation == "bool": record[field] = np.zeros(17, dtype=bool)
    elif mutation == "string": record[field] = values.astype(str)
    elif mutation == "object_numeric": record[field] = values.astype(object)
    else:
        nested = np.empty(17, dtype=object)
        nested[:] = [[value] for value in values]
        record[field] = nested
    with pytest.raises(ValueError):
        write_diagnostic_artifacts(tmp_path, result)
    assert not any(tmp_path.iterdir())

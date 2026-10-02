import copy
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import kcat_rel.modeling_summary as modeling_summary
from kcat_rel.modeling_summary import (
    SummaryInputs,
    SummaryPaths,
    compute_feature_correlations,
    compute_summary_checks,
    load_summary_inputs,
    summary_payload,
    write_summary_outputs,
)
from kcat_rel.registry import ACTIVE_FEATURE_COLUMNS, FROZEN_SOURCE_STUDIES


WORKFLOW_ROOT = Path(__file__).resolve().parents[2] / "workflows/kcat-rel"
DATA_ROOT = Path(__file__).resolve().parents[2] / "data/public/kcat-rel"

_FROZEN_SUMMARY_ARTIFACT_SHA256 = {
    "config": "73e72e6c31da38e6f14d810190709f815e8456a22bbf8bf83be4b105301d113b",
    "features_csv": "6bd0fb2dc2d20eb184bee96cb4c99dad1d16d015bf0725b53b6edb1d0893e372",
    "feature_manifest": "527aec1ffcc2fd21d67211fbb4bcf8a1b3ee4fa9c1e0ec17ddb6802171fa272c",
    "model_table": "bf93b88a548d5393aebd6b1e350f5e99181dda9833beadfe68d902dfdf1cefef",
    "study_oof": "068d7b454d6d1425212e6464fd0e15a05d58c3a126fa3b0ecbc642c44e18bf5e",
    "variant_oof": "7276fa923bac3d26adeb2d3e06139012298ea1798d00b357c337c9118c0fe24d",
    "metrics": "130ae700646e311ea1b146a108d7df3867b598bc4cc71ab3cdcbb45dda15caa7",
    "sensitivity": "0b7bef3e096b48e8acc491b4623aff5d5c1f378399bf4483ae6dde77ba3e47bf",
    "tuning_records": "cdf084be72a25be32e9704cf1fc33b43f6e13ff27b480880d2bd723959446ebf",
    "run_manifest": "a6f3b0dcadcbe03447ff03b05bc48c149a2e8cf1bbd6b73656a4b1323d1aa399",
}


def _copied_summary_paths(tmp_path: Path) -> SummaryPaths:
    """Copy the committed frozen inputs so each mutation is isolated."""
    sources = {
        "config": WORKFLOW_ROOT / "config" / "kcat-rel-v0.json",
        "features": DATA_ROOT / "derived" / "kcat-rel-v0-no-msa-features.csv",
        "feature_manifest": DATA_ROOT / "derived" / "kcat-rel-v0-no-msa-feature-manifest.json",
        "model_table": DATA_ROOT / "derived" / "kcat-rel-v0-no-msa-model-table.csv",
        "study_oof": WORKFLOW_ROOT / "results" / "diagnostic-v0" / "study-oof.csv",
        "variant_oof": WORKFLOW_ROOT / "results" / "diagnostic-v0" / "variant-oof.csv",
        "metrics": WORKFLOW_ROOT / "results" / "diagnostic-v0" / "metrics.json",
        "sensitivity": WORKFLOW_ROOT / "results" / "diagnostic-v0" / "sensitivity.json",
        "tuning_records": WORKFLOW_ROOT / "results" / "diagnostic-v0" / "tuning-records.csv",
        "run_manifest": WORKFLOW_ROOT / "results" / "diagnostic-v0" / "run-manifest.json",
    }
    copies = {}
    for name, source in sources.items():
        destination = tmp_path / source.name
        shutil.copy2(source, destination)
        copies[name] = destination
    return SummaryPaths(**copies)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _refresh_artifact_hash(paths: SummaryPaths, artifact_name: str, source: Path) -> None:
    """Deliberately bypass the byte pin only when testing semantic validation."""
    manifest = json.loads(Path(paths.run_manifest).read_text())
    manifest["artifact_sha256"][artifact_name] = _sha256(source)
    Path(paths.run_manifest).write_text(json.dumps(manifest), encoding="utf-8")


def _repin_test_manifest_root(monkeypatch: pytest.MonkeyPatch, paths: SummaryPaths) -> None:
    """Bypass the root only in tests that need to reach semantic checks."""
    monkeypatch.setattr(modeling_summary, "FROZEN_RUN_MANIFEST_SHA256", _sha256(Path(paths.run_manifest)))


def test_modeling_summary_import_does_not_load_fitting_modules():
    """Catch a read-only summary module that pulls validation or fitting code into memory."""
    environment = {**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[2] / "src")}
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import kcat_rel.modeling_summary; "
            "assert not {'kcat_rel.cli', 'kcat_rel.validation', 'kcat_rel.sensitivity'} & set(sys.modules)",
        ],
        env=environment,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr


def test_compute_feature_correlations_uses_pearson_r_in_registry_order():
    """Catch a non-Pearson statistic, reordered registry, or wrong target vector."""
    features = pd.DataFrame(
        {
            "delta_charge": [-1.0, 0.0, 1.0],
            "grantham_distance": [1.0, 2.0, 3.0],
            "catalytic_motif_min_ca_distance_8act_mean_ab": [3.0, 1.0, 2.0],
        }
    )
    y = np.array([-2.0, 0.0, 2.0])

    got = compute_feature_correlations(features, y)

    assert list(got) == list(ACTIVE_FEATURE_COLUMNS)
    assert got["delta_charge"] == pytest.approx(1.0)
    assert got["grantham_distance"] == pytest.approx(1.0)
    assert got["catalytic_motif_min_ca_distance_8act_mean_ab"] == pytest.approx(-0.5)


@pytest.mark.parametrize(
    "features, target",
    [
        (pd.DataFrame({name: [1.0, 1.0, 1.0] for name in ACTIVE_FEATURE_COLUMNS}), np.array([0.0, 1.0, 2.0])),
        (pd.DataFrame({name: [0.0, 1.0, 2.0] for name in ACTIVE_FEATURE_COLUMNS}), np.array([0.0, np.nan, 2.0])),
    ],
)
def test_compute_feature_correlations_rejects_undefined_or_nonfinite_inputs(features: pd.DataFrame, target: np.ndarray):
    """Catch misleading descriptive correlations from a constant or nonfinite vector."""
    with pytest.raises(ValueError):
        compute_feature_correlations(features, target)


@pytest.mark.parametrize(
    "mutation",
    [
        "reordered_variants",
        "noncanonical_variant",
        "changed_truth",
        "changed_source",
        "missing_oof_row",
        "duplicate_oof_row",
        "renamed_feature",
        "nonfinite_value",
    ],
)
def test_load_summary_inputs_rejects_mutated_frozen_artifact_identity(tmp_path: Path, mutation: str):
    """Catch coordinated-looking inputs that drift from the canonical frozen run."""
    paths = _copied_summary_paths(tmp_path)
    if mutation == "reordered_variants":
        table = pd.read_csv(paths.model_table)
        table = pd.concat([table.iloc[1:], table.iloc[:1]], ignore_index=True)
        table.to_csv(paths.model_table, index=False)
    elif mutation == "noncanonical_variant":
        table = pd.read_csv(paths.model_table)
        table.loc[0, "variant"] = "D382Y"
        table.to_csv(paths.model_table, index=False)
    elif mutation == "changed_truth":
        table = pd.read_csv(paths.study_oof)
        table.loc[0, "truth"] += 1e-6
        table.to_csv(paths.study_oof, index=False)
    elif mutation == "changed_source":
        table = pd.read_csv(paths.variant_oof)
        table.loc[0, "source_id"] = "wrong_study"
        table.to_csv(paths.variant_oof, index=False)
    elif mutation == "missing_oof_row":
        table = pd.read_csv(paths.study_oof)
        table.iloc[1:].to_csv(paths.study_oof, index=False)
    elif mutation == "duplicate_oof_row":
        table = pd.read_csv(paths.variant_oof)
        pd.concat([table, table.iloc[[0]]], ignore_index=True).to_csv(paths.variant_oof, index=False)
    elif mutation == "renamed_feature":
        table = pd.read_csv(paths.model_table)
        table.rename(columns={ACTIVE_FEATURE_COLUMNS[0]: "renamed_feature"}).to_csv(paths.model_table, index=False)
    else:
        table = pd.read_csv(paths.model_table)
        table.loc[0, ACTIVE_FEATURE_COLUMNS[0]] = np.nan
        table.to_csv(paths.model_table, index=False)

    with pytest.raises(ValueError):
        load_summary_inputs(paths)


@pytest.mark.parametrize("mutation", ["truth", "source_and_fold"])
def test_load_summary_inputs_rejects_coordinated_primary_table_and_oof_mutations_by_hash(
    tmp_path: Path, mutation: str
):
    """Catch circular trust when coordinated primary-table/OOF edits retain semantic agreement."""
    paths = _copied_summary_paths(tmp_path)
    model_table = pd.read_csv(paths.model_table)
    study_oof = pd.read_csv(paths.study_oof)
    variant_oof = pd.read_csv(paths.variant_oof)
    if mutation == "truth":
        variant = model_table.loc[0, "variant"]
        model_table.loc[0, "y"] += 0.01
        for table in (study_oof, variant_oof):
            table.loc[table["variant"] == variant, "truth"] += 0.01
    else:
        first, second = model_table.loc[0, "variant"], model_table.loc[1, "variant"]
        source_by_variant = dict(zip(model_table["variant"], model_table["source_id"], strict=True))
        model_table.loc[:1, "source_id"] = [source_by_variant[second], source_by_variant[first]]
        for table in (study_oof, variant_oof):
            table.loc[table["variant"] == first, "source_id"] = source_by_variant[second]
            table.loc[table["variant"] == second, "source_id"] = source_by_variant[first]
        study_oof.loc[study_oof["variant"] == first, "outer_fold"] = FROZEN_SOURCE_STUDIES.index(source_by_variant[second])
        study_oof.loc[study_oof["variant"] == second, "outer_fold"] = FROZEN_SOURCE_STUDIES.index(source_by_variant[first])
    model_table.to_csv(paths.model_table, index=False)
    study_oof.to_csv(paths.study_oof, index=False)
    variant_oof.to_csv(paths.variant_oof, index=False)

    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        load_summary_inputs(paths)


@pytest.mark.parametrize("mutation", ["row_order", "hyperparameters", "scaler"])
def test_load_summary_inputs_rejects_hash_bypassed_oof_protocol_tampering(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mutation: str
):
    """Catch malformed OOF protocol data even if a test fixture updates its artifact hash."""
    paths = _copied_summary_paths(tmp_path)
    table = pd.read_csv(paths.study_oof)
    if mutation == "row_order":
        table.iloc[[0, 1]] = table.iloc[[1, 0]].to_numpy()
    elif mutation == "hyperparameters":
        table.loc[0, "selected_inner_hyperparameters"] = '{"length_scale":0.75,"noise_variance":0.04}'
    else:
        table.loc[0, "scaler_scale"] = "[1.0, 2.0]"
    table.to_csv(paths.study_oof, index=False)
    _refresh_artifact_hash(paths, "study-oof.csv", paths.study_oof)
    _repin_test_manifest_root(monkeypatch, paths)

    with pytest.raises(ValueError, match="row order|hyperparameters|scaler"):
        load_summary_inputs(paths)


@pytest.mark.parametrize("mutation", ["truth", "source_and_fold"])
def test_load_summary_inputs_rejects_repointed_manifest_after_coordinated_primary_mutation(
    tmp_path: Path, mutation: str
):
    """Catch a forged inner manifest whose hashes were repointed to coordinated data edits."""
    paths = _copied_summary_paths(tmp_path)
    model_table = pd.read_csv(paths.model_table)
    study_oof = pd.read_csv(paths.study_oof)
    variant_oof = pd.read_csv(paths.variant_oof)
    if mutation == "truth":
        variant = model_table.loc[0, "variant"]
        model_table.loc[0, "y"] += 0.01
        for table in (study_oof, variant_oof):
            table.loc[table["variant"] == variant, "truth"] += 0.01
    else:
        first, second = model_table.loc[0, "variant"], model_table.loc[1, "variant"]
        source_by_variant = dict(zip(model_table["variant"], model_table["source_id"], strict=True))
        model_table.loc[:1, "source_id"] = [source_by_variant[second], source_by_variant[first]]
        for table in (study_oof, variant_oof):
            table.loc[table["variant"] == first, "source_id"] = source_by_variant[second]
            table.loc[table["variant"] == second, "source_id"] = source_by_variant[first]
        study_oof.loc[study_oof["variant"] == first, "outer_fold"] = FROZEN_SOURCE_STUDIES.index(source_by_variant[second])
        study_oof.loc[study_oof["variant"] == second, "outer_fold"] = FROZEN_SOURCE_STUDIES.index(source_by_variant[first])
    model_table.to_csv(paths.model_table, index=False)
    study_oof.to_csv(paths.study_oof, index=False)
    variant_oof.to_csv(paths.variant_oof, index=False)
    for artifact, path in (("model-table.csv", paths.model_table), ("study-oof.csv", paths.study_oof), ("variant-oof.csv", paths.variant_oof)):
        _refresh_artifact_hash(paths, artifact, path)

    with pytest.raises(ValueError, match="immutable SHA-256"):
        load_summary_inputs(paths)


@pytest.mark.parametrize("mutation", ["alternate_parameter", "scaler"])
def test_load_summary_inputs_reconstructs_tuning_protocol_after_test_only_root_bypass(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mutation: str
):
    """Catch an approved-but-nonminimal parameter or non-training scaler after hash bypass."""
    paths = _copied_summary_paths(tmp_path)
    tuning = pd.read_csv(paths.tuning_records)
    study_oof = pd.read_csv(paths.study_oof)
    first = tuning.iloc[0]
    mask = (study_oof["model_kind"] == first["model_kind"]) & (study_oof["outer_fold"] == first["outer_fold"])
    if mutation == "alternate_parameter":
        replacement = '{"length_scale":0.5,"noise_variance":0.01}'
        tuning.loc[0, "selected_params"] = replacement
        study_oof.loc[mask, "selected_inner_hyperparameters"] = replacement
    else:
        tuning.loc[0, "scaler_mean"] = "[999.0, 999.0, 999.0]"
        study_oof.loc[mask, "scaler_mean"] = "[999.0, 999.0, 999.0]"
    tuning.to_csv(paths.tuning_records, index=False)
    study_oof.to_csv(paths.study_oof, index=False)
    _refresh_artifact_hash(paths, "tuning-records.csv", paths.tuning_records)
    _refresh_artifact_hash(paths, "study-oof.csv", paths.study_oof)
    _repin_test_manifest_root(monkeypatch, paths)

    with pytest.raises(ValueError, match="tuning|scaler|hyperparameter"):
        load_summary_inputs(paths)


def _synthetic_summary_inputs() -> SummaryInputs:
    """Four hand-scored rows: two equally weighted two-row source studies."""
    variants = ["v1", "v2", "v3", "v4"]
    sources = ["study_a", "study_a", "study_b", "study_b"]
    truth = [0.0, 2.0, 0.0, 0.0]
    prediction = [1.0, 1.0, 2.0, 0.0]
    baseline = [0.0, 0.0, 0.0, 0.0]
    model_table = pd.DataFrame(
        {
            "variant": variants,
            "y": truth,
            "y_error": [0.1, 0.1, 0.1, 0.1],
            "source_id": sources,
            "delta_charge": [-1.0, 0.0, 1.0, 2.0],
            "grantham_distance": [1.0, 2.0, 3.0, 4.0],
            "catalytic_motif_min_ca_distance_8act_mean_ab": [4.0, 3.0, 2.0, 1.0],
        }
    )

    def oof(view: str) -> pd.DataFrame:
        rows = []
        for model_kind in ("ridge", "gpr"):
            for index, variant in enumerate(variants):
                rows.append(
                    {
                        "model_kind": model_kind,
                        "variant": variant,
                        "source_id": sources[index],
                        "truth": truth[index],
                        "model_prediction": prediction[index],
                        "baseline_prediction": baseline[index],
                        "outer_fold": index if view == "variant" else index // 2,
                        "selected_inner_hyperparameters": "{}",
                        "scaler_mean": "[]",
                        "scaler_scale": "[]",
                    }
                )
        return pd.DataFrame(rows)

    # Variant metrics: MAE=1 and RMSE=sqrt(1.5).  Equal-study macro metrics:
    # MAE=1 and RMSE=(1 + sqrt(2))/2.  Baseline MAE is 0.5; its study
    # macro RMSE is sqrt(2)/2, while its pooled variant RMSE is 1.
    metrics = {
        "winner": None,
        "models": {
            model: {
                "pass_gate": False,
                "failed_rules": [
                    "study_macro_mae",
                    "variant_mae",
                    "study_macro_rmse",
                    "variant_rmse",
                    "single_study_influence",
                ],
                "mae_improvement": {"study": -1.0, "variant": -1.0},
                "rmse_degradation": {
                    "study": 0.7071067811865475,
                    "variant": 0.22474487139158894,
                },
                "model_views": {
                    "study": {"mae": 1.0, "rmse": 1.2071067811865475},
                    "variant": {"mae": 1.0, "rmse": 1.224744871391589},
                },
                "baseline_views": {
                    "study": {"mae": 0.5, "rmse": 0.7071067811865476},
                    "variant": {"mae": 0.5, "rmse": 1.0},
                },
            }
            for model in ("ridge", "gpr")
        },
    }
    registry = SimpleNamespace(
        source_studies=FROZEN_SOURCE_STUDIES,
        gate_thresholds={
            "minimum_mae_improvement_fraction": 0.05,
            "maximum_rmse_degradation_fraction": 0.05,
        },
        seed=20260916,
        model_grids={
            "ridge_alphas": [0.01],
            "gpr": {"length_scales": [0.5], "noise_variances": [0.04]},
        },
    )

    def parameters(model: str) -> dict[str, float]:
        return {"alpha": 0.01} if model == "ridge" else {"length_scale": 0.5, "noise_variance": 0.04}

    def weighted(model: str, view: str) -> dict:
        folds = 9 if view == "study" else 17
        return {
            "model_kind": model,
            "view": view,
            "prediction": prediction,
            "baseline_prediction": baseline,
            "metrics": metrics["models"][model]["model_views"][view],
            "baseline_metrics": metrics["models"][model]["baseline_views"][view],
            "fold_hyperparameters": [parameters(model)] * folds,
            "fold_sample_weights": [[] for _ in range(folds)],
            "fold_reported_variances": [[] for _ in range(folds)],
        }

    sensitivity = {
        model: {
            "weighted": {view: weighted(model, view) for view in ("study", "variant")},
            "perturbation": {
                "draws": 1000,
                "seed": 20260916,
                "metric_quantiles": {},
                "gate_pass_count": 0,
                "gate_pass_rate": 0.0,
                "fold_hyperparameters": {
                    "study": [parameters(model)] * 9,
                    "variant": [parameters(model)] * 17,
                },
            },
        }
        for model in ("ridge", "gpr")
    }
    paths = SummaryPaths(*(["unused"] * 10))
    return SummaryInputs(
        paths,
        registry,
        model_table.loc[:, ACTIVE_FEATURE_COLUMNS],
        model_table,
        oof("study"),
        oof("variant"),
        None,
        metrics,
        sensitivity,
        {},
    )


def test_compute_summary_checks_reconstructs_hand_scored_variant_and_study_metrics():
    """Catch pooled study scoring or formulas that diverge from frozen OOF definitions."""
    checks = compute_summary_checks(_synthetic_summary_inputs())

    assert checks.model_metrics["ridge"]["variant"] == pytest.approx({"mae": 1.0, "rmse": 1.224744871391589})
    assert checks.model_metrics["ridge"]["study"] == pytest.approx({"mae": 1.0, "rmse": 1.2071067811865475})
    assert checks.baseline_metrics["gpr"]["study"] == pytest.approx({"mae": 0.5, "rmse": 0.7071067811865476})
    assert checks.mae_improvement["ridge"] == pytest.approx({"study": -1.0, "variant": -1.0})
    assert checks.rmse_degradation["ridge"] == pytest.approx(
        {"study": 0.7071067811865475, "variant": 0.22474487139158894}
    )


@pytest.mark.parametrize("mutation", ["metric", "winner", "failed_rule"])
def test_compute_summary_checks_rejects_metrics_not_reconstructed_from_oof(mutation: str):
    """Catch an altered score, unauthorized winner, or omitted gate failure."""
    inputs = _synthetic_summary_inputs()
    if mutation == "metric":
        inputs.metrics["models"]["ridge"]["model_views"]["study"]["mae"] += 1e-9
    elif mutation == "winner":
        inputs.metrics["winner"] = "ridge"
    else:
        inputs.metrics["models"]["gpr"]["failed_rules"].remove("variant_rmse")

    with pytest.raises(ValueError):
        compute_summary_checks(inputs)


@pytest.mark.parametrize("mutation", ["provenance", "metric", "prediction"])
def test_compute_summary_checks_rejects_tampered_weighted_sensitivity(tmp_path: Path, mutation: str):
    """Catch a relabeled, rescored, or changed weighted secondary prediction array."""
    inputs = load_summary_inputs(_copied_summary_paths(tmp_path))
    weighted = inputs.sensitivity["ridge"]["weighted"]["study"]
    if mutation == "provenance":
        weighted["model_kind"] = "gpr"
    elif mutation == "metric":
        weighted["metrics"]["mae"] += 1e-9
    else:
        weighted["prediction"][0] += 1e-6

    with pytest.raises(ValueError, match="weighted sensitivity"):
        compute_summary_checks(inputs)


@pytest.mark.parametrize("mutation", ["ridge_weights", "gpr_variances", "none_weights"])
def test_compute_summary_checks_rejects_invalid_weighted_fold_provenance(
    tmp_path: Path, mutation: str
):
    """Catch non-training uncertainty weights, wrong variances, and non-sequence fold records."""
    inputs = load_summary_inputs(_copied_summary_paths(tmp_path))
    if mutation == "ridge_weights":
        inputs.sensitivity["ridge"]["weighted"]["study"]["fold_sample_weights"][0][0] += 1e-6
    elif mutation == "gpr_variances":
        inputs.sensitivity["gpr"]["weighted"]["variant"]["fold_reported_variances"][0][0] += 1e-6
    else:
        inputs.sensitivity["gpr"]["weighted"]["study"]["fold_sample_weights"][0] = None

    with pytest.raises(ValueError, match="weighted sensitivity"):
        compute_summary_checks(inputs)


@pytest.mark.parametrize("model_kind,view", [(model, view) for model in ("ridge", "gpr") for view in ("study", "variant")])
def test_compute_summary_checks_rejects_perturbation_parameters_that_differ_from_primary_tuning(
    tmp_path: Path, model_kind: str, view: str
):
    """Catch an approved candidate substituted into a perturbation fold sequence."""
    inputs = load_summary_inputs(_copied_summary_paths(tmp_path))
    parameters = inputs.sensitivity[model_kind]["perturbation"]["fold_hyperparameters"][view]
    if model_kind == "ridge":
        parameters[0] = {"alpha": 0.01 if parameters[0]["alpha"] != 0.01 else 0.1}
    else:
        parameters[0] = {
            "length_scale": 1.0 if parameters[0]["length_scale"] != 1.0 else 0.5,
            "noise_variance": parameters[0]["noise_variance"],
        }

    with pytest.raises(ValueError, match="perturbation sensitivity hyperparameters"):
        compute_summary_checks(inputs)


@pytest.fixture(scope="module")
def frozen_parameter_inputs(tmp_path_factory: pytest.TempPathFactory) -> SummaryInputs:
    return load_summary_inputs(_copied_summary_paths(tmp_path_factory.mktemp("parameter-types")))


# Literal frozen folds with integer-valued selections permit same-valued aliases.
_PARAMETER_FOLDS = [
    ("ridge", "study", 2, "alpha", 1.0),
    ("ridge", "variant", 2, "alpha", 1.0),
    ("gpr", "study", 6, "length_scale", 2.0),
    ("gpr", "variant", 1, "length_scale", 1.0),
]


def _mutate_parameter_copy(inputs, source, model_kind, view, fold, key, value):
    """Change one actual artifact parameter while retaining every other fold value."""
    if source in ("oof", "tuning"):
        table = getattr(inputs, f"{view}_oof") if source == "oof" else inputs.tuning_records
        column = "selected_inner_hyperparameters" if source == "oof" else "selected_params"
        mask = (table["model_kind"] == model_kind) & (table["outer_fold"] == fold)
        if source == "tuning":
            mask &= table["view"] == view
        index = table.index[mask][0]
        params = json.loads(table.at[index, column])
        params[key] = value
        # JSON would erase np.float64's type. Preserve NumPy objects in memory;
        # tuning's JSON-text boundary must reject these mappings before parsing.
        if isinstance(value, np.generic):
            table[column] = table[column].astype(object)
            table.at[index, column] = params
        else:
            table.at[index, column] = json.dumps(params)
    elif source == "weighted":
        inputs.sensitivity[model_kind][source][view]["fold_hyperparameters"][fold][key] = value
    else:
        inputs.sensitivity[model_kind][source]["fold_hyperparameters"][view][fold][key] = value


def _check_parameter_copy(inputs, source, view):
    if source == "oof":
        modeling_summary._validate_oof(getattr(inputs, f"{view}_oof"), inputs.model_table, view, inputs.registry)
    elif source == "tuning":
        modeling_summary._validate_tuning_records(
            inputs.tuning_records, inputs.model_table, inputs.study_oof, inputs.variant_oof, inputs.registry
        )
    else:
        return compute_summary_checks(inputs)


@pytest.mark.parametrize("source", ["oof", "tuning", "weighted", "perturbation"])
@pytest.mark.parametrize(
    "model_kind,view,fold,key,original,alias",
    [
        pytest.param(*case, alias, id=f"{case[0]}-{case[1]}-{alias.__name__}")
        for case in _PARAMETER_FOLDS
        for alias in (np.float64, np.int64, np.bool_, bool, str, list)
        if alias not in (np.bool_, bool) or case[-1] == 1.0
    ],
)
def test_parameter_copies_reject_numeric_aliases_with_context(
    frozen_parameter_inputs, source, model_kind, view, fold, key, original, alias
):
    """Catch numeric equality admitting non-built-in parameter values in any copy."""
    inputs = copy.deepcopy(frozen_parameter_inputs)
    primary = inputs.tuning_records
    selected = primary.loc[
        (primary["model_kind"] == model_kind) & (primary["view"] == view) & (primary["outer_fold"] == fold),
        "selected_params",
    ].iloc[0]
    assert json.loads(selected)[key] == original
    replacement = [original] if alias is list else alias(original)
    if alias not in (str, list):
        assert replacement == original
    _mutate_parameter_copy(inputs, source, model_kind, view, fold, key, replacement)
    context = {
        "oof": f"{view} OOF hyperparameters",
        "tuning": "tuning selected hyperparameters",
        "weighted": "weighted sensitivity hyperparameters",
        "perturbation": "perturbation sensitivity hyperparameters",
    }[source]

    with pytest.raises(ValueError, match=context):
        _check_parameter_copy(inputs, source, view)


@pytest.mark.parametrize("source", ["oof", "tuning", "weighted", "perturbation"])
@pytest.mark.parametrize("model_kind,view,fold,key,original", _PARAMETER_FOLDS)
@pytest.mark.parametrize("number_type", [int, float])
def test_parameter_copies_accept_exact_builtin_approved_numbers(
    frozen_parameter_inputs, source, model_kind, view, fold, key, original, number_type
):
    """Catch an over-strict fix rejecting equivalent approved built-in ints or floats."""
    inputs = copy.deepcopy(frozen_parameter_inputs)
    _mutate_parameter_copy(inputs, source, model_kind, view, fold, key, number_type(original))

    checks = _check_parameter_copy(inputs, source, view)

    if source in ("weighted", "perturbation"):
        assert checks.pass_gate == {"ridge": False, "gpr": False}


def test_frozen_inputs_produce_read_only_descriptive_check_payload(tmp_path: Path):
    """Catch a loader that mutates frozen mappings or loses required summary headlines."""
    inputs = load_summary_inputs(_copied_summary_paths(tmp_path))
    metrics_before = copy.deepcopy(inputs.metrics)
    sensitivity_before = copy.deepcopy(inputs.sensitivity)

    checks = compute_summary_checks(inputs)
    payload = summary_payload(checks)

    assert inputs.metrics == metrics_before
    assert inputs.sensitivity == sensitivity_before
    assert list(checks.feature_correlations) == list(ACTIVE_FEATURE_COLUMNS)
    assert len(checks.study_oof) == 34
    assert len(checks.variant_oof) == 34
    assert checks.pass_gate == {"ridge": False, "gpr": False}
    assert checks.perturbation_sensitivity["ridge"] == {
        "draws": 1000,
        "seed": 20260916,
        "gate_pass_count": 1,
        "gate_pass_rate": 0.001,
    }
    assert payload["failed_rules"]["gpr"] == [
        "study_macro_mae",
        "variant_mae",
        "study_macro_rmse",
        "variant_rmse",
        "single_study_influence",
    ]


def test_frozen_checks_serialize_independently_calculated_artifact_hashes(tmp_path: Path):
    """Catch provenance copied from historical fields instead of hashing exact loaded bytes."""
    checks = compute_summary_checks(load_summary_inputs(_copied_summary_paths(tmp_path)))
    payload = summary_payload(checks)

    assert checks.run_manifest_root_sha256 == _FROZEN_SUMMARY_ARTIFACT_SHA256["run_manifest"]
    assert checks.artifact_sha256 == _FROZEN_SUMMARY_ARTIFACT_SHA256
    assert payload["run_manifest_root_sha256"] == _FROZEN_SUMMARY_ARTIFACT_SHA256["run_manifest"]
    assert payload["artifact_sha256"] == _FROZEN_SUMMARY_ARTIFACT_SHA256


def test_load_summary_inputs_rejects_tampered_artifact_before_serializing_its_hash(tmp_path: Path):
    """Catch a summary that serializes a changed frozen artifact instead of rejecting it."""
    paths = _copied_summary_paths(tmp_path)
    with Path(paths.metrics).open("ab") as handle:
        handle.write(b"\n")

    with pytest.raises(ValueError, match="metrics SHA-256 mismatch"):
        load_summary_inputs(paths)


@pytest.fixture(scope="module")
def frozen_checks(tmp_path_factory: pytest.TempPathFactory):
    """Provide the validated frozen summary values for rendering boundaries."""
    paths = _copied_summary_paths(tmp_path_factory.mktemp("summary-figures"))
    return compute_summary_checks(load_summary_inputs(paths))


def test_write_summary_outputs_is_deterministic_and_complete(tmp_path: Path, frozen_checks):
    """Catch nondeterministic JSON, incomplete figures, or partial-success output."""
    first = tmp_path / "first"
    second = tmp_path / "second"

    paths_a = write_summary_outputs(frozen_checks, first)
    paths_b = write_summary_outputs(frozen_checks, second)

    assert (first / "modeling-summary.json").read_bytes() == (second / "modeling-summary.json").read_bytes()
    assert [path.read_bytes() for path in paths_a[1]] == [
        (second / path.name).read_bytes() for path in paths_a[1]
    ]
    assert {path.name for path in paths_a[1]} == {
        "feature-correlations.png",
        "study-oof-observed-vs-predicted.png",
        "held-out-error-vs-baseline.png",
    }
    assert all(path.stat().st_size > 10_000 for path in paths_a[1])
    payload = json.loads(paths_a[0].read_text(encoding="utf-8"))
    assert payload["status"] == "posthoc_descriptive_only"
    assert payload["cohort_size"] == 17
    assert payload["source_study_count"] == 9
    assert payload["feature_names"] == list(ACTIVE_FEATURE_COLUMNS)


def _frozen_summary_cli_args(output_dir: Path) -> list[str]:
    """Pass every frozen input explicitly through the standalone CLI boundary."""
    diagnostic = WORKFLOW_ROOT / "results" / "diagnostic-v0"
    return [
        "--config", str(WORKFLOW_ROOT / "config" / "kcat-rel-v0.json"),
        "--features", str(DATA_ROOT / "derived" / "kcat-rel-v0-no-msa-features.csv"),
        "--feature-manifest", str(DATA_ROOT / "derived" / "kcat-rel-v0-no-msa-feature-manifest.json"),
        "--model-table", str(DATA_ROOT / "derived" / "kcat-rel-v0-no-msa-model-table.csv"),
        "--study-oof", str(diagnostic / "study-oof.csv"),
        "--variant-oof", str(diagnostic / "variant-oof.csv"),
        "--metrics", str(diagnostic / "metrics.json"),
        "--sensitivity", str(diagnostic / "sensitivity.json"),
        "--tuning-records", str(diagnostic / "tuning-records.csv"),
        "--run-manifest", str(diagnostic / "run-manifest.json"),
        "--output-dir", str(output_dir),
    ]


def test_main_renders_exact_outputs_without_mutating_explicit_frozen_inputs(tmp_path: Path):
    """Catch a CLI that omits an output or changes a frozen source artifact."""
    args = _frozen_summary_cli_args(tmp_path / "summary")
    input_paths = [Path(value) for flag, value in zip(args[::2], args[1::2], strict=True) if flag != "--output-dir"]
    hashes_before = {path: _sha256(path) for path in input_paths}

    assert modeling_summary.main(args) == 0

    assert {path.name for path in (tmp_path / "summary").iterdir()} == {
        "modeling-summary.json",
        "feature-correlations.png",
        "study-oof-observed-vs-predicted.png",
        "held-out-error-vs-baseline.png",
    }
    assert {path: _sha256(path) for path in input_paths} == hashes_before


def test_main_rejects_an_occupied_destination_before_rendering(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Catch a CLI that renders into an occupied path despite its atomic contract."""
    occupied = tmp_path / "occupied"
    occupied.mkdir()
    sentinel = occupied / "preserve-me.txt"
    sentinel.write_text("existing artifact\n", encoding="utf-8")

    def rendering_must_not_start(*_args, **_kwargs):
        raise AssertionError("rendering began for an occupied destination")

    monkeypatch.setattr(modeling_summary, "write_summary_outputs", rendering_must_not_start)

    with pytest.raises(ValueError, match="missing or empty"):
        modeling_summary.main(_frozen_summary_cli_args(occupied))

    assert sentinel.read_text(encoding="utf-8") == "existing artifact\n"


def test_main_defaults_provenance_paths_beside_explicit_metrics(tmp_path: Path):
    """Catch a documented CLI invocation that needlessly requires derived provenance flags."""
    args = _frozen_summary_cli_args(tmp_path / "summary")
    for flag in ("--tuning-records", "--run-manifest"):
        index = args.index(flag)
        del args[index : index + 2]

    assert modeling_summary.main(args) == 0


def test_write_summary_outputs_cleans_staging_and_preserves_empty_destination_on_render_failure(
    tmp_path: Path, frozen_checks, monkeypatch: pytest.MonkeyPatch
):
    """Catch a renderer exception that leaves staged files or alters an existing destination."""
    from kcat_rel import summary_figures

    destination = tmp_path / "empty-destination"
    destination.mkdir()

    def fail_after_first_figure(*_args, **_kwargs):
        raise RuntimeError("intentional rendering failure")

    monkeypatch.setattr(summary_figures, "render_study_oof", fail_after_first_figure)

    with pytest.raises(RuntimeError, match="intentional rendering failure"):
        write_summary_outputs(frozen_checks, destination)

    assert list(destination.iterdir()) == []
    assert not list(tmp_path.glob(".empty-destination.staging-*"))


def test_write_summary_outputs_closes_figures_and_cleans_staging_when_savefig_fails(
    tmp_path: Path, frozen_checks, monkeypatch: pytest.MonkeyPatch
):
    """Catch a save failure that leaks a live figure or publishes partial staged output."""
    from kcat_rel import summary_figures

    destination = tmp_path / "empty-destination"
    destination.mkdir()
    initially_live = set(summary_figures.plt.get_fignums())
    failed_figure_numbers: list[int] = []

    def fail_savefig(figure, *_args, **_kwargs):
        failed_figure_numbers.append(figure.number)
        raise OSError("intentional savefig failure")

    monkeypatch.setattr(summary_figures.plt.Figure, "savefig", fail_savefig)

    with pytest.raises(OSError, match="intentional savefig failure"):
        write_summary_outputs(frozen_checks, destination)

    assert failed_figure_numbers
    assert not set(failed_figure_numbers) & set(summary_figures.plt.get_fignums())
    assert set(summary_figures.plt.get_fignums()) == initially_live
    assert list(destination.iterdir()) == []
    assert not list(tmp_path.glob(".empty-destination.staging-*"))


def test_render_study_oof_x_limits_cover_every_uncertainty_endpoint_with_padding(
    tmp_path: Path, frozen_checks, monkeypatch: pytest.MonkeyPatch
):
    """Catch shared study-panel limits that clip horizontal observed-label error bars."""
    from kcat_rel import summary_figures

    rows = frozen_checks.study_oof
    truth = rows["truth"].to_numpy(dtype=float)
    errors = np.asarray([frozen_checks.label_uncertainties[str(variant)] for variant in rows["variant"]])
    independent_values = np.concatenate(
        [
            truth - errors,
            truth + errors,
            rows["model_prediction"].to_numpy(dtype=float),
            rows["baseline_prediction"].to_numpy(dtype=float),
        ]
    )
    independent_minimum = float(independent_values.min())
    independent_maximum = float(independent_values.max())
    padding = max((independent_maximum - independent_minimum) * 0.06, 0.05)
    expected_limits = (independent_minimum - padding, independent_maximum + padding)

    def inspect_limits(figure, _path):
        actual_limits = [axis.get_xlim() for axis in figure.axes]
        assert actual_limits == [pytest.approx(expected_limits), pytest.approx(expected_limits)]
        assert all(
            lower < endpoint < upper
            for lower, upper in actual_limits
            for endpoint in np.concatenate([truth - errors, truth + errors])
        )
        summary_figures.plt.close(figure)

    monkeypatch.setattr(summary_figures, "_save", inspect_limits)

    summary_figures.render_study_oof(frozen_checks, tmp_path / "study-oof.png")


def test_rendered_summary_figures_reserve_caption_and_annotation_space(
    tmp_path: Path, frozen_checks, monkeypatch: pytest.MonkeyPatch
):
    """Catch captions colliding with x labels/ticks or MAE values colliding with its title."""
    from kcat_rel import summary_figures

    original_save = summary_figures._save

    def inspect_layout(figure, path):
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        caption = figure.texts[-1].get_window_extent(renderer)
        label_boxes = []
        for axis in figure.axes:
            label_boxes.append(axis.xaxis.label.get_window_extent(renderer))
            label_boxes.extend(label.get_window_extent(renderer) for label in axis.get_xticklabels())
            title = axis.title.get_window_extent(renderer)
            assert figure.bbox.contains(title.x0, title.y0)
            assert figure.bbox.contains(title.x1, title.y1)
        assert all(not caption.overlaps(box) for box in label_boxes if box.width and box.height)
        if path.name == "held-out-error-vs-baseline.png":
            axis = figure.axes[0]
            title = axis.title.get_window_extent(renderer)
            value_labels = [text.get_window_extent(renderer) for text in axis.texts]
            assert all(not title.overlaps(box) for box in value_labels)
        original_save(figure, path)

    monkeypatch.setattr(summary_figures, "_save", inspect_layout)

    write_summary_outputs(frozen_checks, tmp_path / "summary")

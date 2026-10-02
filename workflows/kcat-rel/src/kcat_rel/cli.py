"""Explicit, label-isolated commands for the frozen ``kcat_rel`` diagnostic."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from kcat_rel.evidence import audit_evidence
from kcat_rel.features import (
    _audit_feature_table,
    _audit_manifest_contract,
    build_feature_table,
    write_feature_artifacts,
)
from kcat_rel.matrix import build_model_table
from kcat_rel.registry import FROZEN_PRIMARY_VARIANTS, load_registry, validate_reference
from kcat_rel.reporting import preflight_diagnostic_output, write_diagnostic_artifacts
from kcat_rel.sensitivity import run_perturbation_sensitivity, run_weighted_sensitivity
from kcat_rel.validation import (
    GateDecision,
    Metrics,
    evaluate_gate,
    run_outer_validation,
    select_winner,
    single_study_influence,
    study_macro_metrics,
    variant_metrics,
)


_P12883_RAW_SHA256 = "45c96586dd50e0ede770a3bb4fb0aab125d103d53b2bc0b0fed89619e07a16ad"


def _sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_pinned_fasta(path: str | Path) -> str:
    fasta = Path(path)
    if not fasta.is_file():
        raise ValueError(f"P12883 FASTA is unavailable: {fasta}")
    if _sha256_file(fasta) != _P12883_RAW_SHA256:
        raise ValueError("P12883 raw FASTA checksum mismatch")
    lines = fasta.read_text(encoding="utf-8").splitlines()
    headers = [line for line in lines if line.startswith(">")]
    if len(headers) != 1 or "P12883" not in headers[0] or "SV=5" not in headers[0]:
        raise ValueError("P12883 FASTA accession or sequence version mismatch")
    return validate_reference("".join(line.strip() for line in lines if not line.startswith(">")))


def _require_pinned_structure(path: str | Path, registry: Any) -> Path:
    """Require the sole reviewed compressed RCSB archive before feature work."""
    structure = Path(path)
    if structure.suffix != ".gz" or not structure.name.lower().endswith(".cif.gz"):
        raise ValueError("features requires the pinned compressed 8ACT .cif.gz archive")
    if not structure.is_file() or _sha256_file(structure) != registry.structure["compressed_sha256"]:
        raise ValueError("8ACT compressed checksum mismatch")
    return structure


def audit_feature_artifacts(features_path: str | Path, manifest_path: str | Path, registry: Any) -> pd.DataFrame:
    """Independently verify the feature-pair schema, output hash, and provenance."""
    csv_path, json_path = Path(features_path), Path(manifest_path)
    if not csv_path.is_file() or not json_path.is_file():
        raise ValueError("audited feature CSV and manifest must both exist")
    try:
        table = pd.read_csv(csv_path)
        manifest = json.loads(json_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, pd.errors.ParserError) as exc:
        raise ValueError("could not read feature artifact pair") from exc
    _audit_feature_table(table, registry)
    _audit_manifest_contract(manifest, table, csv_path.read_bytes(), csv_path)
    structure = manifest.get("structure", {})
    if (
        not isinstance(structure, Mapping)
        or structure.get("compressed_checksum_verified") is not True
        or structure.get("input_sha256") != registry.structure["compressed_sha256"]
        or structure.get("expected_compressed_sha256") != registry.structure["compressed_sha256"]
    ):
        raise ValueError("feature artifact lacks the exact compressed checksum-verified 8ACT provenance")
    output = manifest.get("output")
    if not isinstance(output, Mapping) or output.get("csv_sha256") != _sha256_file(csv_path):
        raise ValueError("feature artifact manifest does not match CSV bytes")
    return table


def _versions() -> dict[str, str]:
    versions = {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__}
    for distribution in ("scikit-learn", "biopython"):
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = "unavailable"
    return versions


def _git_revision() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[4], check=True,
            text=True, capture_output=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout.strip() or None


def _metrics_dict(metrics: Metrics) -> dict[str, float]:
    return {"mae": float(metrics.mae), "rmse": float(metrics.rmse)}


def _gate_dict(gate: GateDecision) -> dict[str, Any]:
    return {
        "pass_gate": gate.pass_gate,
        "failed_rules": list(gate.failed_rules),
        "mae_improvement": dict(gate.mae_improvement),
        "rmse_degradation": dict(gate.rmse_degradation),
        "model_views": {view: _metrics_dict(metric) for view, metric in gate.model_views.items()},
        "baseline_views": {view: _metrics_dict(metric) for view, metric in gate.baseline_views.items()},
    }


def _tuning_frame(results: Mapping[str, Mapping[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for model_kind in sorted(results):
        for view in ("study", "variant"):
            for fold in results[model_kind][view].folds:
                rows.append({
                    "model_kind": model_kind, "view": view, "outer_fold": fold.outer_fold,
                    "train_indices": tuple(int(value) for value in fold.train_indices),
                    "test_indices": tuple(int(value) for value in fold.test_indices),
                    "selected_params": dict(fold.selected_params), "inner_scores": dict(fold.inner_scores),
                    "scaler_mean": fold.scaler_mean, "scaler_scale": fold.scaler_scale,
                    "baseline_mean": fold.baseline_mean,
                })
    return pd.DataFrame(rows)


def _oof_frame(results: Mapping[str, Mapping[str, Any]], view: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for model_kind in sorted(results):
        frame = results[model_kind][view].oof.copy()
        frame.insert(0, "model_kind", model_kind)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _inner_fold_audit(results: Mapping[str, Mapping[str, Any]], groups: np.ndarray, revision: str | None) -> dict[str, Any]:
    """Serialize the inner folds actually retained by each primary outer result."""
    records: list[dict[str, Any]] = []
    for model_kind in sorted(results):
        for view in ("study", "variant"):
            for fold in results[model_kind][view].folds:
                for inner_fold, inner in enumerate(fold.inner_folds):
                    test = [int(value) for value in inner.test_indices]
                    held = tuple(dict.fromkeys(str(groups[index]) for index in test))
                    if len(held) != 1:
                        raise ValueError("inner fold must hold out exactly one source study")
                    records.append({
                        "model_kind": model_kind, "view": view, "outer_fold": int(fold.outer_fold),
                        "inner_fold": inner_fold, "held_out_source_id": held[0],
                        "train_indices": [int(value) for value in inner.train_indices], "test_indices": test,
                        "scaler_mean": [float(value) for value in inner.scaler_mean],
                        "scaler_scale": [float(value) for value in inner.scaler_scale],
                    })
    return {
        "schema_version": 1,
        "provenance": {"method": "recorded_during_execution", "original_execution_revision": revision, "diagnostic_rerun": False},
        "records": records,
    }


def _audit_markdown(winner: str | None, decisions: Mapping[str, GateDecision], evidence: Any) -> str:
    lines = ["# Frozen `kcat_rel` diagnostic audit", "", "| Gate | Result |", "| --- | --- |"]
    lines.append(f"| Evidence audit | PASS ({evidence.canonical_rows} labels, {evidence.canonical_measurements} measurements, {len(evidence.source_studies)} studies) |")
    for model_kind in ("ridge", "gpr"):
        decision = decisions[model_kind]
        outcome = "PASS" if decision.pass_gate else "FAIL: " + ", ".join(decision.failed_rules)
        lines.append(f"| {model_kind} primary gate | {outcome} |")
    lines.extend(["", f"Selected model: `{winner}`." if winner else "Selected model: none; negative diagnostic."])
    return "\n".join(lines) + "\n"


def _diagnostic(args: argparse.Namespace) -> int:
    # Fail before opening inputs or starting the nested fitting/sensitivity work.
    # ``write_diagnostic_artifacts`` repeats this guard just before its atomic write.
    preflight_diagnostic_output(args.output)
    registry = load_registry(args.config)
    # This audit deliberately precedes opening any target-bearing CSV.
    features = audit_feature_artifacts(args.features, args.feature_manifest, registry)
    evidence = audit_evidence(args.canonical_labels, args.measurement_evidence)
    model_table = build_model_table(features, args.canonical_labels, args.measurement_evidence)
    X = model_table.loc[:, ("variant", *registry.feature_columns)]
    y = model_table["y"].to_numpy(dtype=float)
    errors = model_table["y_error"].to_numpy(dtype=float)
    groups = model_table["source_id"].to_numpy(dtype=object)

    results: dict[str, dict[str, Any]] = {}
    decisions: dict[str, GateDecision] = {}
    for model_kind in ("ridge", "gpr"):
        views = {
            view: run_outer_validation(X, y, groups, view, model_kind, registry)
            for view in ("study", "variant")
        }
        results[model_kind] = views
        model_views = {
            "study": study_macro_metrics(y, views["study"].prediction, groups),
            "variant": variant_metrics(y, views["variant"].prediction),
        }
        baseline_views = {
            "study": study_macro_metrics(y, views["study"].baseline_prediction, groups),
            "variant": variant_metrics(y, views["variant"].baseline_prediction),
        }
        decisions[model_kind] = evaluate_gate(
            model_views, baseline_views,
            single_study_influence(y, views["study"].prediction, views["study"].baseline_prediction, groups), registry,
        )
    winner = select_winner(decisions, registry)
    sensitivities: dict[str, Any] = {}
    for model_kind, views in results.items():
        weighted = {view: run_weighted_sensitivity(result, X, y, errors, groups) for view, result in views.items()}
        perturbation = run_perturbation_sensitivity(views, X, y, errors, groups, draws=1000, seed=registry.seed)
        sensitivities[model_kind] = {
            "weighted": {view: asdict(value) for view, value in weighted.items()},
            "perturbation": asdict(perturbation),
        }
    input_hashes = {
        "features_csv": _sha256_file(args.features), "feature_manifest": _sha256_file(args.feature_manifest),
        "canonical_labels": _sha256_file(args.canonical_labels), "measurement_evidence": _sha256_file(args.measurement_evidence),
        "registry": _sha256_file(args.config),
    }
    execution_revision = _git_revision()
    inner_fold_audit = _inner_fold_audit(results, groups, execution_revision)
    run_manifest = {
        "schema_version": 1, "active_profile": registry.active_profile, "fallback_reason": registry.fallback_reason,
        "variants": list(FROZEN_PRIMARY_VARIANTS), "seed": registry.seed, "input_sha256": input_hashes,
        "tool_versions": _versions(), "command": list(sys.argv), "git_revision": execution_revision,
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "inner_fold_audit": dict(inner_fold_audit["provenance"]),
        # Reporting replaces these fixed-shape placeholders after serializing
        # every hashable output (the run manifest and SUCCESS are excluded to
        # avoid self-reference and last-write circularity).
        "artifact_sha256": {
            name: "0" * 64 for name in (
                "AUDIT.md", "metrics.json", "model-table.csv", "run.log", "sensitivity.json",
                "inner-fold-audit.json", "study-oof.csv", "tuning-records.csv", "variant-oof.csv",
            )
        },
    }
    metrics = {"winner": winner, "models": {model: _gate_dict(gate) for model, gate in decisions.items()}}
    report = {
        "run_manifest": run_manifest, "model_table": model_table, "tuning_records": _tuning_frame(results),
        "study_oof": _oof_frame(results, "study"), "variant_oof": _oof_frame(results, "variant"),
        "inner_fold_audit": inner_fold_audit,
        "metrics": metrics, "sensitivity": sensitivities, "audit_markdown": _audit_markdown(winner, decisions, evidence),
        "run_log": "diagnostic completed after evidence and feature audits\n",
    }
    write_diagnostic_artifacts(args.output, report)
    return 0


def _features(args: argparse.Namespace) -> int:
    registry = load_registry(args.config)
    reference = _read_pinned_fasta(args.reference_fasta)
    structure = _require_pinned_structure(args.structure, registry)
    result = build_feature_table(registry.variants, reference, structure, registry)
    output_dir = Path(args.output_dir)
    write_feature_artifacts(
        result, output_dir / f"{registry.active_profile}-features.csv",
        output_dir / f"{registry.active_profile}-feature-manifest.json",
    )
    return 0


def _audit_features(args: argparse.Namespace) -> int:
    registry = load_registry(args.config)
    table = audit_feature_artifacts(args.features, args.feature_manifest, registry)
    print(f"feature audit PASS: {len(table)} rows, {len(table.columns) - 1} features")
    return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="kcat-rel")
    commands = parser.add_subparsers(dest="command", required=True)
    features = commands.add_parser("features", help="build pinned, label-blind features")
    features.add_argument("--reference-fasta", required=True)
    features.add_argument("--structure", required=True)
    features.add_argument("--config", required=True)
    features.add_argument("--output-dir", required=True)
    features.set_defaults(handler=_features)
    audit = commands.add_parser("audit-features", help="verify a feature artifact pair")
    audit.add_argument("--features", required=True)
    audit.add_argument("--feature-manifest", required=True)
    audit.add_argument("--config", required=True)
    audit.set_defaults(handler=_audit_features)
    diagnostic = commands.add_parser("diagnostic", help="run the frozen diagnostic exactly as configured")
    diagnostic.add_argument("--features", required=True)
    diagnostic.add_argument("--feature-manifest", required=True)
    diagnostic.add_argument("--canonical-labels", required=True)
    diagnostic.add_argument("--measurement-evidence", required=True)
    diagnostic.add_argument("--config", required=True)
    diagnostic.add_argument("--output", required=True)
    diagnostic.set_defaults(handler=_diagnostic)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run a command, returning a non-zero status without leaving partial output."""
    try:
        args = _parser().parse_args(argv)
        return int(args.handler(args))
    except (OSError, RuntimeError, ValueError, KeyError, TypeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

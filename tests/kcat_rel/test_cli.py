from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path

import pytest

from kcat_rel.cli import _require_pinned_structure, audit_feature_artifacts, main
from kcat_rel.features import build_feature_table, write_feature_artifacts
from kcat_rel.registry import FROZEN_PRIMARY_VARIANTS, load_registry


WORKFLOW_ROOT = Path(__file__).resolve().parents[2] / "workflows/kcat-rel"
DATA_ROOT = Path(__file__).resolve().parents[2] / "data/public/kcat-rel"
CONFIG = WORKFLOW_ROOT / "config" / "kcat-rel-v0.json"


@pytest.fixture
def feature_pair(tmp_path: Path) -> tuple[Path, Path]:
    source = ast.parse((Path(__file__).with_name("test_chemistry.py")).read_text(encoding="utf-8"))
    reference = next(ast.literal_eval(node.value) for node in source.body if isinstance(node, ast.Assign))
    result = build_feature_table(
        FROZEN_PRIMARY_VARIANTS, reference, Path(__file__).with_name("fixtures") / "8act-mini.cif", load_registry(CONFIG)
    )
    csv_path, manifest_path = tmp_path / "features.csv", tmp_path / "features.json"
    write_feature_artifacts(result, csv_path, manifest_path)
    return csv_path, manifest_path


def test_cli_failure_leaves_no_success_marker(tmp_path: Path) -> None:
    """A missing audited feature artifact must fail before any run is marked complete."""
    missing_features = tmp_path / "missing-features.csv"

    result = main(
        [
            "diagnostic",
            "--features",
            str(missing_features),
            "--feature-manifest",
            str(tmp_path / "missing-features.json"),
            "--canonical-labels",
            str(DATA_ROOT / "canonical-labels.csv"),
            "--measurement-evidence",
            str(DATA_ROOT / "measurement-evidence.csv"),
            "--config",
            str(CONFIG),
            "--output",
            str(tmp_path / "diagnostic"),
        ]
    )

    assert result != 0
    assert not (tmp_path / "diagnostic" / "SUCCESS").exists()


def test_diagnostic_rejects_an_occupied_output_before_fitting(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An occupied destination must fail without beginning the costly diagnostic."""
    output = tmp_path / "diagnostic"
    output.mkdir()
    sentinel = output / "existing-result.txt"
    sentinel.write_text("preserve this result\n", encoding="utf-8")

    def fitting_must_not_run(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("occupied-output diagnostic reached outer fitting")

    monkeypatch.setattr("kcat_rel.cli.run_outer_validation", fitting_must_not_run)

    result = main(
        [
            "diagnostic",
            "--features",
            str(DATA_ROOT / "derived" / "kcat-rel-v0-no-msa-features.csv"),
            "--feature-manifest",
            str(DATA_ROOT / "derived" / "kcat-rel-v0-no-msa-feature-manifest.json"),
            "--canonical-labels",
            str(DATA_ROOT / "canonical-labels.csv"),
            "--measurement-evidence",
            str(DATA_ROOT / "measurement-evidence.csv"),
            "--config",
            str(CONFIG),
            "--output",
            str(output),
        ]
    )

    assert result != 0
    assert sentinel.read_text(encoding="utf-8") == "preserve this result\n"
    assert set(output.iterdir()) == {sentinel}


def test_feature_cli_rejects_an_uncompressed_structure_before_generation(tmp_path: Path) -> None:
    """A local CIF cannot bypass the pinned compressed 8ACT acquisition gate."""
    path = tmp_path / "8ACT.cif"
    path.write_text("data_8ACT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="compressed.*8ACT"):
        _require_pinned_structure(path, load_registry(CONFIG))


def test_feature_audit_rejects_false_compressed_verification_flag(feature_pair: tuple[Path, Path]) -> None:
    """Feature provenance must not trust a manifest that says its archive was unverified."""
    csv_path, manifest_path = feature_pair
    manifest = json.loads(manifest_path.read_text())
    manifest["structure"]["compressed_checksum_verified"] = False
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="compressed checksum"):
        audit_feature_artifacts(csv_path, manifest_path, load_registry(CONFIG))


def test_feature_audit_rejects_arbitrary_64hex_structure_hash(feature_pair: tuple[Path, Path]) -> None:
    """Merely hash-shaped structure provenance cannot substitute for the exact 8ACT pin."""
    csv_path, manifest_path = feature_pair
    manifest = json.loads(manifest_path.read_text())
    manifest["structure"]["compressed_checksum_verified"] = True
    manifest["structure"]["input_sha256"] = "a" * 64
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="exact compressed checksum"):
        audit_feature_artifacts(csv_path, manifest_path, load_registry(CONFIG))

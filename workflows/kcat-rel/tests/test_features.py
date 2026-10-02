from __future__ import annotations

import ast
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import kcat_rel.features as features_module
from kcat_rel.features import FeatureResult, build_feature_table, write_feature_artifacts
from kcat_rel.registry import ACTIVE_FEATURE_COLUMNS, FROZEN_PRIMARY_VARIANTS, load_registry


WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
CONFIG = WORKFLOW_ROOT / "config" / "kcat-rel-v0.json"
FIXTURES = Path(__file__).with_name("fixtures")


@pytest.fixture
def mini_8act() -> Path:
    return FIXTURES / "8act-mini.cif"


@pytest.fixture
def reference_sequence() -> str:
    """Load the independently frozen P12883 test input without label tables."""
    source = ast.parse((Path(__file__).with_name("test_chemistry.py")).read_text(encoding="utf-8"))
    return next(
        ast.literal_eval(node.value)
        for node in source.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "REFERENCE_SEQUENCE" for target in node.targets)
    )


def test_feature_builder_has_no_label_inputs(reference_sequence: str, mini_8act: Path) -> None:
    """Adding a label-derived column or label input must fail this isolation contract."""
    result = build_feature_table(
        variants=FROZEN_PRIMARY_VARIANTS,
        reference_sequence=reference_sequence,
        structure_path=mini_8act,
        registry=load_registry(CONFIG),
    )

    assert tuple(result.table.columns) == ("variant", *ACTIVE_FEATURE_COLUMNS)
    assert set(result.table.columns).isdisjoint({"y", "ln_kcat_rel", "source_id"})
    assert result.manifest["fallback_reason"] == "msa_toolchain_and_snapshot_unavailable"


def test_feature_builder_produces_complete_finite_ordered_artifact(
    reference_sequence: str, mini_8act: Path
) -> None:
    """Dropping/reordering a primary row or emitting non-finite features must fail."""
    result = build_feature_table(
        FROZEN_PRIMARY_VARIANTS, reference_sequence, mini_8act, load_registry(CONFIG)
    )

    assert tuple(result.table.variant) == FROZEN_PRIMARY_VARIANTS
    assert result.table.shape == (17, 4)
    assert np.isfinite(result.table.loc[:, ACTIVE_FEATURE_COLUMNS].to_numpy(dtype=float)).all()
    assert result.manifest["reference"]["normalized_sequence_sha256"] == hashlib.sha256(
        reference_sequence.encode("utf-8")
    ).hexdigest()
    assert result.manifest["structure"]["input_sha256"] == hashlib.sha256(mini_8act.read_bytes()).hexdigest()
    assert result.manifest["feature_columns"] == list(ACTIVE_FEATURE_COLUMNS)
    assert result.manifest["variants"] == list(FROZEN_PRIMARY_VARIANTS)


def test_feature_builder_rejects_a_bad_compressed_8act_checksum(
    reference_sequence: str, mini_8act: Path, tmp_path: Path
) -> None:
    """Accepting an unpinned downloaded 8ACT archive would break input provenance."""
    unpinned = tmp_path / "8ACT.cif.gz"
    unpinned.write_bytes(gzip.compress(mini_8act.read_bytes()))

    with pytest.raises(ValueError, match="8ACT compressed checksum mismatch"):
        build_feature_table(FROZEN_PRIMARY_VARIANTS, reference_sequence, unpinned, load_registry(CONFIG))


def test_artifact_writer_serializes_audited_table_and_manifest_atomically(
    reference_sequence: str, mini_8act: Path, tmp_path: Path
) -> None:
    """A writer that serializes unaudited or mismatched artifact pairs must fail."""
    result = build_feature_table(
        FROZEN_PRIMARY_VARIANTS, reference_sequence, mini_8act, load_registry(CONFIG)
    )
    csv_path = tmp_path / "features.csv"
    manifest_path = tmp_path / "features.json"

    write_feature_artifacts(result, csv_path, manifest_path)

    assert tuple(pd.read_csv(csv_path).columns) == ("variant", *ACTIVE_FEATURE_COLUMNS)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["output"]["csv_sha256"] == hashlib.sha256(csv_path.read_bytes()).hexdigest()
    assert payload["output"]["rows"] == 17


def test_artifact_writer_refuses_invalid_result_without_replacing_outputs(tmp_path: Path) -> None:
    """Replacing valid outputs after an invalid feature audit would corrupt artifacts."""
    csv_path = tmp_path / "features.csv"
    manifest_path = tmp_path / "features.json"
    csv_path.write_text("old-csv\n", encoding="utf-8")
    manifest_path.write_text("old-manifest\n", encoding="utf-8")

    invalid = FeatureResult(
        table=pd.DataFrame({"variant": ["Y115H"], "delta_charge": [0]}),
        manifest={"feature_columns": list(ACTIVE_FEATURE_COLUMNS)},
    )

    with pytest.raises(ValueError, match="schema|frozen cohort|17"):
        write_feature_artifacts(invalid, csv_path, manifest_path)

    assert csv_path.read_text(encoding="utf-8") == "old-csv\n"
    assert manifest_path.read_text(encoding="utf-8") == "old-manifest\n"


def test_artifact_writer_rejects_incomplete_or_tampered_manifest_before_writing(
    reference_sequence: str, mini_8act: Path, tmp_path: Path
) -> None:
    """Serializing a missing or unlinked provenance contract would create unauditable output."""
    result = build_feature_table(
        FROZEN_PRIMARY_VARIANTS, reference_sequence, mini_8act, load_registry(CONFIG)
    )
    csv_path = tmp_path / "features.csv"
    manifest_path = tmp_path / "features.json"
    incomplete = deepcopy(dict(result.manifest))
    del incomplete["reference"]
    tampered = deepcopy(dict(result.manifest))
    tampered["structure"]["expected_compressed_sha256"] = "0" * 64
    unlinked_output = deepcopy(dict(result.manifest))
    unlinked_output["output"] = {
        "csv_filename": csv_path.name,
        "csv_sha256": "0" * 64,
        "rows": 17,
        "columns": ["variant", *ACTIVE_FEATURE_COLUMNS],
    }

    for invalid in (incomplete, tampered, unlinked_output):
        with pytest.raises(ValueError, match="manifest|provenance|hash"):
            write_feature_artifacts(FeatureResult(result.table, invalid), csv_path, manifest_path)

    assert not csv_path.exists()
    assert not manifest_path.exists()


@pytest.mark.parametrize(
    ("path", "replacement"),
    [
        (("reference", "entry_name"), "forged"),
        (("reference", "accession"), "P00000"),
        (("reference", "sequence_version"), True),
        (("reference", "normalized_sequence_sha256"), "0" * 64),
        (("structure", "entry_id"), "forged"),
        (("structure", "expected_compressed_sha256"), "0" * 64),
        (("structure", "author_chains"), ["B", "A"]),
        (("structure", "chain_aggregation"), "forged"),
        (("motifs", "p_loop"), [178]),
        (("active_profile",), "forged"),
        (("fallback_reason",), "forged"),
    ],
)
def test_artifact_writer_rejects_tampered_frozen_identity_provenance(
    reference_sequence: str,
    mini_8act: Path,
    tmp_path: Path,
    path: tuple[str, ...],
    replacement: object,
) -> None:
    """Changing an immutable registry identity field must invalidate the artifact."""
    result = build_feature_table(
        FROZEN_PRIMARY_VARIANTS, reference_sequence, mini_8act, load_registry(CONFIG)
    )
    manifest = deepcopy(dict(result.manifest))
    destination = manifest
    for field in path[:-1]:
        destination = destination[field]
    destination[path[-1]] = replacement

    with pytest.raises(ValueError, match="manifest|provenance|profile|motif"):
        write_feature_artifacts(
            FeatureResult(result.table, manifest), tmp_path / "features.csv", tmp_path / "features.json"
        )


def test_artifact_writer_rolls_back_preexisting_pair_when_second_replace_fails(
    reference_sequence: str, mini_8act: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed manifest replacement must not leave a new CSV paired with an old manifest."""
    result = build_feature_table(
        FROZEN_PRIMARY_VARIANTS, reference_sequence, mini_8act, load_registry(CONFIG)
    )
    csv_path = tmp_path / "features.csv"
    manifest_path = tmp_path / "features.json"
    write_feature_artifacts(result, csv_path, manifest_path)
    previous_csv, previous_manifest = csv_path.read_bytes(), manifest_path.read_bytes()
    revised_table = result.table.copy()
    revised_table.loc[0, "delta_charge"] = 99
    revised = FeatureResult(revised_table, result.manifest)
    real_replace = features_module.os.replace

    def fail_new_manifest_replace(source: str | Path, destination: str | Path) -> None:
        if (
            Path(destination) == manifest_path
            and Path(source).name.startswith(".features.json.")
            and Path(source).suffix == ".tmp"
        ):
            raise OSError("fault-injected manifest replace failure")
        real_replace(source, destination)

    monkeypatch.setattr(features_module.os, "replace", fail_new_manifest_replace)
    with pytest.raises(OSError, match="fault-injected"):
        write_feature_artifacts(revised, csv_path, manifest_path)

    assert csv_path.read_bytes() == previous_csv
    assert manifest_path.read_bytes() == previous_manifest
    assert json.loads(manifest_path.read_text(encoding="utf-8"))["output"]["csv_sha256"] == hashlib.sha256(
        csv_path.read_bytes()
    ).hexdigest()


def test_artifact_writer_removes_first_csv_when_second_replace_fails(
    reference_sequence: str, mini_8act: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed first-write manifest replacement must leave neither final artifact behind."""
    result = build_feature_table(
        FROZEN_PRIMARY_VARIANTS, reference_sequence, mini_8act, load_registry(CONFIG)
    )
    csv_path = tmp_path / "features.csv"
    manifest_path = tmp_path / "features.json"
    real_replace = features_module.os.replace

    def fail_new_manifest_replace(source: str | Path, destination: str | Path) -> None:
        if (
            Path(destination) == manifest_path
            and Path(source).name.startswith(".features.json.")
            and Path(source).suffix == ".tmp"
        ):
            raise OSError("fault-injected manifest replace failure")
        real_replace(source, destination)

    monkeypatch.setattr(features_module.os, "replace", fail_new_manifest_replace)
    with pytest.raises(OSError, match="fault-injected"):
        write_feature_artifacts(result, csv_path, manifest_path)

    assert not csv_path.exists()
    assert not manifest_path.exists()

"""End-to-end feature snapshot generation without outcome leakage."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lsar_na_rel.features import (
    FEATURE_COLUMNS,
    assert_snapshot_matches,
    build_feature_snapshot,
    main,
    write_feature_snapshot,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_ROOT = REPO_ROOT / "workflows" / "lsar-na-rel"
PUBLIC = REPO_ROOT / "data" / "public" / "lsar-na-rel"
PRIMARY = PUBLIC / "na-rel-primary.csv"
CONFIG = WORKFLOW_ROOT / "config" / "lsar-na-rel.json"
FASTA = REPO_ROOT / "tests" / "lsar_na_rel" / "fixtures" / "P12883.fasta"
ASSEMBLY = WORKFLOW_ROOT / "cache" / "8ACT-assembly1.pdb"
APPROVED = PUBLIC / "features.csv"


@pytest.fixture(scope="module")
def snapshot():
    return build_feature_snapshot(
        primary_path=PRIMARY,
        fasta_path=FASTA,
        assembly_pdb_path=ASSEMBLY,
        config_path=CONFIG,
    )


def test_feature_snapshot_matches_approved_values_and_order(snapshot):
    approved = pd.read_csv(APPROVED)
    actual = snapshot.features
    assert list(actual.columns) == list(FEATURE_COLUMNS)
    assert actual["variant"].tolist() == approved["variant"].tolist()
    assert_snapshot_matches(actual, approved)


def test_feature_provenance_is_portable_and_complete(snapshot):
    text = json.dumps(snapshot.provenance, sort_keys=True)
    assert "/Users/" not in text
    assert "research-workhorse" not in text
    assert snapshot.provenance["contact_cutoff_angstrom"] == 4.5
    assert snapshot.provenance["sequence_accession"] == "P12883"
    assert snapshot.provenance["assembly_id"] == "8ACT"
    assert len(snapshot.provenance["fasta_sha256"]) == 64
    assert len(snapshot.provenance["pdb_sha256"]) == 64


def test_write_feature_snapshot_roundtrips_atomically(snapshot, tmp_path: Path):
    features_path = tmp_path / "nested" / "features.csv"
    provenance_path = tmp_path / "nested" / "features-provenance.json"
    write_feature_snapshot(snapshot, features_path, provenance_path)

    loaded = pd.read_csv(features_path)
    assert_snapshot_matches(loaded, snapshot.features)
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["features_sha256"]
    assert not list(features_path.parent.glob("*.partial"))


def test_snapshot_comparison_rejects_numeric_drift(snapshot):
    changed = snapshot.features.copy()
    changed.loc[0, "wt_site_rsa_8act_mean_ab"] += 1e-6
    with pytest.raises(ValueError, match="numeric mismatch"):
        assert_snapshot_matches(snapshot.features, changed)


def test_snapshot_comparison_rejects_schema_drift(snapshot):
    changed = snapshot.features.drop(columns=["coverage_b"])
    with pytest.raises(ValueError, match="schema mismatch"):
        assert_snapshot_matches(snapshot.features, changed)


def test_snapshot_has_no_targets_or_source_groups(snapshot):
    forbidden = {"Na_rel", "ln(Na_rel)", "y", "source_id"}
    assert forbidden.isdisjoint(snapshot.features.columns)
    numeric = snapshot.features.drop(columns=["variant"]).select_dtypes(include=[np.number])
    assert np.isfinite(numeric.to_numpy(dtype=float)).all()


def test_check_mode_rejects_provenance_hash_drift(snapshot, tmp_path: Path, monkeypatch):
    features_path = tmp_path / "features.csv"
    provenance_path = tmp_path / "features-provenance.json"
    write_feature_snapshot(snapshot, features_path, provenance_path)
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["features_sha256"] = "0" * 64
    provenance_path.write_text(json.dumps(provenance), encoding="utf-8")

    class Assembly:
        pdb_path = ASSEMBLY

    monkeypatch.setattr("lsar_na_rel.features.ensure_8act_assembly", lambda path: Assembly())
    monkeypatch.setattr("lsar_na_rel.features.build_feature_snapshot", lambda **kwargs: snapshot)

    with pytest.raises(ValueError, match="provenance features_sha256 mismatch"):
        main(
            [
                "--primary",
                str(PRIMARY),
                "--fasta",
                str(FASTA),
                "--config",
                str(CONFIG),
                "--cache",
                str(tmp_path / "cache"),
                "--features-output",
                str(features_path),
                "--provenance-output",
                str(provenance_path),
                "--check",
            ]
        )


def test_check_mode_rejects_generation_provenance_drift(snapshot, tmp_path: Path, monkeypatch):
    features_path = tmp_path / "features.csv"
    provenance_path = tmp_path / "features-provenance.json"
    write_feature_snapshot(snapshot, features_path, provenance_path)
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["contact_cutoff_angstrom"] = 5.0
    provenance_path.write_text(json.dumps(provenance), encoding="utf-8")

    class Assembly:
        pdb_path = ASSEMBLY

    monkeypatch.setattr("lsar_na_rel.features.ensure_8act_assembly", lambda path: Assembly())
    monkeypatch.setattr("lsar_na_rel.features.build_feature_snapshot", lambda **kwargs: snapshot)

    with pytest.raises(ValueError, match="provenance drift.*contact_cutoff_angstrom"):
        main(
            [
                "--features-output",
                str(features_path),
                "--provenance-output",
                str(provenance_path),
                "--check",
            ]
        )

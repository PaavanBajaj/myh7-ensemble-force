"""Manifest fields for the operator-blocked FoldX gate."""

from __future__ import annotations

from pathlib import Path

from lsar_na_rel.manifests import (
    FOLDX_FALLBACK_REASON,
    FOLDX_GATE,
    build_run_manifest,
    sha256_file,
    write_run_manifest,
)


def test_sha256_known_bytes(tmp_path: Path):
    path = tmp_path / "blob.txt"
    path.write_bytes(b"abc")
    # echo -n abc | shasum -a 256
    assert sha256_file(path) == (
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )


def test_foldx_gate_constants_are_durable_public_status():
    assert FOLDX_GATE == "unavailable_no_values"
    assert FOLDX_FALLBACK_REASON == (
        "No FoldX-derived values are present in this release; the registered "
        "public result uses delta_charge plus checksummed 8ACT features."
    )
    assert "session" not in FOLDX_FALLBACK_REASON.lower()


def test_build_run_manifest_records_blocked_foldx_and_ladder_three_profile(tmp_path: Path):
    primary = tmp_path / "p.csv"
    labels = tmp_path / "l.csv"
    features = tmp_path / "f.csv"
    primary.write_text("a\n", encoding="utf-8")
    labels.write_text("b\n", encoding="utf-8")
    features.write_text("c\n", encoding="utf-8")
    payload = build_run_manifest(
        run_id="ladder-3-nofoldx",
        feature_columns=("delta_charge", "wt_site_rsa_8act_mean_ab", "ihm_interface_flag_8act"),
        seed=20260811,
        splitter_name="LeaveOneGroupOut",
        primary_path=primary,
        labels_path=labels,
        features_path=features,
        repo_root=tmp_path,
        command="python -m lsar_na_rel.cli talk-v0 --run-id ladder-3-nofoldx",
        timestamp="20260816T000000Z",
        n_rows=16,
        n_folds=6,
        git_sha_value="deadbeef",
        package_versions_value={"python": "3.12.13", "scikit-learn": "1.9.0"},
        git_dirty_value=True,
        workflow_content_sha256_value="f" * 64,
    )
    assert payload["foldx_gate"] == "unavailable_no_values"
    assert payload["fallback_reason"] == FOLDX_FALLBACK_REASON
    assert payload["profile"] == "talk-v0-no-foldx"
    assert payload["modeling_family"] == "chemistry_plus_8act_features_csv"
    assert payload["primary_sha256"] == sha256_file(primary)
    assert payload["run_id"] == "ladder-3-nofoldx"
    assert payload["git_dirty"] is True
    assert payload["workflow_content_sha256"] == "f" * 64
    out = tmp_path / "run-manifest.json"
    write_run_manifest(out, payload)
    text = out.read_text(encoding="utf-8")
    assert "unavailable_no_values" in text


def test_diagnostic_runs_use_diagnostic_ladder_profile(tmp_path: Path):
    dummy = tmp_path / "x.csv"
    dummy.write_text("z\n", encoding="utf-8")
    payload = build_run_manifest(
        run_id="ladder-1-rsa",
        feature_columns=("wt_site_rsa_8act_mean_ab",),
        seed=20260811,
        splitter_name="LeaveOneGroupOut",
        primary_path=dummy,
        labels_path=dummy,
        features_path=dummy,
        repo_root=tmp_path,
        command="x",
        timestamp="t",
        n_rows=16,
        n_folds=6,
        git_sha_value="abc",
        package_versions_value={"python": "3.12.13"},
        git_dirty_value=False,
        workflow_content_sha256_value="a" * 64,
    )
    assert payload["profile"] == "diagnostic-ladder"

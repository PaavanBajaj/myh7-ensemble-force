"""Provenance must agree in both copies and retain the frozen execution history."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from kcat_rel.reporting import write_diagnostic_artifacts
from .test_reporting_producers import _producer_result


@pytest.mark.parametrize("side", ["audit", "manifest", "both"])
@pytest.mark.parametrize("mutation", ["method_recorded", "revision", "rerun", "method_unknown"])
def test_frozen_provenance_cannot_be_mismatched_or_relabelled(tmp_path: Path, side: str, mutation: str) -> None:
    result = _producer_result()
    audit = result["inner_fold_audit"]["provenance"]
    manifest = result["run_manifest"]["post_execution_provenance_repair"]["inner_fold_audit"]
    targets = {"audit": [audit], "manifest": [manifest], "both": [audit, manifest]}[side]
    for target in targets:
        if mutation == "method_recorded": target["method"] = "recorded_during_execution"
        elif mutation == "revision": target["original_execution_revision"] = "e06e8f74ea7934a2752f607f969820480b0a6787"
        elif mutation == "rerun": target["diagnostic_rerun"] = True
        else: target["method"] = "recorded_at_execution"
    with pytest.raises(ValueError):
        write_diagnostic_artifacts(tmp_path, result)
    assert not any(tmp_path.iterdir())


@pytest.mark.parametrize("mutation", ["run_revision", "repair_revision", "missing_method_copy", "changed_scientific_files"])
def test_reconstruction_remains_bound_to_original_execution(tmp_path: Path, mutation: str) -> None:
    result = _producer_result()
    manifest = result["run_manifest"]
    repair = manifest["post_execution_provenance_repair"]
    if mutation == "run_revision": manifest["git_revision"] = "e06e8f74ea7934a2752f607f969820480b0a6787"
    elif mutation == "repair_revision": repair["execution_git_head"] = "e06e8f74ea7934a2752f607f969820480b0a6787"
    elif mutation == "missing_method_copy": repair.pop("inner_fold_audit")
    else: repair["inner_fold_audit"]["scientific_result_files_changed"] = True
    with pytest.raises(ValueError):
        write_diagnostic_artifacts(tmp_path, result)
    assert not any(tmp_path.iterdir())


@pytest.mark.parametrize("mutation", ["valid", "missing", "method", "revision", "rerun", "audit_revision", "frozen_relabel"])
def test_recorded_execution_requires_agreeing_manifest_provenance(tmp_path: Path, mutation: str) -> None:
    result = _producer_result()
    manifest = result["run_manifest"]
    manifest.pop("post_execution_provenance_repair")
    revision = "e06e8f74ea7934a2752f607f969820480b0a6787"
    manifest["git_revision"] = revision
    provenance = {"method": "recorded_during_execution", "original_execution_revision": revision, "diagnostic_rerun": False}
    result["inner_fold_audit"]["provenance"] = dict(provenance)
    manifest["inner_fold_audit"] = dict(provenance)
    if mutation == "missing": manifest.pop("inner_fold_audit")
    elif mutation == "method": manifest["inner_fold_audit"]["method"] = "post_execution_reconstruction"
    elif mutation == "revision": manifest["inner_fold_audit"]["original_execution_revision"] = "0" * 40
    elif mutation == "rerun": manifest["inner_fold_audit"]["diagnostic_rerun"] = True
    elif mutation == "audit_revision": result["inner_fold_audit"]["provenance"]["original_execution_revision"] = "0" * 40
    elif mutation == "frozen_relabel":
        frozen_revision = "91affadb0f5bdb3fbf8f09b22209836f38ddf9c5"
        manifest["git_revision"] = frozen_revision
        manifest["inner_fold_audit"]["original_execution_revision"] = frozen_revision
        result["inner_fold_audit"]["provenance"]["original_execution_revision"] = frozen_revision
    if mutation == "valid":
        write_diagnostic_artifacts(tmp_path, result)
        assert (tmp_path / "SUCCESS").read_text() == "completed\n"
        written = json.loads((tmp_path / "run-manifest.json").read_text())
        assert written["inner_fold_audit"]["method"] == "recorded_during_execution"
    else:
        with pytest.raises(ValueError):
            write_diagnostic_artifacts(tmp_path, result)
        assert not any(tmp_path.iterdir())

"""Run manifests: hashes, versions, FoldX blocked gate, command, timestamp."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from importlib import metadata
from pathlib import Path
from typing import Mapping, Sequence

FOLDX_GATE = "unavailable_no_values"
FOLDX_FALLBACK_REASON = (
    "No FoldX-derived values are present in this release; the registered "
    "public result uses delta_charge plus checksummed 8ACT features."
)
MODELING_FAMILY = "chemistry_plus_8act_features_csv"
_PACKAGES = ("scikit-learn", "numpy", "scipy", "pandas", "matplotlib")


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_sha(repo_root: str | Path) -> str:
    output = subprocess.check_output(
        ["git", "rev-parse", "HEAD"],
        cwd=str(repo_root),
        text=True,
    )
    return output.strip()


def git_dirty(repo_root: str | Path) -> bool:
    output = subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=all"],
        cwd=str(repo_root),
        text=True,
    )
    return bool(output.strip())


def workflow_content_sha256(repo_root: str | Path) -> str:
    """Hash the runnable source/config/input state, including untracked promotion files."""
    root = Path(repo_root).resolve()
    candidates = [
        root / "environment.yml",
        root / "environment-osx-arm64.lock",
        root / "workflows" / "lsar-na-rel" / "config" / "lsar-na-rel.json",
        root / "data" / "public" / "lsar-na-rel" / "na-rel-primary.csv",
        root / "data" / "public" / "lsar-na-rel" / "labels.csv",
        root / "data" / "public" / "lsar-na-rel" / "features.csv",
    ]
    candidates.extend(sorted((root / "src" / "lsar_na_rel").glob("*.py")))
    # Private pre-promotion layout.
    candidates.extend(
        sorted((root / "workflows" / "lsar-na-rel" / "src" / "lsar_na_rel").glob("*.py"))
    )
    digest = hashlib.sha256()
    found = False
    for path in sorted({path for path in candidates if path.is_file()}):
        found = True
        try:
            label = path.relative_to(root).as_posix()
        except ValueError:
            label = path.name
        digest.update(label.encode("utf-8"))
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    if not found:
        raise FileNotFoundError(f"no workflow content found under {root}")
    return digest.hexdigest()


def package_versions() -> dict[str, str]:
    versions = {"python": sys.version.split()[0]}
    for name in _PACKAGES:
        versions[name] = metadata.version(name)
    return versions


def _profile_for_run(run_id: str) -> str:
    if run_id == "ladder-3-nofoldx":
        return "talk-v0-no-foldx"
    return "diagnostic-ladder"


def build_run_manifest(
    *,
    run_id: str,
    feature_columns: Sequence[str],
    seed: int,
    splitter_name: str,
    primary_path: str | Path,
    labels_path: str | Path,
    features_path: str | Path,
    repo_root: str | Path,
    command: str,
    timestamp: str,
    n_rows: int,
    n_folds: int,
    git_sha_value: str | None = None,
    package_versions_value: Mapping[str, str] | None = None,
    git_dirty_value: bool | None = None,
    workflow_content_sha256_value: str | None = None,
) -> dict[str, object]:
    return {
        "run_id": run_id,
        "profile": _profile_for_run(run_id),
        "feature_columns": list(feature_columns),
        "seed": int(seed),
        "foldx_gate": FOLDX_GATE,
        "fallback_reason": FOLDX_FALLBACK_REASON,
        "modeling_family": MODELING_FAMILY,
        "splitter_name": splitter_name,
        "primary_sha256": sha256_file(primary_path),
        "labels_sha256": sha256_file(labels_path),
        "features_sha256": sha256_file(features_path),
        "git_sha": git_sha_value if git_sha_value is not None else git_sha(repo_root),
        "git_dirty": (
            bool(git_dirty_value) if git_dirty_value is not None else git_dirty(repo_root)
        ),
        "workflow_content_sha256": (
            workflow_content_sha256_value
            if workflow_content_sha256_value is not None
            else workflow_content_sha256(repo_root)
        ),
        "package_versions": dict(
            package_versions_value if package_versions_value is not None else package_versions()
        ),
        "command": command,
        "timestamp": timestamp,
        "n_rows": int(n_rows),
        "n_folds": int(n_folds),
    }


def write_run_manifest(path: str | Path, payload: Mapping[str, object]) -> None:
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

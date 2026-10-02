"""Label-blind assembly and serialization of the frozen ``kcat_rel`` features."""

from __future__ import annotations

import hashlib
import gzip
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import tempfile
from copy import deepcopy
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from kcat_rel.chemistry import generate_chemistry
from kcat_rel.registry import (
    ACTIVE_FEATURE_COLUMNS,
    ACTIVE_PROFILE,
    FROZEN_PRIMARY_VARIANTS,
    MOTIF_POSITIONS,
    REFERENCE_ACCESSION,
    REFERENCE_SEQUENCE_LENGTH,
    REFERENCE_SEQUENCE_SHA256,
    FrozenRegistry,
    validate_reference,
)
from kcat_rel.structure import generate_catalytic_distances, load_ca_coordinates, validate_8act_structure


_TABLE_COLUMNS = ("variant", *ACTIVE_FEATURE_COLUMNS)
_REFERENCE_ENTRY_NAME = "MYH7_HUMAN"
_REFERENCE_SEQUENCE_VERSION = 5
_STRUCTURE_ENTRY_ID = "8ACT"
_STRUCTURE_COMPRESSED_SHA256 = "69bc81d3b6ee01a09e300a0e9b65f5a02837bbda9394e5566b248bf612a1b5e0"
_STRUCTURE_AUTHOR_CHAINS = ["A", "B"]
_STRUCTURE_CHAIN_AGGREGATION = "arithmetic_mean_ab"


@dataclass(frozen=True)
class FeatureResult:
    """An independently auditable, target-free feature table and provenance record."""

    table: pd.DataFrame
    manifest: Mapping[str, Any]


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_revision() -> str | None:
    """Return the checked-out revision when the source is available in Git."""
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            cwd=Path(__file__).resolve().parents[2],
            text=True,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout.strip() or None


def _package_versions() -> dict[str, str]:
    versions = {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__}
    for distribution, name in (("biopython", "biopython"),):
        try:
            versions[name] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "unavailable"
    return versions


def _audit_feature_table(table: pd.DataFrame, registry: FrozenRegistry) -> pd.DataFrame:
    """Enforce the exact active-profile schema before target data can enter."""
    if registry.variants != FROZEN_PRIMARY_VARIANTS or registry.feature_columns != ACTIVE_FEATURE_COLUMNS:
        raise ValueError("registry does not match the frozen primary cohort and active profile")
    if not isinstance(table, pd.DataFrame) or tuple(table.columns) != _TABLE_COLUMNS:
        raise ValueError("feature table schema does not match the active profile")
    if len(table) != len(FROZEN_PRIMARY_VARIANTS) or tuple(table["variant"]) != FROZEN_PRIMARY_VARIANTS:
        raise ValueError("feature table does not match the frozen cohort and order")
    if table["variant"].duplicated().any():
        raise ValueError("feature table must contain one-to-one frozen cohort rows")
    numeric = table.loc[:, ACTIVE_FEATURE_COLUMNS].to_numpy(dtype=float)
    if numeric.shape != (len(FROZEN_PRIMARY_VARIANTS), len(ACTIVE_FEATURE_COLUMNS)) or not np.isfinite(numeric).all():
        raise ValueError("feature table contains non-finite active-profile values")
    return table.copy()


def _is_sha256(value: object) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(character in "0123456789abcdef" for character in value)


def _manifest_error(message: str) -> ValueError:
    return ValueError(f"feature manifest {message}")


def _require_mapping(payload: Mapping[str, Any], field: str) -> Mapping[str, Any]:
    value = payload.get(field)
    if not isinstance(value, Mapping):
        raise _manifest_error(f"lacks required {field} provenance")
    return value


def _audit_manifest_contract(
    candidate: Mapping[str, Any], table: pd.DataFrame, csv_bytes: bytes, csv_path: Path
) -> dict[str, Any]:
    """Validate provenance and its table/output linkage before any destination changes."""
    if not isinstance(candidate, Mapping):
        raise _manifest_error("must be a mapping")
    manifest = deepcopy(dict(candidate))
    if manifest.get("schema_version") != 1:
        raise _manifest_error("has an unsupported schema version")
    if manifest.get("active_profile") != ACTIVE_PROFILE:
        raise _manifest_error("has an unexpected active profile")
    if manifest.get("fallback_reason") != "msa_toolchain_and_snapshot_unavailable":
        raise _manifest_error("has an unexpected fallback reason")
    if manifest.get("variants") != list(FROZEN_PRIMARY_VARIANTS):
        raise _manifest_error("does not preserve the ordered frozen cohort")
    if manifest.get("feature_columns") != list(ACTIVE_FEATURE_COLUMNS):
        raise _manifest_error("does not preserve the active-profile schema")
    definitions = _require_mapping(manifest, "feature_definitions")
    if set(definitions) != set(ACTIVE_FEATURE_COLUMNS) or any(
        not isinstance(value, str) or not value.strip() for value in definitions.values()
    ):
        raise _manifest_error("has incomplete feature definitions")

    reference = _require_mapping(manifest, "reference")
    if (
        reference.get("accession") != REFERENCE_ACCESSION
        or reference.get("entry_name") != _REFERENCE_ENTRY_NAME
        or type(reference.get("sequence_version")) is not int
        or reference.get("sequence_version") != _REFERENCE_SEQUENCE_VERSION
        or reference.get("normalized_sequence_length") != REFERENCE_SEQUENCE_LENGTH
        or reference.get("normalized_sequence_sha256") != REFERENCE_SEQUENCE_SHA256
        or reference.get("expected_normalized_sequence_sha256") != REFERENCE_SEQUENCE_SHA256
    ):
        raise _manifest_error("has invalid P12883 reference provenance")

    structure = _require_mapping(manifest, "structure")
    if (
        structure.get("entry_id") != _STRUCTURE_ENTRY_ID
        or structure.get("author_chains") != _STRUCTURE_AUTHOR_CHAINS
        or structure.get("chain_aggregation") != _STRUCTURE_CHAIN_AGGREGATION
        or structure.get("expected_compressed_sha256") != _STRUCTURE_COMPRESSED_SHA256
        or not isinstance(structure.get("input_filename"), str)
        or not structure["input_filename"]
        or not _is_sha256(structure.get("input_sha256"))
        or not isinstance(structure.get("compressed_checksum_verified"), bool)
    ):
        raise _manifest_error("has invalid 8ACT structure provenance")
    if manifest.get("motifs") != {name: list(positions) for name, positions in MOTIF_POSITIONS.items()}:
        raise _manifest_error("does not preserve the frozen motif definition")

    versions = _require_mapping(manifest, "tool_versions")
    if any(not isinstance(versions.get(name), str) or not versions[name] for name in ("python", "numpy", "pandas", "biopython")):
        raise _manifest_error("has incomplete tool provenance")
    if not isinstance(manifest.get("command"), list) or not all(isinstance(value, str) for value in manifest["command"]):
        raise _manifest_error("has invalid command provenance")
    if manifest.get("git_revision") is not None and not isinstance(manifest.get("git_revision"), str):
        raise _manifest_error("has invalid Git provenance")
    timestamp = manifest.get("timestamp_utc")
    if not isinstance(timestamp, str):
        raise _manifest_error("lacks timestamp provenance")
    try:
        datetime.fromisoformat(timestamp)
    except ValueError as exc:
        raise _manifest_error("has invalid timestamp provenance") from exc
    audit = _require_mapping(manifest, "audit")
    if audit.get("rows") != len(FROZEN_PRIMARY_VARIANTS) or audit.get("finite") is not True or audit.get("schema") != list(_TABLE_COLUMNS):
        raise _manifest_error("audit does not match the frozen feature table")

    output = manifest.get("output")
    if output is not None:
        if not isinstance(output, Mapping) or (
            output.get("csv_filename") != csv_path.name
            or output.get("csv_sha256") != _sha256_bytes(csv_bytes)
            or output.get("rows") != len(table)
            or output.get("columns") != list(_TABLE_COLUMNS)
        ):
            raise _manifest_error("output hash is not linked to the table bytes")
    return manifest


def _feature_manifest(
    table: pd.DataFrame,
    reference: str,
    structure_path: Path,
    registry: FrozenRegistry,
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "active_profile": registry.active_profile,
        "fallback_reason": registry.fallback_reason,
        "variants": list(FROZEN_PRIMARY_VARIANTS),
        "feature_columns": list(ACTIVE_FEATURE_COLUMNS),
        "feature_definitions": {
            "delta_charge": "q(mutant) - q(WT); D/E=-1, K/R=+1, H=0, other standard residues=0",
            "grantham_distance": "Original symmetric integer Grantham (1974) Table 2 lookup",
            "catalytic_motif_min_ca_distance_8act_mean_ab": (
                "Mean of deposited human 8ACT author chains A and B minimum C-alpha distances "
                "from each WT variant site to P-loop 178-185, Switch I 238-246, and Switch II 461-466"
            ),
        },
        "reference": {
            "accession": registry.reference["accession"],
            "entry_name": registry.reference["entry_name"],
            "sequence_version": registry.reference["sequence_version"],
            "normalized_sequence_length": len(reference),
            "normalized_sequence_sha256": _sha256_bytes(reference.encode("utf-8")),
            "expected_normalized_sequence_sha256": registry.reference["normalized_sequence_sha256"],
        },
        "structure": {
            "entry_id": registry.structure["entry_id"],
            "author_chains": list(registry.structure["author_chains"]),
            "chain_aggregation": registry.structure["chain_aggregation"],
            "input_filename": structure_path.name,
            "input_sha256": _sha256_file(structure_path),
            "expected_compressed_sha256": registry.structure["compressed_sha256"],
            "compressed_checksum_verified": structure_path.suffix == ".gz",
        },
        "motifs": {name: list(positions) for name, positions in registry.motifs.items()},
        "tool_versions": _package_versions(),
        "command": list(sys.argv),
        "git_revision": _git_revision(),
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "audit": {"rows": len(table), "finite": True, "schema": list(_TABLE_COLUMNS)},
    }


def _load_pinned_coordinates(structure_path: Path, registry: FrozenRegistry):
    """Verify and temporarily unpack the checksum-pinned RCSB archive when supplied."""
    if structure_path.suffix != ".gz":
        return load_ca_coordinates(structure_path)
    if _sha256_file(structure_path) != registry.structure["compressed_sha256"]:
        raise ValueError("8ACT compressed checksum mismatch")
    descriptor, temporary_name = tempfile.mkstemp(prefix="8ACT.", suffix=".cif")
    temporary = Path(temporary_name)
    try:
        with gzip.open(structure_path, "rb") as source, os.fdopen(descriptor, "wb") as destination:
            destination.write(source.read())
        return load_ca_coordinates(temporary)
    except OSError as exc:
        raise ValueError(f"could not decompress pinned 8ACT archive: {structure_path}") from exc
    finally:
        temporary.unlink(missing_ok=True)


def build_feature_table(
    variants: Sequence[str],
    reference_sequence: str,
    structure_path: str | Path,
    registry: FrozenRegistry,
) -> FeatureResult:
    """Build and audit the sole 17-row feature artifact without reading evidence."""
    labels = tuple(variants)
    if labels != FROZEN_PRIMARY_VARIANTS:
        raise ValueError("variants do not match the frozen primary cohort")
    if registry.active_profile != "kcat-rel-v0-no-msa":
        raise ValueError("registry does not select the approved active profile")
    reference = validate_reference(reference_sequence)
    cif_path = Path(structure_path)
    if not cif_path.is_file():
        raise ValueError(f"8ACT structure input is unavailable: {cif_path}")

    coordinates = _load_pinned_coordinates(cif_path, registry)
    validate_8act_structure(coordinates, reference, labels)
    chemistry = generate_chemistry(labels, reference)
    geometry = generate_catalytic_distances(labels, coordinates).loc[:, [
        "variant", "catalytic_motif_min_ca_distance_8act_mean_ab"
    ]]
    table = chemistry.merge(geometry, on="variant", how="inner", validate="one_to_one", sort=False)
    table = table.loc[:, _TABLE_COLUMNS]
    audited = _audit_feature_table(table, registry)
    manifest = _feature_manifest(audited, reference, cif_path, registry)
    return FeatureResult(table=audited, manifest=manifest)


def _atomic_write(path: Path, data: bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        return temporary
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _reserve_backup(path: Path) -> Path:
    """Reserve a sibling pathname used only while replacing one artifact pair."""
    descriptor, backup_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".backup", dir=path.parent)
    os.close(descriptor)
    return Path(backup_name)


def _restore_pair(
    csv_path: Path,
    manifest_path: Path,
    csv_backup: Path | None,
    manifest_backup: Path | None,
    csv_replaced: bool,
    manifest_replaced: bool,
) -> None:
    """Restore the former matched pair after a replacement failure."""
    for path, backup, replaced in (
        (manifest_path, manifest_backup, manifest_replaced),
        (csv_path, csv_backup, csv_replaced),
    ):
        if backup is not None:
            path.unlink(missing_ok=True)
            os.replace(backup, path)
        elif replaced:
            path.unlink(missing_ok=True)


def write_feature_artifacts(result: FeatureResult, csv_path: str | Path, manifest_path: str | Path) -> None:
    """Validate and replace an artifact pair, rolling back both paths on failure."""
    if not isinstance(result, FeatureResult):
        raise ValueError("feature result must be an audited FeatureResult")
    # The table has already been audited by the builder; use the frozen registry-shaped
    # contract here again so callers cannot serialize a subsequently mutated DataFrame.
    class _ArtifactRegistry:
        variants = FROZEN_PRIMARY_VARIANTS
        feature_columns = ACTIVE_FEATURE_COLUMNS

    table = _audit_feature_table(result.table, _ArtifactRegistry())
    target_csv, target_manifest = Path(csv_path), Path(manifest_path)
    if target_csv == target_manifest:
        raise ValueError("feature CSV and manifest paths must differ")
    csv_bytes = table.to_csv(index=False, lineterminator="\n").encode("utf-8")
    manifest = _audit_manifest_contract(result.manifest, table, csv_bytes, target_csv)
    manifest["output"] = {
        "csv_filename": target_csv.name,
        "csv_sha256": _sha256_bytes(csv_bytes),
        "rows": len(table),
        "columns": list(_TABLE_COLUMNS),
    }
    manifest = _audit_manifest_contract(manifest, table, csv_bytes, target_csv)
    manifest_bytes = (json.dumps(manifest, sort_keys=True, indent=2) + "\n").encode("utf-8")
    temporary_csv: Path | None = None
    temporary_manifest: Path | None = None
    pending_csv_backup: Path | None = None
    pending_manifest_backup: Path | None = None
    csv_backup: Path | None = None
    manifest_backup: Path | None = None
    csv_replaced = False
    manifest_replaced = False
    try:
        temporary_csv = _atomic_write(target_csv, csv_bytes)
        temporary_manifest = _atomic_write(target_manifest, manifest_bytes)
        if target_csv.exists():
            pending_csv_backup = _reserve_backup(target_csv)
            os.replace(target_csv, pending_csv_backup)
            csv_backup, pending_csv_backup = pending_csv_backup, None
        if target_manifest.exists():
            pending_manifest_backup = _reserve_backup(target_manifest)
            os.replace(target_manifest, pending_manifest_backup)
            manifest_backup, pending_manifest_backup = pending_manifest_backup, None
        os.replace(temporary_csv, target_csv)
        temporary_csv = None
        csv_replaced = True
        os.replace(temporary_manifest, target_manifest)
        temporary_manifest = None
        manifest_replaced = True
    except BaseException:
        _restore_pair(
            target_csv,
            target_manifest,
            csv_backup,
            manifest_backup,
            csv_replaced,
            manifest_replaced,
        )
        csv_backup = None
        manifest_backup = None
        raise
    finally:
        if temporary_csv is not None:
            temporary_csv.unlink(missing_ok=True)
        if temporary_manifest is not None:
            temporary_manifest.unlink(missing_ok=True)
        if pending_csv_backup is not None:
            pending_csv_backup.unlink(missing_ok=True)
        if pending_manifest_backup is not None:
            pending_manifest_backup.unlink(missing_ok=True)
        if csv_backup is not None:
            csv_backup.unlink(missing_ok=True)
        if manifest_backup is not None:
            manifest_backup.unlink(missing_ok=True)

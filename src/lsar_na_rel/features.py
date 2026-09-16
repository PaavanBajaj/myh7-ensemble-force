"""Generate the public, outcome-free LSAR Na-rel feature snapshot."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

from lsar_na_rel.chemistry import generate_delta_charge_features, load_reference_sequence
from lsar_na_rel.config import load_talk_v0_config
from lsar_na_rel.paths import DEFAULT_PATHS
from lsar_na_rel.structure import (
    ASSEMBLY_URL,
    BH_CHAIN,
    DECOMPRESSED_PDB_SHA256,
    FH_CHAIN,
    GZIP_SHA256,
    MOTOR_LEVER_RANGES,
    PROXIMAL_S2_RANGES,
    ensure_8act_assembly,
    generate_structure_features,
)


FEATURE_COLUMNS: tuple[str, ...] = (
    "variant",
    "delta_charge",
    "wt_site_rsa_8act_a",
    "wt_site_rsa_8act_b",
    "wt_site_rsa_8act_mean_ab",
    "coverage_a",
    "coverage_b",
    "wt_site_rsa_8act_bh",
    "wt_site_rsa_8act_fh",
    "ihm_contact_bh_fh",
    "ihm_contact_bh_s2",
    "ihm_interface_flag_8act",
)


@dataclass(frozen=True)
class FeatureSnapshot:
    features: pd.DataFrame
    provenance: Mapping[str, object]


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: str | Path) -> str:
    return _sha256_bytes(Path(path).read_bytes())


def _feature_csv_bytes(features: pd.DataFrame) -> bytes:
    ordered = features.loc[:, list(FEATURE_COLUMNS)]
    return ordered.to_csv(index=False, float_format="%.12g").encode("utf-8")


def build_feature_snapshot(
    *,
    primary_path: str | Path,
    fasta_path: str | Path,
    assembly_pdb_path: str | Path,
    config_path: str | Path,
    mkdssp: str = "mkdssp",
) -> FeatureSnapshot:
    """Build features using variant identities only; never load outcome columns."""
    primary_path = Path(primary_path)
    fasta_path = Path(fasta_path)
    config_path = Path(config_path)
    variants = pd.read_csv(primary_path, usecols=["variant"])["variant"].astype(str).tolist()
    config = load_talk_v0_config(config_path)
    sequence = load_reference_sequence(fasta_path)
    chemistry = generate_delta_charge_features(variants, reference_sequence=sequence)
    structure = generate_structure_features(
        variants,
        assembly_pdb_path=assembly_pdb_path,
        mkdssp=mkdssp,
        contact_cutoff=config.contact_cutoff_angstrom,
    )

    site = structure.site_provenance.copy()
    model = structure.features.copy()
    merged = chemistry.features.merge(site, on="variant", how="left", validate="one_to_one")
    merged = merged.merge(
        model.loc[:, ["variant", "wt_site_rsa_8act_mean_ab", "ihm_interface_flag_8act"]],
        on="variant",
        how="left",
        validate="one_to_one",
    )
    merged["coverage_a"] = merged["coverage_a"].astype(bool)
    merged["coverage_b"] = merged["coverage_b"].astype(bool)
    features = merged.loc[:, list(FEATURE_COLUMNS)].reset_index(drop=True)
    if features["variant"].tolist() != variants:
        raise ValueError("generated feature order differs from primary variant order")
    if features.isna().any().any():
        raise ValueError("generated features contain missing values")

    provenance: dict[str, object] = {
        "schema_version": 1,
        "assembly_id": "8ACT",
        "assembly_url": ASSEMBLY_URL,
        "pdb_sha256": DECOMPRESSED_PDB_SHA256,
        "gzip_sha256": GZIP_SHA256,
        "sequence_accession": "P12883",
        "sequence_refseq": "NP_000248.2",
        "sequence_length": len(sequence),
        "fasta_sha256": _sha256_file(fasta_path),
        "primary_sha256": _sha256_file(primary_path),
        "config_sha256": _sha256_file(config_path),
        "contact_cutoff_angstrom": config.contact_cutoff_angstrom,
        "bh_chain": BH_CHAIN,
        "fh_chain": FH_CHAIN,
        "motor_lever_ranges": MOTOR_LEVER_RANGES,
        "proximal_s2_ranges": PROXIMAL_S2_RANGES,
        "dssp_executable": Path(mkdssp).name,
        "dssp_acc_array": "Sander",
        "feature_columns": list(FEATURE_COLUMNS),
        "feature_generation_uses_outcomes": False,
    }
    return FeatureSnapshot(features=features, provenance=provenance)


def assert_snapshot_matches(
    actual: pd.DataFrame,
    expected: pd.DataFrame,
    *,
    atol: float = 1e-12,
) -> None:
    if list(actual.columns) != list(expected.columns):
        raise ValueError(
            f"feature schema mismatch: {list(actual.columns)!r} != {list(expected.columns)!r}"
        )
    if actual["variant"].astype(str).tolist() != expected["variant"].astype(str).tolist():
        raise ValueError("feature variant order mismatch")
    for column in actual.columns[1:]:
        left = actual[column]
        right = expected[column]
        if pd.api.types.is_bool_dtype(left) or pd.api.types.is_bool_dtype(right):
            if left.astype(bool).tolist() != right.astype(bool).tolist():
                raise ValueError(f"feature value mismatch in {column}")
            continue
        try:
            left_values = pd.to_numeric(left, errors="raise").to_numpy(dtype=float)
            right_values = pd.to_numeric(right, errors="raise").to_numpy(dtype=float)
        except (TypeError, ValueError):
            if left.astype(str).tolist() != right.astype(str).tolist():
                raise ValueError(f"feature value mismatch in {column}")
            continue
        if not np.allclose(left_values, right_values, rtol=0.0, atol=atol):
            delta = float(np.max(np.abs(left_values - right_values)))
            raise ValueError(f"feature numeric mismatch in {column}; max |delta|={delta}")


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".partial", dir=str(path.parent)
    )
    temp_path = Path(temporary)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp_path, path)
    except Exception:
        temp_path.unlink(missing_ok=True)
        raise


def write_feature_snapshot(
    snapshot: FeatureSnapshot,
    features_path: str | Path,
    provenance_path: str | Path,
) -> None:
    features_path = Path(features_path)
    provenance_path = Path(provenance_path)
    csv_bytes = _feature_csv_bytes(snapshot.features)
    provenance = dict(snapshot.provenance)
    provenance["features_sha256"] = _sha256_bytes(csv_bytes)
    json_bytes = (json.dumps(provenance, indent=2, sort_keys=True) + "\n").encode("utf-8")
    _atomic_write(features_path, csv_bytes)
    _atomic_write(provenance_path, json_bytes)


def verify_feature_provenance(
    features_path: str | Path,
    provenance_path: str | Path,
    expected_provenance: Mapping[str, object] | None = None,
) -> None:
    features_path = Path(features_path)
    provenance_path = Path(provenance_path)
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    expected = provenance.get("features_sha256")
    observed = _sha256_file(features_path)
    if expected != observed:
        raise ValueError(
            f"provenance features_sha256 mismatch: {expected!r} != {observed!r}"
        )
    if expected_provenance is not None:
        normalized = json.loads(json.dumps(dict(expected_provenance), sort_keys=True))
        for key, expected_value in normalized.items():
            observed_value = provenance.get(key)
            if observed_value != expected_value:
                raise ValueError(
                    f"provenance drift for {key}: {observed_value!r} != {expected_value!r}"
                )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="lsar_na_rel.features")
    parser.add_argument("--primary", type=Path, default=DEFAULT_PATHS.primary_path)
    parser.add_argument("--fasta", type=Path, default=DEFAULT_PATHS.reference_fasta)
    parser.add_argument("--config", type=Path, default=DEFAULT_PATHS.config_path)
    parser.add_argument("--cache", type=Path, default=DEFAULT_PATHS.cache_dir)
    parser.add_argument("--features-output", type=Path, default=DEFAULT_PATHS.features_path)
    parser.add_argument(
        "--provenance-output", type=Path, default=DEFAULT_PATHS.features_provenance_path
    )
    parser.add_argument("--mkdssp", default="mkdssp")
    parser.add_argument("--check", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)
    assembly = ensure_8act_assembly(args.cache)
    snapshot = build_feature_snapshot(
        primary_path=args.primary,
        fasta_path=args.fasta,
        assembly_pdb_path=assembly.pdb_path,
        config_path=args.config,
        mkdssp=args.mkdssp,
    )
    if args.check:
        committed = pd.read_csv(args.features_output)
        assert_snapshot_matches(snapshot.features, committed)
        verify_feature_provenance(
            args.features_output,
            args.provenance_output,
            expected_provenance=snapshot.provenance,
        )
        print(f"features match {args.features_output} rows={len(committed)}")
        return 0
    write_feature_snapshot(snapshot, args.features_output, args.provenance_output)
    print(f"wrote {args.features_output} rows={len(snapshot.features)}")
    print(f"wrote {args.provenance_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

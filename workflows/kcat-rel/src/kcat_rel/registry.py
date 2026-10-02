"""Immutable, label-blind inputs for the MYH7 ``kcat_rel`` diagnostic."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping


STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
FROZEN_PRIMARY_VARIANTS = (
    "Y115H",
    "D239N",
    "R249Q",
    "H251N",
    "R403Q",
    "R453C",
    "I457T",
    "E497D",
    "R663H",
    "R719W",
    "G741R",
    "G768R",
    "D778V",
    "L781P",
    "S782N",
    "A797T",
    "F834L",
)
FROZEN_SOURCE_STUDIES = (
    "adhikari2016",
    "adhikari2019",
    "kawana2017",
    "morck2022",
    "nag2015",
    "nandwani2025",
    "pathak2026",
    "sarkar2020",
    "sommese2013",
)
ACTIVE_PROFILE = "kcat-rel-v0-no-msa"
ACTIVE_FEATURE_COLUMNS = (
    "delta_charge",
    "grantham_distance",
    "catalytic_motif_min_ca_distance_8act_mean_ab",
)
REFERENCE_ACCESSION = "P12883"
REFERENCE_SEQUENCE_LENGTH = 1935
REFERENCE_SEQUENCE_SHA256 = "e98e5d01a359820bcd1a00c3422025d9fedf1fdbc46ab586e3f5b731f4585be6"
MOTIF_POSITIONS = {
    "p_loop": (178, 179, 180, 181, 182, 183, 184, 185),
    "switch_i": (238, 239, 240, 241, 242, 243, 244, 245, 246),
    "switch_ii": (461, 462, 463, 464, 465, 466),
}


@dataclass(frozen=True)
class ParsedVariant:
    label: str
    wt: str
    position: int
    mutant: str


@dataclass(frozen=True)
class FrozenRegistry:
    variants: tuple[str, ...]
    source_studies: tuple[str, ...]
    active_profile: str
    fallback_reason: str
    feature_columns: tuple[str, ...]
    reference: Mapping[str, Any]
    structure: Mapping[str, Any]
    motifs: Mapping[str, tuple[int, ...]]
    model_grids: Mapping[str, Any]
    seed: int
    gate_thresholds: Mapping[str, Any]


def _frozen_payload() -> dict[str, Any]:
    """Return the complete configuration contract; a file may not amend it."""
    return {
        "name": "kcat-rel-v0",
        "variants": list(FROZEN_PRIMARY_VARIANTS),
        "source_studies": list(FROZEN_SOURCE_STUDIES),
        "active_profile": ACTIVE_PROFILE,
        "fallback_reason": "msa_toolchain_and_snapshot_unavailable",
        "feature_columns": list(ACTIVE_FEATURE_COLUMNS),
        "reference": {
            "accession": REFERENCE_ACCESSION,
            "entry_name": "MYH7_HUMAN",
            "sequence_version": 5,
            "normalized_sequence_length": REFERENCE_SEQUENCE_LENGTH,
            "normalized_sequence_sha256": REFERENCE_SEQUENCE_SHA256,
            "fasta_url": "https://rest.uniprot.org/uniprotkb/P12883.fasta",
        },
        "structure": {
            "entry_id": "8ACT",
            "url": "https://files.rcsb.org/download/8ACT.cif.gz",
            "compressed_sha256": "69bc81d3b6ee01a09e300a0e9b65f5a02837bbda9394e5566b248bf612a1b5e0",
            "author_chains": ["A", "B"],
            "chain_aggregation": "arithmetic_mean_ab",
        },
        "motifs": {name: list(positions) for name, positions in MOTIF_POSITIONS.items()},
        "model_grids": {
            "ridge_alphas": [0.01, 0.1, 1.0, 10.0, 100.0],
            "gpr": {
                "signal_amplitude": 1.0,
                "kernel": "matern_3_2",
                "length_scales": [0.5, 1.0, 2.0],
                "noise_variances": [0.0025, 0.01, 0.04],
                "optimizer": None,
                "normalize_y": False,
            },
        },
        "seed": 20260916,
        "gate_thresholds": {
            "minimum_mae_improvement_fraction": 0.05,
            "maximum_rmse_degradation_fraction": 0.05,
            "winner_effective_tie_fraction": 0.02,
            "single_study_influence_gain_strictly_positive": True,
        },
    }


def parse_variant(label: str) -> ParsedVariant:
    """Parse one canonical, positive-position amino-acid substitution."""
    if not isinstance(label, str):
        raise ValueError(f"invalid variant: {label!r}")
    match = re.fullmatch(r"([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])", label)
    if match is None:
        raise ValueError(f"invalid variant: {label!r}")
    return ParsedVariant(label, match[1], int(match[2]), match[3])


def validate_reference(sequence: str) -> str:
    """Normalize and checksum-validate the sole permitted MYH7 reference."""
    if not isinstance(sequence, str):
        raise ValueError("P12883 reference sequence must be text")
    normalized = "".join(sequence.split()).upper()
    if len(normalized) != REFERENCE_SEQUENCE_LENGTH:
        raise ValueError("P12883 reference sequence has an unexpected length")
    if any(residue not in STANDARD_AMINO_ACIDS for residue in normalized):
        raise ValueError("P12883 reference sequence contains non-standard residues")
    digest = hashlib.sha256(normalized.encode("utf-8")).hexdigest()
    if digest != REFERENCE_SEQUENCE_SHA256:
        raise ValueError("P12883 reference sequence checksum mismatch")
    return normalized


def load_registry(path: str | Path) -> FrozenRegistry:
    """Load only the exact reviewed registry; reject any local amendment."""
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid frozen registry: {path}") from exc
    expected = _frozen_payload()
    if not isinstance(payload, dict):
        raise ValueError("invalid frozen registry")
    if payload.get("variants") != expected["variants"]:
        raise ValueError("registry does not match frozen primary cohort")
    if payload != expected:
        raise ValueError("registry differs from the frozen kcat-rel-v0 contract")
    return FrozenRegistry(
        variants=FROZEN_PRIMARY_VARIANTS,
        source_studies=FROZEN_SOURCE_STUDIES,
        active_profile=ACTIVE_PROFILE,
        fallback_reason=expected["fallback_reason"],
        feature_columns=ACTIVE_FEATURE_COLUMNS,
        reference=expected["reference"],
        structure=expected["structure"],
        motifs=MOTIF_POSITIONS,
        model_grids=expected["model_grids"],
        seed=expected["seed"],
        gate_thresholds=expected["gate_thresholds"],
    )

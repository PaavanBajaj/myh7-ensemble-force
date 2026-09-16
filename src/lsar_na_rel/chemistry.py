"""Frozen formal-charge convention and ``delta_charge`` (Talk-v0 / #28).

Charge convention (issue #26 approved-with-edits; do not retune vs labels):

- Asp/Glu (D/E) = -1
- Lys/Arg (K/R) = +1
- His (H) = 0 (formal neutral for the no-FoldX profile)
- all other standard residues = 0

Feature definition: ``delta_charge = q_mut - q_WT``.

Identity: every variant is validated against UniProt P12883 / RefSeq
NP_000248.2 by sequence lookup and WT match at the stated position.
Mismatches and non-standard residue codes are hard errors (no imputation).

Volume and hydropathy columns are intentionally not generated in Batch 1.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd

from lsar_na_rel.inputs import parse_variant

# Frozen formal charges at the approved Talk-v0 physiological convention.
# His is deliberately neutral (0), not +1 / +0.5.
CHARGE_CONVENTION: dict[str, int] = {
    "D": -1,
    "E": -1,
    "K": 1,
    "R": 1,
    "H": 0,
}

STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

CHARGE_CONVENTION_CITATION = (
    "Talk-v0 formal-charge freeze (GitHub issue #26 approved-with-edits; "
    "plan §6.3 / research §3.1): D/E=-1, K/R=+1, H=0, else 0; "
    "delta_charge = q_mut - q_WT. His remains formally neutral for the "
    "global no-FoldX profile and must not be tuned against labels."
)


@dataclass(frozen=True)
class ChemistryResult:
    """Chemistry feature table plus provenance for the no-FoldX profile."""

    features: pd.DataFrame
    provenance: Mapping[str, object]


def formal_charge(residue: str) -> int:
    code = str(residue).strip().upper()
    if len(code) != 1 or code not in STANDARD_AMINO_ACIDS:
        raise ValueError(f"Non-standard amino-acid residue code: {residue!r}")
    return int(CHARGE_CONVENTION.get(code, 0))


def load_reference_sequence(fasta_path: str | Path) -> str:
    path = Path(fasta_path)
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Expected FASTA header in {path}")
    sequence = "".join(line.strip() for line in lines[1:] if line.strip())
    if not sequence:
        raise ValueError(f"Empty sequence in {path}")
    if any(ch not in STANDARD_AMINO_ACIDS for ch in sequence.upper()):
        bad = sorted({ch for ch in sequence.upper() if ch not in STANDARD_AMINO_ACIDS})
        raise ValueError(f"Non-standard residue(s) in reference sequence: {bad}")
    return sequence.upper()


def _assert_variant_identity(label: str, reference_sequence: str):
    parsed = parse_variant(label)
    if parsed.wt not in STANDARD_AMINO_ACIDS or parsed.mut not in STANDARD_AMINO_ACIDS:
        raise ValueError(
            f"Non-standard amino-acid residue in variant {label!r} "
            f"(wt={parsed.wt!r}, mut={parsed.mut!r}); do not impute"
        )
    if parsed.position < 1 or parsed.position > len(reference_sequence):
        raise ValueError(
            f"Variant {label!r} position {parsed.position} is out of range for "
            f"reference sequence of length {len(reference_sequence)} (identity)"
        )
    observed = reference_sequence[parsed.position - 1]
    if observed != parsed.wt:
        raise ValueError(
            f"WT identity mismatch for {label!r}: claimed WT {parsed.wt!r} but "
            f"reference has {observed!r} at position {parsed.position}"
        )
    return parsed


def generate_delta_charge_features(
    variant_labels: Sequence[str],
    *,
    reference_sequence: str,
    sequence_accession: str = "P12883",
    sequence_refseq: str = "NP_000248.2",
) -> ChemistryResult:
    """Compute ``delta_charge`` for each variant without loading outcomes.

    Does not read labels/targets. Does not emit volume or hydropathy columns.
    """
    if not reference_sequence:
        raise ValueError("reference_sequence must be non-empty")

    rows: list[dict[str, object]] = []
    for label in variant_labels:
        parsed = _assert_variant_identity(str(label), reference_sequence)
        delta = formal_charge(parsed.mut) - formal_charge(parsed.wt)
        rows.append({"variant": parsed.label, "delta_charge": int(delta)})

    features = pd.DataFrame(rows, columns=["variant", "delta_charge"])
    provenance: dict[str, object] = {
        "sequence_accession": sequence_accession,
        "sequence_refseq": sequence_refseq,
        "sequence_length": len(reference_sequence),
        "charge_convention": dict(CHARGE_CONVENTION),
        "default_charge_other_residues": 0,
        "delta_charge_definition": "q_mut - q_WT",
        "citation": CHARGE_CONVENTION_CITATION,
        "feature_columns": ("delta_charge",),
        "excluded_model_columns": ("delta_volume", "delta_hydropathy"),
    }
    return ChemistryResult(features=features, provenance=provenance)

"""Frozen formal-charge and Grantham Table-2 substitution descriptors."""

from __future__ import annotations

from numbers import Integral
from typing import Sequence

import pandas as pd

from kcat_rel.registry import (
    FROZEN_PRIMARY_VARIANTS,
    STANDARD_AMINO_ACIDS,
    parse_variant,
    validate_reference,
)


FORMAL_CHARGES = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0}

# Grantham (1974), Table 2 original integer distances.  This published order
# is retained with the literal matrix so its transcription can be audited.
GRANTHAM_AMINO_ACIDS = ("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
GRANTHAM_TABLE = (
    (0, 112, 111, 126, 195, 91, 107, 60, 86, 94, 96, 106, 84, 113, 27, 99, 58, 148, 112, 64),
    (112, 0, 86, 96, 180, 43, 54, 125, 29, 97, 102, 26, 91, 97, 103, 110, 71, 101, 77, 96),
    (111, 86, 0, 23, 139, 46, 42, 80, 68, 149, 153, 94, 142, 158, 91, 46, 65, 174, 143, 133),
    (126, 96, 23, 0, 154, 61, 45, 94, 81, 168, 172, 101, 160, 177, 108, 65, 85, 181, 160, 152),
    (195, 180, 139, 154, 0, 154, 170, 159, 174, 198, 198, 202, 196, 205, 169, 112, 149, 215, 194, 192),
    (91, 43, 46, 61, 154, 0, 29, 87, 24, 109, 113, 53, 101, 116, 76, 68, 42, 130, 99, 96),
    (107, 54, 42, 45, 170, 29, 0, 98, 40, 134, 138, 56, 126, 140, 93, 80, 65, 152, 122, 121),
    (60, 125, 80, 94, 159, 87, 98, 0, 98, 135, 138, 127, 127, 153, 42, 56, 59, 184, 147, 109),
    (86, 29, 68, 81, 174, 24, 40, 98, 0, 94, 99, 32, 87, 100, 77, 89, 47, 115, 83, 84),
    (94, 97, 149, 168, 198, 109, 134, 135, 94, 0, 5, 102, 10, 21, 95, 142, 89, 61, 33, 29),
    (96, 102, 153, 172, 198, 113, 138, 138, 99, 5, 0, 107, 15, 22, 98, 145, 92, 61, 36, 32),
    (106, 26, 94, 101, 202, 53, 56, 127, 32, 102, 107, 0, 95, 102, 103, 121, 78, 110, 85, 97),
    (84, 91, 142, 160, 196, 101, 126, 127, 87, 10, 15, 95, 0, 28, 87, 127, 81, 67, 36, 21),
    (113, 97, 158, 177, 205, 116, 140, 153, 100, 21, 22, 102, 28, 0, 114, 155, 103, 40, 22, 50),
    (27, 103, 91, 108, 169, 76, 93, 42, 77, 95, 98, 103, 87, 114, 0, 74, 38, 147, 110, 68),
    (99, 110, 46, 65, 112, 68, 80, 56, 89, 142, 145, 121, 127, 155, 74, 0, 58, 177, 144, 124),
    (58, 71, 65, 85, 149, 42, 65, 59, 47, 89, 92, 78, 81, 103, 38, 58, 0, 128, 92, 69),
    (148, 101, 174, 181, 215, 130, 152, 184, 115, 61, 61, 110, 67, 40, 147, 177, 128, 0, 37, 88),
    (112, 77, 143, 160, 194, 99, 122, 147, 83, 33, 36, 85, 36, 22, 110, 144, 92, 37, 0, 55),
    (64, 96, 133, 152, 192, 96, 121, 109, 84, 29, 32, 97, 21, 50, 68, 124, 69, 88, 55, 0),
)


def _validate_grantham_table() -> None:
    size = len(GRANTHAM_AMINO_ACIDS)
    if len(set(GRANTHAM_AMINO_ACIDS)) != size or set(GRANTHAM_AMINO_ACIDS) != STANDARD_AMINO_ACIDS:
        raise RuntimeError("Grantham matrix residue order is invalid")
    if len(GRANTHAM_TABLE) != size or any(len(row) != size for row in GRANTHAM_TABLE):
        raise RuntimeError("Grantham matrix must be 20 by 20")
    for i, row in enumerate(GRANTHAM_TABLE):
        for j, value in enumerate(row):
            if not isinstance(value, Integral):
                raise RuntimeError("Grantham matrix values must be integers")
            if value < 0 or value != GRANTHAM_TABLE[j][i]:
                raise RuntimeError("Grantham matrix must be symmetric and non-negative")
        if row[i] != 0:
            raise RuntimeError("Grantham matrix diagonal must be zero")
    index = {residue: i for i, residue in enumerate(GRANTHAM_AMINO_ACIDS)}
    if GRANTHAM_TABLE[index["I"]][index["L"]] != 5:
        raise RuntimeError("Grantham I/L anchor is invalid")
    if GRANTHAM_TABLE[index["C"]][index["W"]] != 215:
        raise RuntimeError("Grantham C/W anchor is invalid")


_validate_grantham_table()
_GRANTHAM_INDEX = {residue: i for i, residue in enumerate(GRANTHAM_AMINO_ACIDS)}


def _canonical_residue(residue: str) -> str:
    if not isinstance(residue, str):
        raise ValueError(f"non-standard amino-acid residue: {residue!r}")
    code = residue.strip().upper()
    if len(code) != 1 or code not in STANDARD_AMINO_ACIDS:
        raise ValueError(f"non-standard amino-acid residue: {residue!r}")
    return code


def formal_charge(residue: str) -> int:
    """Return the frozen formal-residue charge (histidine is neutral)."""
    return FORMAL_CHARGES.get(_canonical_residue(residue), 0)


def grantham_distance(wt: str, mutant: str) -> int:
    """Look up one original Grantham (1974) Table-2 distance."""
    return int(GRANTHAM_TABLE[_GRANTHAM_INDEX[_canonical_residue(wt)]][_GRANTHAM_INDEX[_canonical_residue(mutant)]])


def generate_chemistry(variants: Sequence[str], reference_sequence: str) -> pd.DataFrame:
    """Build the complete, ordered, label-blind chemistry feature table."""
    labels = tuple(variants)
    if labels != FROZEN_PRIMARY_VARIANTS:
        raise ValueError("variants do not match the frozen primary cohort")
    reference = validate_reference(reference_sequence)
    rows: list[dict[str, int | str]] = []
    for label in labels:
        variant = parse_variant(label)
        observed = reference[variant.position - 1]
        if observed != variant.wt:
            raise ValueError(
                f"WT identity mismatch for {variant.label}: expected {variant.wt}, "
                f"P12883 has {observed} at position {variant.position}"
            )
        rows.append(
            {
                "variant": variant.label,
                "delta_charge": formal_charge(variant.mutant) - formal_charge(variant.wt),
                "grantham_distance": grantham_distance(variant.wt, variant.mutant),
            }
        )
    return pd.DataFrame(rows, columns=["variant", "delta_charge", "grantham_distance"])

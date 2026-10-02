"""Strict, label-blind C-alpha geometry from deposited human 8ACT."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from kcat_rel.registry import FROZEN_PRIMARY_VARIANTS, MOTIF_POSITIONS, parse_variant, validate_reference


AUTHOR_CHAINS = ("A", "B")
_THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
_MOTIF_POSITIONS = tuple(position for motif in MOTIF_POSITIONS.values() for position in motif)
_REQUIRED_POSITIONS = frozenset(
    [*(parse_variant(label).position for label in FROZEN_PRIMARY_VARIANTS), *_MOTIF_POSITIONS]
)


@dataclass(frozen=True)
class ResidueCoordinate:
    """One unambiguous, author-numbered C-alpha record."""

    residue: str
    coordinate: np.ndarray


@dataclass(frozen=True)
class CoordinateSet:
    """Selected author-chain C-alpha coordinates keyed by P12883 position."""

    chains: Mapping[str, Mapping[int, ResidueCoordinate]]


def _required_column(data: Mapping[str, list[str]], name: str) -> list[str]:
    try:
        return data[name]
    except KeyError as exc:
        raise ValueError(f"8ACT mmCIF lacks required column {name}") from exc


def _parse_float(value: str, *, description: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {description}: {value!r}") from exc
    if not np.isfinite(parsed):
        raise ValueError(f"invalid {description}: {value!r}")
    return parsed


def load_ca_coordinates(cif_path: str | Path) -> CoordinateSet:
    """Load unambiguous C-alphas for the frozen 8ACT site and motif positions.

    The loader deliberately uses author chain and residue identifiers, because
    P12883 is the numbering authority for the downstream identity gate.
    """
    try:
        data = MMCIF2Dict(str(cif_path))
    except (OSError, ValueError) as exc:
        raise ValueError(f"could not read 8ACT mmCIF: {cif_path}") from exc

    columns = {
        name: _required_column(data, name)
        for name in (
            "_atom_site.group_PDB",
            "_atom_site.label_atom_id",
            "_atom_site.auth_asym_id",
            "_atom_site.auth_seq_id",
            "_atom_site.auth_comp_id",
            "_atom_site.label_alt_id",
            "_atom_site.occupancy",
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.pdbx_PDB_model_num",
        )
    }
    lengths = {len(values) for values in columns.values()}
    if len(lengths) != 1:
        raise ValueError("8ACT mmCIF atom-site columns have inconsistent lengths")

    candidates: dict[tuple[str, int], list[tuple[float, str, np.ndarray]]] = {}
    row_count = next(iter(lengths), 0)
    for index in range(row_count):
        if (
            columns["_atom_site.pdbx_PDB_model_num"][index] != "1"
            or columns["_atom_site.group_PDB"][index] != "ATOM"
            or columns["_atom_site.label_atom_id"][index] != "CA"
        ):
            continue
        chain = columns["_atom_site.auth_asym_id"][index]
        if chain not in AUTHOR_CHAINS:
            continue
        try:
            position = int(columns["_atom_site.auth_seq_id"][index])
        except ValueError as exc:
            raise ValueError(f"invalid author residue number in chain {chain}") from exc
        if position not in _REQUIRED_POSITIONS:
            continue
        residue = _THREE_TO_ONE.get(columns["_atom_site.auth_comp_id"][index])
        if residue is None:
            raise ValueError(f"non-canonical residue in chain {chain} at position {position}")
        coordinate = np.asarray(
            [
                _parse_float(columns["_atom_site.Cartn_x"][index], description="C-alpha x coordinate"),
                _parse_float(columns["_atom_site.Cartn_y"][index], description="C-alpha y coordinate"),
                _parse_float(columns["_atom_site.Cartn_z"][index], description="C-alpha z coordinate"),
            ],
            dtype=float,
        )
        occupancy = _parse_float(columns["_atom_site.occupancy"][index], description="C-alpha occupancy")
        if occupancy < 0:
            raise ValueError(f"invalid C-alpha occupancy for chain {chain} at position {position}")
        candidates.setdefault((chain, position), []).append((occupancy, residue, coordinate))

    chains: dict[str, dict[int, ResidueCoordinate]] = {chain: {} for chain in AUTHOR_CHAINS}
    for (chain, position), options in candidates.items():
        identities = {residue for _, residue, _ in options}
        if len(identities) != 1:
            raise ValueError(f"ambiguous residue identity in chain {chain} at position {position}")
        maximum = max(occupancy for occupancy, _, _ in options)
        best = [(residue, coordinate) for occupancy, residue, coordinate in options if occupancy == maximum]
        if len(best) != 1:
            raise ValueError(f"ambiguous maximum-occupancy C-alpha in chain {chain} at position {position}")
        residue, coordinate = best[0]
        chains[chain][position] = ResidueCoordinate(residue=residue, coordinate=coordinate)
    return CoordinateSet(chains=chains)


def _parse_requested_variants(variants: Sequence[str]) -> tuple:
    labels = tuple(variants)
    if not labels:
        raise ValueError("at least one frozen primary variant is required")
    if len(set(labels)) != len(labels) or any(label not in FROZEN_PRIMARY_VARIANTS for label in labels):
        raise ValueError("variants do not match the frozen primary cohort")
    return tuple(parse_variant(label) for label in labels)


def _coordinate_for(coordinates: CoordinateSet, chain: str, position: int, label: str) -> ResidueCoordinate:
    try:
        return coordinates.chains[chain][position]
    except KeyError as exc:
        raise ValueError(f"missing C-alpha coordinate for chain {chain} at {label}") from exc


def validate_8act_structure(coordinates: CoordinateSet, reference_sequence: str, variants: tuple[str, ...]) -> None:
    """Require complete, identity-matching 8ACT coordinates before feature use."""
    reference = validate_reference(reference_sequence)
    if variants != FROZEN_PRIMARY_VARIANTS:
        raise ValueError("variants do not match the frozen primary cohort")
    parsed_variants = tuple(parse_variant(label) for label in FROZEN_PRIMARY_VARIANTS)
    for chain in AUTHOR_CHAINS:
        for variant in parsed_variants:
            record = _coordinate_for(coordinates, chain, variant.position, variant.label)
            expected = reference[variant.position - 1]
            if expected != variant.wt:
                raise ValueError(
                    f"WT identity mismatch for {variant.label}: expected {variant.wt}, P12883 has {expected}"
                )
            if record.residue != expected:
                raise ValueError(
                    f"residue identity mismatch for chain {chain} at {variant.label}: "
                    f"expected {expected}, 8ACT has {record.residue}"
                )
        for position in _MOTIF_POSITIONS:
            record = _coordinate_for(coordinates, chain, position, f"motif position {position}")
            expected = reference[position - 1]
            if record.residue != expected:
                raise ValueError(
                    f"residue identity mismatch for chain {chain} at motif position {position}: "
                    f"expected {expected}, 8ACT has {record.residue}"
                )


def minimum_ca_distance(site: np.ndarray, motif: np.ndarray) -> float:
    """Return the finite minimum Euclidean distance from one site to a motif."""
    distances = np.linalg.norm(motif - site, axis=1)
    value = float(distances.min())
    if not np.isfinite(value) or value < 0:
        raise ValueError("invalid C-alpha distance")
    return value


def generate_catalytic_distances(variants: tuple[str, ...], coordinates: CoordinateSet) -> pd.DataFrame:
    """Calculate the chain-specific and mean catalytic-motif C-alpha distances."""
    parsed_variants = _parse_requested_variants(variants)
    rows: list[dict[str, float | str]] = []
    for variant in parsed_variants:
        row: dict[str, float | str] = {"variant": variant.label}
        values: list[float] = []
        for chain in AUTHOR_CHAINS:
            site = _coordinate_for(coordinates, chain, variant.position, variant.label).coordinate
            motif = np.asarray(
                [_coordinate_for(coordinates, chain, position, f"motif position {position}").coordinate for position in _MOTIF_POSITIONS],
                dtype=float,
            )
            value = minimum_ca_distance(site, motif)
            row[f"catalytic_motif_min_ca_distance_8act_{chain}"] = value
            values.append(value)
        row["catalytic_motif_min_ca_distance_8act_mean_ab"] = float(np.mean(values))
        rows.append(row)
    return pd.DataFrame(
        rows,
        columns=[
            "variant",
            "catalytic_motif_min_ca_distance_8act_A",
            "catalytic_motif_min_ca_distance_8act_B",
            "catalytic_motif_min_ca_distance_8act_mean_ab",
        ],
    )

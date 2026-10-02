from __future__ import annotations

import ast
from pathlib import Path

import pytest

from kcat_rel.registry import FROZEN_PRIMARY_VARIANTS
from kcat_rel.structure import generate_catalytic_distances, load_ca_coordinates, validate_8act_structure


FIXTURES = Path(__file__).with_name("fixtures")


@pytest.fixture
def mini_8act() -> Path:
    return FIXTURES / "8act-mini.cif"


@pytest.fixture
def mini_8act_missing_b(mini_8act: Path, tmp_path: Path) -> Path:
    missing = tmp_path / "8act-mini-missing-b.cif"
    missing.write_text(
        mini_8act.read_text(encoding="utf-8").replace(
            "ATOM 25 C CA . TYR B 1 115 ? 0.0 0.0 0.0 1.00 0.0 ? 115 TYR B CA 1",
            "ATOM 25 C CB . TYR B 1 115 ? 0.0 0.0 0.0 1.00 0.0 ? 115 TYR B CB 1",
        ),
        encoding="utf-8",
    )
    return missing


@pytest.fixture
def reference_sequence() -> str:
    """Load the independently frozen P12883 test input without label tables."""
    source = ast.parse((Path(__file__).with_name("test_chemistry.py")).read_text(encoding="utf-8"))
    return next(
        ast.literal_eval(node.value)
        for node in source.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "REFERENCE_SEQUENCE" for target in node.targets)
    )


def test_minimum_distance_and_chain_mean_are_exact(mini_8act: Path) -> None:
    """Changing the chain-specific minimum or mean calculation must fail this test."""
    table = generate_catalytic_distances(("Y115H",), load_ca_coordinates(mini_8act))

    row = table.iloc[0]
    assert row.catalytic_motif_min_ca_distance_8act_A == pytest.approx(3.0)
    assert row.catalytic_motif_min_ca_distance_8act_B == pytest.approx(5.0)
    assert row.catalytic_motif_min_ca_distance_8act_mean_ab == pytest.approx(4.0)


def test_structure_rejects_missing_site_in_either_chain(mini_8act_missing_b: Path) -> None:
    """Omitting a required author-chain C-alpha must invalidate geometry generation."""
    with pytest.raises(ValueError, match=r"chain B.*Y115"):
        generate_catalytic_distances(("Y115H",), load_ca_coordinates(mini_8act_missing_b))


def test_validation_accepts_the_complete_frozen_cohort(mini_8act: Path, reference_sequence: str) -> None:
    """Removing any frozen site, motif, or identity gate must fail this audit."""
    validate_8act_structure(load_ca_coordinates(mini_8act), reference_sequence, FROZEN_PRIMARY_VARIANTS)


def test_validation_rejects_a_proper_primary_cohort_subset(mini_8act: Path, reference_sequence: str) -> None:
    """Permitting a subset would weaken the required 17-site structure gate."""
    with pytest.raises(ValueError, match="frozen primary cohort"):
        validate_8act_structure(load_ca_coordinates(mini_8act), reference_sequence, ("Y115H",))


def test_validation_rejects_a_wt_identity_mismatch(mini_8act: Path, reference_sequence: str, tmp_path: Path) -> None:
    """Changing a deposited WT identity must invalidate the whole structure."""
    mutated = tmp_path / "8act-mini-wt-identity-mismatch.cif"
    mutated.write_text(
        mini_8act.read_text(encoding="utf-8")
        .replace("? 115 TYR A CA 1", "? 115 ALA A CA 1"),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match=r"chain A.*Y115"):
        validate_8act_structure(load_ca_coordinates(mutated), reference_sequence, FROZEN_PRIMARY_VARIANTS)


def test_validation_rejects_a_motif_identity_mismatch(mini_8act: Path, reference_sequence: str, tmp_path: Path) -> None:
    """Changing a deposited motif identity must invalidate the whole structure."""
    mutated = tmp_path / "8act-mini-motif-identity-mismatch.cif"
    mutated.write_text(
        mini_8act.read_text(encoding="utf-8").replace("? 178 GLY A CA 1", "? 178 ALA A CA 1"),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match=r"chain A.*motif position 178"):
        validate_8act_structure(load_ca_coordinates(mutated), reference_sequence, FROZEN_PRIMARY_VARIANTS)

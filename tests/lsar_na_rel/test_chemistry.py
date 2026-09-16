"""Failing-first contract tests for substitution chemistry (#28 Task 2)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from lsar_na_rel.chemistry import (
    CHARGE_CONVENTION,
    ChemistryResult,
    formal_charge,
    generate_delta_charge_features,
    load_reference_sequence,
)
from lsar_na_rel.inputs import parse_variant

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_ROOT = REPO_ROOT / "workflows" / "lsar-na-rel"
PUBLIC_DIR = REPO_ROOT / "data" / "public" / "lsar-na-rel"
PRIMARY_CSV = PUBLIC_DIR / "na-rel-primary.csv"
REFERENCE_FASTA = Path(__file__).resolve().parent / "fixtures" / "P12883.fasta"

# Frozen formal-charge table (issue #26): D/E=-1, K/R=+1, H=0, else 0.
EXPECTED_CHARGES = {
    "D": -1,
    "E": -1,
    "K": 1,
    "R": 1,
    "H": 0,
    "A": 0,
    "N": 0,
    "Y": 0,
    "Q": 0,
    "G": 0,
    "P": 0,
    "W": 0,
    "I": 0,
    "T": 0,
    "S": 0,
    "F": 0,
    "L": 0,
    "V": 0,
}


@pytest.fixture(scope="module")
def p12883_sequence() -> str:
    return load_reference_sequence(REFERENCE_FASTA)


@pytest.fixture(scope="module")
def primary_variants() -> list[str]:
    return pd.read_csv(PRIMARY_CSV)["variant"].astype(str).tolist()


def test_formal_charge_convention_matches_approved_freeze():
    for residue, charge in EXPECTED_CHARGES.items():
        assert formal_charge(residue) == charge
    assert CHARGE_CONVENTION["D"] == -1
    assert CHARGE_CONVENTION["E"] == -1
    assert CHARGE_CONVENTION["K"] == 1
    assert CHARGE_CONVENTION["R"] == 1
    assert CHARGE_CONVENTION["H"] == 0


def test_y115h_and_h251n_are_zero_under_neutral_histidine(p12883_sequence, primary_variants):
    assert "Y115H" in primary_variants
    assert "H251N" in primary_variants
    result = generate_delta_charge_features(
        ["Y115H", "H251N"],
        reference_sequence=p12883_sequence,
    )
    by_variant = result.features.set_index("variant")["delta_charge"]
    assert by_variant.loc["Y115H"] == 0  # Y=0, H=0
    assert by_variant.loc["H251N"] == 0  # H=0, N=0


def test_sign_convention_is_mutant_minus_wt(p12883_sequence):
    # R(+1)->Q(0) => -1; G(0)->R(+1) => +1; D(-1)->Y(0) => +1
    result = generate_delta_charge_features(
        ["R249Q", "G768R", "D382Y"],
        reference_sequence=p12883_sequence,
    )
    by_variant = result.features.set_index("variant")["delta_charge"]
    assert by_variant.loc["R249Q"] == formal_charge("Q") - formal_charge("R")
    assert by_variant.loc["R249Q"] == -1
    assert by_variant.loc["G768R"] == formal_charge("R") - formal_charge("G")
    assert by_variant.loc["G768R"] == 1
    assert by_variant.loc["D382Y"] == formal_charge("Y") - formal_charge("D")
    assert by_variant.loc["D382Y"] == 1


def test_representative_charge_changes_across_primary_cohort(p12883_sequence, primary_variants):
    result = generate_delta_charge_features(
        primary_variants,
        reference_sequence=p12883_sequence,
    )
    by_variant = result.features.set_index("variant")["delta_charge"].to_dict()
    expected = {
        "Y115H": 0,
        "E497D": 0,  # E=-1, D=-1
        "R249Q": -1,
        "H251N": 0,
        "D382Y": 1,
        "I457T": 0,
        "R719W": -1,
        "P710R": 1,
        "D778V": 1,
        "L781P": 0,
        "S782N": 0,
        "A797T": 0,
        "F834L": 0,
        "R403Q": -1,
        "R663H": -1,  # R=+1, H=0 (not +1)
        "G768R": 1,
    }
    assert by_variant == expected


def test_identity_mismatch_is_hard_error(p12883_sequence):
    with pytest.raises(ValueError, match="WT mismatch|identity"):
        generate_delta_charge_features(
            ["A115H"],  # position 115 is Y in P12883
            reference_sequence=p12883_sequence,
        )


def test_out_of_range_position_is_hard_error(p12883_sequence):
    with pytest.raises(ValueError, match="position|out of range|identity"):
        generate_delta_charge_features(
            ["Y9999H"],
            reference_sequence=p12883_sequence,
        )


def test_invalid_variant_label_is_hard_error(p12883_sequence):
    with pytest.raises(ValueError, match="Invalid variant"):
        generate_delta_charge_features(
            ["not-a-variant"],
            reference_sequence=p12883_sequence,
        )


def test_nonstandard_residue_code_is_hard_error(p12883_sequence):
    # Syntactically parseable, but B is not a standard residue; do not impute.
    parsed = parse_variant("Y115B")
    assert parsed.mut == "B"
    with pytest.raises(ValueError, match="[Nn]on-standard|residue|amino"):
        generate_delta_charge_features(
            ["Y115B"],
            reference_sequence=p12883_sequence,
        )


def test_no_volume_or_hydropathy_columns_are_generated(p12883_sequence, primary_variants):
    result = generate_delta_charge_features(
        primary_variants,
        reference_sequence=p12883_sequence,
    )
    columns = set(result.features.columns)
    forbidden = {
        "delta_volume",
        "volume",
        "delta_hydropathy",
        "hydropathy",
        "delta_kyte_doolittle",
        "kyte_doolittle",
    }
    assert columns.isdisjoint(forbidden)
    assert "delta_charge" in columns
    assert "variant" in columns


def test_chemistry_does_not_load_target_or_labels(p12883_sequence, primary_variants):
    result = generate_delta_charge_features(
        primary_variants,
        reference_sequence=p12883_sequence,
    )
    assert isinstance(result, ChemistryResult)
    assert "y" not in result.features.columns
    assert "Na_rel" not in result.features.columns
    assert "ln(Na_rel)" not in result.features.columns
    assert "source_id" not in result.features.columns
    # Labels file is never required for chemistry generation.
    assert not hasattr(result, "target")


def test_provenance_metadata_records_charge_convention_and_sequence(p12883_sequence, primary_variants):
    result = generate_delta_charge_features(
        primary_variants,
        reference_sequence=p12883_sequence,
        sequence_accession="P12883",
        sequence_refseq="NP_000248.2",
    )
    provenance = result.provenance
    assert provenance["sequence_accession"] == "P12883"
    assert provenance["sequence_refseq"] == "NP_000248.2"
    assert provenance["sequence_length"] == len(p12883_sequence) == 1935
    assert provenance["charge_convention"]["H"] == 0
    assert provenance["delta_charge_definition"] == "q_mut - q_WT"
    assert "citation" in provenance
    assert provenance["feature_columns"] == ("delta_charge",)


def test_reference_fasta_loads_p12883_sequence():
    sequence = load_reference_sequence(REFERENCE_FASTA)
    assert len(sequence) == 1935
    assert sequence[114] == "Y"
    assert sequence[250] == "H"
    assert sequence[248] == "R"

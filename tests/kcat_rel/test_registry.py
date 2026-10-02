from __future__ import annotations

import json
from pathlib import Path

import pytest

from kcat_rel.registry import (
    FROZEN_PRIMARY_VARIANTS,
    load_registry,
    parse_variant,
    validate_reference,
)


WORKFLOW_ROOT = Path(__file__).resolve().parents[2] / "workflows/kcat-rel"
CONFIG = WORKFLOW_ROOT / "config" / "kcat-rel-v0.json"
EXPECTED_PRIMARY_VARIANTS = (
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


def write_registry(tmp_path: Path, *, variants: list[str]) -> Path:
    payload = json.loads(CONFIG.read_text(encoding="utf-8"))
    payload["variants"] = variants
    path = tmp_path / "kcat-rel-v0.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_registry_freezes_exact_active_profile() -> None:
    registry = load_registry(CONFIG)

    assert FROZEN_PRIMARY_VARIANTS == EXPECTED_PRIMARY_VARIANTS
    assert registry.variants == EXPECTED_PRIMARY_VARIANTS
    assert "D382Y" not in FROZEN_PRIMARY_VARIANTS
    assert "P710R" not in FROZEN_PRIMARY_VARIANTS
    assert registry.active_profile == "kcat-rel-v0-no-msa"
    assert registry.feature_columns == (
        "delta_charge",
        "grantham_distance",
        "catalytic_motif_min_ca_distance_8act_mean_ab",
    )


def test_registry_rejects_primary_cohort_drift(tmp_path: Path) -> None:
    altered = write_registry(tmp_path, variants=[*FROZEN_PRIMARY_VARIANTS, "D382Y"])

    with pytest.raises(ValueError, match="frozen primary cohort"):
        load_registry(altered)


@pytest.mark.parametrize(
    ("label", "wt", "position", "mutant"),
    [("R403Q", "R", 403, "Q"), ("Y115H", "Y", 115, "H")],
)
def test_parse_variant_preserves_the_strict_canonical_substitution(
    label: str, wt: str, position: int, mutant: str
) -> None:
    parsed = parse_variant(label)

    assert (parsed.label, parsed.wt, parsed.position, parsed.mutant) == (
        label,
        wt,
        position,
        mutant,
    )


@pytest.mark.parametrize("label", ["r403q", "R0Q", "R-1Q", "B403Q", "R403B"])
def test_parse_variant_rejects_noncanonical_or_nonpositive_labels(label: str) -> None:
    with pytest.raises(ValueError, match="invalid variant"):
        parse_variant(label)


def test_validate_reference_rejects_a_sequence_outside_the_pinned_p12883_identity() -> None:
    with pytest.raises(ValueError, match="P12883"):
        validate_reference("ACDEFGHIKLMNPQRSTVWY")

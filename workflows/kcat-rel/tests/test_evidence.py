from __future__ import annotations

import csv
import math
from pathlib import Path

import pytest

from kcat_rel.evidence import EvidenceAuditError, audit_evidence, main


WORKFLOW_ROOT = Path(__file__).resolve().parents[1]
CANONICAL_PATH = WORKFLOW_ROOT / "data" / "canonical-labels.csv"
MEASUREMENTS_PATH = WORKFLOW_ROOT / "data" / "measurement-evidence.csv"

EXPECTED_PRIMARY_VARIANTS = {
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
}


def _copy_csv_with_change(
    source: Path, destination: Path, *, variant: str, column: str, value: str
) -> None:
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fieldnames = list(rows[0])
    for row in rows:
        if row["variant"] == variant:
            row[column] = value
            break
    else:
        raise AssertionError(f"fixture variant not found: {variant}")
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_frozen_evidence_package_passes_all_gates() -> None:
    report = audit_evidence(CANONICAL_PATH, MEASUREMENTS_PATH)

    assert report.canonical_rows == 17
    assert report.variants == EXPECTED_PRIMARY_VARIANTS
    assert report.source_studies == {
        "sommese2013",
        "nag2015",
        "adhikari2016",
        "kawana2017",
        "adhikari2019",
        "sarkar2020",
        "morck2022",
        "nandwani2025",
        "pathak2026",
    }
    assert report.canonical_measurements == 22


def test_exact_log_label_is_recomputed_from_raw_pair() -> None:
    report = audit_evidence(CANONICAL_PATH, MEASUREMENTS_PATH)
    r403q = report.labels_by_variant["R403Q"]

    assert r403q.kcat_rel == pytest.approx(7.6 / 6.0, abs=1e-12)
    assert r403q.ln_kcat_rel == pytest.approx(0.23638877806423034, abs=1e-12)
    assert r403q.ln_kcat_rel == pytest.approx(math.log(r403q.kcat_rel), abs=1e-12)


def test_audit_rejects_a_tampered_log_label(tmp_path: Path) -> None:
    tampered = tmp_path / "canonical-labels.csv"
    _copy_csv_with_change(
        CANONICAL_PATH,
        tampered,
        variant="R403Q",
        column="ln_kcat_rel",
        value="0.25",
    )

    with pytest.raises(EvidenceAuditError, match="R403Q.*ln_kcat_rel"):
        audit_evidence(tampered, MEASUREMENTS_PATH)


def test_audit_rejects_a_tampered_published_replicate_aggregate(tmp_path: Path) -> None:
    tampered = tmp_path / "canonical-labels.csv"
    _copy_csv_with_change(
        CANONICAL_PATH,
        tampered,
        variant="L781P",
        column="kcat_rel",
        value="1.10",
    )

    with pytest.raises(EvidenceAuditError, match="L781P.*printed replicates"):
        audit_evidence(tampered, MEASUREMENTS_PATH)


def test_audit_rejects_sensitivity_rows_in_primary_cohort(tmp_path: Path) -> None:
    expanded = tmp_path / "canonical-labels.csv"
    expanded.write_text(
        CANONICAL_PATH.read_text(encoding="utf-8").replace("Y115H", "D382Y", 1),
        encoding="utf-8",
    )

    with pytest.raises(EvidenceAuditError, match="frozen primary cohort"):
        audit_evidence(expanded, MEASUREMENTS_PATH)


def test_audit_rejects_unmatched_wild_type_context(tmp_path: Path) -> None:
    tampered = tmp_path / "measurement-evidence.csv"
    _copy_csv_with_change(
        MEASUREMENTS_PATH,
        tampered,
        variant="Y115H",
        column="wt_context_match",
        value="false",
    )

    with pytest.raises(EvidenceAuditError, match="Y115H.*matched WT"):
        audit_evidence(CANONICAL_PATH, tampered)


def test_audit_rejects_nonpreferred_construct_selection(tmp_path: Path) -> None:
    tampered = tmp_path / "canonical-labels.csv"
    _copy_csv_with_change(
        CANONICAL_PATH,
        tampered,
        variant="H251N",
        column="construct_family",
        value="2-hep",
    )

    with pytest.raises(EvidenceAuditError, match="H251N.*construct"):
        audit_evidence(tampered, MEASUREMENTS_PATH)


def test_audit_rejects_duplicate_canonical_variants(tmp_path: Path) -> None:
    duplicated = tmp_path / "canonical-labels.csv"
    lines = CANONICAL_PATH.read_text(encoding="utf-8").splitlines()
    duplicated.write_text("\n".join([*lines, lines[1], ""]), encoding="utf-8")

    with pytest.raises(EvidenceAuditError, match="one canonical label per variant"):
        audit_evidence(duplicated, MEASUREMENTS_PATH)


def test_cli_reports_a_compact_success_summary(capsys: pytest.CaptureFixture[str]) -> None:
    assert main([str(CANONICAL_PATH), str(MEASUREMENTS_PATH)]) == 0
    assert capsys.readouterr().out == (
        "PASS: 17 canonical labels; 22 linked measurements; 9 source studies\n"
    )

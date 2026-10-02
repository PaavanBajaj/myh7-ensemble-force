"""Validate canonical ``kcat_rel`` labels against measurement evidence."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


class EvidenceAuditError(ValueError):
    """Raised when the evidence package violates a frozen-cohort gate."""


@dataclass(frozen=True)
class CanonicalLabel:
    variant: str
    kcat_rel: float
    ln_kcat_rel: float


@dataclass(frozen=True)
class AuditReport:
    canonical_rows: int
    canonical_measurements: int
    variants: frozenset[str]
    source_studies: frozenset[str]
    labels_by_variant: dict[str, CanonicalLabel]


FROZEN_PRIMARY_VARIANTS = frozenset(
    {
        "Y115H", "D239N", "R249Q", "H251N", "R403Q", "R453C",
        "I457T", "E497D", "R663H", "R719W", "G741R", "G768R",
        "D778V", "L781P", "S782N", "A797T", "F834L",
    }
)

EXPECTED_CONSTRUCT_FAMILY = {
    variant: "sS1/S1"
    for variant in (
        "Y115H", "D239N", "H251N", "R403Q", "R453C", "E497D",
        "R663H", "R719W", "G741R", "G768R",
    )
}
EXPECTED_CONSTRUCT_FAMILY.update(
    {
        variant: "2-hep HMM"
        for variant in (
            "R249Q", "I457T", "D778V", "L781P", "S782N", "A797T", "F834L",
        )
    }
)

CANONICAL_COLUMNS = frozenset(
    {
        "variant", "gene", "source_id", "canonical_measurement_ids",
        "construct_family", "aggregation_method", "wt_kcat_s_1",
        "wt_kcat_error", "mutant_kcat_s_1", "mutant_kcat_error",
        "kcat_rel", "kcat_rel_error", "ln_kcat_rel", "ln_kcat_rel_error",
    }
)

MEASUREMENT_COLUMNS = frozenset(
    {
        "measurement_id", "variant", "gene", "source_id",
        "biological_replicate", "canonical_evidence", "construct_family",
        "actin_system", "assay_regime", "detection_chemistry",
        "wt_kcat_s_1", "wt_kcat_error", "mutant_kcat_s_1",
        "mutant_kcat_error", "reported_kcat_rel", "reported_kcat_rel_error",
        "wt_context_id", "mutant_context_id", "wt_context_match",
        "source_location", "doi",
    }
)

ABS_TOLERANCE = 1e-12
PRINTED_REPLICATE_ROUNDING_TOLERANCE = 0.0050000001


def _read_csv(path: str | Path, required: frozenset[str], name: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise EvidenceAuditError(f"{name} missing columns: {sorted(missing)}")
        return list(reader)


def _float(row: dict[str, str], column: str, *, context: str) -> float:
    try:
        value = float(row[column])
    except (TypeError, ValueError) as exc:
        raise EvidenceAuditError(f"{context}: invalid {column}") from exc
    if not math.isfinite(value):
        raise EvidenceAuditError(f"{context}: non-finite {column}")
    return value


def _assert_close(actual: float, expected: float, *, message: str) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=ABS_TOLERANCE):
        raise EvidenceAuditError(message)


def _display_values(row: dict[str, str], column: str, *, context: str) -> list[float]:
    try:
        return [float(value) for value in row[column].split("|")]
    except ValueError as exc:
        raise EvidenceAuditError(f"{context}: invalid displayed {column}") from exc


def audit_evidence(canonical_path: str | Path, measurements_path: str | Path) -> AuditReport:
    """Audit cohort identity, provenance, label math, and replicate aggregation."""
    canonical = _read_csv(canonical_path, CANONICAL_COLUMNS, "canonical table")
    measurements = _read_csv(measurements_path, MEASUREMENT_COLUMNS, "measurement table")

    variants = [row["variant"] for row in canonical]
    if len(variants) != len(set(variants)):
        raise EvidenceAuditError("expected one canonical label per variant")
    if set(variants) != FROZEN_PRIMARY_VARIANTS:
        raise EvidenceAuditError("canonical table does not match frozen primary cohort")

    by_id: dict[str, dict[str, str]] = {}
    for row in measurements:
        measurement_id = row["measurement_id"]
        if not measurement_id or measurement_id in by_id:
            raise EvidenceAuditError("measurement_id values must be unique and nonempty")
        by_id[measurement_id] = row

    used_ids: list[str] = []
    labels: dict[str, CanonicalLabel] = {}
    studies: set[str] = set()

    for row in canonical:
        variant = row["variant"]
        if row["gene"] != "MYH7":
            raise EvidenceAuditError(f"{variant}: gene must be MYH7")
        expected_construct = EXPECTED_CONSTRUCT_FAMILY[variant]
        if row["construct_family"] != expected_construct:
            raise EvidenceAuditError(f"{variant}: construct violates sS1/S1-first selection rule")

        ids = [value for value in row["canonical_measurement_ids"].split("|") if value]
        if not ids:
            raise EvidenceAuditError(f"{variant}: missing canonical measurement IDs")
        try:
            selected = [by_id[measurement_id] for measurement_id in ids]
        except KeyError as exc:
            raise EvidenceAuditError(
                f"{variant}: unknown canonical measurement ID {exc.args[0]}"
            ) from exc
        used_ids.extend(ids)

        for measurement in selected:
            if measurement["variant"] != variant:
                raise EvidenceAuditError(f"{variant}: linked measurement has wrong variant")
            if measurement["gene"] != "MYH7":
                raise EvidenceAuditError(f"{variant}: linked measurement gene is not MYH7")
            if measurement["source_id"] != row["source_id"]:
                raise EvidenceAuditError(f"{variant}: source study mismatch")
            if measurement["construct_family"] != expected_construct:
                raise EvidenceAuditError(f"{variant}: linked measurement has wrong construct")
            if measurement["canonical_evidence"].lower() != "true":
                raise EvidenceAuditError(f"{variant}: linked measurement is not canonical evidence")
            if measurement["actin_system"] != "pure_actin":
                raise EvidenceAuditError(f"{variant}: assay is not pure-actin")
            if measurement["assay_regime"] != "steady_state_actin_activated_kcat":
                raise EvidenceAuditError(f"{variant}: assay is not eligible actin-activated kcat")
            if not measurement["detection_chemistry"]:
                raise EvidenceAuditError(f"{variant}: detection chemistry metadata missing")
            matched = measurement["wt_context_match"].lower() == "true"
            same_context = measurement["wt_context_id"] == measurement["mutant_context_id"]
            if not matched or not same_context:
                raise EvidenceAuditError(f"{variant}: canonical evidence lacks matched WT context")
            if not measurement["source_location"] or not measurement["doi"]:
                raise EvidenceAuditError(f"{variant}: incomplete provenance")

        displayed_wt = _display_values(row, "wt_kcat_s_1", context=variant)
        displayed_mutant = _display_values(row, "mutant_kcat_s_1", context=variant)
        measurement_wt = [
            _float(measurement, "wt_kcat_s_1", context=variant) for measurement in selected
        ]
        measurement_mutant = [
            _float(measurement, "mutant_kcat_s_1", context=variant)
            for measurement in selected
        ]
        if displayed_wt != measurement_wt or displayed_mutant != measurement_mutant:
            raise EvidenceAuditError(f"{variant}: displayed kcat values disagree with evidence")

        stored_ratio = _float(row, "kcat_rel", context=variant)
        stored_ratio_error = _float(row, "kcat_rel_error", context=variant)
        stored_log = _float(row, "ln_kcat_rel", context=variant)
        stored_log_error = _float(row, "ln_kcat_rel_error", context=variant)

        method = row["aggregation_method"]
        if method == "ratio_of_summary_means":
            if len(selected) != 1 or selected[0]["biological_replicate"] != "summary":
                raise EvidenceAuditError(f"{variant}: summary label must link one summary measurement")
            measurement = selected[0]
            wt = measurement_wt[0]
            mutant = measurement_mutant[0]
            wt_error = _float(measurement, "wt_kcat_error", context=variant)
            mutant_error = _float(measurement, "mutant_kcat_error", context=variant)
            expected_ratio = mutant / wt
            expected_ln_error = math.sqrt((mutant_error / mutant) ** 2 + (wt_error / wt) ** 2)
            expected_ratio_error = expected_ratio * expected_ln_error
        elif method == "mean_of_reported_paired_ratios":
            if len(selected) < 2 or any(
                measurement["biological_replicate"] == "summary" for measurement in selected
            ):
                raise EvidenceAuditError(f"{variant}: paired aggregation requires separate replicates")
            ratios = [
                _float(measurement, "reported_kcat_rel", context=variant)
                for measurement in selected
            ]
            ratio_errors = [
                _float(measurement, "reported_kcat_rel_error", context=variant)
                for measurement in selected
            ]
            printed_ratio_mean = sum(ratios) / len(ratios)
            printed_error_mean = sum(ratio_errors) / len(ratio_errors)
            if abs(stored_ratio - printed_ratio_mean) > PRINTED_REPLICATE_ROUNDING_TOLERANCE:
                raise EvidenceAuditError(
                    f"{variant}: published aggregate ratio is inconsistent with printed replicates"
                )
            if abs(stored_ratio_error - printed_error_mean) > PRINTED_REPLICATE_ROUNDING_TOLERANCE:
                raise EvidenceAuditError(
                    f"{variant}: published aggregate error is inconsistent with printed replicates"
                )
            expected_ratio = stored_ratio
            expected_ratio_error = stored_ratio_error
            expected_ln_error = expected_ratio_error / expected_ratio
        else:
            raise EvidenceAuditError(f"{variant}: unknown aggregation_method {method}")

        _assert_close(
            stored_ratio, expected_ratio,
            message=f"{variant}: kcat_rel disagrees with measurement evidence",
        )
        _assert_close(
            stored_ratio_error, expected_ratio_error,
            message=f"{variant}: kcat_rel_error disagrees with propagated uncertainty",
        )
        _assert_close(
            stored_log, math.log(expected_ratio),
            message=f"{variant}: ln_kcat_rel disagrees with exact log ratio",
        )
        _assert_close(
            stored_log_error, expected_ln_error,
            message=f"{variant}: ln_kcat_rel_error disagrees with propagated uncertainty",
        )

        studies.add(row["source_id"])
        labels[variant] = CanonicalLabel(variant, stored_ratio, stored_log)

    if len(used_ids) != len(set(used_ids)):
        raise EvidenceAuditError("a measurement is reused by multiple canonical labels")
    canonical_measurement_ids = {
        row["measurement_id"]
        for row in measurements
        if row["canonical_evidence"].lower() == "true"
    }
    if set(used_ids) != canonical_measurement_ids:
        raise EvidenceAuditError("canonical measurement links are incomplete")

    return AuditReport(
        canonical_rows=len(canonical),
        canonical_measurements=len(used_ids),
        variants=frozenset(variants),
        source_studies=frozenset(studies),
        labels_by_variant=labels,
    )


def main(argv: list[str] | None = None) -> int:
    """Run the evidence audit from the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("canonical_labels", type=Path)
    parser.add_argument("measurement_evidence", type=Path)
    args = parser.parse_args(argv)
    report = audit_evidence(args.canonical_labels, args.measurement_evidence)
    print(
        f"PASS: {report.canonical_rows} canonical labels; "
        f"{report.canonical_measurements} linked measurements; "
        f"{len(report.source_studies)} source studies"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

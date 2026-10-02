"""Assemble the three measured MYH7 ratios with their source contexts."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, stdev


FIELDS = (
    "variant", "outcome", "cohort", "value", "ln_value", "uncertainty",
    "error_type", "uncertainty_method", "source_id", "doi", "source_location",
    "source_row_id", "measurement_ids", "construct_family", "construct",
    "assay", "filament", "metric", "temperature_c", "wt_value", "mutant_value",
    "unit", "wt_error", "mutant_error", "wt_context_id", "aggregation_method",
    "table1_number", "table1_direction", "table1_agreement", "modeling_ready",
    "evidence_epoch", "notes",
)
NA_SOURCE_DOI = {
    "adhikari2019": "10.1038/s41467-019-10555-9",
    "morck2022": "10.7554/eLife.76805",
    "nandwani2025": "10.1038/s41467-025-63816-1",
    "pathak2026": "10.64898/2026.06.02.729673",
    "sarkar2020": "10.1126/sciadv.aax0069",
    "vanderroest2021": "10.1073/pnas.2025030118",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _numeric(value: str, *, description: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {description}: {value!r}") from exc
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"invalid {description}: {value!r}")
    return result


def _ratio(mutant: str, wt: str) -> float:
    return _numeric(mutant, description="mutant absolute") / _numeric(wt, description="WT absolute")


def _propagated_ratio_error(ratio: float, mutant: float, wt: float, mut_error: str, wt_error: str) -> str:
    if not mut_error or not wt_error:
        return ""
    mut_sd = float(mut_error)
    wt_sd = float(wt_error)
    if mut_sd < 0 or wt_sd < 0:
        raise ValueError("negative reported velocity error")
    return f"{ratio * math.hypot(mut_sd / mutant, wt_sd / wt):.12g}"


def _base(variant: str, outcome: str, cohort: str, value: float) -> dict[str, str]:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"nonpositive outcome for {variant} {outcome}")
    return {"variant": variant, "outcome": outcome, "cohort": cohort,
            "value": f"{value:.12g}", "ln_value": f"{math.log(value):.12g}"}


def assemble_outcomes(*, kcat_labels: Path, kcat_evidence: Path, na_primary: Path,
                      na_labels: Path, na_measurements: Path, velocity_summary: Path,
                      morck_pairs: Path, output_dir: Path) -> None:
    paths = {
        "kcat_labels": kcat_labels, "kcat_evidence": kcat_evidence,
        "na_primary": na_primary, "na_labels": na_labels,
        "na_measurements": na_measurements, "velocity_summary": velocity_summary,
        "morck_pairs": morck_pairs,
    }
    rows: list[dict[str, str]] = []
    kcat_evidence_by_id = {row["measurement_id"]: row for row in read_csv(kcat_evidence)}
    for line, label in enumerate(read_csv(kcat_labels), 2):
        ids = label["canonical_measurement_ids"].split("|")
        evidence = []
        for measurement_id in ids:
            record = kcat_evidence_by_id.get(measurement_id)
            if record is None or record["variant"] != label["variant"] or record["source_id"] != label["source_id"]:
                raise ValueError(f"kcat evidence mismatch for {label['variant']} {measurement_id}")
            evidence.append(record)
        value = _numeric(label["kcat_rel"], description="canonical kcat_rel")
        row = _base(label["variant"], "kcat_rel", "frozen_primary", value)
        row.update({
            "uncertainty": label["kcat_rel_error"], "error_type": label["error_type"],
            "uncertainty_method": label["uncertainty_method"],
            "source_id": label["source_id"],
            "doi": "|".join(dict.fromkeys(record["doi"] for record in evidence)),
            "source_location": "|".join(dict.fromkeys(record["source_location"] for record in evidence)),
            "source_row_id": f"workflows/kcat-rel/data/canonical-labels.csv:{line}",
            "measurement_ids": label["canonical_measurement_ids"],
            "construct_family": label["construct_family"], "construct": label["construct"],
            "assay": "|".join(dict.fromkeys(f"{record['assay_regime']}/{record['detection_chemistry']}" for record in evidence)),
            "filament": "|".join(dict.fromkeys(record["actin_system"] for record in evidence)),
            "temperature_c": "|".join(dict.fromkeys(record["temperature_c"] for record in evidence)),
            "wt_value": label["wt_kcat_s_1"], "mutant_value": label["mutant_kcat_s_1"],
            "unit": "s^-1", "wt_error": label["wt_kcat_error"],
            "mutant_error": label["mutant_kcat_error"],
            "wt_context_id": "|".join(dict.fromkeys(record["wt_context_id"] for record in evidence)),
            "aggregation_method": label["aggregation_method"], "modeling_ready": "true",
            "evidence_epoch": "|".join(dict.fromkeys(record["evidence_epoch"] for record in evidence)),
            "notes": label["selection_basis"],
        })
        rows.append(row)
    labels_by_variant = {row["variant"]: row for row in read_csv(na_labels)}
    na_measurements_by_variant: dict[str, list[dict[str, str]]] = defaultdict(list)
    for source_line, record in enumerate(read_csv(na_measurements), 2):
        record = {**record, "_source_line": str(source_line)}
        na_measurements_by_variant[record["variant"]].append(record)
    component_rows: list[dict[str, str]] = []
    for line, primary in enumerate(read_csv(na_primary), 2):
        variant = primary["variant"]
        label = labels_by_variant.get(variant)
        if label is None or label["source_id"] == "":
            raise ValueError(f"Na primary lacks label context for {variant}")
        measurements = [record for record in na_measurements_by_variant[variant]
                        if record["source_id"] == label["source_id"]]
        if not measurements:
            raise ValueError(f"Na primary lacks measurement links for {variant}")
        recomputed = _ratio(primary["LSAR_mut"], primary["LSAR_WT"])
        value = _numeric(primary["Na_rel"], description="canonical Na_rel")
        if abs(recomputed - value) > 1e-9 or abs(float(label["na_rel"]) - value) > 1e-9:
            raise ValueError(f"Na primary/within-study WT disagreement for {variant}")
        row = _base(variant, "Na_rel", "public_primary", value)
        row.update({
            "uncertainty": label["na_rel_error"],
            "error_type": "|".join(dict.fromkeys(record["error_type"] for record in measurements)),
            "uncertainty_method": "public_release_canonical",
            "source_id": label["source_id"],
            "doi": NA_SOURCE_DOI[label["source_id"]],
            "source_location": "|".join(dict.fromkeys(record["source_location"] for record in measurements)),
            "source_row_id": f"data/public/lsar-na-rel/na-rel-primary.csv:{line}",
            "measurement_ids": "|".join(record["measurement_id"] for record in measurements),
            "construct_family": "short/long_HMM",
            "construct": f"{label['short_construct']} / {label['long_construct']}",
            "assay": label["assay_method"], "filament": "actin_activated_ATPase",
            "metric": "LSAR_mut_over_study_WT_LSAR",
            "wt_value": primary["LSAR_WT"], "mutant_value": primary["LSAR_mut"],
            "unit": "LSAR", "wt_error": label["study_wt_lsar_error"],
            "mutant_error": label["reported_lsar_error"],
            "wt_context_id": f"{label['source_id']}-WT-LSAR",
            "aggregation_method": label["canonical_lsar_basis"],
            "modeling_ready": label["modeling_ready"],
            "evidence_epoch": label["evidence_epoch"],
            "notes": label["notes"],
        })
        rows.append(row)
        for measurement in measurements:
            component_rows.append({
                "measurement_id": measurement["measurement_id"],
                "variant": variant,
                "source_id": measurement["source_id"],
                "short_construct": measurement["short_construct"],
                "long_construct": measurement["long_construct"],
                "short_kcat_s_1": measurement["short_kcat_s_1"],
                "short_kcat_error": measurement["short_kcat_error"],
                "long_kcat_s_1": measurement["long_kcat_s_1"],
                "long_kcat_error": measurement["long_kcat_error"],
                "error_type": measurement["error_type"],
                "reported_lsar": measurement["reported_lsar"],
                "recomputed_lsar": measurement["recomputed_lsar"],
                "study_wt_lsar": primary["LSAR_WT"],
                "source_location": measurement["source_location"],
                "source_row_id": f"data/public/lsar-na-rel/measurements.csv:{measurement['_source_line']}",
            })
    summary_rows = read_csv(velocity_summary)
    if len({row["variant"] for row in summary_rows}) != len(summary_rows):
        raise ValueError("duplicate velocity summary variants")
    for line, evidence in enumerate(summary_rows, 2):
        variant = evidence["variant"]
        value = _ratio(evidence["mutant_velocity"], evidence["wt_velocity"])
        row = _base(variant, "v_rel", evidence["cohort"], value)
        typed_error = evidence["error_type"] in ("SEM", "SD")
        distributional_sd = variant == "M493I" and evidence["error_type"] == "SD"
        row.update({
            "uncertainty": _propagated_ratio_error(value, float(evidence["mutant_velocity"]),
                                                     float(evidence["wt_velocity"]),
                                                     evidence["mutant_error"], evidence["wt_error"]) if typed_error and not distributional_sd else "",
            "error_type": evidence["error_type"],
            "uncertainty_method": "not_available_distributional_SD" if distributional_sd else (
                "first_order_independent_errors" if typed_error and evidence["wt_error"] and evidence["mutant_error"]
                else "not_available_or_untyped_source_error"),
            "source_id": evidence["source_id"], "doi": evidence["doi"],
            "source_location": evidence["source_location"],
            "source_row_id": f"workflows/myh7-correlation-screen/data/velocity-summary-evidence-v1.csv:{line};{evidence['source_row_id']}",
            "measurement_ids": f"VS{line - 1:03d}",
            "construct_family": evidence["construct_family"], "construct": evidence["construct"],
            "assay": "unloaded_in_vitro_motility", "filament": evidence["filament"],
            "metric": evidence["metric"], "temperature_c": evidence["temperature_c"],
            "wt_value": evidence["wt_velocity"], "mutant_value": evidence["mutant_velocity"],
            "unit": evidence["velocity_unit"], "wt_error": evidence["wt_error"],
            "mutant_error": evidence["mutant_error"],
            "wt_context_id": evidence["wt_context_id"],
            "aggregation_method": "ratio_of_summary_means",
            "table1_number": evidence["table1_number"],
            "table1_direction": evidence["table1_direction"],
            "table1_agreement": evidence["table1_agreement"],
            "evidence_epoch": "preprint" if evidence["source_id"] in ("pathak2026", "lehman2023") else "published_primary",
            "notes": evidence["notes"],
        })
        rows.append(row)
    grouped_pairs: dict[str, list[dict[str, str]]] = defaultdict(list)
    for pair in read_csv(morck_pairs):
        if pair["source_id"] != "morck2022":
            raise ValueError("Morck pair source mismatch")
        exact = _ratio(pair["mutant_velocity_nm_s"], pair["wt_velocity_nm_s"])
        if abs(exact - float(pair["v_rel"])) > 1e-10:
            raise ValueError(f"Morck pair {pair['measurement_id']} ratio mismatch")
        grouped_pairs[pair["variant"]].append(pair)
    expected_counts = {"D778V": 5, "L781P": 5, "S782N": 9, "A797T": 7, "F834L": 6}
    if {variant: len(group) for variant, group in grouped_pairs.items()} != expected_counts:
        raise ValueError("Morck pair cohort/count mismatch")
    for variant, group in grouped_pairs.items():
        ratios = [float(pair["mutant_velocity_nm_s"]) / float(pair["wt_velocity_nm_s"]) for pair in group]
        value = mean(ratios)
        row = _base(variant, "v_rel", "table1_primary", value)
        row.update({
            "uncertainty": f"{stdev(ratios):.12g}", "error_type": "slide_ratio_sample_SD",
            "uncertainty_method": "sample_SD_of_paired_slide_ratios_not_biological_SEM",
            "source_id": "morck2022", "doi": "10.7554/eLife.76805",
            "source_location": "Figure 3F publisher source workbook",
            "source_row_id": "|".join(pair["summary_cell"] for pair in group),
            "measurement_ids": "|".join(pair["measurement_id"] for pair in group),
            "construct_family": "2-hep_HMM", "construct": "2-hep HMM",
            "assay": "unloaded_in_vitro_motility", "filament": "pure_F_actin",
            "metric": "filtered_MVEL", "temperature_c": "21-23",
            "wt_value": "|".join(pair["wt_velocity_nm_s"] for pair in group),
            "mutant_value": "|".join(pair["mutant_velocity_nm_s"] for pair in group),
            "unit": "nm/s", "wt_context_id": "|".join(dict.fromkeys(pair["wt_context_id"] for pair in group)),
            "aggregation_method": "mean_of_paired_slide_ratios",
            "table1_number": {"D778V": "46", "L781P": "30", "S782N": "0", "A797T": "5", "F834L": "0"}[variant],
            "table1_direction": {"D778V": "increase", "L781P": "decrease", "S782N": "none", "A797T": "decrease", "F834L": "none"}[variant],
            "table1_agreement": "aligned_no_change" if variant in ("S782N", "F834L") else "aligned",
            "evidence_epoch": "published_primary",
            "notes": f"{len(group)} selected slide/channel pairs; biological preparations fewer; shared WT controls across variants",
        })
        rows.append(row)
    keys = [(row["variant"], row["outcome"]) for row in rows]
    if len(keys) != len(set(keys)) or len(rows) != 51:
        raise ValueError(f"outcome identity/count mismatch: {len(rows)} rows")
    rows.sort(key=lambda row: (row["variant"], row["outcome"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "outcomes-v1.csv"
    with table_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    component_path = output_dir / "na-components-v1.csv"
    with component_path.open("w", newline="") as stream:
        component_fields = ("measurement_id", "variant", "source_id", "short_construct", "long_construct",
                            "short_kcat_s_1", "short_kcat_error", "long_kcat_s_1", "long_kcat_error",
                            "error_type", "reported_lsar", "recomputed_lsar", "study_wt_lsar",
                            "source_location", "source_row_id")
        writer = csv.DictWriter(stream, fieldnames=component_fields)
        writer.writeheader()
        writer.writerows(component_rows)
    manifest = {
        "schema_version": 1, "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": {name: sha256(path) for name, path in paths.items()},
        "output_sha256": {table_path.name: sha256(table_path), component_path.name: sha256(component_path)},
        "cohort_counts": {"kcat_rel": 17, "Na_rel": 16, "v_rel_table1_primary": 15, "v_rel_expanded_sensitivity": 3},
        "primary_scale": "unlogged_mutant_over_matched_study_WT_ratio",
        "na_definition": "canonical LSAR_mut / LSAR_WT_study from public release",
        "velocity_definition": "mutant_unloaded_velocity / matched_WT; Morck mean of within-slide ratios",
        "source_release_read_only": True,
    }
    (output_dir / "outcome-manifest-v1.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

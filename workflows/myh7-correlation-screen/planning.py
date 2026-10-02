"""Freeze the eligible comparison set before computing any associations."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _by_variant(path: Path) -> dict[str, dict[str, str]]:
    rows = read_csv(path)
    indexed = {row["variant"]: row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError(f"duplicate feature variant in {path}")
    return indexed


def _csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _study_counts(rows: list[dict[str, str]], variable: str) -> str:
    if variable not in ("kcat_rel", "Na_rel", "v_rel"):
        return ""
    counts = Counter(row[f"{variable}_study"] for row in rows if row[f"{variable}_study"])
    return "|".join(f"{study}:{counts[study]}" for study in sorted(counts))


def freeze_plan(*, outcomes_path: Path, new_features_path: Path, kcat_features_path: Path,
                na_features_path: Path, config_path: Path, output_dir: Path) -> None:
    config = json.loads(config_path.read_text())
    if config.get("schema_version") != 1 or config.get("primary_velocity_cohort") != "table1_primary":
        raise ValueError("unsupported or non-primary comparison config")
    outcomes = read_csv(outcomes_path)
    outcome_index = {(row["variant"], row["outcome"]): row for row in outcomes}
    if len(outcome_index) != len(outcomes):
        raise ValueError("duplicate outcome variant/target")
    new = _by_variant(new_features_path)
    kcat = _by_variant(kcat_features_path)
    na = _by_variant(na_features_path)
    features = [entry["name"] for entry in config["features"]]
    if len(features) != len(set(features)) or set(config["outcomes"]) != {"kcat_rel", "Na_rel", "v_rel"}:
        raise ValueError("feature or outcome config mismatch")
    table = []
    for variant in sorted(new):
        row = {"variant": variant}
        study_context = []
        for target in config["outcomes"]:
            evidence = outcome_index.get((variant, target))
            if evidence and (target != "v_rel" or evidence["cohort"] == "table1_primary"):
                row[target] = evidence["value"]
                row[f"{target}_study"] = evidence["source_id"]
                study_context.append(f"{target}:{evidence['source_id']}")
            else:
                row[target] = ""
                row[f"{target}_study"] = ""
        expanded = outcome_index.get((variant, "v_rel"))
        row["v_rel_expanded"] = expanded["value"] if expanded and expanded["cohort"] == "expanded_sensitivity" else ""
        row["v_rel_expanded_study"] = expanded["source_id"] if row["v_rel_expanded"] else ""
        if row["v_rel_expanded"]:
            study_context.append(f"v_expanded:{expanded['source_id']}")
        row["study_context"] = "|".join(study_context)
        for feature in features:
            candidates = [source[variant][feature] for source in (new, kcat, na)
                          if variant in source and feature in source[variant] and source[variant][feature] != ""]
            if len(candidates) > 1 and any(abs(float(candidate) - float(candidates[0])) > 1e-9 for candidate in candidates[1:]):
                raise ValueError(f"conflicting duplicate {feature} for {variant}")
            row[feature] = candidates[0] if candidates else ""
            if row[feature] and not math.isfinite(float(row[feature])):
                raise ValueError(f"nonfinite {feature} for {variant}")
        table.append(row)
    if len(table) != 22 or sum(bool(row["v_rel"]) for row in table) != 15 or sum(bool(row["v_rel_expanded"]) for row in table) != 3:
        raise ValueError("primary/expanded cohort boundary mismatch")
    variables = [*config["outcomes"], *features]
    kinds = {target: "continuous" for target in config["outcomes"]}
    kinds.update({entry["name"]: entry["kind"] for entry in config["features"]})
    units = {target: "ratio" for target in config["outcomes"]}
    units.update({entry["name"]: entry["unit"] for entry in config["features"]})
    requested_pairs = [
        *(("outcome_outcome", x, y) for x, y in itertools.combinations(config["outcomes"], 2)),
        *(("outcome_feature", x, y) for x in config["outcomes"] for y in features),
        *(("feature_feature", x, y) for x, y in itertools.combinations(features, 2)),
    ]
    inventory = []
    for number, (pair_type, x, y) in enumerate(requested_pairs, 1):
        paired = [row for row in table if row[x] and row[y]]
        inventory.append({
            "comparison_id": f"C{number:03d}", "pair_type": pair_type,
            "x": x, "y": y, "x_kind": kinds[x], "y_kind": kinds[y],
            "x_unit": units[x], "y_unit": units[y],
            "n_universe": len(table), "n_x_available": sum(bool(row[x]) for row in table),
            "n_y_available": sum(bool(row[y]) for row in table), "n_paired": len(paired),
            "n_missing_x": sum(not bool(row[x]) for row in table),
            "n_missing_y": sum(not bool(row[y]) for row in table),
            "n_missing_either": len(table) - len(paired),
            "paired_variants": "|".join(row["variant"] for row in paired),
            "x_study_composition": _study_counts(paired, x),
            "y_study_composition": _study_counts(paired, y),
            "missing_x_variants": "|".join(row["variant"] for row in table if not row[x]),
            "missing_y_variants": "|".join(row["variant"] for row in table if not row[y]),
        })
    if len(inventory) != 78:
        raise ValueError("comparison plan count mismatch")
    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "analysis-table-v1.csv"
    inventory_path = output_dir / "comparison-inventory-v1.csv"
    table_fields = ["variant", *config["outcomes"], "v_rel_expanded", *features,
                    *(f"{target}_study" for target in config["outcomes"]),
                    "v_rel_expanded_study", "study_context"]
    inventory_fields = list(inventory[0])
    _csv(table_path, table, table_fields)
    _csv(inventory_path, inventory, inventory_fields)
    manifest = {
        "schema_version": 1, "planned_at_utc": datetime.now(timezone.utc).isoformat(),
        "association_estimates_computed": False,
        "input_sha256": {"outcomes": sha256(outcomes_path), "new_features": sha256(new_features_path),
                         "kcat_features": sha256(kcat_features_path), "na_features": sha256(na_features_path),
                         "config": sha256(config_path)},
        "output_sha256": {table_path.name: sha256(table_path), inventory_path.name: sha256(inventory_path)},
        "comparison_counts": {"outcome_outcome": 3, "outcome_feature": 30, "feature_feature": 45},
        "primary_scale": config["primary_outcome_scale"],
        "missingness_universe": len(table),
    }
    (output_dir / "comparison-plan-manifest-v1.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

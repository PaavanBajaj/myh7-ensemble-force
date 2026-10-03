"""Frozen descriptive nonlinear dependence screen. No actual-cohort inference."""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import platform
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import rankdata
import scipy

from .statistics import distance_scores, fmt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "data/public/myh7-correlation-screen"
HISTORICAL = ROOT / "results/myh7-correlation-screen/results/screen-v1"
CONFIG = ROOT / "results/myh7-nonlinear-dependence/config/dependence-config-v1.json"
PLAN = ROOT / "data/public/myh7-nonlinear-dependence/dependence-plan-v1.json"
OUTPUT = ROOT / "results/myh7-nonlinear-dependence/results/screen-v1"
OUTCOMES = {"kcat_rel", "Na_rel", "v_rel"}
LABELS = {
    "kcat_rel": "relative catalytic turnover", "Na_rel": "relative available heads", "v_rel": "relative unloaded velocity",
    "adp_rigor_ca_shift_A": "site movement from ADP to rigor", "adp_adenine_ca_distance_8efe_A": "distance to ADP adenine",
    "rigor_actin_min_heavy_distance_8efi_A": "distance to actin", "relay_landmark_min_ca_distance_8act_mean_ab": "distance to the relay landmark",
    "elc_min_heavy_distance_8act_mean_ab": "distance to the essential light chain", "delta_charge": "mutation charge change",
    "grantham_distance": "amino-acid chemistry difference", "catalytic_motif_min_ca_distance_8act_mean_ab": "distance to the catalytic motif",
    "wt_site_rsa_8act_mean_ab": "wild-type site solvent exposure", "ihm_interface_flag_8act": "wild-type IHM-interface flag",
    "chemistry": "chemistry", "adp_geometry": "ADP geometry", "motor_geometry": "motor geometry", "all_features": "all ten features"
}
FILES = ["analysis-table-v1.csv", "comparison-inventory-v1.csv",
         "comparison-results-v1.csv", "comparison-plan-manifest-v1.json",
         "outcomes-v1.csv", "dependence-ledger-v1.csv", "velocity-dependence-v1.csv",
         "sensitivity-results-v1.csv", "feature-coverage-v1.csv"]
LIMIT = ("Descriptive only; heterogeneous assays, shared measurements and static WT annotations; "
         "no p-value, mechanism, interaction attribution or predictive accuracy; deletion ranges are not confidence intervals")
SCORE_FIELDS = ["distance_correlation", "u_centered_squared_distance_correlation", "status",
                "u_status", "x_constant_columns", "y_constant_columns", "x_varying_dimensions",
                "y_varying_dimensions", "x_scale_sample_sd_json", "y_scale_sample_sd_json"]
COMBINED_FIELDS = ["comparison_id", "pair_type", "analysis_role", "x", "y", "x_columns", "y_columns",
                   "x_dimensions", "y_dimensions", "original_method", "original_status", "original_effect",
                   "original_effect_unit", "spearman_rho", "spearman_blank_reason", "group_0_n", "group_1_n",
                   "original_loo_variant_min", "original_loo_variant_max", "original_loo_study_min", "original_loo_study_max",
                   "n_paired", "paired_variants", "n_universe", "n_missing_either", "missing_variants",
                   "missing_columns_by_variant_json", "feature_missing_reasons_json", "n_studies", "study_ids",
                   "x_study_composition", "y_study_composition"] + SCORE_FIELDS + [
                   "rank_distance_correlation", "rank_u_centered_squared_distance_correlation", "rank_status",
                   "log_outcome_distance_correlation", "log_outcome_u_centered_squared_distance_correlation",
                   "log_outcome_status", "loo_variant_distance_min", "loo_variant_distance_max",
                   "loo_variant_u_min", "loo_variant_u_max", "loo_variant_estimable", "loo_study_distance_min",
                   "loo_study_distance_max", "loo_study_u_min", "loo_study_u_max", "loo_study_estimable",
                   "interpretation_limit", "description"]
SENS_FIELDS = ["comparison_id", "scenario", "omitted_unit", "n_paired", "paired_variants", "excluded_variants",
               "study_ids", "original_spearman_rho", "note"] + SCORE_FIELDS


def input_path(source: Path, name: str) -> Path:
    """Resolve the published split layout or an explicit nine-file input pack."""
    if source.resolve() == SOURCE.resolve() and name in {"comparison-results-v1.csv", "sensitivity-results-v1.csv"}:
        return HISTORICAL / name
    return source / name


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(value: object) -> str:
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"missing or duplicate columns: {path.name}")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"malformed CSV: {path.name}")
    return rows


def write(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def unique(rows: list[dict], key: str, label: str) -> dict[str, dict]:
    mapping = {row[key]: row for row in rows}
    if len(mapping) != len(rows) or "" in mapping:
        raise ValueError(f"duplicate or missing {key}: {label}")
    return mapping


def complete(table: list[dict], spec: dict) -> list[dict]:
    columns = spec["x_columns"] + spec["y_columns"]
    return [row for row in table if all(row[column] != "" for column in columns)]


def studies(row: dict, spec: dict) -> set[str]:
    return {row[f"{column}_study"] for column in spec["x_columns"] + spec["y_columns"]
            if column in OUTCOMES and row[f"{column}_study"]}


def composition(rows: list[dict], columns: list[str]) -> str:
    counts = Counter(row[f"{c}_study"] for row in rows for c in columns if c in OUTCOMES)
    return "|".join(f"{key}:{counts[key]}" for key in sorted(counts))



def score(rows: list[dict], spec: dict, transform: str = "raw") -> dict:
    scales: dict[str, float] = {}
    constants: dict[str, list[str]] = {}
    arrays = []
    for side in ("x", "y"):
        columns = spec[f"{side}_columns"]
        array = np.asarray([[float(row[c]) for c in columns] for row in rows], dtype=float).reshape(len(rows), len(columns))
        if transform == "rank" and len(rows):
            array = rankdata(array, axis=0, method="average")
        elif transform == "log":
            for index, column in enumerate(columns):
                if column in OUTCOMES:
                    array[:, index] = np.log(array[:, index])
        varying = np.ptp(array, axis=0) > 0 if len(rows) else np.zeros(len(columns), dtype=bool)
        constants[side] = [c for c, varies in zip(columns, varying, strict=True) if not varies]
        for index, column in enumerate(columns):
            if varying[index] and column not in OUTCOMES:
                sd = float(np.std(array[:, index], ddof=1))
                scales[column] = sd
                array[:, index] = (array[:, index] - array[:, index].mean()) / sd
        arrays.append(array[:, varying])
    return distance_scores(*arrays) | {
        "x_constant_columns": "|".join(constants["x"]), "y_constant_columns": "|".join(constants["y"]),
        "x_varying_dimensions": str(arrays[0].shape[1]), "y_varying_dimensions": str(arrays[1].shape[1]),
        "x_scale_sample_sd_json": canonical({c: scales[c] for c in spec["x_columns"] if c in scales}),
        "y_scale_sample_sd_json": canonical({c: scales[c] for c in spec["y_columns"] if c in scales})}


def load(source: Path, config_path: Path) -> tuple[list[dict], list[dict], dict, dict, list[dict]]:
    table = read(input_path(source, FILES[0]))
    unique(table, "variant", "analysis table")
    if len(table) != 22 or [r["variant"] for r in table] != sorted(r["variant"] for r in table):
        raise ValueError("expected frozen sorted 22-variant analysis universe")
    config = json.loads(config_path.read_text())
    if config != json.loads(CONFIG.read_text()):
        raise ValueError("configuration differs from fixed implementation contract")
    variables = list(config["outcomes"]) + config["all_feature_sensitivity"]
    for row in table:
        for column in variables:
            value = row[column]
            if value and not math.isfinite(float(value)):
                raise ValueError(f"nonfinite value in {row['variant']}/{column}")
            if value and column in OUTCOMES and (float(value) <= 0 or not row[f"{column}_study"]):
                raise ValueError("outcome must be positive and have a source study")
            if value and column == "ihm_interface_flag_8act" and float(value) not in (0, 1):
                raise ValueError("nonbinary IHM annotation")
    inventory = read(input_path(source, FILES[1]))
    original = unique(read(input_path(source, FILES[2])), "comparison_id", "historical results")
    unique(inventory, "comparison_id", "inventory")
    expected_ids = [f"C{i:03d}" for i in range(1, 79)]
    if [r["comparison_id"] for r in inventory] != expected_ids or set(original) != set(expected_ids):
        raise ValueError("expected all 78 historical comparisons")
    if {frozenset((r["x"], r["y"])) for r in inventory} != {frozenset(p) for p in itertools.combinations(variables, 2)}:
        raise ValueError("historical inventory pair set differs")
    old_plan = json.loads((input_path(source, FILES[3])).read_text())
    for name in FILES[:2]:
        if sha(input_path(source, name)) != old_plan["output_sha256"][name]:
            raise ValueError(f"historical frozen plan hash mismatch: {name}")
    if sha(input_path(source, "outcomes-v1.csv")) != old_plan["input_sha256"]["outcomes"]:
        raise ValueError("historical outcome hash mismatch")
    specs = []
    for inv in inventory:
        spec = {"comparison_id": inv["comparison_id"], "pair_type": inv["pair_type"], "analysis_role": "historical_pair",
                "x": inv["x"], "y": inv["y"], "x_columns": [inv["x"]], "y_columns": [inv["y"]]}
        rows = complete(table, spec)
        prior = original[inv["comparison_id"]]
        variants = "|".join(row["variant"] for row in rows)
        if any(str(len(rows)) != r["n_paired"] or variants != r["paired_variants"] for r in (inv, prior)):
            raise ValueError("historical paired count/variant mismatch")
        if prior["x"] != inv["x"] or prior["y"] != inv["y"] or prior["pair_type"] != inv["pair_type"]:
            raise ValueError("historical comparison identity mismatch")
        binary = "ihm_interface_flag_8act" in (inv["x"], inv["y"])
        expected_method = "median_flag1_minus_flag0" if binary else "spearman_tie_aware"
        if prior["method"] != expected_method:
            raise ValueError("historical method mismatch")
        xs, ys = [float(r[inv["x"]]) for r in rows], [float(r[inv["y"]]) for r in rows]
        if binary:
            bc = "ihm_interface_flag_8act"
            other = inv["y"] if inv["x"] == bc else inv["x"]
            groups = [[float(r[other]) for r in rows if float(r[bc]) == flag] for flag in (0, 1)]
            if any(str(len(groups[i])) != prior[f"group_{i}_n"] for i in (0, 1)):
                raise ValueError("historical binary group mismatch")
            expected_effect = float(np.median(groups[1]) - np.median(groups[0])) if all(groups) else None
        else:
            expected_effect = float(np.corrcoef(rankdata(xs), rankdata(ys))[0, 1]) if np.ptp(xs) and np.ptp(ys) and len(rows) >= 3 else None
        if (prior["effect"] == "") != (expected_effect is None) or expected_effect is not None and not math.isclose(float(prior["effect"]), expected_effect, abs_tol=1e-10):
            raise ValueError("historical effect differs from paired data")
        specs.append(spec)
    counter = 1
    for block, columns in list(config["primary_blocks"].items()) + [("all_features", config["all_feature_sensitivity"])]:
        for outcome in config["outcomes"]:
            specs.append({"comparison_id": f"B{counter:03d}", "pair_type": "block_outcome", "analysis_role": "all_feature_sensitivity" if block == "all_features" else "primary_block",
                          "x": block, "y": outcome, "x_columns": columns, "y_columns": [outcome]})
            counter += 1
    deps = unique(read(input_path(source, "dependence-ledger-v1.csv")), "variant", "dependence ledger")
    for row in deps.values():
        if row["same_mutant_short_S"] not in {"yes", "no"} or row["same_WT_short_S"] not in {"yes", "no", "uncertain"}:
            raise ValueError("invalid dependence ledger flags")
    return table, specs, original, deps, read(input_path(source, "velocity-dependence-v1.csv"))


def snapshot(source: Path, config: Path, table: list[dict], specs: list[dict]) -> dict:
    return {"schema_version": 1, "association_estimates_computed": False,
            "post_hoc_relative_to_historical_spearman": True,
            "input_sha256": {name: sha(input_path(source, name)) for name in FILES},
            "config_sha256": sha(config), "implementation_sha256": {name: sha(HERE / name) for name in ("__init__.py", "analysis.py", "statistics.py", "screen.py")},
            "method_settings": json.loads(config.read_text()), "comparisons": specs,
            "variant_ids": [r["variant"] for r in table], "comparison_count": len(specs)}


def plan(source: Path, config: Path, output: Path) -> None:
    table, specs, _, _, _ = load(source, config)
    frozen = snapshot(source, config, table, specs)
    frozen["plan_content_sha256"] = digest(frozen)
    output.mkdir(parents=True, exist_ok=True)
    (output / "dependence-plan-v1.json").write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")


def analyze(source: Path, config: Path, plan_path: Path, output: Path) -> None:
    frozen = json.loads(plan_path.read_text())
    saved_digest = frozen.pop("plan_content_sha256")
    if digest(frozen) != saved_digest:
        raise ValueError("plan content tampering detected")
    table, specs, original, dependencies, velocity = load(source, config)
    if frozen != snapshot(source, config, table, specs):
        raise ValueError("inputs, configuration or implementation differ from frozen plan")
    coverage = {(r["variant"], r["feature"]): r["missing_reason"] for r in read(input_path(source, "feature-coverage-v1.csv")) if r["available"] == "0"}
    shared_s = {r["variant"] for r in dependencies.values() if r["same_mutant_short_S"] == "yes"}
    shared_wt = {r["variant"] for r in dependencies.values() if r["same_WT_short_S"] == "yes"}
    morck = {r["variant"] for r in velocity if r["cohort"] == "table1_primary" and r["shared_context_group"] == "morck2022_shared_slides"}
    historical_sensitivity = read(input_path(source, "sensitivity-results-v1.csv"))
    combined, sensitivity, point_rows = [], [], []
    for spec in specs:
        cid = spec["comparison_id"]
        rows = complete(table, spec)
        ids = {r["variant"] for r in rows}
        prior = original.get(cid, {})
        raw = score(rows, spec)
        rank = score(rows, spec, "rank")
        log = score(rows, spec, "log") if set(spec["x_columns"] + spec["y_columns"]) & OUTCOMES else None
        study_set = sorted(set().union(*(studies(r, spec) for r in rows))) if rows else []
        row_result = spec | raw | {"x_columns": "|".join(spec["x_columns"]), "y_columns": "|".join(spec["y_columns"]),
            "x_dimensions": len(spec["x_columns"]), "y_dimensions": len(spec["y_columns"]),
            "original_method": prior.get("method", "not_applicable_multivariate_block"),
            "original_status": prior.get("status", ""), "original_effect": prior.get("effect", ""),
            "original_effect_unit": prior.get("effect_unit", ""), "spearman_rho": prior.get("effect", "") if prior.get("method") == "spearman_tie_aware" else "",
            "spearman_blank_reason": "" if prior.get("method") == "spearman_tie_aware" else "binary_original_is_median_difference" if prior else "multivariate_vector_has_no_single_rho",
            "group_0_n": prior.get("group_0_n", ""), "group_1_n": prior.get("group_1_n", ""),
            "original_loo_variant_min": prior.get("loo_variant_min", ""), "original_loo_variant_max": prior.get("loo_variant_max", ""),
            "original_loo_study_min": prior.get("loo_study_min", ""), "original_loo_study_max": prior.get("loo_study_max", ""),
            "n_paired": len(rows), "paired_variants": "|".join(r["variant"] for r in rows),
            "n_universe": len(table), "n_missing_either": len(table) - len(rows),
            "missing_variants": "|".join(r["variant"] for r in table if r["variant"] not in ids),
            "missing_columns_by_variant_json": canonical({r["variant"]: [c for c in spec["x_columns"] + spec["y_columns"] if not r[c]] for r in table if r["variant"] not in ids}),
            "feature_missing_reasons_json": canonical({r["variant"]: {c: coverage.get((r["variant"], c), "missing_in_frozen_analysis_table") for c in spec["x_columns"] + spec["y_columns"] if c not in OUTCOMES and not r[c]} for r in table if any(not r[c] and c not in OUTCOMES for c in spec["x_columns"] + spec["y_columns"])}),
            "n_studies": len(study_set), "study_ids": "|".join(study_set),
            "x_study_composition": composition(rows, spec["x_columns"]), "y_study_composition": composition(rows, spec["y_columns"]),
            "rank_distance_correlation": rank["distance_correlation"], "rank_u_centered_squared_distance_correlation": rank["u_centered_squared_distance_correlation"],
            "rank_status": rank["u_status"], "log_outcome_distance_correlation": log["distance_correlation"] if log else "",
            "log_outcome_u_centered_squared_distance_correlation": log["u_centered_squared_distance_correlation"] if log else "",
            "log_outcome_status": log["u_status"] if log else "not_applicable_feature_pair", "interpretation_limit": LIMIT}
        def add_sensitivity(subset: list[dict], scenario: str, unit: str = "", note: str = "", old_rho: str = "") -> dict:
            result = {"comparison_id": cid, "scenario": scenario, "omitted_unit": unit, "n_paired": len(subset),
                      "paired_variants": "|".join(r["variant"] for r in subset),
                      "excluded_variants": "|".join(sorted(ids - {r["variant"] for r in subset})),
                      "study_ids": "|".join(sorted(set().union(*(studies(r, spec) for r in subset)))) if subset else "",
                      "original_spearman_rho": old_rho, "note": note} | score(subset, spec)
            sensitivity.append(result)
            return result
        for mode, units in (("variant", sorted(ids)), ("study", study_set)):
            deletions = []
            for unit in units:
                subset = [r for r in rows if (r["variant"] != unit if mode == "variant" else unit not in studies(r, spec))]
                deletions.append(add_sensitivity(subset, f"leave_one_{mode}_out", unit,
                                                "Scaling recomputed in retained cohort; study deletion uses union of outcome sources."))
            for field, short in (("distance_correlation", "distance"), ("u_centered_squared_distance_correlation", "u")):
                values = [float(r[field]) for r in deletions if r[field]]
                row_result[f"loo_{mode}_{short}_min"] = fmt(min(values)) if values else ""
                row_result[f"loo_{mode}_{short}_max"] = fmt(max(values)) if values else ""
            row_result[f"loo_{mode}_estimable"] = sum(bool(r["u_centered_squared_distance_correlation"]) for r in deletions) if units else "not_applicable_feature_pair" if mode == "study" else 0
        if cid == "C001":
            add_sensitivity([r for r in rows if r["variant"] not in shared_s], "exclude_shared_short_S", note="Exact reused mutant short-head ATPase components from dependence ledger; remaining rows are not independent replication.")
            add_sensitivity([r for r in rows if r["variant"] not in shared_s | shared_wt], "exclude_shared_short_S_or_WT", note="Union of exact reused mutant and WT short-head components; remaining study controls may still be shared.")
            for old in historical_sensitivity:
                if old["comparison_id"] == cid and old["scenario"] == "exclude_shared_short_S_and_risky_Na":
                    selected = set(old["paired_variants"].split("|"))
                    subset = [r for r in rows if r["variant"] in selected]
                    if len(subset) != int(old["n_paired"]):
                        raise ValueError("historical shared/risky sensitivity cohort mismatch")
                    add_sensitivity(subset, old["scenario"], note=old["note"], old_rho=old["effect"])
        if "v_rel" in spec["x_columns"] + spec["y_columns"]:
            add_sensitivity([r for r in rows if r["variant"] not in morck], "exclude_morck_shared_WT", note="Drops primary Morck variants with reused slide/channel WT records; no independence claim for remainder.")
            for old in historical_sensitivity:
                if old["comparison_id"] == cid and old["scenario"] == "strict_short_head_MVEL20":
                    selected = set(old["paired_variants"].split("|"))
                    subset = [r for r in rows if r["variant"] in selected]
                    if len(subset) != int(old["n_paired"]):
                        raise ValueError("historical strict velocity cohort mismatch")
                    add_sensitivity(subset, old["scenario"], note=old["note"], old_rho=old["effect"] if old["method"] == "spearman_tie_aware" else "")
        x_label, y_label = LABELS[spec['x']], LABELS[spec['y']]
        if prior.get("method") == "median_flag1_minus_flag0":
            binary_limit = "A singleton group makes the U score undefined, and shared study context limits interpretation." if raw["u_status"] == "zero_u_distance_variance" else "Sparse binary groups and shared study context limit interpretation."
            description = f"This row compares {x_label} with {y_label} across {len(rows)} measured variants. The earlier score is a binary-group median difference, while distance scores describe dependence without an increasing or decreasing direction. {binary_limit}"
        elif prior:
            description = f"This row compares {x_label} with {y_label} across {len(rows)} complete variants. Spearman summarizes increasing or decreasing rank patterns, while the distance scores also describe curved patterns. These exploratory scores do not establish a biological mechanism."
        else:
            description = f"This row compares the {x_label} feature block with {y_label} across {len(rows)} complete variants. The distance scores describe the block jointly and do not identify an individual feature or prove an interaction. Small cohorts, missing measurements and shared studies limit interpretation."
        row_result["description"] = description
        combined.append(row_result)
        for row in rows:
            for side in ("x", "y"):
                for c in spec[f"{side}_columns"]:
                    point_rows.append({"comparison_id": cid, "variant": row["variant"], "side": side,
                                       "column": c, "raw_value": row[c], "source_study": row.get(f"{c}_study", "")})
    output.mkdir(parents=True, exist_ok=True)
    write(output / "combined-correlations-v1.csv", combined, COMBINED_FIELDS)
    write(output / "dependence-sensitivities-v1.csv", sensitivity, SENS_FIELDS)
    write(output / "dependence-points-v1.csv", point_rows, ["comparison_id", "variant", "side", "column", "raw_value", "source_study"])
    manifests = {"schema_version": 1, "plan_sha256": sha(plan_path), "plan_content_sha256": saved_digest,
                 "input_sha256": frozen["input_sha256"], "config_sha256": frozen["config_sha256"],
                 "implementation_sha256": frozen["implementation_sha256"],
                 "output_sha256": {name: sha(output / name) for name in ("combined-correlations-v1.csv", "dependence-sensitivities-v1.csv", "dependence-points-v1.csv")},
                 "runtime_versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
                 "comparison_count": len(combined), "historical_pair_count": len(original),
                 "sensitivity_count": len(sensitivity), "association_estimates_computed": True,
                 "inferential_p_values_computed": False}
    (output / "dependence-run-manifest-v1.json").write_text(json.dumps(manifests, indent=2, sort_keys=True) + "\n")

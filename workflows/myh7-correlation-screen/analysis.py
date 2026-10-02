"""Compute descriptive associations from the previously frozen comparison plan."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from statistics import median


OUTCOMES = {"kcat_rel", "Na_rel", "v_rel"}
RESULT_FIELDS = [
    "comparison_id", "pair_type", "x", "y", "n_paired", "n_studies",
    "study_ids", "paired_variants", "status", "method", "effect", "effect_unit",
    "group_0_n", "group_1_n", "secondary_log_pearson", "loo_variant_min",
    "loo_variant_max", "loo_study_min", "loo_study_max", "loo_study_estimable",
    "n_missing_either", "interpretation_limit",
]
SENSITIVITY_FIELDS = [
    "comparison_id", "scenario", "n_paired", "paired_variants", "study_ids",
    "status", "method", "effect", "effect_unit", "note",
]
POINT_FIELDS = [
    "comparison_id", "variant", "x_value", "y_value", "x_study", "y_study",
    "study_context",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    position = 0
    while position < len(values):
        end = position + 1
        while end < len(values) and values[order[end]] == values[order[position]]:
            end += 1
        average_rank = (position + 1 + end) / 2
        for index in order[position:end]:
            result[index] = average_rank
        position = end
    return result


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3:
        return None
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    centered = [(x - x_mean, y - y_mean) for x, y in zip(xs, ys, strict=True)]
    x_variance = sum(x * x for x, _ in centered)
    y_variance = sum(y * y for _, y in centered)
    if x_variance <= 0 or y_variance <= 0:
        return None
    result = sum(x * y for x, y in centered) / math.sqrt(x_variance * y_variance)
    return max(-1.0, min(1.0, result))


def study_ids(row: dict[str, str], spec: dict[str, str]) -> set[str]:
    return {row[f"{variable}_study"] for variable in (spec["x"], spec["y"])
            if variable in OUTCOMES and row.get(f"{variable}_study")}


def paired(table: list[dict[str, str]], spec: dict[str, str]) -> list[dict[str, str]]:
    return [row for row in table if row.get(spec["x"]) and row.get(spec["y"])]


def estimate(rows: list[dict[str, str]], spec: dict[str, str]) -> dict[str, str]:
    x, y = spec["x"], spec["y"]
    xs = [float(row[x]) for row in rows]
    ys = [float(row[y]) for row in rows]
    binary = x if spec["x_kind"] == "binary" else y if spec["y_kind"] == "binary" else None
    if binary:
        other = y if binary == x else x
        groups = {flag: [float(row[other]) for row in rows if float(row[binary]) == flag]
                  for flag in (0.0, 1.0)}
        if any(float(row[binary]) not in (0.0, 1.0) for row in rows):
            raise ValueError(f"nonbinary value in {binary}")
        n0, n1 = len(groups[0.0]), len(groups[1.0])
        effect = median(groups[1.0]) - median(groups[0.0]) if n0 and n1 else None
        return {
            "status": "estimated" if n0 >= 2 and n1 >= 2 else "fragile_group_n_lt2" if effect is not None else "group_absent",
            "method": "median_flag1_minus_flag0", "effect": fmt(effect),
            "effect_unit": spec["y_unit"] if binary == x else spec["x_unit"],
            "group_0_n": str(n0), "group_1_n": str(n1),
            "secondary_log_pearson": "",
        }
    rank_correlation = pearson(ranks(xs), ranks(ys))
    transformed_x = [math.log(value) if x in OUTCOMES else value for value in xs]
    transformed_y = [math.log(value) if y in OUTCOMES else value for value in ys]
    secondary = pearson(transformed_x, transformed_y) if x in OUTCOMES or y in OUTCOMES else None
    return {
        "status": "estimated" if rank_correlation is not None else "n_lt3_or_constant",
        "method": "spearman_tie_aware", "effect": fmt(rank_correlation),
        "effect_unit": "rho", "group_0_n": "", "group_1_n": "",
        "secondary_log_pearson": fmt(secondary),
    }


def fmt(value: float | None) -> str:
    return "" if value is None else f"{value:.12g}"


def effect_number(rows: list[dict[str, str]], spec: dict[str, str]) -> float | None:
    result = estimate(rows, spec)["effect"]
    return float(result) if result else None


def leave_one_out(rows: list[dict[str, str]], spec: dict[str, str], *, by_study: bool) -> tuple[str, str, str]:
    if by_study:
        if not ({spec["x"], spec["y"]} & OUTCOMES):
            return "", "", "not_applicable_feature_pair"
        omitted = sorted(set().union(*(study_ids(row, spec) for row in rows))) if rows else []
        subsets = [[row for row in rows if study not in study_ids(row, spec)] for study in omitted]
    else:
        subsets = [[item for item in rows if item["variant"] != row["variant"]] for row in rows]
    effects = [value for subset in subsets if (value := effect_number(subset, spec)) is not None]
    if not effects:
        return "", "", "0"
    return fmt(min(effects)), fmt(max(effects)), str(len(effects))


def _svg_plot(path: Path, spec: dict[str, str], rows: list[dict[str, str]], effect: dict[str, str]) -> None:
    width, height = 1050, max(550, 160 + 25 * len(rows))
    left, right, top, bottom = 85, 620, 80, height - 75
    xs = [float(row[spec["x"]]) for row in rows]
    ys = [float(row[spec["y"]]) for row in rows]
    def limits(values: list[float], binary: bool) -> tuple[float, float]:
        if binary:
            return -0.2, 1.2
        if not values:
            return 0, 1
        low, high = min(values), max(values)
        spread = high - low or max(abs(low) * 0.15, 1.0)
        return low - spread * 0.1, high + spread * 0.1
    x0, x1 = limits(xs, spec["x_kind"] == "binary")
    y0, y1 = limits(ys, spec["y_kind"] == "binary")
    esc = html.escape
    title = f"{spec['comparison_id']}  {spec['x']} vs {spec['y']}"
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{esc(title)}">',
        '<rect width="100%" height="100%" fill="#000"/>',
        f'<text x="50" y="37" fill="#fff" font-family="sans-serif" font-size="20">{esc(title)}</text>',
        f'<text x="50" y="60" fill="#aaa" font-family="sans-serif" font-size="14">n={len(rows)} · {esc(effect["method"])}={esc(effect["effect"] or "undefined")} · descriptive</text>',
        f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#777"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#777"/>',
        f'<text x="{(left + right) / 2:.0f}" y="{height - 20}" fill="#fff" text-anchor="middle" font-family="sans-serif" font-size="14">{esc(spec["x"])} ({esc(spec["x_unit"])})</text>',
        f'<text x="18" y="{(top + bottom) / 2:.0f}" fill="#fff" text-anchor="middle" transform="rotate(-90 18 {(top + bottom) / 2:.0f})" font-family="sans-serif" font-size="14">{esc(spec["y"])} ({esc(spec["y_unit"])})</text>',
    ]
    for index in range(5):
        fraction = index / 4
        x_tick = left + fraction * (right - left)
        y_tick = bottom - fraction * (bottom - top)
        parts.append(f'<text x="{x_tick:.1f}" y="{bottom + 21}" fill="#aaa" text-anchor="middle" font-family="sans-serif" font-size="12">{x0 + fraction * (x1 - x0):.3g}</text>')
        parts.append(f'<text x="{left - 11}" y="{y_tick + 4:.1f}" fill="#aaa" text-anchor="end" font-family="sans-serif" font-size="12">{y0 + fraction * (y1 - y0):.3g}</text>')
    parts.append('<text x="660" y="90" fill="#fff" font-family="sans-serif" font-size="14">Variant · source study</text>')
    for index, row in enumerate(rows, 1):
        x_value, y_value = float(row[spec["x"]]), float(row[spec["y"]])
        px = left + (x_value - x0) / (x1 - x0) * (right - left)
        py = bottom - (y_value - y0) / (y1 - y0) * (bottom - top)
        study = ", ".join(sorted(study_ids(row, spec))) or "intrinsic features"
        parts.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="6" fill="#fff" stroke="#000" stroke-width="1"><title>{esc(row["variant"])} · {esc(study)} · ({x_value:.5g}, {y_value:.5g})</title></circle>')
        parts.append(f'<text x="{px + 9:.2f}" y="{py - 9:.2f}" fill="#fff" font-family="sans-serif" font-size="11">{index}</text>')
        parts.append(f'<text x="660" y="{100 + 25 * index}" fill="#ddd" font-family="sans-serif" font-size="13">{index}. {esc(row["variant"])} · {esc(study)}</text>')
    parts.append('</svg>')
    path.write_text("\n".join(parts) + "\n")


def analyze(*, table_path: Path, inventory_path: Path, plan_manifest_path: Path,
            outcomes_path: Path, dependence_ledger_path: Path,
            velocity_alternatives_path: Path, output_dir: Path) -> None:
    plan = json.loads(plan_manifest_path.read_text())
    expected = plan["output_sha256"]
    if plan["association_estimates_computed"] or sha256(table_path) != expected[table_path.name] or sha256(inventory_path) != expected[inventory_path.name]:
        raise ValueError("analysis table or inventory differs from frozen pre-analysis plan")
    if sha256(outcomes_path) != plan["input_sha256"]["outcomes"]:
        raise ValueError("outcomes differ from frozen pre-analysis plan")
    table = read_csv(table_path)
    inventory = read_csv(inventory_path)
    outcome_rows = {(row["variant"], row["outcome"]): row for row in read_csv(outcomes_path)}
    dependencies = {row["variant"]: row for row in read_csv(dependence_ledger_path)}
    alternatives = read_csv(velocity_alternatives_path)
    if len(inventory) != 78 or len(table) != 22:
        raise ValueError("frozen plan count mismatch")
    output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    results: list[dict[str, str]] = []
    points: list[dict[str, str]] = []
    sensitivity: list[dict[str, str]] = []

    def add_sensitivity(spec: dict[str, str], scenario: str, candidates: list[dict[str, str]], note: str) -> None:
        subset = paired(candidates, spec)
        stat = estimate(subset, spec)
        sensitivity.append({
            "comparison_id": spec["comparison_id"], "scenario": scenario,
            "n_paired": len(subset), "paired_variants": "|".join(row["variant"] for row in subset),
            "study_ids": "|".join(sorted(set().union(*(study_ids(row, spec) for row in subset)))) if subset else "",
            "status": stat["status"], "method": stat["method"], "effect": stat["effect"],
            "effect_unit": stat["effect_unit"], "note": note,
        })

    for spec in inventory:
        selected = paired(table, spec)
        if len(selected) != int(spec["n_paired"]) or "|".join(row["variant"] for row in selected) != spec["paired_variants"]:
            raise ValueError(f"inventory mismatch for {spec['comparison_id']}")
        stat = estimate(selected, spec)
        variant_min, variant_max, _ = leave_one_out(selected, spec, by_study=False)
        study_min, study_max, study_n = leave_one_out(selected, spec, by_study=True)
        studies = sorted(set().union(*(study_ids(row, spec) for row in selected))) if selected else []
        results.append({
            **{field: spec[field] for field in ("comparison_id", "pair_type", "x", "y", "n_paired", "paired_variants", "n_missing_either")},
            **stat, "n_studies": len(studies), "study_ids": "|".join(studies),
            "loo_variant_min": variant_min, "loo_variant_max": variant_max,
            "loo_study_min": study_min, "loo_study_max": study_max,
            "loo_study_estimable": study_n,
            "interpretation_limit": "exploratory; correlated measures and heterogeneous assays; no inferential p-value",
        })
        for row in selected:
            points.append({
                "comparison_id": spec["comparison_id"], "variant": row["variant"],
                "x_value": row[spec["x"]], "y_value": row[spec["y"]],
                "x_study": row.get(f"{spec['x']}_study", ""),
                "y_study": row.get(f"{spec['y']}_study", ""),
                "study_context": row["study_context"],
            })
        _svg_plot(plots_dir / f"{spec['comparison_id']}.svg", spec, selected, stat)
        if spec["pair_type"] == "outcome_outcome":
            add_sensitivity(spec, "same_source_study", [row for row in table if row[f"{spec['x']}_study"] and row[f"{spec['x']}_study"] == row[f"{spec['y']}_study"]], "Both outcomes drawn from the same publication; construct and measurement independence still require the ledger.")
        if spec["comparison_id"] == "C001":
            independent = [row for row in table if row["variant"] not in dependencies or dependencies[row["variant"]]["same_mutant_short_S"] != "yes"]
            add_sensitivity(spec, "exclude_shared_short_S", independent, "Drops seven variants whose Na and kcat reuse the exact mutant short-head ATPase value.")
            independent_low_risk = [row for row in independent if row["variant"] not in dependencies or dependencies[row["variant"]]["na_risky"] != "yes"]
            add_sensitivity(spec, "exclude_shared_short_S_and_risky_Na", independent_low_risk, "Also drops three public Na values flagged for low n or epoch/uncertainty concerns.")
        if "v_rel" in (spec["x"], spec["y"]):
            strict = {variant for (variant, target), evidence in outcome_rows.items()
                      if target == "v_rel" and evidence["cohort"] == "table1_primary"
                      and evidence["construct_family"] == "sS1/S1"
                      and evidence["metric"] == "MVEL20"
                      and evidence["filament"] == "pure_F_actin"}
            add_sensitivity(spec, "strict_short_head_MVEL20", [row for row in table if row["variant"] in strict], "Four selected pure-actin short-head MVEL20 variants; source conditions still differ.")
            add_sensitivity(spec, "exclude_morck_shared_WT", [row for row in table if row["v_rel_study"] != "morck2022"], "Drops the five Morck variants whose slide/channel WT records are reused across mutants.")
            expanded = []
            for row in table:
                changed = dict(row)
                if row["v_rel_expanded"]:
                    changed["v_rel"] = row["v_rel_expanded"]
                    changed["v_rel_study"] = row["v_rel_expanded_study"]
                expanded.append(changed)
            add_sensitivity(spec, "expanded_human_velocity", expanded, "Adds three outside-Table-1 human MYH7 pairs; different constructs and temperatures remain heterogeneous.")
            for alternative in alternatives:
                if alternative["filament"] != "pure_F_actin":
                    continue
                variant = alternative["variant"]
                candidate = expanded if variant == "R671C" else table
                changed = []
                for row in candidate:
                    item = dict(row)
                    if row["variant"] == variant:
                        item["v_rel"] = fmt(float(alternative["mutant_velocity"]) / float(alternative["wt_velocity"]))
                        item["v_rel_study"] = alternative["source_id"]
                    changed.append(item)
                scenario = f"alternative_{variant}_{alternative['metric']}_{alternative['temperature_c']}C_{alternative['status']}"
                add_sensitivity(spec, scenario, changed, alternative["reason"])
    write_csv(output_dir / "comparison-results-v1.csv", results, RESULT_FIELDS)
    write_csv(output_dir / "comparison-points-v1.csv", points, POINT_FIELDS)
    write_csv(output_dir / "sensitivity-results-v1.csv", sensitivity, SENSITIVITY_FIELDS)
    manifest = {
        "schema_version": 1, "analyzed_at_utc": datetime.now(timezone.utc).isoformat(),
        "plan_sha256": sha256(plan_manifest_path),
        "input_sha256": {path.name: sha256(path) for path in (table_path, inventory_path, outcomes_path, dependence_ledger_path, velocity_alternatives_path)},
        "comparison_count": len(results), "point_count": len(points),
        "sensitivity_count": len(sensitivity), "plot_count": len(list(plots_dir.glob("C*.svg"))),
        "output_sha256": {name: sha256(output_dir / name) for name in ("comparison-results-v1.csv", "comparison-points-v1.csv", "sensitivity-results-v1.csv")},
        "policy": "descriptive exploratory estimates; no p-values, confidence intervals, predictive claims, or causal interpretation",
    }
    (output_dir / "analysis-manifest-v1.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

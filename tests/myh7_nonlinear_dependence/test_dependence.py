"""External artifact contract plus independent estimator capability fixtures."""
from __future__ import annotations

import csv
import hashlib
import os
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import spearmanr

from myh7_nonlinear_dependence import analysis as screen

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / "src/myh7_nonlinear_dependence"
SOURCE = screen.SOURCE
FROZEN = screen.OUTPUT


def copy_inputs(destination: Path) -> None:
    destination.mkdir()
    for name in screen.FILES:
        shutil.copyfile(screen.input_path(SOURCE, name), destination / name)


def rows(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def cli(command: str, output: Path, source: Path = SOURCE, **kwargs: str) -> subprocess.CompletedProcess:
    args = [sys.executable, "-m", "myh7_nonlinear_dependence.screen", command, "--source-dir", str(source), "--output-dir", str(output)]
    for key, value in kwargs.items():
        args += ["--" + key.replace("_", "-"), str(value)]
    return subprocess.run(args, text=True, capture_output=True, check=False)


def assert_ok(completed: subprocess.CompletedProcess) -> None:
    assert completed.returncode == 0, completed.stderr


def reference(xs: np.ndarray, ys: np.ndarray, u: bool) -> float:
    """Scalar loop reference, independently expressed using paper's sums."""
    n = len(xs)
    matrices = []
    for points in (xs, ys):
        d = [[math.dist(points[i], points[j]) for j in range(n)] for i in range(n)]
        row_sums = [sum(row) for row in d]
        total = sum(row_sums)
        if u:
            matrices.append([[0.0 if i == j else d[i][j] - row_sums[i] / (n - 2) - row_sums[j] / (n - 2) + total / ((n - 1) * (n - 2)) for j in range(n)] for i in range(n)])
        else:
            matrices.append([[d[i][j] - row_sums[i] / n - row_sums[j] / n + total / n**2 for j in range(n)] for i in range(n)])
    a, b = matrices
    product = sum(a[i][j] * b[i][j] for i in range(n) for j in range(n))
    norm_a = sum(a[i][j] ** 2 for i in range(n) for j in range(n))
    norm_b = sum(b[i][j] ** 2 for i in range(n) for j in range(n))
    ratio = product / math.sqrt(norm_a * norm_b)
    return ratio if u else math.sqrt(max(0, ratio))


def test_estimator_matches_independent_reference_symmetry_and_units() -> None:
    xs = np.array([[0, 1], [1, 2], [2, 4], [4, 5], [8, 9]], dtype=float)
    ys = np.array([[5], [1], [6], [0], [9]], dtype=float)
    actual = screen.distance_scores(xs, ys)
    for key, u in (("distance_correlation", False), ("u_centered_squared_distance_correlation", True)):
        assert float(actual[key]) == pytest.approx(reference(xs, ys, u), abs=1e-11)
        assert float(actual[key]) == pytest.approx(float(screen.distance_scores(ys, xs)[key]), abs=1e-11)
        assert float(actual[key]) == pytest.approx(float(screen.distance_scores(xs * 7 + 99, ys * 0.2 - 2)[key]), abs=1e-11)


def test_negative_u_is_preserved_and_small_n_is_explicit() -> None:
    actual = screen.distance_scores([0, 1, 2, 3], [0, 2, 1, 3])
    assert actual["u_centered_squared_distance_correlation"] == "-0.5"
    assert float(actual["distance_correlation"]) == pytest.approx(math.sqrt(9 / 13))
    small = screen.distance_scores([1, 2, 3], [2, 5, 1])
    assert small["status"] == "estimated"
    assert small["u_status"] == "n_lt4"
    assert small["u_centered_squared_distance_correlation"] == ""
    constant = screen.distance_scores([1, 1, 1, 1], [0, 2, 4, 6])
    assert constant["status"] == constant["u_status"] == "constant_side"
    assert constant["distance_correlation"] == ""
    singleton = screen.distance_scores(np.arange(14), [0] * 13 + [1])
    assert singleton["status"] == "estimated"
    assert singleton["u_status"] == "zero_u_distance_variance"
    assert singleton["u_centered_squared_distance_correlation"] == ""
    assert screen.distance_scores([0] * 13 + [1], np.arange(14))["u_status"] == "zero_u_distance_variance"
    with pytest.raises(ValueError, match="nonfinite"):
        screen.distance_scores([1, 2, np.nan, 4], [1, 2, 3, 4])


def test_u_shape_and_joint_xor_capability_not_real_data_significance() -> None:
    x = np.linspace(-1, 1, 81)
    y = x**2
    assert abs(spearmanr(x, y).statistic) < 0.02
    assert float(screen.distance_scores(x, y)["u_centered_squared_distance_correlation"]) > 0.2
    # All four binary patterns repeated evenly: neither marginal carries label dependence.
    xy = np.tile(np.array([[0, 0], [0, 1], [1, 0], [1, 1]]), (20, 1))
    outcome = (xy[:, 0] != xy[:, 1]).astype(float)
    assert all(abs(spearmanr(xy[:, i], outcome).statistic) < 1e-12 for i in range(2))
    assert float(screen.distance_scores(xy, outcome)["u_centered_squared_distance_correlation"]) > 0.1
    assert all(float(screen.distance_scores(xy[:, i], outcome)["u_centered_squared_distance_correlation"]) < 0 for i in range(2))


def test_constant_dimensions_exposed_and_scaling_recomputed() -> None:
    spec = {"x_columns": ["f", "constant"], "y_columns": ["kcat_rel"]}
    table = [{"f": str(i), "constant": "4", "kcat_rel": str(i * i + 1)} for i in range(5)]
    result = screen.score(table, spec)
    assert result["x_constant_columns"] == "constant"
    assert result["x_varying_dimensions"] == "1"
    assert json.loads(result["x_scale_sample_sd_json"])["f"] == pytest.approx(np.std(np.arange(5), ddof=1))
    deleted = screen.score(table[:-1], spec)
    assert json.loads(deleted["x_scale_sample_sd_json"])["f"] == pytest.approx(np.std(np.arange(4), ddof=1))


def test_cli_preserves_pairs_binary_and_fixed_blocks_reproducibly(tmp_path: Path) -> None:
    assert_ok(cli("plan", tmp_path))
    plan = json.loads((tmp_path / "dependence-plan-v1.json").read_text())
    assert plan["association_estimates_computed"] is False
    assert len(plan["comparisons"]) == 90
    assert_ok(cli("analyze", tmp_path))
    for name in ("combined-correlations-v1.csv", "dependence-sensitivities-v1.csv", "dependence-points-v1.csv"):
        assert (tmp_path / name).read_bytes() == (FROZEN / name).read_bytes()
    combined = rows(tmp_path / "combined-correlations-v1.csv")
    prior = {r["comparison_id"]: r for r in rows(screen.input_path(SOURCE, "comparison-results-v1.csv"))}
    assert len(combined) == 90
    assert [r["comparison_id"] for r in combined[:78]] == [f"C{i:03d}" for i in range(1, 79)]
    assert sum(r["original_method"] == "spearman_tie_aware" for r in combined) == 66
    assert sum(r["original_method"] == "median_flag1_minus_flag0" for r in combined) == 12
    for result in combined[:78]:
        old = prior[result["comparison_id"]]
        assert result["n_paired"] == old["n_paired"]
        assert result["paired_variants"] == old["paired_variants"]
        assert result["original_effect"] == old["effect"]
        assert result["spearman_rho"] == (old["effect"] if old["method"] == "spearman_tie_aware" else "")
    blocks = combined[78:]
    assert sum(r["analysis_role"] == "primary_block" for r in blocks) == 9
    assert sum(r["analysis_role"] == "all_feature_sensitivity" for r in blocks) == 3
    assert all(r["spearman_rho"] == "" for r in blocks)
    assert [r["n_paired"] for r in blocks[-3:]] == ["11", "11", "9"]
    assert all(r["description"].count(".") == 3 for r in combined)
    manifest = json.loads((tmp_path / "dependence-run-manifest-v1.json").read_text())
    assert manifest["inferential_p_values_computed"] is False
    assert "/Users/" not in json.dumps(manifest) + json.dumps(plan)
    initial = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in tmp_path.iterdir()}
    assert_ok(cli("analyze", tmp_path))
    assert initial == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in tmp_path.iterdir()}


def test_cli_exact_shared_components_and_union_study_deletion(tmp_path: Path) -> None:
    assert_ok(cli("plan", tmp_path))
    assert_ok(cli("analyze", tmp_path))
    sensitivities = rows(tmp_path / "dependence-sensitivities-v1.csv")
    shared = next(r for r in sensitivities if r["comparison_id"] == "C001" and r["scenario"] == "exclude_shared_short_S")
    assert shared["n_paired"] == "7"
    assert set(shared["excluded_variants"].split("|")) == {"R249Q", "I457T", "D778V", "L781P", "S782N", "A797T", "F834L"}
    strict = next(r for r in sensitivities if r["comparison_id"] == "C002" and r["scenario"] == "strict_short_head_MVEL20")
    assert strict["n_paired"] == "4"
    assert strict["original_spearman_rho"] == "0.8"
    morck = next(r for r in sensitivities if r["comparison_id"] == "C002" and r["scenario"] == "exclude_morck_shared_WT")
    assert morck["n_paired"] == "9"
    deleted = next(r for r in sensitivities if r["comparison_id"] == "C001" and r["scenario"] == "leave_one_study_out" and r["omitted_unit"] == "adhikari2019")
    table = rows(SOURCE / "analysis-table-v1.csv")
    expected = {r["variant"] for r in table if r["kcat_rel"] and r["Na_rel"] and "adhikari2019" not in {r["kcat_rel_study"], r["Na_rel_study"]}}
    assert set(deleted["paired_variants"].split("|")) == expected
    assert not any(r["comparison_id"] == "C034" and r["scenario"] == "leave_one_study_out" for r in sensitivities)


def test_cli_rejects_plan_and_input_tampering_before_outputs(tmp_path: Path) -> None:
    source = tmp_path / "source"
    copy_inputs(source)
    output = tmp_path / "output"
    assert_ok(cli("plan", output, source))
    plan_path = output / "dependence-plan-v1.json"
    plan = json.loads(plan_path.read_text())
    plan["comparisons"][78]["x_columns"].append("wt_site_rsa_8act_mean_ab")
    plan_path.write_text(json.dumps(plan))
    failure = cli("analyze", output, source)
    assert failure.returncode == 2 and "tampering" in failure.stderr
    assert not (output / "combined-correlations-v1.csv").exists()
    assert_ok(cli("plan", output, source))
    ledger = source / "dependence-ledger-v1.csv"
    ledger.write_text(ledger.read_text().replace("same_paper_different_construct,", "new_note,"))
    failure = cli("analyze", output, source)
    assert failure.returncode == 2 and "frozen plan" in failure.stderr
    assert not (output / "combined-correlations-v1.csv").exists()


@pytest.mark.parametrize("mutation, message", [("duplicate", "duplicate"), ("nonfinite", "nonfinite"), ("malformed", "malformed"), ("wrong_count", "hash mismatch")])
def test_cli_rejects_invalid_tables(tmp_path: Path, mutation: str, message: str) -> None:
    source = tmp_path / "source"
    copy_inputs(source)
    path = source / "analysis-table-v1.csv"
    lines = path.read_text().splitlines()
    if mutation == "duplicate":
        lines.append(lines[1])
    elif mutation == "nonfinite":
        lines[1] = lines[1].replace("0.93,", "nan,", 1)
    elif mutation == "malformed":
        lines[1] += ",surplus"
    else:
        inv = source / "comparison-inventory-v1.csv"
        contents = inv.read_text().replace(",22,17,16,14,", ",22,17,16,99,", 1)
        inv.write_text(contents)
    path.write_text("\n".join(lines) + "\n")
    failure = cli("plan", tmp_path / "output", source)
    assert failure.returncode == 2 and message in failure.stderr
    assert not (tmp_path / "output/dependence-plan-v1.json").exists()


def test_committed_public_plan_reproduces_release_outside_repo(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "myh7_nonlinear_dependence.screen", "analyze",
         "--plan", str(screen.PLAN), "--output-dir", str(tmp_path)],
        cwd=tmp_path, env=os.environ | {"PYTHONPATH": str(ROOT / "src")},
        text=True, capture_output=True, check=False,
    )
    assert_ok(result)
    for name in ("combined-correlations-v1.csv", "dependence-sensitivities-v1.csv",
                 "dependence-points-v1.csv", "dependence-run-manifest-v1.json"):
        assert (tmp_path / name).read_bytes() == (FROZEN / name).read_bytes()

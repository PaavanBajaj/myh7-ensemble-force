"""The frozen CLI plan produces a complete, inspectable descriptive screen."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
WORKFLOW = ROOT / "workflows/myh7-correlation-screen"


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def test_analyze_cli_covers_frozen_pairs_and_dependence(tmp_path: Path) -> None:
    command = [
        sys.executable, str(WORKFLOW / "screen.py"), "analyze",
        "--table", str(WORKFLOW / "data/analysis-table-v1.csv"),
        "--inventory", str(WORKFLOW / "data/comparison-inventory-v1.csv"),
        "--plan-manifest", str(WORKFLOW / "data/comparison-plan-manifest-v1.json"),
        "--outcomes", str(WORKFLOW / "data/outcomes-v1.csv"),
        "--dependence-ledger", str(WORKFLOW / "data/dependence-ledger-v1.csv"),
        "--velocity-alternatives", str(WORKFLOW / "data/velocity-alternatives-v1.csv"),
        "--output-dir", str(tmp_path),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    assert completed.returncode == 0, completed.stderr
    results = _rows(tmp_path / "comparison-results-v1.csv")
    points = _rows(tmp_path / "comparison-points-v1.csv")
    sensitivity = _rows(tmp_path / "sensitivity-results-v1.csv")
    assert len(results) == 78
    assert {r["comparison_id"] for r in results} == {f"C{i:03d}" for i in range(1, 79)}
    assert len(list((tmp_path / "plots").glob("C*.svg"))) == 78
    assert len([p for p in points if p["comparison_id"] == "C001"]) == 14
    assert next(r for r in results if r["comparison_id"] == "C001")["n_paired"] == "14"
    assert next(r for r in sensitivity if r["comparison_id"] == "C001" and r["scenario"] == "exclude_shared_short_S")["n_paired"] == "7"
    assert next(r for r in sensitivity if r["comparison_id"] == "C001" and r["scenario"] == "exclude_shared_short_S_and_risky_Na")["n_paired"] == "4"
    assert any(r["scenario"] == "expanded_human_velocity" for r in sensitivity)
    assert "Y115H" in (tmp_path / "plots/C001.svg").read_text()
    manifest = json.loads((tmp_path / "analysis-manifest-v1.json").read_text())
    assert manifest["plan_sha256"]
    assert manifest["comparison_count"] == 78

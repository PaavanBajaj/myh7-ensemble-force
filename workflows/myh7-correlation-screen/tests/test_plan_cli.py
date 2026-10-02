"""Planning freezes all eligible pairs before effect estimates are produced."""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
WORKFLOW = ROOT / "workflows/myh7-correlation-screen"
PUBLIC = ROOT / "data/public/lsar-na-rel"


def test_plan_cli_inventories_complete_pair_set_before_analysis(tmp_path: Path) -> None:
    command = [
        sys.executable, str(WORKFLOW / "screen.py"), "plan",
        "--outcomes", str(WORKFLOW / "data/outcomes-v1.csv"),
        "--new-features", str(WORKFLOW / "data/features-v1.csv"),
        "--kcat-features", str(ROOT / "workflows/kcat-rel/data/derived/kcat-rel-v0-no-msa-features.csv"),
        "--na-features", str(PUBLIC / "features.csv"),
        "--config", str(WORKFLOW / "comparison-config-v1.json"),
        "--output-dir", str(tmp_path),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    assert completed.returncode == 0, completed.stderr
    with (tmp_path / "analysis-table-v1.csv").open(newline="") as stream:
        table = list(csv.DictReader(stream))
    with (tmp_path / "comparison-inventory-v1.csv").open(newline="") as stream:
        inventory = list(csv.DictReader(stream))
    assert len(table) == 22
    assert len(inventory) == 78  # 3 outcome pairs + 30 outcome-feature + 45 feature pairs
    assert {row["pair_type"] for row in inventory} == {"outcome_outcome", "outcome_feature", "feature_feature"}
    kcat_na = next(row for row in inventory if row["x"] == "kcat_rel" and row["y"] == "Na_rel")
    assert int(kcat_na["n_paired"]) == 14
    v_elc = next(row for row in inventory if row["x"] == "v_rel" and row["y"] == "elc_min_heavy_distance_8act_mean_ab")
    assert int(v_elc["n_paired"]) == 15
    assert "A797T" in v_elc["paired_variants"]
    assert not (tmp_path / "comparison-results-v1.csv").exists()

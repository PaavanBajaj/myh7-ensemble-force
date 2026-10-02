"""The outcome command preserves published pairings and canonical ratios."""

from __future__ import annotations

import csv
import os
import math
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / "workflows" / "myh7-correlation-screen"
DATA = ROOT / "data/public/myh7-correlation-screen"
PUBLIC = ROOT / "data/public/lsar-na-rel"


def test_outcomes_cli_keeps_absolute_pairs_and_within_study_na(tmp_path: Path) -> None:
    output = tmp_path / "out"
    command = [
        sys.executable, "-m", "myh7_correlation_screen.screen", "outcomes",
        "--kcat-labels", str(ROOT / "data/public/kcat-rel/canonical-labels.csv"),
        "--kcat-evidence", str(ROOT / "data/public/kcat-rel/measurement-evidence.csv"),
        "--na-primary", str(PUBLIC / "na-rel-primary.csv"),
        "--na-labels", str(PUBLIC / "labels.csv"),
        "--na-measurements", str(PUBLIC / "measurements.csv"),
        "--velocity-summary", str(DATA / "velocity-summary-evidence-v1.csv"),
        "--morck-pairs", str(DATA / "morck-velocity-pairs-v1.csv"),
        "--output-dir", str(output),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False, cwd=ROOT, env={**os.environ, "PYTHONPATH": str(ROOT / "src")})
    assert completed.returncode == 0, completed.stderr
    with (output / "outcomes-v1.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 51  # 17 kcat + 16 Na + 15 Table 1 velocity + 3 expanded
    by_key = {(row["variant"], row["outcome"]): row for row in rows}
    assert math.isclose(float(by_key["Y115H", "v_rel"]["value"]), 406 / 903, rel_tol=1e-11)
    assert by_key["Y115H", "v_rel"]["wt_value"] == "903"
    assert by_key["Y115H", "v_rel"]["error_type"] == "SEM"
    assert math.isclose(float(by_key["Y115H", "Na_rel"]["value"]), 0.9 / 0.63, rel_tol=1e-11)
    assert math.isclose(float(by_key["D778V", "v_rel"]["value"]), 1.458718421, rel_tol=1e-8)
    assert by_key["D778V", "v_rel"]["aggregation_method"] == "mean_of_paired_slide_ratios"
    assert by_key["M493I", "v_rel"]["cohort"] == "expanded_sensitivity"
    assert by_key["M493I", "v_rel"]["uncertainty"] == ""
    assert by_key["M493I", "v_rel"]["uncertainty_method"] == "not_available_distributional_SD"
    with (output / "na-components-v1.csv").open(newline="") as stream:
        components = list(csv.DictReader(stream))
    assert len(components) == 21
    assert not any(row["variant"] == "G256E" for row in components)
    r249q = next(row for row in components if row["variant"] == "R249Q")
    assert r249q["short_kcat_s_1"] == "1.69"
    assert r249q["long_kcat_s_1"] == "1.63"

"""The feature command is the public, label-blind seam."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
WORKFLOW = ROOT / "workflows" / "myh7-correlation-screen"


def test_feature_cli_keeps_deposited_missingness_and_full_8act_coverage(tmp_path: Path) -> None:
    variants = tmp_path / "variants.csv"
    variants.write_text("variant\nY115H\nS782N\nA797T\nF834L\n")
    output = tmp_path / "out"
    command = [
        sys.executable,
        str(WORKFLOW / "screen.py"),
        "features",
        "--variants",
        str(variants),
        "--fasta",
        str(ROOT / "tests/lsar_na_rel/fixtures/P12883.fasta"),
        "--assembly-8act",
        str(ROOT / "workflows/lsar-na-rel/cache/8ACT-assembly1.pdb"),
        "--cif-8efe",
        str(WORKFLOW / "cache/8EFE.cif"),
        "--cif-8efd",
        str(WORKFLOW / "cache/8EFD.cif"),
        "--cif-8efi",
        str(WORKFLOW / "cache/8EFI.cif"),
        "--output-dir",
        str(output),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    assert completed.returncode == 0, completed.stderr

    with (output / "features-v1.csv").open(newline="") as handle:
        rows = {row["variant"]: row for row in csv.DictReader(handle)}
    assert set(rows) == {"Y115H", "S782N", "A797T", "F834L"}
    assert rows["Y115H"]["adp_rigor_ca_shift_A"]
    assert rows["S782N"]["adp_rigor_ca_shift_A"]
    assert rows["S782N"]["rigor_actin_min_heavy_distance_8efi_A"] == ""
    assert rows["A797T"]["adp_rigor_ca_shift_A"] == ""
    assert rows["F834L"]["adp_adenine_ca_distance_8efe_A"] == ""
    assert rows["F834L"]["relay_landmark_min_ca_distance_8act_mean_ab"]
    assert rows["F834L"]["elc_min_heavy_distance_8act_mean_ab"]
    # A797 lies in the light-chain-binding lever; a 50+ Å result indicates
    # that the deposited ELC chains were paired to the wrong myosin heads.
    assert float(rows["A797T"]["elc_min_heavy_distance_8act_mean_ab"]) < 10.0

    manifest = json.loads((output / "feature-manifest-v1.json").read_text())
    assert manifest["feature_generation_uses_outcomes"] is False
    assert len(manifest["features"]) == 5
    assert manifest["structures"]["8EFD"]["sha256"] == "fed8cd752b452399b8409ee72dfa8785625f3a6cd9fb5919f8e961aa94e12e58"

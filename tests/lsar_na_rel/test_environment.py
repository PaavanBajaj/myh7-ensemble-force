"""Exact osx-arm64 environment release contract."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_ROOT = REPO_ROOT / "workflows" / "lsar-na-rel"
LOCK = REPO_ROOT / "environment-osx-arm64.lock"


def test_explicit_lock_pins_osx_arm64_runtime():
    text = LOCK.read_text(encoding="utf-8")
    assert text.startswith("# This file may be used to create an environment using:\n")
    assert "# platform: osx-arm64" in text
    assert "@EXPLICIT" in text
    for package in (
        "python-3.12.13-",
        "numpy-2.5.2-",
        "pandas-3.0.5-",
        "scipy-1.18.0-",
        "scikit-learn-1.9.0-",
        "matplotlib-3.11.1-",
        "biopython-1.88-",
        "dssp-4.6.1-",
        "freesasa-2.2.1-",
        "pytest-9.1.1-",
    ):
        assert package in text
    package_lines = [line for line in text.splitlines() if line.startswith("https://")]
    assert package_lines
    assert all("#" in line for line in package_lines)

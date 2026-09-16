"""Portable paths for the public LSAR Na-rel repository layout."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


def discover_repo_root(source_file: str | Path | None = None) -> Path:
    """Resolve the checkout root from public or private workflow source layout."""
    source = Path(__file__ if source_file is None else source_file).resolve()
    package_dir = source.parent
    src_dir = package_dir.parent
    candidate = src_dir.parent
    if candidate.name == "lsar-na-rel" and candidate.parent.name == "workflows":
        return candidate.parent.parent.resolve()
    return candidate.resolve()


@dataclass(frozen=True)
class RepoPaths:
    """All default workflow paths rooted in one checkout."""

    repo_root: Path

    @classmethod
    def from_repo_root(cls, root: str | Path) -> "RepoPaths":
        return cls(repo_root=Path(root).resolve())

    @property
    def workflow_root(self) -> Path:
        return self.repo_root / "workflows" / "lsar-na-rel"

    @property
    def public_data(self) -> Path:
        return self.repo_root / "data" / "public" / "lsar-na-rel"

    @property
    def config_path(self) -> Path:
        return self.workflow_root / "config" / "lsar-na-rel.json"

    @property
    def cache_dir(self) -> Path:
        return self.workflow_root / "cache"

    @property
    def results_dir(self) -> Path:
        return self.workflow_root / "results"

    @property
    def reference_fasta(self) -> Path:
        return self.repo_root / "tests" / "lsar_na_rel" / "fixtures" / "P12883.fasta"

    @property
    def features_path(self) -> Path:
        return self.public_data / "features.csv"

    @property
    def features_provenance_path(self) -> Path:
        return self.public_data / "features-provenance.json"

    @property
    def primary_path(self) -> Path:
        return self.public_data / "na-rel-primary.csv"

    @property
    def labels_path(self) -> Path:
        return self.public_data / "labels.csv"


DEFAULT_PATHS = RepoPaths.from_repo_root(discover_repo_root())

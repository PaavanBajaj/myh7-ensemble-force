"""Repository-local path resolution for a portable public checkout."""

from pathlib import Path

from lsar_na_rel.paths import RepoPaths, discover_repo_root


def test_repo_paths_are_local_to_arbitrarily_named_checkout(tmp_path: Path):
    root = tmp_path / "renamed-checkout"
    paths = RepoPaths.from_repo_root(root)

    assert paths.repo_root == root.resolve()
    assert paths.workflow_root == root.resolve() / "workflows" / "lsar-na-rel"
    assert paths.public_data == root.resolve() / "data" / "public" / "lsar-na-rel"
    assert paths.config_path == paths.workflow_root / "config" / "lsar-na-rel.json"
    assert paths.cache_dir == paths.workflow_root / "cache"
    assert paths.results_dir == paths.workflow_root / "results"
    assert paths.reference_fasta == root.resolve() / "tests" / "lsar_na_rel" / "fixtures" / "P12883.fasta"


def test_discover_repo_root_supports_public_src_layout(tmp_path: Path):
    root = tmp_path / "any-name"
    source = root / "src" / "lsar_na_rel" / "paths.py"
    source.parent.mkdir(parents=True)
    source.write_text("", encoding="utf-8")

    assert discover_repo_root(source) == root.resolve()


def test_discover_repo_root_supports_private_workflow_src_layout(tmp_path: Path):
    root = tmp_path / "private-any-name"
    source = root / "workflows" / "lsar-na-rel" / "src" / "lsar_na_rel" / "paths.py"
    source.parent.mkdir(parents=True)
    source.write_text("", encoding="utf-8")

    assert discover_repo_root(source) == root.resolve()

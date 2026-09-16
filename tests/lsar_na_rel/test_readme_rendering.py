"""README math must use GitHub-supported dollar delimiters."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]


def test_readmes_do_not_use_unsupported_backslash_math_delimiters():
    offenders: list[str] = []
    for path in REPO_ROOT.rglob("README.md"):
        if any(part.startswith(".") for part in path.relative_to(REPO_ROOT).parts):
            continue
        text = path.read_text(encoding="utf-8")
        if any(token in text for token in (r"\(", r"\)", r"\[", r"\]")):
            offenders.append(path.relative_to(REPO_ROOT).as_posix())
    assert offenders == []

"""Static checks for the hand-maintained documentation navigation.

These tests need neither Sphinx nor native SPR-KKR binaries.
"""

from pathlib import Path
import re
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "doc_src"
DOC_ROLE = re.compile(r":doc:`([^`]+)`")
REF_ROLE = re.compile(r":ref:`([^`]+)`")
LABEL = re.compile(r"^\.\. _([\w-]+):", re.MULTILINE)
LOCAL_RST = re.compile(r"<((?:doc_src/)?[\w-]+\.rst)(#[\w-]+)?>")
README_LINK = re.compile(r"\]\((doc_src/[^)]+)\)")
SPECIAL_REFS = {"genindex", "modindex", "search"}


def _target(value):
    return value.split("<")[-1].removesuffix(">").strip()


def test_rst_pages_and_cross_references_exist():
    pages = {path.stem: path for path in DOCS.glob("*.rst")}
    labels = set()
    for path in pages.values():
        labels.update(LABEL.findall(path.read_text(encoding="utf-8")))

    for path in pages.values():
        source = path.read_text(encoding="utf-8")
        for reference in DOC_ROLE.findall(source):
            target = _target(reference).removesuffix(".rst")
            assert target in pages, (path.name, "missing document", target)
        for reference in REF_ROLE.findall(source):
            target = _target(reference)
            assert target in labels or target in SPECIAL_REFS, (
                path.name, "missing reference", target
            )
        for reference, _fragment in LOCAL_RST.findall(source):
            target = ROOT / reference if reference.startswith("doc_src/") else DOCS / reference
            assert target.is_file(), (path.name, "missing local link", reference)


def test_readme_file_links_and_anchors():
    source = (ROOT / "README.md").read_text(encoding="utf-8")
    for link in README_LINK.findall(source):
        filename, _, anchor = link.partition("#")
        target = ROOT / filename
        assert target.is_file(), link
        if anchor:
            lines = target.read_text(encoding="utf-8").splitlines()
            headings = (
                line.lower().strip().replace(" ", "-")
                for index, line in enumerate(lines[:-1])
                if re.fullmatch(r"-{3,}", lines[index + 1])
            )
            assert anchor in headings, (link, "missing heading anchor")


def test_rst_heading_underlines_are_not_shorter_than_titles():
    for path in DOCS.glob("*.rst"):
        lines = path.read_text(encoding="utf-8").splitlines()
        for index, line in enumerate(lines[1:], 1):
            if re.fullmatch(r"([=~-])\1{2,}", line) and lines[index - 1].strip():
                assert len(line) >= len(lines[index - 1]), (
                    path.name, index + 1, lines[index - 1]
                )


@pytest.mark.skipif(shutil.which("pandoc") is None, reason="Pandoc is not installed")
def test_readme_is_synchronized_with_rst_source():
    """Catch drift between the canonical RST and the GitHub-facing Markdown."""
    rendered = subprocess.run(
        [
            "pandoc", "-f", "rst", "-t", "gfm", "--wrap=none",
            "doc_src/README.rst",
        ],
        cwd=ROOT,
        capture_output=True,
        check=True,
        text=True,
    ).stdout
    assert rendered == (ROOT / "README.md").read_text(encoding="utf-8")

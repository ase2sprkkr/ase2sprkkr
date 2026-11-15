"""Guard the public documentation's introductory links and positioning."""

from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "path",
    (
        "doc_src/introduction.rst",
        "doc_src/README.rst",
    ),
)
def test_public_introduction_uses_clickable_rst_links(path):
    text = (ROOT / path).read_text(encoding="utf-8")
    assert "modern, comprehensive replacement" in text
    assert "`SPR-KKR <https://sprkkr.org/>`_" in text
    assert "`XBAND <https://sprkkr.org/xband/>`_" in text
    assert "`Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/>`_" in text


def test_markdown_readme_links_match_rst():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "[SPR-KKR](https://sprkkr.org/)" in text
    assert "[XBAND](https://sprkkr.org/xband/)" in text
    assert "[Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/)" in text

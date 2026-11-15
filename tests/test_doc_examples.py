"""Tests for the source-driven documentation examples gallery."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "doc_src" / "examples.py"
SPEC = importlib.util.spec_from_file_location("docs_example_generator", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
generator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(generator)


def test_gallery_discovers_sources_without_importing_them(tmp_path):
    source = tmp_path / "examples"
    later = source / "A03_DOS"
    earlier = source / "A01_Fe"
    later.mkdir(parents=True)
    earlier.mkdir(parents=True)
    (later / "dos.py").write_text(
        '"""DOS calculation.\n\nMore information."""\n'
        "import intentionally_nonexistent_package\n", encoding="utf-8"
    )
    (earlier / "scf.py").write_text(
        '"""SCF example.\nA second line."""\n', encoding="utf-8"
    )
    (earlier / "__init__.py").write_text('"""Package init."""\n', encoding="utf-8")
    (source / "not_an_example").mkdir()

    rst = generator.render_gallery(source)
    assert rst.index("A01_Fe") < rst.index("A03_DOS")
    assert "SCF example. A second line." in rst
    assert "DOS calculation." in rst
    assert "More information." not in rst
    assert "intentionally_nonexistent_package" not in rst
    assert rst.count(":download:") == 2
    assert rst.count(".. literalinclude::") == 2
    assert ".. literalinclude:: ../src/ase2sprkkr/examples/A03_DOS/dos.py" in rst
    assert "../src/ase2sprkkr/examples/A03_DOS/dos.py" in rst
    assert "ase2sprkkr example 3 --copy" in rst


def test_gallery_generation_and_staleness_check(tmp_path):
    source = tmp_path / "examples"
    script = source / "A02_Sample" / "test_case.py"
    script.parent.mkdir(parents=True)
    script.write_text('"""An example."""\n', encoding="utf-8")
    output = tmp_path / "gallery.inc"

    with pytest.raises(RuntimeError, match="missing or outdated"):
        generator.generate_gallery(source, output, check=True)

    generator.generate_gallery(source, output)
    generator.generate_gallery(source, output, check=True)
    assert output.read_text(encoding="utf-8") == generator.render_gallery(source)

    script.write_text('"""Changed example."""\n', encoding="utf-8")
    with pytest.raises(RuntimeError, match="missing or outdated"):
        generator.generate_gallery(source, output, check=True)


def test_real_examples_are_all_indexed():
    source = ROOT / "src" / "ase2sprkkr" / "examples"
    expected = sum(
        bool(generator.EXAMPLE_NAME.fullmatch(folder.name))
        for folder in source.iterdir()
        if folder.is_dir()
        for script in folder.glob("*.py")
        if not script.name.startswith("_")
    )
    rst = generator.render_gallery(source)
    assert expected > 0
    assert rst.count(":download:") == expected
    assert rst.count(".. literalinclude::") == expected

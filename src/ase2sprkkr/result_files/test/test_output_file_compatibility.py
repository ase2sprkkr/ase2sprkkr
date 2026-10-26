import importlib
import re
from pathlib import Path

import pytest


def test_top_level_output_file_alias_warns():
    import ase2sprkkr

    with pytest.warns(
        DeprecationWarning,
        match=r"ase2sprkkr\.OutputFile is deprecated",
    ):
        output_file = ase2sprkkr.OutputFile
    assert output_file is ase2sprkkr.ResultFile


def test_internal_output_file_alias_is_result_file():
    from ase2sprkkr import ResultFile
    from ase2sprkkr.result_files.output_files import OutputFile

    assert OutputFile is ResultFile

def test_old_output_files_package_is_deprecated():
    with pytest.warns(DeprecationWarning):
        importlib.reload(
            importlib.import_module("ase2sprkkr.output_files")
        )

def test_canonical_result_file_api_does_not_use_legacy_names():
    package_root = Path(__file__).resolve().parents[2]
    allowed = {
        package_root / "__init__.py",
        package_root / "result_files" / "output_files.py",
        Path(__file__).resolve(),
    }
    failures = []
    for path in package_root.rglob("*.py"):
        if path in allowed:
            continue
        text = path.read_text()
        if "ase2sprkkr.output_files." in text:
            failures.append((path, "ase2sprkkr.output_files"))
        if re.search(
            r"from\s+[.\w]*output_files(?:\s+import|\.)", text
        ):
            failures.append((path, "legacy output_files import"))
        if re.search(r"\bOutputFile\b", text):
            failures.append((path, "OutputFile"))
    assert not failures, failures

import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    ("statement", "expected_warning"),
    [
        (
            "from ase2sprkkr import OutputFile, ResultFile",
            "ase2sprkkr.OutputFile is deprecated; use ase2sprkkr.ResultFile instead.",
        ),
        (
            "from ase2sprkkr.output_files import OutputFile; from ase2sprkkr import ResultFile",
            "ase2sprkkr.output_files is deprecated; use ase2sprkkr.result_files instead.",
        ),
        (
            "from ase2sprkkr.output_files.output_files import OutputFile; from ase2sprkkr import ResultFile",
            "ase2sprkkr.output_files is deprecated; use ase2sprkkr.result_files instead.",
        ),
    ],
)
def test_deprecated_output_file_imports_warn_and_alias_result_file(statement, expected_warning):
    code = f"""
import warnings
warnings.simplefilter('always')
with warnings.catch_warnings(record=True) as caught:
    {statement}
assert OutputFile is ResultFile
matching = [item for item in caught if str(item.message) == {expected_warning!r}]
assert len(matching) == 1, [(item.category.__name__, str(item.message)) for item in caught]
assert issubclass(matching[0].category, DeprecationWarning)
"""
    subprocess.run([sys.executable, "-c", code], check=True)


def test_canonical_result_file_api_does_not_use_deprecated_package():
    code = """
import sys
import warnings
warnings.simplefilter('always')
with warnings.catch_warnings(record=True) as caught:
    import ase2sprkkr
    from ase2sprkkr import ResultFile
    from ase2sprkkr.result_files import ResultFile as PackageResultFile
    assert ResultFile is PackageResultFile
    assert ResultFile.definitions
assert 'ase2sprkkr.output_files' not in sys.modules
unexpected = [
    item for item in caught
    if 'OutputFile is deprecated' in str(item.message)
    or 'ase2sprkkr.output_files is deprecated' in str(item.message)
]
assert not unexpected, [(item.category.__name__, str(item.message)) for item in unexpected]
"""
    subprocess.run([sys.executable, "-c", code], check=True)

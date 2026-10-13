"""Public terminology for persistent result-file metadata."""

from ase2sprkkr.execution import (
    ArtifactRecord,
    ArtifactState,
    ResultFilesRecord,
    ResultFileState,
)


def test_result_file_names_are_canonical_compatible_aliases():
    assert ResultFilesRecord is ArtifactRecord
    assert ResultFileState is ArtifactState
    assert ResultFileState.CURRENT.value == "current"
    assert ResultFileState.OVERWRITTEN.value == "overwritten"

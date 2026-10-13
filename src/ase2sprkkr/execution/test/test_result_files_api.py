"""Public terminology for persistent result-file metadata."""

from ase2sprkkr.execution import ResultFilesRecord, ResultFilesState


def test_result_file_names_use_canonical_public_api():
    record = ResultFilesRecord()

    assert record.state is ResultFilesState.CURRENT
    assert ResultFilesState.CURRENT.value == "current"
    assert ResultFilesState.OVERWRITTEN.value == "overwritten"

"""Public terminology and capabilities for persistent result files."""

from ase2sprkkr.execution import Executor, ResultFilesRecord, ResultFilesState


def test_result_file_names_use_canonical_public_api():
    record = ResultFilesRecord()

    assert record.state is ResultFilesState.CURRENT
    assert ResultFilesState.CURRENT.value == "current"
    assert ResultFilesState.OVERWRITTEN.value == "overwritten"


def test_result_file_retrieval_capability_is_inferred_from_override():
    class PlainExecutor(Executor):
        pass

    class RetrievingExecutor(Executor):
        def fetch_result_files(self, reference, files):
            pass

    class DerivedRetrievingExecutor(RetrievingExecutor):
        pass

    assert PlainExecutor.supports_result_file_retrieval is False
    assert PlainExecutor().supports_result_file_retrieval is False
    assert RetrievingExecutor.supports_result_file_retrieval is True
    assert RetrievingExecutor().supports_result_file_retrieval is True
    assert DerivedRetrievingExecutor.supports_result_file_retrieval is True
    assert DerivedRetrievingExecutor().supports_result_file_retrieval is True

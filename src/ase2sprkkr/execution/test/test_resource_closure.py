"""Shared close semantics across output files and owned resources."""

import pytest

from ase2sprkkr.execution.base import ExecutionResources, close_output_files


class _Closable:
    def __init__(self, name, calls, error=None):
        self.name = name
        self.calls = calls
        self.error = error
        self.closed = False

    def close(self):
        self.calls.append(self.name)
        self.closed = True
        if self.error is not None:
            raise self.error


def test_close_output_files_attempts_all_after_error_and_deduplicates():
    calls = []
    error = OSError("first output refused close")
    first = _Closable("stdout", calls, error)
    second = _Closable("stderr", calls)
    with pytest.raises(OSError, match="first output refused close") as caught:
        close_output_files(first, first, second, None, "/tmp/a-path")
    assert caught.value is error
    assert calls == ["stdout", "stderr"]
    assert first.closed and second.closed


def test_execution_resources_close_attempts_all_and_preserves_first_failure():
    calls = []
    error = ValueError("broken close")
    first = _Closable("first", calls, error)
    second = _Closable("second", calls)
    resources = ExecutionResources(first, first, second)
    with pytest.raises(ValueError, match="broken close") as caught:
        resources.close()
    assert caught.value is error
    assert resources.closed
    resources.close()
    assert calls == ["first", "second"]


def test_execution_resources_can_suppress_close_errors_and_skip_preclosed():
    calls = []
    already = _Closable("preclosed", calls)
    already.closed = True
    failure = _Closable("failure", calls, RuntimeError("ignored"))
    final = _Closable("final", calls)
    resources = ExecutionResources(already, failure, final)
    resources.close(suppress_errors=True)
    assert calls == ["failure", "final"]
    assert resources.closed

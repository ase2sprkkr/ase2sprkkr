from types import SimpleNamespace

import pytest

import ase2sprkkr.execution.local_persistent as local_persistent
from ase2sprkkr.execution import ExecutionSpec, LocalExecutor


class _Reader:
    print_output = False

    def create_result(self, output_file, input_file):
        return SimpleNamespace(output_file=output_file, input_file=input_file)


def test_local_persistent_programming_error_is_not_submission_failure(
    tmp_path, monkeypatch
):
    runner = LocalExecutor().create_persistent_run(
        _Reader(),
        ExecutionSpec(
            argv=("unused",),
            input_file=None,
            output_file=str(tmp_path / "calculation.out"),
            directory=str(tmp_path),
        ),
        run_id="programming-error",
    )

    def fail_with_programming_error(request):
        raise ValueError("invalid detached request")

    monkeypatch.setattr(
        local_persistent,
        "submit_detached",
        fail_with_programming_error,
    )

    with pytest.raises(ValueError, match="invalid detached request"):
        runner.submit()

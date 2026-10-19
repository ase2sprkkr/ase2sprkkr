"""Explicit stream persistence must not depend on the other stream."""

import sys

import pytest

from ..base import ExecutionSpec, RunState
from ..local import LocalExecutor


@pytest.mark.parametrize(
    "store_stdout, store_stderr",
    ((True, True), (True, False), (False, True), (False, False)),
)
def test_local_stream_targets_are_independent(tmp_path, store_stdout, store_stderr):
    stdout_path = tmp_path / "stdout.txt"
    stderr_path = tmp_path / "stderr.txt"
    stdout = stdout_path.open("wb") if store_stdout else None
    stderr = stderr_path.open("wb") if store_stderr else None

    runner = LocalExecutor().create_run(ExecutionSpec(
        argv=(sys.executable, "-c",
              "import sys; print('from stdout'); print('from stderr', file=sys.stderr)"),
        input_file=None,
        output_file=stdout,
        error_file=stderr,
        directory=str(tmp_path),
    ))
    observed = {"out": [], "err": []}
    runner.add_output_listener(lambda data, kind: observed[kind].append(data))
    assert runner.run().return_code == 0
    assert runner.state == RunState.FINISHED

    assert "".join(observed["out"]) == "from stdout\n"
    assert "".join(observed["err"]) == "from stderr\n"
    assert stdout_path.exists() is store_stdout
    assert stderr_path.exists() is store_stderr
    if store_stdout:
        assert stdout_path.read_bytes() == b"from stdout\n"
        assert stdout.closed
    if store_stderr:
        assert stderr_path.read_bytes() == b"from stderr\n"
        assert stderr.closed

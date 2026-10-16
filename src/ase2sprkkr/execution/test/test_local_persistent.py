import os
import sys
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import ExecutionSpec, LocalExecutor, RunState


pytestmark = pytest.mark.skipif(
    os.name == "nt",
    reason="the detached helper currently uses POSIX file locking",
)


class _Reader:
    print_output = False

    def create_result(self, output_file, input_file):
        return SimpleNamespace(
            output_file=getattr(output_file, "name", output_file),
            input_file=getattr(input_file, "name", input_file),
            stdout=None,
            stderr=None,
            return_code=None,
        )

    async def parse_streams(self, stdout, stderr, result):
        return await stdout.read(), await stderr.read()

    def complete_result(self, result, output, error, return_code):
        result.stdout = output
        result.stderr = error
        result.return_code = return_code
        return result


def _spec(tmp_path):
    return ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            "import time; print('before', flush=True); "
            "time.sleep(0.1); print('after', flush=True)",
        ),
        input_file=None,
        output_file=str(tmp_path / "calculation.out"),
        directory=str(tmp_path),
    )


def test_persistent_local_run_can_be_reattached_and_replayed(tmp_path):
    executor = LocalExecutor()
    runner = executor.create_persistent_run(
        _Reader(), _spec(tmp_path), run_id="local-persistent"
    )
    runner.submit()
    reference = runner.reference
    assert reference.native_id == "local-persistent"

    result = runner.wait()
    assert result.return_code == 0
    assert result.stdout == b"before\nafter\n"
    assert executor.status(reference)[0].state == RunState.FINISHED

    # Reattach as a new client. Replaying the helper-owned stream must rebuild
    # even a stale client-side output artifact and deliver all output to the
    # listener installed atomically during attach.
    (tmp_path / "calculation.out").write_text("stale\n")
    replay = []
    attached = LocalExecutor().attach(
        reference,
        _Reader(),
        output_listener=lambda text, kind: replay.append((text, kind)),
        start=False,
    )
    assert replay == []

    attached.start()
    attached_result = attached.wait()

    assert attached_result.return_code == 0
    assert attached_result.stdout == b"before\nafter\n"
    assert (tmp_path / "calculation.out").read_bytes() == b"before\nafter\n"
    assert "".join(text for text, kind in replay if kind == "out") == "before\nafter\n"


def test_persistent_local_submission_can_be_recovered_by_run_id(tmp_path):
    executor = LocalExecutor()
    runner = executor.create_persistent_run(
        _Reader(), _spec(tmp_path), run_id="recover-local"
    )
    runner.submit()

    recovered = LocalExecutor().recover_submission(
        "recover-local", SimpleNamespace(directory=str(tmp_path))
    )

    assert recovered is not None
    assert recovered.reference.executor == "local"
    assert recovered.reference.native_id == "recover-local"
    assert recovered.reference.output_file == str(tmp_path / "calculation.out")
    assert recovered.state in {
        RunState.SUBMITTING,
        RunState.RUNNING,
        RunState.FINISHED,
    }
    runner.wait()

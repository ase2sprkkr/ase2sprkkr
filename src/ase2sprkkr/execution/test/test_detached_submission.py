import os
import sys
from types import SimpleNamespace

import pytest

import ase2sprkkr.execution.detached as detached
import ase2sprkkr.execution.detached_local as detached_local
import ase2sprkkr.execution.local_persistent as local_persistent
from ase2sprkkr.execution import (
    ExecutionManager,
    ExecutionSpec,
    FailurePhase,
    LocalExecutor,
    RunState,
    RunStore,
)
from ase2sprkkr.execution.base import SubmissionNotStarted
from ase2sprkkr.execution.detached import (
    run_detached_worker,
    status_detached,
    submit_detached,
)


pytestmark = pytest.mark.skipif(
    os.name == "nt",
    reason="the detached helper currently uses POSIX file locking",
)


class _Reader:
    print_output = False

    def __init__(self):
        self.parse_calls = 0
        self.complete_calls = 0

    def create_result(self, output_file, input_file):
        return SimpleNamespace(
            output_file=getattr(output_file, "name", output_file),
            input_file=getattr(input_file, "name", input_file),
            return_code=None,
        )

    async def parse_streams(self, stdout, stderr, result):
        self.parse_calls += 1
        return await stdout.read(), await stderr.read()

    def complete_result(self, result, output, error, return_code):
        self.complete_calls += 1
        result.return_code = return_code
        return result


def _popen_failure(*args, **kwargs):
    raise OSError("worker could not be launched")


async def _process_start_failure(*args, **kwargs):
    raise OSError("calculation could not be launched")


def test_detached_helper_launch_failure_is_terminal_without_output(tmp_path, monkeypatch):
    output = tmp_path / "calculation.out"
    monkeypatch.setattr(detached, "_spawn_detached_worker", _popen_failure)

    state = submit_detached(
        {
            "id": "launch-failure",
            "executor": "local",
            "directory": str(tmp_path),
            "argv": ["unused"],
            "output_file": str(output),
        }
    )

    assert state["state"] == "failed"
    assert state["failure_phase"] == "launch"
    assert state["finished_at"] is not None
    assert "worker could not be launched" in state["error"]
    assert not output.exists()
    assert status_detached(tmp_path, "launch-failure") == state


def test_detached_helper_exit_before_handoff_is_terminal_launch_failure(
    tmp_path, monkeypatch
):
    class ExitedWorker:
        def poll(self):
            return 7

    monkeypatch.setattr(
        detached,
        "_spawn_detached_worker",
        lambda *args: ExitedWorker(),
    )

    state = submit_detached(
        {
            "id": "helper-exit-before-handoff",
            "executor": "local",
            "directory": str(tmp_path),
            "argv": ["unused"],
        }
    )

    assert state["state"] == "failed"
    assert state["failure_phase"] == "launch"
    assert state["started_at"] is None
    assert "return code 7" in state["error"]


def test_calculation_process_failure_does_not_publish_or_overwrite_output(
    tmp_path, monkeypatch
):
    output = tmp_path / "calculation.out"
    output.write_text("previous result\n")

    monkeypatch.setattr(detached, "_spawn_detached_worker", lambda *args: object())
    submit_detached(
        {
            "id": "calculation-launch-failure",
            "executor": "local",
            "directory": str(tmp_path),
            "argv": ["unused"],
            "output_file": str(output),
        }
    )

    monkeypatch.setattr(
        detached_local,
        "start_local_process",
        _process_start_failure,
    )
    with pytest.raises(OSError, match="calculation could not be launched"):
        run_detached_worker(tmp_path, "calculation-launch-failure")

    state = status_detached(tmp_path, "calculation-launch-failure")
    assert state["state"] == "failed"
    assert state["failure_phase"] == "launch"
    assert output.read_text() == "previous result\n"


def test_running_state_persistence_failure_stops_process_and_preserves_output(
    tmp_path, monkeypatch
):
    output = tmp_path / "calculation.out"
    output.write_text("previous result\n")
    run_id = "running-state-failure"
    monkeypatch.setattr(detached, "_spawn_detached_worker", lambda *args: object())
    submit_detached(
        {
            "id": run_id,
            "executor": "local",
            "directory": str(tmp_path),
            "argv": [sys.executable, "-c", "import time; time.sleep(30)"],
            "output_file": str(output),
        }
    )

    processes = []
    real_start = detached_local.start_local_process
    real_write = detached._write_state_unlocked

    async def record_start(*args, **kwargs):
        process = await real_start(*args, **kwargs)
        processes.append(process)
        return process

    def fail_running_state(workspace, current_run_id, **changes):
        if changes.get("state") == "running":
            raise OSError("running state write failed")
        return real_write(workspace, current_run_id, **changes)

    monkeypatch.setattr(detached_local, "start_local_process", record_start)
    monkeypatch.setattr(detached, "_write_state_unlocked", fail_running_state)

    with pytest.raises(OSError, match="running state write failed"):
        run_detached_worker(tmp_path, run_id)

    assert len(processes) == 1
    assert processes[0].returncode is not None
    assert processes[0].stdout is None
    assert processes[0].stderr is None
    assert output.read_text() == "previous result\n"

    state = status_detached(tmp_path, run_id)
    assert state["state"] == "failed"
    assert state["failure_phase"] == "launch"
    assert state["pid"] is None
    assert state["started_at"] is None


def test_finalization_sync_failure_preserves_private_stream(tmp_path, monkeypatch):
    run_id = "finalization-sync-failure"
    output = tmp_path / "calculation.out"
    monkeypatch.setattr(detached, "_spawn_detached_worker", lambda *args: object())
    submit_detached(
        {
            "id": run_id,
            "executor": "local",
            "directory": str(tmp_path),
            "argv": [sys.executable, "-c", "print('final output')"],
            "output_file": str(output),
        }
    )

    def fail_sync(_stream):
        raise OSError("stream fsync failed")

    monkeypatch.setattr(detached_local, "sync_stream", fail_sync)

    with pytest.raises(OSError, match="stream fsync failed"):
        run_detached_worker(tmp_path, run_id)

    state = status_detached(tmp_path, run_id)
    assert state["state"] == "failed"
    assert state["failure_phase"] == FailurePhase.FINALIZATION.value
    stream = detached.detached_stream_file(state, "output_file")
    assert stream == state["stream_output_file"]
    assert os.path.isfile(stream)
    with open(stream) as file:
        assert file.read() == "final output\n"
    assert not output.exists()


def test_partial_publish_failure_keeps_each_stream_replayable(tmp_path, monkeypatch):
    run_id = "partial-publish-failure"
    output = tmp_path / "calculation.out"
    error = tmp_path / "calculation.err"
    monkeypatch.setattr(detached, "_spawn_detached_worker", lambda *args: object())
    submit_detached(
        {
            "id": run_id,
            "executor": "local",
            "directory": str(tmp_path),
            "argv": [
                sys.executable,
                "-c",
                "import sys; print('stdout'); print('stderr', file=sys.stderr)",
            ],
            "output_file": str(output),
            "error_file": str(error),
        }
    )

    real_replace = detached_local.durable_replace
    replacements = []

    def fail_second_replace(source, destination):
        replacements.append((source, destination))
        if len(replacements) == 2:
            raise OSError("stderr publish failed")
        return real_replace(source, destination)

    monkeypatch.setattr(detached_local, "durable_replace", fail_second_replace)

    with pytest.raises(OSError, match="stderr publish failed"):
        run_detached_worker(tmp_path, run_id)

    state = status_detached(tmp_path, run_id)
    assert state["state"] == "failed"
    assert state["failure_phase"] == FailurePhase.FINALIZATION.value

    replay_output = detached.detached_stream_file(state, "output_file")
    replay_error = detached.detached_stream_file(state, "error_file")
    assert replay_output == str(output)
    assert replay_error == state["stream_error_file"]
    assert os.path.isfile(replay_output)
    assert os.path.isfile(replay_error)
    with open(replay_output) as file:
        assert file.read() == "stdout\n"
    with open(replay_error) as file:
        assert file.read() == "stderr\n"


def test_local_persistent_launch_failure_is_submission_not_started(tmp_path, monkeypatch):
    output = tmp_path / "calculation.out"
    monkeypatch.setattr(detached, "_spawn_detached_worker", _popen_failure)
    reader = _Reader()
    runner = LocalExecutor().create_persistent_run(
        reader,
        ExecutionSpec(
            argv=("unused",),
            input_file=None,
            output_file=str(output),
            directory=str(tmp_path),
        ),
        run_id="persistent-launch-failure",
    )

    with pytest.raises(SubmissionNotStarted, match="Detached worker launch failed"):
        runner.submit()

    assert runner.state == RunState.FAILED
    assert runner.failure_phase == FailurePhase.LAUNCH
    assert runner.backend_id is None
    assert reader.parse_calls == 0
    assert reader.complete_calls == 0
    assert not output.exists()


def test_local_persistent_metadata_failure_is_submission_failure(tmp_path, monkeypatch):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    run = manager.create_run(
        LocalExecutor(),
        _Reader(),
        ExecutionSpec(
            argv=("unused",),
            input_file=None,
            output_file=str(tmp_path / "calculation.out"),
            directory=str(tmp_path),
        ),
        task="SCF",
    )

    def fail_before_helper_start(request):
        raise OSError("state metadata could not be written")

    monkeypatch.setattr(
        local_persistent,
        "submit_detached",
        fail_before_helper_start,
    )

    with pytest.raises(SubmissionNotStarted, match="state metadata could not be written"):
        run.submit()

    record = manager.store.get(run.id)
    assert record.state == RunState.EXECUTION_FAILED
    assert record.failure_phase == FailurePhase.SUBMISSION
    assert record.execution.backend_id is None
    assert "state metadata could not be written" in record.message


def test_calculation_launch_failure_is_reported_by_submit(tmp_path):
    output = tmp_path / "calculation.out"
    reader = _Reader()
    runner = LocalExecutor().create_persistent_run(
        reader,
        ExecutionSpec(
            argv=(str(tmp_path / "missing-executable"),),
            input_file=None,
            output_file=str(output),
            directory=str(tmp_path),
        ),
        run_id="calculation-popen-failure",
    )

    with pytest.raises(SubmissionNotStarted, match="Detached worker launch failed"):
        runner.submit()

    assert runner.state == RunState.FAILED
    assert runner.failure_phase == FailurePhase.LAUNCH
    assert runner.backend_id is None
    assert reader.parse_calls == 0
    assert reader.complete_calls == 0
    assert not output.exists()
    state = status_detached(tmp_path, "calculation-popen-failure")
    assert state["failure_phase"] == "launch"


def test_launch_failure_releases_persistent_submission_namespace(tmp_path, monkeypatch):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = LocalExecutor()
    spec = ExecutionSpec(
        argv=(sys.executable, "-c", "print('ok')"),
        input_file=None,
        output_file=str(tmp_path / "calculation.out"),
        directory=str(tmp_path),
    )
    failed = manager.create_run(
        executor,
        _Reader(),
        spec,
        task="SCF",
        dataset="same-dataset",
    )

    monkeypatch.setattr(detached, "_spawn_detached_worker", _popen_failure)
    with pytest.raises(SubmissionNotStarted, match="Detached worker launch failed"):
        failed.submit()

    failed_record = manager.store.get(failed.id)
    assert failed_record.state == RunState.FAILED
    assert failed_record.failure_phase == FailurePhase.LAUNCH
    assert failed_record.execution.backend_id is None
    assert failed_record.files.captured_at is None
    assert failed_record.files.output_versions == ()
    assert "worker could not be launched" in failed_record.message

    monkeypatch.undo()
    replacement = manager.create_run(
        executor,
        _Reader(),
        spec,
        task="SCF",
        dataset="same-dataset",
    )
    replacement.submit()
    assert replacement.wait().return_code == 0


def test_calculation_launch_failure_never_enters_finalizing(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = LocalExecutor()
    reader = _Reader()
    run = manager.create_run(
        executor,
        reader,
        ExecutionSpec(
            argv=(str(tmp_path / "missing-executable"),),
            input_file=None,
            output_file=str(tmp_path / "calculation.out"),
            directory=str(tmp_path),
        ),
        task="SCF",
    )

    with pytest.raises(SubmissionNotStarted, match="Detached worker launch failed"):
        run.submit()

    record = manager.store.get(run.id)
    assert record.state == RunState.EXECUTION_FAILED
    assert record.failure_phase == FailurePhase.LAUNCH
    assert record.execution.backend_id is None
    assert record.files.captured_at is None
    assert reader.parse_calls == 0
    assert reader.complete_calls == 0


def test_psij_detached_submission_is_submitting_before_scheduler_acceptance(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(detached, "_spawn_detached_worker", lambda *args: object())

    state = submit_detached(
        {
            "id": "psij-submitting",
            "executor": "psij",
            "directory": str(tmp_path),
            "argv": ["unused"],
        }
    )

    assert state["state"] == RunState.SUBMITTING.value
    assert state["submitted_at"] is None
    assert state["backend_reference"] is None
    assert state["failure_phase"] is None

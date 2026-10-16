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


def test_running_state_persistence_failure_stops_spawned_process(tmp_path, monkeypatch):
    output = tmp_path / "calculation.out"
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

    state = status_detached(tmp_path, run_id)
    assert state["state"] == "failed"
    assert state["failure_phase"] == "launch"
    assert state["pid"] is None
    assert state["started_at"] is None


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
    assert runner.reference is None
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
    assert record.reference is None
    assert "state metadata could not be written" in record.message


def test_calculation_launch_failure_skips_parser_and_output_finalization(tmp_path):
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

    runner.submit()
    result = runner.wait()

    assert runner.state == RunState.FAILED
    assert runner.failure_phase == FailurePhase.LAUNCH
    assert result.return_code is None
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
    assert failed_record.reference is None
    assert failed_record.result_files.captured_at is None
    assert failed_record.result_files.output_versions == ()
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


def test_calculation_launch_failure_with_reference_never_enters_finalizing(tmp_path):
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

    run.submit()
    reference = run.reference
    assert reference is not None

    result = run.wait()
    assert result.return_code is None
    record = manager.store.get(run.id)
    assert record.state == RunState.EXECUTION_FAILED
    assert record.failure_phase == FailurePhase.LAUNCH
    assert record.result_files.captured_at is None
    assert reader.parse_calls == 0
    assert reader.complete_calls == 0

    # Simulate a fresh process refreshing the already-referenced failed run.
    recovered_manager = ExecutionManager(manager.store)
    recovered_manager.register_executor(executor, record.endpoint, reference)
    refreshed = recovered_manager.refresh((record,))[0]
    assert refreshed.state == RunState.EXECUTION_FAILED
    assert refreshed.failure_phase == FailurePhase.LAUNCH
    assert refreshed.result_files.captured_at is None


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

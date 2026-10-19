import asyncio
import os
import sys

import pytest

from ase2sprkkr.execution.base import (
    ExecutionSpec,
    FailurePhase,
    RunState,
    SubmissionNotStarted,
)
from ase2sprkkr.execution.endpoint import endpoint_for_executor
from ase2sprkkr.execution.local import LocalExecutor, LocalRunner
from ase2sprkkr.execution.manager import (
    ArtifactAvailabilityState,
    ExecutionManager,
    ProjectRunner,
    RunStore,
)
from ase2sprkkr.execution.postprocessing import PostprocessingRunner
from ase2sprkkr.execution.run_model import FileManifest, FileEntry
from ase2sprkkr.execution.transport import LocalTransport


def _local_spec(tmp_path, code="print('ok')"):
    return ExecutionSpec(
        argv=(sys.executable, "-c", code),
        input_file=None,
        output_file=(tmp_path / "out.txt").open("w+b"),
        error_file=(tmp_path / "err.txt").open("w+b"),
        directory=str(tmp_path),
    )


def test_output_observer_failure_is_reported_but_does_not_fail_process(tmp_path):
    runner = LocalRunner(_local_spec(tmp_path))
    observer_errors = []

    def broken_observer(*_):
        raise ValueError("observer failed")

    runner.add_output_listener(broken_observer)
    runner.add_output_error_listener(
        lambda callback, exc: observer_errors.append((callback, exc))
    )

    result = runner.run()

    assert runner.state == RunState.FINISHED
    assert result.return_code == 0
    assert len(observer_errors) == 1
    assert observer_errors[0][0] is broken_observer
    assert str(observer_errors[0][1]) == "observer failed"


def test_postprocessing_failure_does_not_change_raw_execution_state(tmp_path):
    raw = LocalRunner(_local_spec(tmp_path))

    def fail(_outcome):
        raise ValueError("parser failed")

    runner = PostprocessingRunner(raw, fail)
    runner.submit()

    with pytest.raises(ValueError, match="parser failed"):
        runner.wait()

    assert raw.state == RunState.FINISHED
    assert runner.state == RunState.FAILED
    assert runner.failure_phase == FailurePhase.POSTPROCESSING


def test_callback_failure_does_not_change_local_execution_state(tmp_path):
    runner = LocalRunner(
        _local_spec(tmp_path),
        callback=lambda result: (_ for _ in ()).throw(ValueError("callback failed")),
    )

    runner.submit()
    with pytest.raises(ValueError, match="callback failed"):
        runner.wait()

    assert runner.state == RunState.FINISHED


class _UnknownRunner:
    result = object()
    backend_id = None
    state = RunState.NEW
    exception = None
    failure_phase = None
    failure_message = None

    def submit(self):
        self.state = RunState.UNKNOWN

    def wait(self):
        return self.result

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass


def test_unknown_raw_state_is_preserved_and_not_postprocessed():
    raw = _UnknownRunner()
    processed = []
    runner = PostprocessingRunner(raw, lambda outcome: processed.append(outcome))

    runner.submit()
    with pytest.raises(RuntimeError, match="unresolved backend state"):
        runner.wait()

    assert processed == []
    assert runner.result is raw.result
    assert runner.state == RunState.UNKNOWN
    assert runner.failure_phase is None


class _NotStartedRunner:
    backend_id = None
    result = None
    state = RunState.NEW
    exception = None
    failure_phase = None
    failure_message = None

    def submit(self):
        raise SubmissionNotStarted("input staging failed")

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass


class _UncertainRunner(_NotStartedRunner):
    def submit(self):
        raise OSError("connection lost during submit")


def _project_wrapper(tmp_path, runner):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = LocalExecutor()
    endpoint = endpoint_for_executor(executor)
    record = manager.new_record(
        task="SCF",
        directory=tmp_path,
        endpoint=endpoint,
    )
    return manager, ProjectRunner(manager, executor, runner, record)


def test_definite_pre_submit_failure_becomes_submission_failure(tmp_path):
    manager, runner = _project_wrapper(tmp_path, _NotStartedRunner())

    with pytest.raises(SubmissionNotStarted, match="input staging failed"):
        runner.submit()

    record = manager.store.get(runner.id)
    assert record.state == RunState.FAILED
    assert record.failure_phase == FailurePhase.SUBMISSION
    assert record.execution.backend_id is None
    assert record.message == "input staging failed"


def test_uncertain_submit_without_reference_stays_submitting(tmp_path):
    manager, runner = _project_wrapper(tmp_path, _UncertainRunner())

    with pytest.raises(OSError, match="connection lost during submit"):
        runner.submit()

    record = manager.store.get(runner.id)
    assert record.state == RunState.SUBMITTING
    assert record.execution.backend_id is None


def test_artifact_availability_failure_does_not_change_finished_state(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    record = manager.new_record(
        task="SCF",
        directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(tmp_path / "out"))),
    )
    record = manager.store.save(record)
    record = manager.transition(record, state=RunState.FINISHED)

    record = manager.update_artifact_availability(
        record, str(tmp_path / "out"), ArtifactAvailabilityState.FAILED, "ssh disconnected"
    )

    assert record.state == RunState.FINISHED
    assert (
        record.files.availability_state
        == ArtifactAvailabilityState.FAILED
    )
    assert record.files.availability_message == "ssh disconnected"


def test_local_transport_failed_copy_preserves_destination(tmp_path):
    destination = tmp_path / "output"
    destination.write_bytes(b"old")

    with pytest.raises(FileNotFoundError):
        asyncio.run(LocalTransport().put_file(tmp_path / "missing", destination))

    assert destination.read_bytes() == b"old"



import asyncio
import os
import sys

import pytest

from ase2sprkkr.execution.base import (
    file_path,
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
    assert runner.wait() is raw.result

    # UNKNOWN is an active, recoverable outcome, not a fabricated failure.
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




def test_file_path_resolves_named_files_and_tolerates_unnamed(tmp_path):
    output = tmp_path / "output.dat"
    output.write_text("a")
    with output.open("rb") as stream:
        assert file_path(stream, tmp_path) == str(output)
    assert file_path("relative.dat", tmp_path) == str(tmp_path / "relative.dat")
    assert file_path(None, tmp_path) is None
    assert file_path(1, tmp_path) is None
    assert file_path(object(), tmp_path) is None


def test_file_path_preserves_pathlike_directory(tmp_path):
    from pathlib import Path

    absolute = tmp_path / "different" / "result.dat"
    assert file_path(absolute, tmp_path / "work") == str(absolute)
    assert file_path(Path("subdir/result.dat"), tmp_path) == str(
        tmp_path / "subdir" / "result.dat"
    )


def test_file_path_flush_is_explicit(tmp_path):
    class FileLike:
        name = "output.dat"
        def __init__(self):
            self.flush_calls = 0
        def flush(self):
            self.flush_calls += 1

    value = FileLike()
    assert file_path(value, tmp_path) == str(tmp_path / value.name)
    assert value.flush_calls == 0
    assert file_path(value, tmp_path, flush=True) == str(tmp_path / value.name)
    assert value.flush_calls == 1


def test_file_path_strict_contract_preserves_optional_none(tmp_path):
    from ase2sprkkr.execution.base import file_path

    assert file_path(None, tmp_path, required=True) is None
    assert file_path("case.out", tmp_path, required=True) == str(tmp_path / "case.out")
    with pytest.raises(ValueError, match="named file"):
        file_path(1, tmp_path, required=True)
    with pytest.raises(ValueError, match="named file"):
        file_path(object(), tmp_path, required=True)


def test_strict_runner_path_errors_remain_backend_specific(tmp_path):
    from ase2sprkkr.execution.psij import _path as psij_path
    from ase2sprkkr.execution.remote import RemoteRunner

    spec = ExecutionSpec(
        argv=(), input_file=None, output_file=None, directory=str(tmp_path),
    )
    runner = RemoteRunner(object(), spec)
    with pytest.raises(ValueError, match="Remote execution requires named files"):
        runner._local_path(object())
    with pytest.raises(ValueError, match="Remote execution requires named files"):
        runner._local_path(1)
    with pytest.raises(ValueError, match="PSI/J execution requires named"):
        psij_path(object(), tmp_path)
    assert runner._local_path(None) is None
    assert psij_path(None, tmp_path) is None



def test_common_command_request_preserves_mpi_and_launch_contract(tmp_path):
    from ase2sprkkr.execution import ExecutionSpec, LaunchSpec

    spec = ExecutionSpec(
        argv=("serial", "input.inp"),
        mpi_argv=("parallel", "input.inp"),
        mpi_uses_stdin=False,
        input_file=None, output_file=None,
        directory=str(tmp_path),
        launch=LaunchSpec.from_mpi(2),
    )
    request = spec.command_request()
    assert request["argv"] == ["serial", "input.inp"]
    assert request["mpi_argv"] == ["parallel", "input.inp"]
    assert request["mpi_uses_stdin"] is False

    transformed = spec.command_request(
        lambda command: [arg.upper() for arg in command]
    )
    assert transformed["argv"] == ["SERIAL", "INPUT.INP"]
    assert transformed["mpi_argv"] == ["PARALLEL", "INPUT.INP"]
    assert transformed["launch"] == request["launch"]
    assert transformed["mpi_uses_stdin"] is False
    assert spec.argv == ("serial", "input.inp")

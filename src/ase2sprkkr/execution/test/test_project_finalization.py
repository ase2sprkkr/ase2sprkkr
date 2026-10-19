from threading import Event
from types import SimpleNamespace

import pytest

import ase2sprkkr.outputs.project_results as project_output
from ase2sprkkr.execution import (
    ExecutionOutcome,
    Executor,
    FailurePhase,
    RunState,
    RunStore,
)
from ase2sprkkr.execution.manager import ExecutionManager
from ase2sprkkr.outputs.project_results import (
    UnrecoverableProjectResultError,
    attach_project_run,
)


class _RawRunner:
    def __init__(self):
        self.state = RunState.NEW
        self.backend_id = None
        self.result = ExecutionOutcome(return_code=0)
        self.exception = None
        self.failure_phase = None
        self.failure_message = None

    @property
    def return_code(self):
        return self.result.return_code

    def submit(self):
        self.backend_id = "job-1"
        self.state = RunState.FINISHED

    def wait(self):
        return self.result

    def cancel(self):
        self.state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass


class _Executor(Executor):
    name = "test"

    def create_project_run(self, spec, *, run_id, callback=None):
        self.runner = _RawRunner()
        return self.runner


def test_finished_run_is_persisted_before_callback_without_wait(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    finished = Event()
    run = manager.create_run(
        _Executor(),
        SimpleNamespace(directory=str(tmp_path)),
        process_result=lambda outcome: outcome,
        task="DOS",
        callback=lambda result: finished.set(),
    )
    run.submit()
    assert finished.wait(timeout=5), "completion callback was not called"
    record = store.get(run.id)
    assert record.state == RunState.FINISHED
    assert record.finished_at is not None


def test_result_processor_failure_is_postprocessing_failure(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _Executor()

    def fail(_outcome):
        raise ValueError("result parsing failed")

    run = manager.create_run(
        executor,
        SimpleNamespace(directory=str(tmp_path)),
        process_result=fail,
        task="DOS",
    )
    run.submit()

    with pytest.raises(ValueError, match="result parsing failed"):
        run.wait()

    record = store.get(run.id)
    assert record.state == RunState.FAILED
    assert record.failure_phase == FailurePhase.POSTPROCESSING
    assert record.message == "result parsing failed"
    assert record.finished_at is not None


def test_unrecoverable_processor_reconstruction_fails_postprocessing(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="NO_SUCH_TASK", directory=tmp_path)
    record = manager.transition(record, state=RunState.POSTPROCESSING)

    with pytest.raises(UnrecoverableProjectResultError):
        attach_project_run(manager, record)

    failed = store.get(record.id)
    assert failed.state == RunState.FAILED
    assert failed.failure_phase == FailurePhase.POSTPROCESSING
    assert failed.finished_at is not None
    assert "unknown persisted task" in failed.message


def test_reader_constructor_programming_error_is_not_reclassified(
    tmp_path, monkeypatch
):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="DOS", directory=tmp_path)
    record = manager.transition(record, state=RunState.POSTPROCESSING)
    error = RuntimeError("reader constructor failed")

    class BrokenReader:
        def __init__(self, *_args, **_kwargs):
            raise error

    monkeypatch.setattr(
        project_output.reader_registry,
        "class_for_task",
        lambda _task_definition: BrokenReader,
    )

    with pytest.raises(RuntimeError, match="reader constructor failed") as exc_info:
        attach_project_run(manager, record)

    assert exc_info.value is error
    current = store.get(record.id)
    assert current.state == RunState.POSTPROCESSING
    assert current.message is None
    assert current.finished_at is None


def test_transient_attach_failure_keeps_postprocessing_active(tmp_path, monkeypatch):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="DOS", directory=tmp_path)
    record = manager.transition(
        record,
        state=RunState.POSTPROCESSING,
        execution=record.execution.submitted("job-1"),
    )

    monkeypatch.setattr(
        project_output,
        "_result_processor_for_record",
        lambda _manager, _record: lambda outcome: outcome,
    )

    def fail_attach(*_args, **_kwargs):
        raise OSError("temporary ssh failure")

    monkeypatch.setattr(manager, "attach", fail_attach)

    with pytest.raises(OSError, match="temporary ssh failure"):
        attach_project_run(manager, record)

    current = store.get(record.id)
    assert current.state == RunState.POSTPROCESSING
    assert current.finished_at is None


def test_attaching_failed_download_reopens_transfer_without_resubmitting(
    tmp_path, monkeypatch
):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    from ase2sprkkr.execution import FileEntry, FileManifest
    from ase2sprkkr.outputs.project_results import attach_project_run

    output = tmp_path / "run.out"
    output.write_bytes(b"partial output")
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(
            local_path=str(output), availability_state="failed",
            availability_message="connection dropped",
        )),
    )
    record = manager.transition(
        record, state=RunState.DOWNLOADING,
        execution=record.execution.submitted("existing-backend-id"),
    )
    manager.complete(
        record, state=RunState.FAILED, return_code=0,
        failure_phase=FailurePhase.DOWNLOAD,
        message="connection dropped", artifacts_available=False,
    )
    import ase2sprkkr.outputs.project_results as project_results

    monkeypatch.setattr(
        project_results, "_result_processor_for_record",
        lambda *_: lambda result: result,
    )
    attached = []
    def fake_attach(record, **kwargs):
        attached.append(record)
        return record

    monkeypatch.setattr(manager, "attach", fake_attach)
    resumed = attach_project_run(manager, record.id)
    assert resumed.state == RunState.DOWNLOADING
    assert resumed.execution.backend_id == "existing-backend-id"
    assert manager.store.get(record.id).finished_at is None
    assert output.read_bytes() == b"partial output"
    assert attached == [resumed]


def test_failed_download_reader_construction_does_not_reopen_run(tmp_path, monkeypatch):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(
            local_path=str(tmp_path / "out"),
            availability_state="failed",
            availability_message="original download failure",
        )),
    )
    record = manager.transition(
        record, state=RunState.DOWNLOADING,
        execution=record.execution.submitted("existing-backend"),
    )
    failed = manager.complete(
        record, state=RunState.FAILED, return_code=0,
        failure_phase=FailurePhase.DOWNLOAD, artifacts_available=False,
    )

    def broken_reader(*args):
        raise RuntimeError("Reader constructor is broken")

    monkeypatch.setattr(
        project_output, "_result_processor_for_record", broken_reader
    )
    with pytest.raises(RuntimeError, match="Reader constructor is broken"):
        attach_project_run(manager, failed.id)

    current = manager.store.get(failed.id)
    assert current.state == RunState.FAILED
    assert current.failure_phase == FailurePhase.DOWNLOAD
    assert current.files == failed.files


def test_wait_does_not_finalize_a_second_time(tmp_path, monkeypatch):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    completed = []
    original_complete = manager.complete

    def complete(*args, **kwargs):
        completed.append(kwargs["state"])
        return original_complete(*args, **kwargs)

    monkeypatch.setattr(manager, "complete", complete)
    run = manager.create_run(
        _Executor(),
        SimpleNamespace(directory=str(tmp_path)),
        process_result=lambda outcome: outcome,
        task="DOS",
    )
    run.submit()
    run.wait()
    run.wait()

    assert completed == [RunState.FINISHED]
    assert store.get(run.id).state == RunState.FINISHED


def test_state_property_does_not_query_backend_or_persist(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    run = manager.create_run(
        _Executor(),
        SimpleNamespace(directory=str(tmp_path)),
        process_result=lambda outcome: outcome,
        task="DOS",
    )

    class StateQueryForbidden:
        @property
        def state(self):
            raise AssertionError("reading ProjectRunner.state queried the backend")

    run._runner = StateQueryForbidden()
    assert run.state == RunState.NEW
    assert store.get(run.id).state == RunState.NEW

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

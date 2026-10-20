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


def test_postprocessing_retry_uses_saved_files_without_backend(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    output = tmp_path / "run.out"
    output.write_text("original result")
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(output))),
    )
    record = manager.transition(
        record, state=RunState.POSTPROCESSING,
        execution=record.execution.submitted("existing-backend"),
    )
    failed = manager.complete(
        record, state=RunState.FAILED, return_code=0,
        failure_phase=FailurePhase.POSTPROCESSING,
        message="transient parse error", artifacts_available=True,
    )
    assert failed.files.captured_at is not None
    calls = []

    def parser(outcome):
        calls.append(outcome)
        return SimpleNamespace(value=output.read_text(), run_record=None)

    result = manager.retry_postprocessing(failed.id, parser)
    assert result.value == "original result"
    assert len(calls) == 1
    latest = manager.store.get(failed.id)
    assert latest.state == RunState.FINISHED
    assert latest.execution.backend_id == "existing-backend"
    assert latest.files.input_versions == failed.files.input_versions
    assert latest.files.output_versions == failed.files.output_versions
    with pytest.raises(ValueError, match="Only postprocessing failures"):
        manager.retry_postprocessing(failed.id, parser)


def test_postprocessing_retry_failure_is_retryable_and_does_not_resubmit(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    output = tmp_path / "run.out"
    output.write_text("original")
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(output))),
    )
    record = manager.transition(
        record, state=RunState.POSTPROCESSING,
        execution=record.execution.submitted("existing-backend"),
    )
    failed = manager.complete(
        record, state=RunState.FAILED, return_code=0,
        failure_phase=FailurePhase.POSTPROCESSING, artifacts_available=True,
    )
    def broken(_outcome):
        raise ValueError("parse failed again")
    with pytest.raises(ValueError, match="parse failed again"):
        manager.retry_postprocessing(failed.id, broken)
    current = manager.store.get(failed.id)
    assert current.state == RunState.FAILED
    assert current.failure_phase == FailurePhase.POSTPROCESSING
    assert current.execution.backend_id == failed.execution.backend_id
    manager.retry_postprocessing(failed.id, lambda outcome: outcome)
    assert manager.store.get(failed.id).state == RunState.FINISHED


def test_postprocessing_retry_refuses_modified_or_uncaptured_artifacts(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    output = tmp_path / "run.out"
    output.write_text("original")
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(output))),
    )
    record = manager.transition(
        record, state=RunState.POSTPROCESSING,
        execution=record.execution.submitted("existing-backend"),
    )
    failed = manager.complete(
        record, state=RunState.FAILED, return_code=0,
        failure_phase=FailurePhase.POSTPROCESSING, artifacts_available=True,
    )
    output.write_text("overwritten")
    with pytest.raises(ValueError, match="unchanged historical artifacts"):
        manager.retry_postprocessing(failed.id, lambda outcome: outcome)
    assert manager.store.get(failed.id).state == RunState.FAILED


def test_postprocessing_retry_reader_constructor_error_does_not_reopen(tmp_path, monkeypatch):
    from ase2sprkkr.execution import FileEntry, FileManifest
    from ase2sprkkr.outputs import retry_project_postprocessing
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    output = tmp_path / "run.out"
    output.write_text("original")
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(output))),
    )
    record = manager.transition(
        record, state=RunState.POSTPROCESSING,
        execution=record.execution.submitted("existing-backend"),
    )
    failed = manager.complete(
        record, state=RunState.FAILED, return_code=0,
        failure_phase=FailurePhase.POSTPROCESSING, artifacts_available=True,
    )
    monkeypatch.setattr(
        project_output, "_result_processor_for_record",
        lambda *_args: (_ for _ in ()).throw(RuntimeError("broken reader")),
    )
    with pytest.raises(RuntimeError, match="broken reader"):
        retry_project_postprocessing(manager, failed.id)
    assert manager.store.get(failed.id).state == RunState.FAILED


def test_project_retry_attaches_run_provenance(tmp_path, monkeypatch):
    from ase2sprkkr.execution import FileEntry, FileManifest
    from ase2sprkkr.outputs import retry_project_postprocessing

    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    output = tmp_path / "run.out"
    output.write_text("saved data")
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(output))),
    )
    record = manager.transition(
        record, state=RunState.POSTPROCESSING,
        execution=record.execution.submitted("job-existing"),
    )
    failed = manager.complete(
        record, state=RunState.FAILED, return_code=0,
        failure_phase=FailurePhase.POSTPROCESSING, artifacts_available=True,
    )
    monkeypatch.setattr(
        project_output, "_result_processor_for_record",
        lambda _manager, _record: lambda outcome: SimpleNamespace(value="parsed"),
    )
    result = retry_project_postprocessing(manager, failed.id)
    assert result.value == "parsed"
    assert result.run_record.id == failed.id
    assert result.run_record.state == RunState.FINISHED


def test_postprocessing_retry_only_one_caller_can_claim_run(tmp_path):
    from threading import Thread
    from ase2sprkkr.execution import FileEntry, FileManifest

    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    output = tmp_path / "run.out"
    output.write_text("saved data")
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(output))),
    )
    record = manager.transition(
        record, state=RunState.POSTPROCESSING,
        execution=record.execution.submitted("job-existing"),
    )
    failed = manager.complete(
        record, state=RunState.FAILED, return_code=0,
        failure_phase=FailurePhase.POSTPROCESSING, artifacts_available=True,
    )

    processing = Event()
    release = Event()
    errors = []

    def parse(outcome):
        processing.set()
        if not release.wait(5):
            raise TimeoutError("test retry parser was not released")
        return outcome

    def first_attempt():
        try:
            manager.retry_postprocessing(failed.id, parse)
        except BaseException as exc:
            errors.append(exc)

    worker = Thread(target=first_attempt)
    worker.start()
    try:
        assert processing.wait(5), "first retry did not acquire the run"
        with pytest.raises(ValueError, match="Only postprocessing failures"):
            manager.retry_postprocessing(failed.id, lambda outcome: outcome)
        assert manager.store.get(failed.id).state == RunState.POSTPROCESSING
    finally:
        release.set()
        worker.join(5)
    assert not worker.is_alive()
    assert not errors
    assert manager.store.get(failed.id).state == RunState.FINISHED


def test_retry_accepts_plain_processing_result(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    output = tmp_path / "run.out"
    output.write_text("contents")
    record = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(output))),
    )
    record = manager.transition(
        record, state=RunState.POSTPROCESSING,
        execution=record.execution.submitted("job-existing"),
    )
    failed = manager.complete(
        record, state=RunState.FAILED, return_code=0,
        failure_phase=FailurePhase.POSTPROCESSING, artifacts_available=True,
    )
    assert manager.retry_postprocessing(failed.id, lambda _: 42) == 42
    assert manager.store.get(failed.id).state == RunState.FINISHED

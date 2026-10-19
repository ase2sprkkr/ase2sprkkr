from types import SimpleNamespace
from threading import Event

import pytest

from ase2sprkkr.execution import (
    DuplicateSubmissionError,
    Executor,
    RunState,
    RunStore,
)

from ase2sprkkr.execution.manager import ExecutionManager


class _Runner:
    def __init__(self, callback=None):
        self.callback = callback
        self.state = RunState.NEW
        self.backend_id = None
        self._done = Event()
        self.result = SimpleNamespace(return_code=None)
        self.exception = None

    def submit(self):
        self.backend_id = "job"
        self.state = RunState.QUEUED

    def wait(self):
        self._done.wait(2)
        return self.result

    def run(self):
        self.submit()
        return self.wait()

    def cancel(self):
        self.state = RunState.CANCELLED
        self._done.set()

    def add_output_listener(self, callback):
        pass

    def finish(self):
        self.result.return_code = 0
        self.state = RunState.FINISHED
        self._done.set()
        if self.callback:
            self.callback(self.result)


class _Executor(Executor):
    name = "test"

    def __init__(self):
        self.runners = []

    def create_run(self, spec, callback=None):
        runner = _Runner(callback)
        self.runners.append(runner)
        return runner


def _run(manager, executor, directory, *, task="DOS", dataset="case"):
    return manager.create_run(
        executor,
        spec=SimpleNamespace(directory=str(directory)),
        process_result=lambda outcome: outcome,
        task=task,
        dataset=dataset,
    )


def test_duplicate_task_and_dataset_is_rejected_across_managers(tmp_path):
    directory = tmp_path / "work"
    store_path = tmp_path / "runs"
    first_manager = ExecutionManager(RunStore(store_path))
    second_manager = ExecutionManager(RunStore(store_path))
    first = _run(first_manager, _Executor(), directory)
    second = _run(second_manager, _Executor(), directory)

    first.submit()

    with pytest.raises(DuplicateSubmissionError) as exc:
        second.submit()

    assert exc.value.conflicts[0].id == first.id
    assert second.record.state == RunState.NEW


def test_different_dataset_can_share_non_scf_directory(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    first = _run(manager, executor, tmp_path, dataset="one")
    second = _run(manager, executor, tmp_path, dataset="two")

    first.submit()
    second.submit()

    assert first.record.state == RunState.QUEUED
    assert second.record.state == RunState.QUEUED


def test_scf_conflicts_with_other_tasks_in_same_namespace(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    scf = _run(manager, executor, tmp_path, task="SCF", dataset="scf")
    dos = _run(manager, executor, tmp_path, task="DOS", dataset="dos")

    scf.submit()

    with pytest.raises(DuplicateSubmissionError):
        dos.submit()


def test_final_run_releases_submission_namespace(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    first = _run(manager, executor, tmp_path)
    first.submit()
    executor.runners[0].finish()
    # The raw runner is finished, but the project's asynchronous postprocessor
    # must complete its durable callback before the namespace can be released.
    first.wait()
    assert manager.store.get(first.id).state == RunState.FINISHED

    second = _run(manager, executor, tmp_path)
    second.submit()

    assert second.record.state == RunState.QUEUED


def test_submitting_record_contains_submission_key_before_backend_submit(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _Executor()
    run = _run(manager, executor, tmp_path)
    backend = executor.runners[0]
    seen = []

    original_submit = backend.submit

    def submit():
        record = store.get(run.id)
        seen.append(record)
        original_submit()

    backend.submit = submit
    run.submit()

    assert seen[0].state == RunState.SUBMITTING
    assert seen[0].task == "DOS"
    assert seen[0].dataset == "case"
    assert seen[0].execution.directory.local == str(tmp_path.resolve())

from types import SimpleNamespace
from threading import Event

import pytest

from ase2sprkkr.execution import (
    DuplicateSubmissionError,
    Executor,
    FileEntry,
    FileManifest,
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


@pytest.mark.parametrize("overlap", (True, False), ids=("same-file", "separate-files"))
def test_declared_exact_outputs_conflict_only_if_they_share_a_target(
    tmp_path, overlap,
):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    first_path = str(tmp_path / "one.out")
    second_path = first_path if overlap else str(tmp_path / "two.out")
    first = manager.new_record(
        task="DOS", directory=tmp_path, dataset="one",
        files=FileManifest(stdout=FileEntry(local_path=first_path)),
    )
    second = manager.new_record(
        task="BSF" if overlap else "DOS", directory=tmp_path,
        dataset="two" if overlap else "one",
        files=(
            FileManifest(outputs=(FileEntry(local_path=second_path),))
            if overlap else
            FileManifest(stdout=FileEntry(local_path=second_path))
        ),
    )

    manager.reserve(first)
    if overlap:
        with pytest.raises(DuplicateSubmissionError):
            manager.reserve(second)
    else:
        assert manager.reserve(second).state == RunState.SUBMITTING


@pytest.mark.parametrize(
    "first_access, same_directory, conflict",
    [
        pytest.param("output", True, True, id="writer-reader-same-directory"),
        pytest.param("output", False, True, id="writer-reader-different-directories"),
        pytest.param("input", False, False, id="two-readers-different-directories"),
    ],
)
def test_declared_read_write_conflicts_respect_access_and_absolute_paths(
    tmp_path, first_access, same_directory, conflict,
):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    shared = tmp_path / ("shared.pot" if same_directory or first_access == "input" else "shared.out")
    if same_directory:
        shared.write_text("density")
    first_directory = tmp_path if same_directory else tmp_path / "one"
    second_directory = tmp_path if same_directory else tmp_path / "two"

    entry = FileEntry(local_path=str(shared))
    first_files = (
        FileManifest(outputs=(entry,))
        if first_access == "output" else FileManifest(inputs=(entry,))
    )
    second_files = FileManifest(
        inputs=(FileEntry(local_path=str(shared)),),
        stdout=FileEntry(local_path=str(tmp_path / "other.out")) if same_directory else None,
    )
    first = manager.new_record(
        task="DOS", directory=first_directory, dataset="one", files=first_files,
    )
    second = manager.new_record(
        task="BSF", directory=second_directory, dataset="two", files=second_files,
    )

    manager.reserve(first)
    if conflict:
        with pytest.raises(DuplicateSubmissionError):
            manager.reserve(second)
    else:
        assert manager.reserve(second).state == RunState.SUBMITTING


@pytest.mark.parametrize(
    "first_pattern, second_pattern, second_is_stdout, conflict",
    [
        pytest.param("result-*.dat", "result-1.dat", True, True,
                     id="wildcard-versus-exact-stdout"),
        pytest.param("spectrum-*.dat", "spectrum-?.dat", False, True,
                     id="overlapping-wildcard-patterns"),
        pytest.param("bsf-*.dat", "dos-*.dat", False, False,
                     id="disjoint-wildcard-patterns"),
    ],
)
def test_wildcard_output_collisions_follow_pattern_overlap(
    tmp_path, first_pattern, second_pattern, second_is_stdout, conflict,
):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    first = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(outputs=(
            FileEntry(local_path=str(tmp_path / first_pattern)),)),
    )
    second_entry = FileEntry(local_path=str(tmp_path / second_pattern))
    second = manager.new_record(
        task="BSF", directory=tmp_path,
        files=(
            FileManifest(stdout=second_entry)
            if second_is_stdout else FileManifest(outputs=(second_entry,))
        ),
    )

    manager.reserve(first)
    if conflict:
        with pytest.raises(DuplicateSubmissionError):
            manager.reserve(second)
    else:
        assert manager.reserve(second).state == RunState.SUBMITTING



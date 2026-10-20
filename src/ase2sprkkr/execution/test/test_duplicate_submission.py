from types import SimpleNamespace
from threading import Event

import pytest

from ase2sprkkr.execution import (
    DuplicateSubmissionError,
    ExecutionSpec,
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
        spec=ExecutionSpec(
            argv=(), input_file=None, output_file=None, directory=str(directory),
        ),
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


def test_different_tasks_conflict_on_shared_declared_output(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    target = str(tmp_path / "same.out")
    first = manager.new_record(
        task="DOS", directory=tmp_path, dataset="one",
        files=FileManifest(stdout=FileEntry(local_path=target)),
    )
    second = manager.new_record(
        task="BSF", directory=tmp_path, dataset="two",
        files=FileManifest(outputs=(FileEntry(local_path=target),)),
    )
    manager.reserve(first)
    with pytest.raises(DuplicateSubmissionError):
        manager.reserve(second)


def test_declared_input_read_conflicts_with_other_task_output(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    input_path = tmp_path / "shared.pot"
    input_path.write_text("density")
    writer = manager.new_record(
        task="DOS", directory=tmp_path, dataset="one",
        files=FileManifest(outputs=(FileEntry(local_path=str(input_path)),)),
    )
    reader = manager.new_record(
        task="BSF", directory=tmp_path, dataset="two",
        files=FileManifest(
            inputs=(FileEntry(local_path=str(input_path)),),
            stdout=FileEntry(local_path=str(tmp_path / "other.out")),
        ),
    )
    manager.reserve(writer)
    with pytest.raises(DuplicateSubmissionError):
        manager.reserve(reader)


def test_nonoverlapping_declared_outputs_can_run_in_same_directory(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    first = manager.new_record(
        task="DOS", directory=tmp_path, dataset="one",
        files=FileManifest(stdout=FileEntry(local_path=str(tmp_path / "one.out"))),
    )
    second = manager.new_record(
        task="DOS", directory=tmp_path, dataset="one",
        files=FileManifest(stdout=FileEntry(local_path=str(tmp_path / "two.out"))),
    )
    manager.reserve(first)
    assert manager.reserve(second).state == RunState.SUBMITTING


def test_wildcard_declared_output_conflicts_with_concrete_output(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    first = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(outputs=(FileEntry(local_path=str(tmp_path / "result-*.dat")),)),
    )
    second = manager.new_record(
        task="BSF", directory=tmp_path,
        files=FileManifest(stdout=FileEntry(local_path=str(tmp_path / "result-1.dat"))),
    )
    manager.reserve(first)
    with pytest.raises(DuplicateSubmissionError):
        manager.reserve(second)


def test_shared_absolute_output_conflicts_across_working_directories(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    shared = str(tmp_path / "shared.out")
    first = manager.new_record(
        task="DOS", directory=tmp_path / "one",
        files=FileManifest(outputs=(FileEntry(local_path=shared),)),
    )
    second = manager.new_record(
        task="BSF", directory=tmp_path / "two",
        files=FileManifest(inputs=(FileEntry(local_path=shared),)),
    )
    manager.reserve(first)
    with pytest.raises(DuplicateSubmissionError):
        manager.reserve(second)


def test_wildcard_to_wildcard_conflicts_even_if_patterns_differ(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    first = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(outputs=(
            FileEntry(local_path=str(tmp_path / "spectrum-*.dat")),)),
    )
    second = manager.new_record(
        task="BSF", directory=tmp_path,
        files=FileManifest(outputs=(
            FileEntry(local_path=str(tmp_path / "spectrum-?.dat")),)),
    )
    manager.reserve(first)
    with pytest.raises(DuplicateSubmissionError):
        manager.reserve(second)


def test_disjoint_wildcard_prefixes_do_not_conflict(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    first = manager.new_record(
        task="DOS", directory=tmp_path,
        files=FileManifest(outputs=(
            FileEntry(local_path=str(tmp_path / "bsf-*.dat")),)),
    )
    second = manager.new_record(
        task="BSF", directory=tmp_path,
        files=FileManifest(outputs=(
            FileEntry(local_path=str(tmp_path / "dos-*.dat")),)),
    )
    manager.reserve(first)
    assert manager.reserve(second).state == RunState.SUBMITTING


def test_shared_input_is_not_a_conflict(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    shared = str(tmp_path / "shared.pot")
    first = manager.new_record(
        task="DOS", directory=tmp_path / "one",
        files=FileManifest(inputs=(FileEntry(local_path=shared),)),
    )
    second = manager.new_record(
        task="BSF", directory=tmp_path / "two",
        files=FileManifest(inputs=(FileEntry(local_path=shared),)),
    )
    manager.reserve(first)
    assert manager.reserve(second).state == RunState.SUBMITTING


@pytest.mark.parametrize("cache_exists", [False, True])
def test_shared_cache_blocks_concurrent_writers_with_distinct_outputs(
    tmp_path, cache_exists,
):
    from ase2sprkkr.execution import CacheEntry, FileEntry, FileManifest

    cache = tmp_path / "shared.tau.dir"
    if cache_exists:
        cache.write_text("existing cache")
    manager = ExecutionManager(RunStore(tmp_path / "runs"))

    def record(task):
        return manager.new_record(
            task=task, dataset=task, directory=tmp_path,
            files=FileManifest(
                outputs=(FileEntry(local_path=str(tmp_path / f"{task}.out")),),
                caches=(CacheEntry(FileEntry(local_path=str(cache))),),
            ),
        )

    first = manager.reserve(record("DOS"))
    assert first.state == RunState.SUBMITTING
    with pytest.raises(DuplicateSubmissionError):
        manager.reserve(record("BSF"))


def test_missing_cache_producer_still_blocks_later_reader(tmp_path):
    from ase2sprkkr.execution import CacheEntry, FileEntry, FileManifest

    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    cache = tmp_path / "shared.tau.dir"

    def record(task):
        return manager.new_record(
            task=task, directory=tmp_path,
            files=FileManifest(
                outputs=(FileEntry(local_path=str(tmp_path / f"{task}.out")),),
                caches=(CacheEntry(FileEntry(local_path=str(cache))),),
            ),
        )

    manager.reserve(record("DOS"))
    cache.write_text("partial data")
    with pytest.raises(DuplicateSubmissionError):
        manager.reserve(record("BSF"))


def test_remote_file_collision_uses_target_identity(tmp_path):
    from ase2sprkkr.execution import ExecutorEndpoint, FileEntry, FileManifest

    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    endpoint = ExecutorEndpoint(executor="ssh", host="remote-host")
    first = manager.new_record(
        task="DOS", directory=tmp_path / "one", endpoint=endpoint,
        remote_directory="/remote/job",
        files=FileManifest(outputs=(
            FileEntry(local_path=str(tmp_path / "one" / "a"),
                      remote_path="/remote/job/result"),)),
    )
    second = manager.new_record(
        task="BSF", directory=tmp_path / "two", endpoint=endpoint,
        remote_directory="/remote/job",
        files=FileManifest(inputs=(
            FileEntry(local_path=str(tmp_path / "two" / "b"),
                      remote_path="/remote/job/result"),)),
    )
    manager.reserve(first)
    with pytest.raises(DuplicateSubmissionError):
        manager.reserve(second)



def test_stdin_read_conflicts_with_output_writer(tmp_path):
    from ase2sprkkr.execution import FileEntry, FileManifest

    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    shared = tmp_path / "shared.in"
    reader = manager.new_record(
        task="DOS", directory=tmp_path, dataset="reader",
        files=FileManifest.from_spec(ExecutionSpec(
            argv=("cat",), input_file=None, output_file="reader.out",
            directory=str(tmp_path), stdin=str(shared),
        )),
    )
    writer = manager.new_record(
        task="BSF", directory=tmp_path, dataset="writer",
        files=FileManifest(outputs=(FileEntry(local_path=str(shared)),)),
    )
    assert reader.files.input_file is None
    manager.reserve(reader)
    with pytest.raises(DuplicateSubmissionError):
        manager.reserve(writer)


def test_two_stdin_readers_can_run_concurrently(tmp_path):
    from ase2sprkkr.execution import FileManifest

    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    for task in ("DOS", "BSF"):
        record = manager.new_record(
            task=task, directory=tmp_path, dataset=task,
            files=FileManifest.from_spec(ExecutionSpec(
                argv=("cat",), input_file=None, output_file=f"{task}.out",
                directory=str(tmp_path), stdin="shared.in",
            )),
        )
        assert manager.reserve(record).state == RunState.SUBMITTING

from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    DuplicateSubmissionError,
    ExecutionManager,
    Executor,
    RunReference,
    RunState,
    RunStore,
)


class _Runner:
    def __init__(self, callback=None):
        self.callback = callback
        self.state = RunState.NEW
        self.reference = None
        self.result = SimpleNamespace(return_code=None)
        self.exception = None

    def submit(self):
        self.reference = RunReference(executor="test", native_id="job")
        self.state = RunState.QUEUED

    def wait(self):
        return self.result

    def run(self):
        self.submit()
        return self.wait()

    def cancel(self):
        self.state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass

    def finish(self):
        self.result.return_code = 0
        self.state = RunState.FINISHED
        if self.callback:
            self.callback(self.result)


class _Executor(Executor):
    name = "test"

    def __init__(self):
        self.runners = []

    def create_run(self, reader, spec, callback=None):
        runner = _Runner(callback)
        self.runners.append(runner)
        return runner


def _spec(directory, *, inputs=(), outputs=(), caches=()):
    return SimpleNamespace(
        directory=str(directory),
        input_file=None,
        output_file=None,
        error_file=None,
        input_files=tuple(inputs),
        output_files=tuple(outputs),
        cache_files=tuple(caches),
    )


def _run(
    manager,
    executor,
    directory,
    *,
    inputs=(),
    outputs=(),
    caches=(),
    task="DOS",
    dataset="case",
):
    return manager.create_run(
        executor,
        reader=object(),
        spec=_spec(
            directory,
            inputs=inputs,
            outputs=outputs,
            caches=caches,
        ),
        task=task,
        dataset=dataset,
    )


def test_two_readers_can_share_an_input_across_managers(tmp_path):
    shared = tmp_path / "shared.dat"
    store_path = tmp_path / "runs"
    first = _run(
        ExecutionManager(RunStore(store_path)),
        _Executor(),
        tmp_path,
        inputs=(shared,),
    )
    second = _run(
        ExecutionManager(RunStore(store_path)),
        _Executor(),
        tmp_path,
        inputs=(shared,),
    )

    first.submit()
    second.submit()

    assert first.record.state == RunState.QUEUED
    assert second.record.state == RunState.QUEUED


@pytest.mark.parametrize(
    ("first_files", "second_files"),
    [
        ({"outputs": ("shared.dat",)}, {"outputs": ("shared.dat",)}),
        ({"outputs": ("shared.dat",)}, {"inputs": ("shared.dat",)}),
        ({"inputs": ("shared.dat",)}, {"outputs": ("shared.dat",)}),
    ],
)
def test_writer_conflicts_with_reader_or_writer(tmp_path, first_files, second_files):
    store_path = tmp_path / "runs"
    first = _run(
        ExecutionManager(RunStore(store_path)),
        _Executor(),
        tmp_path,
        **first_files,
    )
    second = _run(
        ExecutionManager(RunStore(store_path)),
        _Executor(),
        tmp_path,
        **second_files,
    )

    first.submit()

    with pytest.raises(DuplicateSubmissionError) as exc:
        second.submit()

    assert exc.value.conflicts[0].id == first.id
    assert second.record.state == RunState.NEW


def test_task_and_dataset_do_not_create_a_collision_without_shared_files(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    first = _run(manager, executor, tmp_path, outputs=("one.dat",))
    second = _run(manager, executor, tmp_path, outputs=("two.dat",))

    first.submit()
    second.submit()

    assert first.record.state == RunState.QUEUED
    assert second.record.state == RunState.QUEUED


def test_existing_cache_is_materialized_as_shared_input(tmp_path):
    cache = tmp_path / "shared.tau"
    cache.write_text("ready")
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    first = _run(manager, executor, tmp_path, caches=(cache,))
    second = _run(manager, executor, tmp_path, caches=(cache,))

    first.submit()
    second.submit()

    assert first.record.files.input == (str(cache),)
    assert first.record.files.output == ()
    assert first.record.files.cache == (str(cache),)
    assert second.record.files.input == (str(cache),)


def test_missing_cache_is_materialized_as_output_and_conflicts(tmp_path):
    cache = tmp_path / "shared.tau"
    store_path = tmp_path / "runs"
    first = _run(
        ExecutionManager(RunStore(store_path)),
        _Executor(),
        tmp_path,
        caches=(cache,),
    )
    second = _run(
        ExecutionManager(RunStore(store_path)),
        _Executor(),
        tmp_path,
        caches=(cache,),
    )

    first.submit()

    assert first.record.files.input == ()
    assert first.record.files.output == (str(cache),)
    assert first.record.files.cache == (str(cache),)
    with pytest.raises(DuplicateSubmissionError):
        second.submit()


def test_cache_created_by_active_run_stays_write_locked(tmp_path):
    cache = tmp_path / "shared.tau"
    store_path = tmp_path / "runs"
    first = _run(
        ExecutionManager(RunStore(store_path)),
        _Executor(),
        tmp_path,
        caches=(cache,),
    )
    first.submit()

    # Simulate the backend creating the cache while it is still running. A new
    # reservation now sees the physical file, but must still conflict with the
    # persisted write role of the first run.
    cache.write_text("partial")
    second = _run(
        ExecutionManager(RunStore(store_path)),
        _Executor(),
        tmp_path,
        caches=(cache,),
    )

    with pytest.raises(DuplicateSubmissionError):
        second.submit()

    assert first.record.files.output == (str(cache),)


def test_file_pattern_conflicts_with_matching_concrete_output(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    first = _run(manager, executor, tmp_path, outputs=(tmp_path / "case*.agr",))
    second = _run(manager, executor, tmp_path, outputs=(tmp_path / "case1.agr",))

    first.submit()

    with pytest.raises(DuplicateSubmissionError):
        second.submit()


def test_final_run_releases_file_locks(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    first = _run(manager, executor, tmp_path, outputs=("shared.dat",))
    first.submit()
    executor.runners[0].finish()

    second = _run(manager, executor, tmp_path, outputs=("shared.dat",))
    second.submit()

    assert second.record.state == RunState.QUEUED


def test_submitting_record_contains_materialized_files_before_backend_submit(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _Executor()
    cache = tmp_path / "cache.dat"
    run = _run(
        manager,
        executor,
        tmp_path,
        inputs=("read.dat",),
        outputs=("write.dat",),
        caches=(cache,),
    )
    backend = executor.runners[0]
    seen = []

    original_submit = backend.submit

    def submit():
        seen.append(store.get(run.id))
        original_submit()

    backend.submit = submit
    run.submit()

    record = seen[0]
    assert record.state == RunState.SUBMITTING
    assert record.files.input == (str(tmp_path / "read.dat"),)
    assert record.files.output == (
        str(tmp_path / "write.dat"),
        str(cache),
    )
    assert record.files.cache == (str(cache),)
    assert record.submission_directory == str(tmp_path.resolve())

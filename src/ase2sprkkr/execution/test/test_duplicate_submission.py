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


def _run(
    manager,
    executor,
    directory,
    *,
    task="DOS",
    dataset="case",
    inputs=(),
    outputs=(),
    caches=(),
):
    spec = SimpleNamespace(
        directory=str(directory),
        input_file=None,
        output_file=None,
        input_files=tuple(inputs),
        output_files=tuple(outputs),
        cache_files=tuple(caches),
    )
    return manager.create_run(
        executor,
        reader=object(),
        spec=spec,
        task=task,
        dataset=dataset,
    )


def test_same_output_is_rejected_across_managers(tmp_path):
    store_path = tmp_path / "runs"
    output = tmp_path / "work" / "shared.dat"
    first_manager = ExecutionManager(RunStore(store_path))
    second_manager = ExecutionManager(RunStore(store_path))
    first = _run(first_manager, _Executor(), tmp_path / "work", outputs=(output,))
    second = _run(
        second_manager,
        _Executor(),
        tmp_path / "work",
        task="JXC",
        dataset="different",
        outputs=(output,),
    )

    first.submit()

    with pytest.raises(DuplicateSubmissionError) as exc:
        second.submit()

    assert exc.value.conflicts[0].id == first.id
    assert second.record.state == RunState.NEW


def test_disjoint_outputs_can_share_directory_even_for_same_task(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    first = _run(manager, executor, tmp_path, outputs=(tmp_path / "one.dat",))
    second = _run(manager, executor, tmp_path, outputs=(tmp_path / "two.dat",))

    first.submit()
    second.submit()

    assert first.record.state == RunState.QUEUED
    assert second.record.state == RunState.QUEUED


def test_input_input_is_shared_but_input_output_conflicts(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    shared = tmp_path / "shared.dat"
    shared.write_text("ready")

    reader1 = _run(manager, executor, tmp_path, inputs=(shared,))
    reader2 = _run(manager, executor, tmp_path, inputs=(shared,))
    writer = _run(manager, executor, tmp_path, outputs=(shared,))

    reader1.submit()
    reader2.submit()
    with pytest.raises(DuplicateSubmissionError):
        writer.submit()


def test_existing_cache_materializes_as_input_and_can_be_shared(tmp_path):
    cache = tmp_path / "shared.tau"
    cache.write_text("ready")
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    first = _run(manager, executor, tmp_path, caches=(cache,))
    second = _run(manager, executor, tmp_path, caches=(cache,))

    first.submit()
    second.submit()

    assert str(cache) in first.record.result_files.input_files
    assert str(cache) in second.record.result_files.input_files
    assert str(cache) in first.record.result_files.cache_files
    assert str(cache) not in first.record.result_files.output_files


def test_missing_cache_materializes_as_output_and_conflicts(tmp_path):
    cache = tmp_path / "shared.tau"
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    first = _run(manager, executor, tmp_path, caches=(cache,))
    second = _run(manager, executor, tmp_path, caches=(cache,))

    first.submit()
    with pytest.raises(DuplicateSubmissionError):
        second.submit()

    assert str(cache) in first.record.result_files.output_files
    assert str(cache) in first.record.result_files.cache_files


def test_cache_created_by_active_producer_still_conflicts_with_new_reader(tmp_path):
    cache = tmp_path / "shared.tau"
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    producer = _run(manager, executor, tmp_path, caches=(cache,))
    producer.submit()

    # Simulate the backend having created the file while still writing it. The
    # second run now classifies the cache as input, but must conflict with the
    # producer's persisted output role.
    cache.write_text("partial")
    consumer = _run(manager, executor, tmp_path, caches=(cache,))

    with pytest.raises(DuplicateSubmissionError):
        consumer.submit()


def test_cache_existence_is_sampled_at_reservation_time(tmp_path):
    cache = tmp_path / "shared.tau"
    cache.write_text("prepared")
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    run = _run(manager, executor, tmp_path, caches=(cache,))

    # create_run() must not freeze a stale cache role.
    cache.unlink()
    run.submit()

    assert str(cache) in run.record.result_files.output_files
    assert str(cache) not in run.record.result_files.input_files


def test_final_run_releases_file_locks(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _Executor()
    output = tmp_path / "result.dat"
    first = _run(manager, executor, tmp_path, outputs=(output,))
    first.submit()
    executor.runners[0].finish()

    second = _run(manager, executor, tmp_path, outputs=(output,))
    second.submit()

    assert second.record.state == RunState.QUEUED


def test_submitting_record_contains_materialized_files_before_backend_submit(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _Executor()
    cache = tmp_path / "cache.dat"
    output = tmp_path / "result.dat"
    run = _run(manager, executor, tmp_path, outputs=(output,), caches=(cache,))
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
    assert record.submission.directory == str(tmp_path.resolve())
    assert not hasattr(record.submission, "task")
    assert not hasattr(record.submission, "dataset")
    assert str(output) in record.result_files.output_files
    assert str(cache) in record.result_files.output_files
    assert str(cache) in record.result_files.cache_files

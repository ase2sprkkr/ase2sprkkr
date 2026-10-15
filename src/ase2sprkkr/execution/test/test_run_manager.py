import json
import os
import warnings
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    ExecutionManager,
    Executor,
    ResultFilesState,
    RunReference,
    RunState,
    RunStatus,
    RunStore,
)
from ase2sprkkr.outputs.persistent import load_persistent_result


class _FakeRunner:
    def __init__(self, callback=None):
        self.callback = callback
        self.state = RunState.NEW
        self.reference = None
        self.result = SimpleNamespace(return_code=None)
        self.exception = None
        self.on_submit = None

    def submit(self):
        if self.on_submit:
            self.on_submit()
        self.reference = RunReference(executor="test", native_id="native-1")
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

    def finish(self, return_code=0):
        self.result.return_code = return_code
        self.state = RunState.FINISHED if return_code == 0 else RunState.FAILED
        if self.callback:
            self.callback(self.result)


class _FakeExecutor(Executor):
    name = "test"

    def __init__(self):
        self.runners = []
        self.status_calls = []

    def create_run(self, reader, spec, callback=None):
        runner = _FakeRunner(callback=callback)
        self.runners.append(runner)
        return runner

    def status(self, *references):
        self.status_calls.append(references)
        return tuple(
            RunStatus(reference, RunState.RUNNING) for reference in references
        )


def _spec(tmp_path):
    return SimpleNamespace(directory=str(tmp_path))


def _result_files_spec(tmp_path):
    input_file = tmp_path / "run.inp"
    output_file = tmp_path / "run.out"
    generated = tmp_path / "case.dat"
    cache = tmp_path / "cache.dat"
    input_file.write_text("input")
    cache.write_text("cache-before")
    return SimpleNamespace(
        directory=str(tmp_path),
        input_file=str(input_file),
        output_file=str(output_file),
        input_files=(str(input_file),),
        output_files=(str(output_file), str(generated)),
        cache_files=(str(cache),),
    )


def test_persistent_run_is_first_saved_as_submitting(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    run = manager.create_run(
        executor,
        reader=object(),
        spec=_spec(tmp_path),
        task="ARPES",
        dataset="case",
    )
    backend = executor.runners[0]
    seen = []
    backend.on_submit = lambda: seen.append(store.get(run.id).state)

    run.submit()

    assert seen == [RunState.SUBMITTING]
    stored = store.get(run.id)
    assert stored.state == RunState.QUEUED
    assert stored.finished_at is None
    assert stored.reference == backend.reference
    assert stored.task == "ARPES"
    assert stored.dataset == "case"


def test_run_id_exists_before_submit_and_double_submit_is_rejected(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    executor = _FakeExecutor()
    run = manager.create_run(
        executor,
        reader=object(),
        spec=_spec(tmp_path),
        task="SCF",
    )

    assert run.id
    assert run.record.state == RunState.NEW
    assert run.record.finished_at is None

    run.submit()
    with pytest.raises(RuntimeError, match="already queued"):
        run.submit()


def test_completion_callback_persists_final_result(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    run = manager.create_run(
        executor,
        reader=object(),
        spec=_spec(tmp_path),
        task="DOS",
    )

    run.submit()
    executor.runners[0].finish(return_code=0)

    stored = store.get(run.id)
    assert stored.state == RunState.FINISHED
    assert stored.return_code == 0
    assert stored.finished_at == stored.updated_at


@pytest.mark.parametrize(
    "state", (RunState.FINISHED, RunState.FAILED, RunState.CANCELLED)
)
def test_first_terminal_transition_sets_immutable_finished_at(tmp_path, state):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    record = manager.new_record(task="DOS", directory=tmp_path)
    record = manager.transition(record, state=RunState.RUNNING)
    assert record.finished_at is None

    record = manager.transition(record, state=RunState.FINISHING)
    assert record.finished_at is None

    record = manager.transition(record, state=state)
    finished_at = record.finished_at
    assert finished_at == record.updated_at

    updated = manager.transition(record, message="later metadata update")
    assert updated.finished_at == finished_at


def test_completion_captures_inputs_and_outputs_separately(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _result_files_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS", dataset="case")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    (tmp_path / "case.dat").write_text("data")
    executor.runners[0].finish()

    stored = store.get(run.id)
    assert stored.state == RunState.FINISHED
    assert stored.result_files.state == ResultFilesState.CURRENT
    assert stored.result_files.captured_at is not None
    assert {item.path for item in stored.result_files.input_fingerprints} == {
        str(tmp_path / "run.inp")
    }
    assert {item.path for item in stored.result_files.fingerprints} == {
        str(tmp_path / "run.out"),
        str(tmp_path / "case.dat"),
    }
    assert str(tmp_path / "cache.dat") not in stored.result_files.manifests


def test_missing_captured_glob_member_counts_as_missing(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    input_file = tmp_path / "run.inp"
    output_file = tmp_path / "run.out"
    pattern = tmp_path / "case-*.dat"
    input_file.write_text("input")
    spec = SimpleNamespace(
        directory=str(tmp_path),
        input_file=str(input_file),
        output_file=str(output_file),
        input_files=(str(input_file),),
        output_files=(str(output_file), str(pattern)),
        cache_files=(),
    )
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    output_file.write_text("output")
    first = tmp_path / "case-1.dat"
    second = tmp_path / "case-2.dat"
    first.write_text("one")
    second.write_text("two")
    executor.runners[0].finish()
    stored = store.get(run.id)
    second.unlink()

    assert manager._missing_result_files(stored.result_files) == (str(second),)


def test_new_glob_member_does_not_change_historical_snapshot(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    input_file = tmp_path / "run.inp"
    output_file = tmp_path / "run.out"
    pattern = tmp_path / "case-*.dat"
    input_file.write_text("input")
    spec = SimpleNamespace(
        directory=str(tmp_path),
        input_file=str(input_file),
        output_file=str(output_file),
        input_files=(str(input_file),),
        output_files=(str(output_file), str(pattern)),
        cache_files=(),
    )
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    output_file.write_text("output")
    historical = tmp_path / "case-1.dat"
    historical.write_text("one")
    executor.runners[0].finish()
    later = tmp_path / "case-2.dat"
    later.write_text("from another calculation")

    updated = manager.check_result_files(run.id)

    assert updated.result_files.state == ResultFilesState.CURRENT
    assert {item.path for item in updated.result_files.fingerprints} == {
        str(output_file),
        str(historical),
    }


def test_changed_result_file_marks_overwritten_without_changing_outcome(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _result_files_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    (tmp_path / "case.dat").write_text("first")
    executor.runners[0].finish()
    finished = store.get(run.id)
    (tmp_path / "case.dat").write_text("changed contents")

    updated = manager.check_result_files(run.id)

    assert updated.state == RunState.FINISHED
    assert updated.return_code == 0
    assert updated.result_files.state == ResultFilesState.OVERWRITTEN
    assert updated.finished_at == finished.finished_at


def test_same_size_content_change_with_new_mtime_marks_result_overwritten(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _result_files_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    generated = tmp_path / "case.dat"
    generated.write_text("first")
    executor.runners[0].finish()
    original = generated.stat()
    generated.write_text("other")
    os.utime(
        generated,
        ns=(original.st_atime_ns, original.st_mtime_ns + 1_000_000),
    )

    updated = manager.check_result_files(run.id)

    assert updated.result_files.state == ResultFilesState.OVERWRITTEN


@pytest.mark.parametrize("filename", ["run.out", "case.dat"])
def test_output_metadata_change_has_uniform_overwritten_semantics(tmp_path, filename):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _result_files_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    (tmp_path / "case.dat").write_text("data")
    executor.runners[0].finish()
    path = tmp_path / filename
    stat = path.stat()
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))

    updated = manager.check_result_files(run.id)

    assert updated.result_files.state == ResultFilesState.OVERWRITTEN


def test_input_change_does_not_invalidate_current_output(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _result_files_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    (tmp_path / "case.dat").write_text("data")
    executor.runners[0].finish()
    input_file = tmp_path / "run.inp"
    input_file.write_text("changed input")

    updated = manager.check_result_files(run.id)

    assert updated.result_files.state == ResultFilesState.CURRENT
    assert manager.input_file_is_current(updated) is False


def test_changed_cache_does_not_mark_result_overwritten(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _result_files_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    (tmp_path / "case.dat").write_text("data")
    executor.runners[0].finish()
    (tmp_path / "cache.dat").write_text("cache-after-with-different-size")

    updated = manager.check_result_files(run.id)

    assert updated.state == RunState.FINISHED
    assert updated.result_files.state == ResultFilesState.CURRENT


def test_load_persistent_result_parses_output_without_backend_refresh(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _result_files_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    (tmp_path / "case.dat").write_text("data")
    executor.runners[0].finish()

    calls = []

    class Reader:
        def read_from_file(self, output, **kwargs):
            calls.append((output, kwargs))
            return SimpleNamespace(return_code=kwargs["return_code"])

    result = load_persistent_result(manager, run.id, reader=Reader())

    assert calls == [
        (
            str(tmp_path / "run.out"),
            {
                "return_code": 0,
                "input_file": str(tmp_path / "run.inp"),
                "directory": str(tmp_path),
            },
        )
    ]
    assert result.return_code == 0
    assert result.run_record.state == RunState.FINISHED
    assert executor.status_calls == []


def test_load_persistent_result_omits_stale_input_context(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _result_files_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    (tmp_path / "case.dat").write_text("data")
    executor.runners[0].finish()
    (tmp_path / "run.inp").write_text("changed input")
    calls = []

    class Reader:
        def read_from_file(self, output, **kwargs):
            calls.append(kwargs)
            return SimpleNamespace(return_code=kwargs["return_code"])

    result = load_persistent_result(manager, run.id, reader=Reader())

    assert calls[0]["input_file"] is None
    assert result.run_record.result_files.state == ResultFilesState.CURRENT


def test_refresh_batches_referenced_runs(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    representative = RunReference(executor="test", native_id="registered")
    manager.register_executor(executor, reference=representative)

    records = []
    for native_id in ("one", "two"):
        record = manager.new_record(task="DOS", directory=tmp_path)
        record = manager.transition(
            record,
            state=RunState.QUEUED,
            reference=RunReference(executor="test", native_id=native_id),
        )
        records.append(record)

    refreshed = manager.refresh(records)

    assert len(executor.status_calls) == 1
    assert executor.status_calls[0] == tuple(record.reference for record in records)
    assert [record.state for record in refreshed] == [
        RunState.RUNNING,
        RunState.RUNNING,
    ]


def test_run_store_shards_records_and_marks_schema(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="ARPES", directory=tmp_path)
    record = manager.transition(record, state=RunState.SUBMITTING, reference=None)

    path = tmp_path / "runs" / record.created_at[:10] / f"{record.id}.json"
    data = json.loads(path.read_text())

    assert path.is_file()
    assert data["schema_version"] == 3
    assert data["task"] == "ARPES"
    assert data["input_definition"] is None
    assert data["requested_task"] is None
    assert "task_definition" not in data
    assert data["finished_at"] is None


def test_run_store_migrates_schema_one_records(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    spec = _result_files_spec(tmp_path)
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=manager._result_files_for_spec(spec),
    )
    store.save(record)
    record = manager.finalize(record, state=RunState.FINISHED)
    path = tmp_path / "runs" / record.created_at[:10] / f"{record.id}.json"
    data = json.loads(path.read_text())
    data["schema_version"] = 1
    result_files = data["result_files"]
    result_files["fingerprints"] = [
        *result_files.pop("input_fingerprints"),
        *result_files["fingerprints"],
    ]
    path.write_text(json.dumps(data))

    migrated = store.get(record.id)

    assert migrated.schema_version == 3
    assert migrated.task == "DOS"
    assert migrated.input_definition is None
    assert migrated.requested_task is None
    assert {item.path for item in migrated.result_files.input_fingerprints} == {
        str(tmp_path / "run.inp")
    }
    assert str(tmp_path / "run.inp") not in {
        item.path for item in migrated.result_files.fingerprints
    }


def test_development_dual_task_schema_preserves_effective_task(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    record = manager.new_record(task="ARPES", directory=tmp_path)
    store.save(record)
    path = tmp_path / "runs" / record.created_at[:10] / f"{record.id}.json"
    data = json.loads(path.read_text())
    data["task"] = "SPLEED"
    data["task_definition"] = "ARPES"
    path.write_text(json.dumps(data))

    migrated = store.get(record.id)

    assert migrated.task == "SPLEED"
    assert migrated.requested_task == "ARPES"
    assert migrated.input_definition is None
    assert not hasattr(migrated, "task_definition")


def test_schema_two_development_digest_is_ignored(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    spec = _result_files_spec(tmp_path)
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=manager._result_files_for_spec(spec),
    )
    store.save(record)
    path = tmp_path / "runs" / record.created_at[:10] / f"{record.id}.json"
    data = json.loads(path.read_text())
    data["schema_version"] = 2
    data["result_files"]["output_digest"] = "obsolete-development-field"
    path.write_text(json.dumps(data))

    loaded = store.get(record.id)

    assert loaded.schema_version == 3
    assert not hasattr(loaded.result_files, "output_digest")


def test_run_store_listing_warns_once_and_skips_corrupt_record(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    valid = manager.new_record(task="SCF", directory=tmp_path)
    valid = manager.transition(valid, state=RunState.SUBMITTING, reference=None)
    broken = tmp_path / "runs" / valid.created_at[:10] / "broken.json"
    broken.write_text("{")

    with pytest.warns(RuntimeWarning, match="Corrupt run record"):
        assert store.runs() == (valid,)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert store.runs() == (valid,)
    assert not caught


def test_weaker_observation_does_not_degrade_final_record(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    reference = RunReference(executor="test", native_id="one")
    record = manager.new_record(task="DOS", directory=tmp_path)
    record = manager.transition(
        record,
        state=RunState.FINISHED,
        reference=reference,
        return_code=0,
        message="completed",
    )
    finished_at = record.finished_at

    updated = manager.update(record, RunStatus(reference, RunState.RUNNING))

    assert updated.state == RunState.FINISHED
    assert updated.return_code == 0
    assert updated.message == "completed"
    assert updated.finished_at == finished_at
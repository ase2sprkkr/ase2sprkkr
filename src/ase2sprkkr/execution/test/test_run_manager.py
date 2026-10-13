import json
import warnings
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    ArtifactState,
    ExecutionManager,
    Executor,
    RunReference,
    RunState,
    RunStatus,
    RunStore,
)


class _FakeRunner:
    def __init__(self, callback=None):
        self.callback = callback
        self.state = RunState.NEW
        self.reference = None
        self.result = SimpleNamespace(return_code=None)
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


def _artifact_spec(tmp_path):
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
        task="SCF",
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
    assert stored.task == "SCF"
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


def test_completion_captures_artifacts_but_not_cache(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _artifact_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS", dataset="case")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    (tmp_path / "case.dat").write_text("data")
    executor.runners[0].finish()

    stored = store.get(run.id)
    assert stored.state == RunState.FINISHED
    assert stored.artifacts.state == ArtifactState.CURRENT
    assert stored.artifacts.captured_at is not None
    assert {fingerprint.path for fingerprint in stored.artifacts.fingerprints} == {
        str(tmp_path / "run.inp"),
        str(tmp_path / "run.out"),
        str(tmp_path / "case.dat"),
    }
    assert str(tmp_path / "cache.dat") not in stored.artifacts.manifests


def test_changed_artifact_marks_overwritten_without_changing_outcome(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _artifact_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    (tmp_path / "case.dat").write_text("first")
    executor.runners[0].finish()
    finished = store.get(run.id)
    (tmp_path / "case.dat").write_text("changed contents")

    updated = manager.check_artifacts(run.id)

    assert updated.state == RunState.FINISHED
    assert updated.return_code == 0
    assert updated.artifacts.state == ArtifactState.OVERWRITTEN
    assert updated.finished_at == finished.finished_at


def test_changed_cache_does_not_mark_result_overwritten(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _artifact_spec(tmp_path)
    run = manager.create_run(executor, object(), spec, task="DOS")

    run.submit()
    (tmp_path / "run.out").write_text("output")
    (tmp_path / "case.dat").write_text("data")
    executor.runners[0].finish()
    (tmp_path / "cache.dat").write_text("cache-after-with-different-size")

    updated = manager.check_artifacts(run.id)

    assert updated.state == RunState.FINISHED
    assert updated.artifacts.state == ArtifactState.CURRENT


def test_load_result_parses_persisted_output_without_backend_refresh(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    spec = _artifact_spec(tmp_path)
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

    result = manager.load_result(run.id, reader=Reader())

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


def test_refresh_batches_referenced_runs(tmp_path):
    store = RunStore(tmp_path / "runs")
    manager = ExecutionManager(store)
    executor = _FakeExecutor()
    manager.register_executor(executor)

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
    record = manager.new_record(task="SCF", directory=tmp_path)
    record = manager.transition(record, state=RunState.SUBMITTING, reference=None)

    path = tmp_path / "runs" / record.created_at[:10] / f"{record.id}.json"
    data = json.loads(path.read_text())

    assert path.is_file()
    assert data["schema_version"] == 1
    assert data["finished_at"] is None


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

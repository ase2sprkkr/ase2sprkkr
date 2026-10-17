import json

import pytest

from ase2sprkkr.execution import CorruptRunRecord, RunRecord, RunState, RunStore, Execution, WorkingDirectory, ExecutorEndpoint, FileManifest


@pytest.mark.parametrize("schema_version", [0, 1, 2, 4, 5])
def test_only_current_run_record_schema_is_supported(tmp_path, schema_version):
    store = RunStore(tmp_path / "runs")
    created_at = "2026-10-17T12:34:56+00:00"
    directory = store.directory / "2026-10-17"
    directory.mkdir(parents=True)
    path = directory / "foreign.json"
    path.write_text(
        json.dumps(
            {
                "id": "foreign",
                "state": "finished",
                "created_at": created_at,
                "updated_at": created_at,
                "task": "DOS",
                "directory": str(tmp_path),
                "schema_version": schema_version,
            }
        )
    )

    with pytest.raises(CorruptRunRecord, match="unsupported schema version"):
        store.get("foreign")


def test_current_schema_round_trips_without_legacy_fields(tmp_path):
    store = RunStore(tmp_path / "runs")
    created_at = "2026-10-17T12:34:56+00:00"
    record = RunRecord(
        id="record",
        state=RunState.RUNNING,
        created_at=created_at,
        updated_at=created_at,
        task="DOS",
        execution=Execution(ExecutorEndpoint(), WorkingDirectory(str(tmp_path))),
        files=FileManifest(),
    )

    store.save(record)
    loaded = store.get(record.id)

    assert loaded == record
    assert loaded.schema_version == 3
    assert not hasattr(loaded, "backend_state")
    assert not hasattr(loaded, "requested_task")


def test_run_record_has_only_one_execution_identity(tmp_path):
    record = RunRecord(
        id="native", task="DOS", state=RunState.RUNNING,
        created_at="2026-10-17T12:34:56+00:00",
        updated_at="2026-10-17T12:34:56+00:00",
        execution=Execution(
            ExecutorEndpoint(), WorkingDirectory(str(tmp_path)), "backend-42"
        ),
        files=FileManifest(),
    )
    assert record.execution.backend_id == "backend-42"
    for obsolete in ("reference", "submission", "endpoint", "directory", "result_files"):
        assert not hasattr(record, obsolete)


def test_persisted_execution_identity_cannot_be_replaced(tmp_path):
    from ase2sprkkr.execution.manager import ExecutionManager

    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    record = manager.new_record(task="DOS", directory=tmp_path)
    first = manager.transition(
        record, state=RunState.RUNNING,
        execution=record.execution.submitted("job-42"),
    )
    assert manager.store.get(record.id).execution.backend_id == "job-42"
    with pytest.raises(ValueError, match="Backend identity"):
        manager.transition(
            first, execution=first.execution.submitted("job-42").__class__(
                first.execution.endpoint, first.execution.directory, "job-43"
            ),
        )
    assert manager.store.get(record.id).execution.backend_id == "job-42"

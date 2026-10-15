import json

from ase2sprkkr.execution import RunRecord, RunState, RunStore


def test_schema_v1_record_without_finished_at_remains_readable(tmp_path):
    store = RunStore(tmp_path / "runs")
    created_at = "2026-10-12T12:34:56+00:00"
    directory = store.directory / "2026-10-12"
    directory.mkdir(parents=True)
    path = directory / "legacy.json"
    path.write_text(
        json.dumps(
            {
                "id": "legacy",
                "state": "finished",
                "created_at": created_at,
                "updated_at": created_at,
                "task": "DOS",
                "directory": str(tmp_path),
                "schema_version": 1,
            }
        )
    )

    record = store.get("legacy")

    assert record.state == RunState.FINISHED
    assert record.finished_at is None
    assert record.input_definition is None
    assert record.requested_task is None


def test_legacy_requested_task_does_not_replace_effective_task(tmp_path):
    store = RunStore(tmp_path / "runs")
    created_at = "2026-10-12T12:34:56+00:00"
    directory = store.directory / "2026-10-12"
    directory.mkdir(parents=True)
    path = directory / "legacy-task.json"
    path.write_text(
        json.dumps(
            {
                "id": "legacy-task",
                "state": "finished",
                "created_at": created_at,
                "updated_at": created_at,
                "task": "SPLEED",
                "task_definition": "ARPES",
                "directory": str(tmp_path),
                "schema_version": 2,
            }
        )
    )

    record = store.get("legacy-task")

    assert record.task == "SPLEED"
    assert record.requested_task == "ARPES"
    assert record.input_definition is None


def test_run_record_constructor_does_not_require_finished_at(tmp_path):
    created_at = "2026-10-12T12:34:56+00:00"

    record = RunRecord(
        "legacy",
        RunState.RUNNING,
        created_at,
        created_at,
        "DOS",
        str(tmp_path),
    )

    assert record.task == "DOS"
    assert record.directory == str(tmp_path)
    assert record.finished_at is None

import json

from ase2sprkkr.execution import RunState, RunStore


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

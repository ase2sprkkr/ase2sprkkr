import json

import pytest

from ase2sprkkr.execution import CorruptRunRecord, RunRecord, RunState, RunStore


@pytest.mark.parametrize("schema_version", [0, 2, 3, 4, 5])
def test_other_run_record_schemas_are_not_supported(tmp_path, schema_version):
    store = RunStore(tmp_path / "runs")
    created_at = "2026-10-12T12:34:56+00:00"
    directory = store.directory / "2026-10-12"
    directory.mkdir(parents=True)
    path = directory / "record.json"
    path.write_text(
        json.dumps(
            {
                "id": "record",
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
        store.get("record")


def test_unknown_run_record_fields_are_rejected(tmp_path):
    store = RunStore(tmp_path / "runs")
    created_at = "2026-10-12T12:34:56+00:00"
    directory = store.directory / "2026-10-12"
    directory.mkdir(parents=True)
    path = directory / "record.json"
    path.write_text(
        json.dumps(
            {
                "id": "record",
                "state": "finished",
                "created_at": created_at,
                "updated_at": created_at,
                "task": "SPLEED",
                "input_definition": "ARPES",
                "requested_task": "ARPES",
                "directory": str(tmp_path),
                "schema_version": 1,
            }
        )
    )

    with pytest.raises(CorruptRunRecord, match="requested_task"):
        store.get("record")


def test_run_record_constructor_does_not_require_finished_at(tmp_path):
    created_at = "2026-10-12T12:34:56+00:00"

    record = RunRecord(
        "record",
        RunState.RUNNING,
        created_at,
        created_at,
        "DOS",
        str(tmp_path),
    )

    assert record.task == "DOS"
    assert record.directory == str(tmp_path)
    assert record.finished_at is None
    assert record.schema_version == 1
    assert not hasattr(record, "backend_state")

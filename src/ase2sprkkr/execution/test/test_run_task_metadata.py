from ase2sprkkr.execution import RunStore

from ase2sprkkr.execution.manager import ExecutionManager


def test_run_record_persists_task_and_schema(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    record = manager.new_record(
        task="spleed",
        input_definition="arpes",
        directory=tmp_path,
    )
    manager.store.save(record)

    loaded = manager.store.get(record.id)

    assert loaded.task == "SPLEED"
    assert loaded.input_definition == "ARPES"
    assert not hasattr(loaded, "requested_task")

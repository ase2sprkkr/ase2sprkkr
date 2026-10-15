from ase2sprkkr.execution import ExecutionManager, RunStore


def test_run_record_persists_effective_task_schema_and_requested_task(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    record = manager.new_record(
        task="spleed",
        input_definition="arpes",
        requested_task="arpes",
        directory=tmp_path,
    )
    manager.store.save(record)

    loaded = manager.store.get(record.id)

    assert loaded.task == "SPLEED"
    assert loaded.input_definition == "ARPES"
    assert loaded.requested_task == "ARPES"

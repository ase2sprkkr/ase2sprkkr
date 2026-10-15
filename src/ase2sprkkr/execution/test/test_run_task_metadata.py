from types import SimpleNamespace

from ase2sprkkr.execution import ExecutionManager, RunState, RunStore
from ase2sprkkr.outputs.persistent import load_persistent_result


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


def test_legacy_record_without_schema_does_not_guess_input_parameters(tmp_path):
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    input_file = tmp_path / "run.inp"
    output_file = tmp_path / "run.out"
    input_file.write_text("legacy input without persisted schema")
    output_file.write_text("legacy output")
    spec = SimpleNamespace(
        directory=str(tmp_path),
        input_file=str(input_file),
        output_file=str(output_file),
        input_files=(str(input_file),),
        output_files=(str(output_file),),
        cache_files=(),
    )
    record = manager.new_record(
        task="DOS",
        directory=tmp_path,
        result_files=manager._result_files_for_spec(spec),
    )
    manager.store.save(record)
    record = manager.finalize(record, state=RunState.FINISHED, return_code=0)

    seen = {}

    class Reader:
        input_parameters = None

        def read_from_file(self, output, **kwargs):
            seen["input_file"] = kwargs["input_file"]
            return SimpleNamespace(return_code=kwargs["return_code"])

    result = load_persistent_result(manager, record.id, reader=Reader())

    assert seen["input_file"] == str(input_file)
    assert result.run_record.id == record.id

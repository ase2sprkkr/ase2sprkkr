from types import SimpleNamespace

import pytest

from ase2sprkkr.common.module_registry import RegistryEntryNotFound
from ase2sprkkr.execution import (
    FileManifest,
    FileEntry,
    Execution,
    WorkingDirectory,
    ExecutorEndpoint,
    RunRecord,
    RunState,
)
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.outputs.project import (
    ProjectResultUnavailableError,
    attach_project_run,
    _reader_for_record,
)
from ase2sprkkr.outputs.readers.bsf import BsfOutputReader
from ase2sprkkr.outputs.readers.dos import DosOutputReader


class _Store:
    def __init__(self, record):
        self.record = record

    def get(self, run_id):
        assert run_id == self.record.id
        return self.record


class _Manager:
    def __init__(self, record, *, input_current=True):
        self.store = _Store(record)
        self.input_current = input_current

    def input_file_is_current(self, record):
        return self.input_current

    def attach(
        self,
        record,
        *,
        process_result,
        executor=None,
        callback=None,
        output_listener=None,
    ):
        return SimpleNamespace(record=record, process_result=process_result, executor=executor)


def _record(
    tmp_path,
    task,
    *,
    input_definition=None,
    input_file=None,
):
    files = (
        FileManifest()
        if input_file is None else FileManifest(
            inputs=(FileEntry(local_path=str(input_file)),)
        )
    )
    return RunRecord(
        id="run", task=task, input_definition=input_definition,
        state=RunState.RUNNING,
        created_at="2026-01-01T00:00:00+00:00",
        updated_at="2026-01-01T00:00:00+00:00",
        execution=Execution(
            ExecutorEndpoint(), WorkingDirectory(str(tmp_path)), "1",
        ),
        files=files,
    )


def _input_file(tmp_path, task):
    parameters = InputParameters.create(task)
    parameters.CONTROL.POTFIL = "Fe.pot"
    if task == "BSFEK":
        parameters.TASK.KPATH = 3
    path = tmp_path / "run.inp"
    path.write_text(parameters.to_string(validate=True))
    return path


def test_project_attach_uses_persisted_task(tmp_path):
    record = _record(
        tmp_path,
        "BSFEK",
        input_definition="bsf",
    )

    reader = _reader_for_record(_Manager(record), record)

    assert isinstance(reader, BsfOutputReader)
    assert reader.task_name == "BSFEK"


def test_project_attach_uses_registered_reader(tmp_path):
    record = _record(tmp_path, "DOS", input_definition="dos")

    reader = _reader_for_record(_Manager(record), record)

    assert isinstance(reader, DosOutputReader)
    assert reader.input_parameters is None
    assert reader.task_name == "DOS"


def test_project_attach_reconstructs_exact_input_definition(tmp_path):
    input_file = _input_file(tmp_path, "BSFEK")
    record = _record(
        tmp_path,
        "BSFEK",
        input_definition="bsf",
        input_file=input_file,
    )

    reader = _reader_for_record(_Manager(record), record)
    parameters = reader.input_parameters

    assert isinstance(reader, BsfOutputReader)
    assert parameters._definition is InputParameters.input_parameters_registry["BSF"]
    assert parameters.definition_name == "BSF"
    assert parameters.requested_task_name == "BSFEK"
    assert parameters.task_name == "BSFEK"


def test_project_attach_reconstructs_standard_input_definition(tmp_path):
    input_file = _input_file(tmp_path, "DOS")
    record = _record(
        tmp_path,
        "DOS",
        input_definition="dos",
        input_file=input_file,
    )

    reader = _reader_for_record(_Manager(record), record)
    parameters = reader.input_parameters

    assert isinstance(reader, DosOutputReader)
    assert parameters._definition is InputParameters.input_parameters_registry["DOS"]
    assert parameters.task_name == "DOS"


def test_project_attach_does_not_parse_stale_input(tmp_path):
    input_file = tmp_path / "run.inp"
    input_file.write_text("not valid SPR-KKR input")
    record = _record(
        tmp_path,
        "DOS",
        input_definition="dos",
        input_file=input_file,
    )

    reader = _reader_for_record(_Manager(record, input_current=False), record)

    assert isinstance(reader, DosOutputReader)
    assert reader.input_parameters is None


def test_project_attach_rejects_unknown_task(tmp_path):
    record = _record(tmp_path, "UNKNOWN_TASK", input_definition="dos")

    with pytest.raises(
        ProjectResultUnavailableError,
        match="unknown persisted task 'UNKNOWN_TASK'",
    ):
        attach_project_run(_Manager(record), record.id)


def test_project_attach_rejects_unknown_input_definition(tmp_path):
    input_file = _input_file(tmp_path, "DOS")
    record = _record(
        tmp_path,
        "DOS",
        input_definition="UNKNOWN_SCHEMA",
        input_file=input_file,
    )

    with pytest.raises(
        ProjectResultUnavailableError,
        match="unknown persisted input definition 'UNKNOWN_SCHEMA'",
    ):
        attach_project_run(_Manager(record), record.id)


def test_project_attach_rejects_input_with_different_task(tmp_path):
    input_file = _input_file(tmp_path, "DOS")
    record = _record(
        tmp_path,
        "SCF",
        input_definition="dos",
        input_file=input_file,
    )

    with pytest.raises(
        ProjectResultUnavailableError,
        match="persisted task does not match historical input",
    ):
        attach_project_run(_Manager(record), record.id)


def test_project_attach_does_not_hide_task_registry_programming_errors(
    monkeypatch, tmp_path
):
    class BrokenRegistry:
        def __getitem__(self, name):
            raise KeyError(f"broken task definition for {name}")

    monkeypatch.setattr(InputParameters, "task_registry", BrokenRegistry())
    record = _record(tmp_path, "DOS", input_definition="dos")

    with pytest.raises(KeyError, match="broken task definition"):
        attach_project_run(_Manager(record), record.id)


def test_project_attach_only_translates_missing_task_entries(monkeypatch, tmp_path):
    class MissingRegistry:
        def __getitem__(self, name):
            raise RegistryEntryNotFound(name)

    monkeypatch.setattr(InputParameters, "task_registry", MissingRegistry())
    record = _record(tmp_path, "DOS", input_definition="dos")

    with pytest.raises(ProjectResultUnavailableError):
        attach_project_run(_Manager(record), record.id)


def test_project_attach_passes_executor_to_manager(tmp_path):
    record = _record(tmp_path, "BSFEK", input_definition="bsf")
    executor = object()
    attached = attach_project_run(_Manager(record), record.id, executor=executor)
    assert attached.executor is executor
    assert callable(attached.process_result)

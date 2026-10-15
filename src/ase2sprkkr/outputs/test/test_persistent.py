from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import RunReference, RunState
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.outputs.persistent import attach_persistent_run
from ase2sprkkr.outputs.readers.bsf import BsfOutputReader


class _Store:
    def __init__(self, record):
        self.record = record

    def get(self, run_id):
        assert run_id == self.record.id
        return self.record


class _Manager:
    def __init__(self, record):
        self.store = _Store(record)

    def attach(
        self,
        record,
        *,
        reader,
        executor=None,
        callback=None,
        output_listener=None,
    ):
        return SimpleNamespace(record=record, reader=reader)


def _record(tmp_path, task, *, task_definition=None, input_file=None):
    return SimpleNamespace(
        id="run",
        task=task,
        task_definition=task_definition,
        directory=str(tmp_path),
        reference=RunReference(executor="test", native_id="1"),
        state=RunState.RUNNING,
        result_files=(
            None
            if input_file is None
            else SimpleNamespace(input_file=str(input_file))
        ),
    )


def test_persistent_attach_uses_legacy_task_alias_when_definition_not_persisted(tmp_path):
    record = _record(tmp_path, "BSFEK")

    attached = attach_persistent_run(_Manager(record), record.id)

    assert isinstance(attached.reader, BsfOutputReader)
    assert attached.reader.task_name == "BSFEK"


def test_persistent_attach_uses_bound_task_definition_when_native_task_differs(
    monkeypatch, tmp_path
):
    class SelectedReader:
        def __init__(self, input_parameters, calculator, directory, task_name=None):
            self.input_parameters = input_parameters
            self.task_name = task_name

    task = InputParameters.task_registry["DOS"]
    monkeypatch.setattr(task, "result_reader", SelectedReader)
    record = _record(tmp_path, "SPLEED", task_definition="DOS")

    attached = attach_persistent_run(_Manager(record), record.id)

    assert isinstance(attached.reader, SelectedReader)
    assert attached.reader.input_parameters is None
    assert attached.reader.task_name == "SPLEED"


def test_persistent_attach_reconstructs_input_parameters_when_available(
    monkeypatch, tmp_path
):
    class SelectedReader:
        def __init__(self, input_parameters, calculator, directory, task_name=None):
            self.input_parameters = input_parameters
            self.task_name = task_name

    task = InputParameters.task_registry["DOS"]
    parameters = SimpleNamespace(task_definition=task)
    input_file = tmp_path / "run.inp"
    input_file.write_text("historical input")
    monkeypatch.setattr(task, "result_reader", SelectedReader)
    monkeypatch.setattr(task, "read_input_parameters", lambda file: parameters)
    record = _record(
        tmp_path,
        "DOS",
        task_definition="DOS",
        input_file=input_file,
    )

    attached = attach_persistent_run(_Manager(record), record.id)

    assert isinstance(attached.reader, SelectedReader)
    assert attached.reader.input_parameters is parameters


def test_persistent_attach_does_not_hide_task_registry_programming_errors(
    monkeypatch, tmp_path
):
    class BrokenRegistry:
        def __getitem__(self, name):
            raise KeyError(f"broken task definition for {name}")

    monkeypatch.setattr(InputParameters, "task_registry", BrokenRegistry())
    record = _record(tmp_path, "DOS", task_definition="DOS")

    with pytest.raises(KeyError, match="broken task definition"):
        attach_persistent_run(_Manager(record), record.id)


def test_persistent_attach_respects_explicit_reader(tmp_path):
    record = _record(tmp_path, "BSFEK")
    reader = object()

    attached = attach_persistent_run(_Manager(record), record.id, reader=reader)

    assert attached.reader is reader

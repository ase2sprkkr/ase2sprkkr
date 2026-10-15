from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import RunReference, RunState
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.outputs.persistent import (
    PersistentResultUnavailableError,
    attach_persistent_run,
)
from ase2sprkkr.outputs.readers.bsf import BsfOutputReader


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
        reader,
        executor=None,
        callback=None,
        output_listener=None,
    ):
        return SimpleNamespace(record=record, reader=reader)


def _record(
    tmp_path,
    task,
    *,
    input_file=None,
    input_definition=None,
    requested_task=None,
):
    return SimpleNamespace(
        id="run",
        task=task,
        input_definition=input_definition,
        requested_task=requested_task,
        directory=str(tmp_path),
        reference=RunReference(executor="test", native_id="1"),
        state=RunState.RUNNING,
        result_files=(
            None
            if input_file is None
            else SimpleNamespace(input_file=str(input_file))
        ),
    )


def test_persistent_attach_uses_effective_bsf_task_reader(tmp_path):
    record = _record(tmp_path, "BSFEK")

    attached = attach_persistent_run(_Manager(record), record.id)

    assert isinstance(attached.reader, BsfOutputReader)
    assert attached.reader.task_name == "BSFEK"


def test_persistent_attach_uses_effective_task_definition_reader(monkeypatch, tmp_path):
    class SelectedReader:
        def __init__(self, input_parameters, calculator, directory, task_name=None):
            self.input_parameters = input_parameters
            self.task_name = task_name

    task = InputParameters.task_registry["DOS"]
    monkeypatch.setattr(task, "output_reader", SelectedReader)
    record = _record(tmp_path, "DOS")

    attached = attach_persistent_run(_Manager(record), record.id)

    assert isinstance(attached.reader, SelectedReader)
    assert attached.reader.input_parameters is None
    assert attached.reader.task_name == "DOS"


def test_persistent_attach_reconstructs_exact_input_schema_when_available(
    monkeypatch, tmp_path
):
    class SelectedReader:
        def __init__(self, input_parameters, calculator, directory, task_name=None):
            self.input_parameters = input_parameters
            self.task_name = task_name

    task = InputParameters.task_registry["DOS"]
    schema = InputParameters.input_parameters_registry["DOS"]
    parameters = SimpleNamespace(task_definition=task)
    input_file = tmp_path / "run.inp"
    input_file.write_text("historical input")
    monkeypatch.setattr(task, "output_reader", SelectedReader)

    seen = {}

    def read_from_file(file, *, task_name=None, **kwargs):
        seen["file"] = file
        seen["task_name"] = task_name
        return parameters

    monkeypatch.setattr(schema, "read_from_file", read_from_file)
    record = _record(
        tmp_path,
        "DOS",
        input_file=input_file,
        input_definition="DOS",
        requested_task="DOS",
    )

    attached = attach_persistent_run(_Manager(record), record.id)

    assert isinstance(attached.reader, SelectedReader)
    assert attached.reader.input_parameters is parameters
    assert seen == {"file": str(input_file), "task_name": "DOS"}


def test_persistent_attach_does_not_parse_stale_input(monkeypatch, tmp_path):
    class SelectedReader:
        def __init__(self, input_parameters, calculator, directory, task_name=None):
            self.input_parameters = input_parameters
            self.task_name = task_name

    task = InputParameters.task_registry["DOS"]
    schema = InputParameters.input_parameters_registry["DOS"]
    input_file = tmp_path / "run.inp"
    input_file.write_text("changed historical input")
    monkeypatch.setattr(task, "output_reader", SelectedReader)
    monkeypatch.setattr(
        schema,
        "read_from_file",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("stale input must not be parsed")
        ),
    )
    record = _record(
        tmp_path,
        "DOS",
        input_file=input_file,
        input_definition="DOS",
        requested_task="DOS",
    )

    attached = attach_persistent_run(
        _Manager(record, input_current=False),
        record.id,
    )

    assert isinstance(attached.reader, SelectedReader)
    assert attached.reader.input_parameters is None


def test_persistent_attach_rejects_unknown_effective_task(tmp_path):
    record = _record(tmp_path, "UNKNOWN_TASK")

    with pytest.raises(
        PersistentResultUnavailableError,
        match="unknown persisted effective task",
    ):
        attach_persistent_run(_Manager(record), record.id)


def test_persistent_attach_does_not_hide_task_registry_programming_errors(
    monkeypatch, tmp_path
):
    class BrokenRegistry:
        def __getitem__(self, name):
            raise KeyError(f"broken task definition for {name}")

    monkeypatch.setattr(InputParameters, "task_registry", BrokenRegistry())
    record = _record(tmp_path, "DOS")

    with pytest.raises(KeyError, match="broken task definition"):
        attach_persistent_run(_Manager(record), record.id)


def test_persistent_attach_respects_explicit_reader(tmp_path):
    record = _record(tmp_path, "BSFEK")
    reader = object()

    attached = attach_persistent_run(_Manager(record), record.id, reader=reader)

    assert attached.reader is reader

from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import RunRecord, RunReference, RunState
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.outputs.persistent import attach_persistent_run


class _Store:
    def __init__(self, record):
        self.record = record

    def get(self, run_id):
        assert run_id == self.record.id
        return self.record


class _Manager:
    def __init__(self, record):
        self.store = _Store(record)
        self.attached_reader = "unset"

    def attach(
        self,
        record,
        *,
        reader,
        executor=None,
        callback=None,
        output_listener=None,
    ):
        self.attached_reader = reader
        return SimpleNamespace(record=record, reader=reader)


def _record(tmp_path, task="DOS", *, state=RunState.RUNNING):
    return RunRecord(
        id="run",
        task=task,
        input_definition="UNKNOWN_SCHEMA",
        directory=str(tmp_path),
        state=state,
        created_at="2026-01-01T00:00:00+00:00",
        updated_at="2026-01-01T00:00:00+00:00",
        reference=RunReference(executor="test", native_id="1"),
    )


def test_persistent_attach_does_not_reconstruct_reader_or_input(monkeypatch, tmp_path):
    class BrokenRegistry:
        def __getitem__(self, name):
            raise AssertionError(f"registry must not be used while attaching: {name}")

    monkeypatch.setattr(InputParameters, "task_registry", BrokenRegistry())
    monkeypatch.setattr(InputParameters, "input_parameters_registry", BrokenRegistry())
    record = _record(tmp_path, task="UNKNOWN_TASK")
    manager = _Manager(record)

    attached = attach_persistent_run(manager, record.id)

    assert attached.reader is None
    assert manager.attached_reader is None


def test_persistent_attach_respects_explicit_live_reader(tmp_path):
    record = _record(tmp_path)
    manager = _Manager(record)
    reader = object()

    attached = attach_persistent_run(manager, record.id, reader=reader)

    assert attached.reader is reader
    assert manager.attached_reader is reader


@pytest.mark.parametrize(
    "state",
    [
        RunState.AVAILABLE,
        RunState.CHECKED,
        RunState.FINISHED,
        RunState.EXECUTION_FAILED,
        RunState.FAILED,
        RunState.CANCELLED,
    ],
)
def test_persistent_attach_rejects_terminal_runs(tmp_path, state):
    record = _record(tmp_path, state=state)

    with pytest.raises(RuntimeError, match="use result loading/checking"):
        attach_persistent_run(_Manager(record), record.id)

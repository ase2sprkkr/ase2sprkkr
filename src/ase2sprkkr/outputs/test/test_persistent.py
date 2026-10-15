from types import SimpleNamespace

import pytest

from ase2sprkkr.common.module_registry import RegistryEntryNotFound
from ase2sprkkr.execution import (
    ResultFilesRecord,
    RunRecord,
    RunReference,
    RunState,
)
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.outputs.persistent import (
    PersistentResultUnavailableError,
    attach_persistent_run,
)
from ase2sprkkr.outputs.readers.bsf import BsfOutputReader
from ase2sprkkr.outputs.readers.dos import DosOutputReader


_DOS_INPUT = """CONTROL  DATASET=Fe
         ADSI=Dirac
         POTFIL=Fe.pot
         KRWS=1
TAU      BZINT=POINTS
ENERGY   GRID={3}
         NE=500
         ImE=0.001
         EMIN=-0.2
TASK     DOS
"""

_BSFEK_INPUT = """CONTROL
    DATASET=Fe_BSFEK
    ADSI=BLOCHSF
    POTFIL=Fe.pot_new
    KRWS=1
    PRINT=0
TAU
    BZINT=POINTS
    NKTAB=250
TASK
    BSF
    NK=300
    KPATH=3
ENERGY
    GRID={3}
    NE={200}
    ImE=0.001
    EMIN=-0.2
    EMAX=0.2
SITES
    NL={3}
"""


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
    input_definition=None,
    requested_task=None,
    input_file=None,
):
    result_files = (
        None
        if input_file is None
        else ResultFilesRecord(input_file=str(input_file))
    )
    return RunRecord(
        id="run",
        task=task,
        requested_task=requested_task or task,
        input_definition=input_definition,
        directory=str(tmp_path),
        state=RunState.RUNNING,
        created_at="2026-01-01T00:00:00+00:00",
        updated_at="2026-01-01T00:00:00+00:00",
        reference=RunReference(executor="test", native_id="1"),
        result_files=result_files,
    )


def test_persistent_attach_uses_persisted_effective_task(tmp_path):
    record = _record(
        tmp_path,
        "BSFEK",
        requested_task="BSFKK",
        input_definition="BSF",
    )

    attached = attach_persistent_run(_Manager(record), record.id)

    assert isinstance(attached.reader, BsfOutputReader)
    assert attached.reader.task_name == "BSFEK"


def test_persistent_attach_uses_registered_reader(tmp_path):
    record = _record(tmp_path, "DOS", input_definition="DOS")

    attached = attach_persistent_run(_Manager(record), record.id)

    assert isinstance(attached.reader, DosOutputReader)
    assert attached.reader.input_parameters is None
    assert attached.reader.task_name == "DOS"


def test_persistent_attach_reconstructs_exact_input_definition(tmp_path):
    input_file = tmp_path / "run.inp"
    input_file.write_text(_BSFEK_INPUT)
    record = _record(
        tmp_path,
        "BSFEK",
        requested_task="BSFEK",
        input_definition="BSF",
        input_file=input_file,
    )

    attached = attach_persistent_run(_Manager(record), record.id)
    parameters = attached.reader.input_parameters

    assert isinstance(attached.reader, BsfOutputReader)
    assert parameters._definition is InputParameters.input_parameters_registry["BSF"]
    assert parameters.task_name == "BSFEK"
    assert parameters.task_definition.name == "BSFEK"


def test_persistent_attach_reconstructs_standard_input_definition(tmp_path):
    input_file = tmp_path / "run.inp"
    input_file.write_text(_DOS_INPUT)
    record = _record(
        tmp_path,
        "DOS",
        input_definition="DOS",
        input_file=input_file,
    )

    attached = attach_persistent_run(_Manager(record), record.id)
    parameters = attached.reader.input_parameters

    assert isinstance(attached.reader, DosOutputReader)
    assert parameters._definition is InputParameters.input_parameters_registry["DOS"]
    assert parameters.task_definition.name == "DOS"


def test_persistent_attach_does_not_parse_stale_input(tmp_path):
    input_file = tmp_path / "run.inp"
    input_file.write_text("not valid SPR-KKR input")
    record = _record(
        tmp_path,
        "DOS",
        input_definition="DOS",
        input_file=input_file,
    )

    attached = attach_persistent_run(
        _Manager(record, input_current=False),
        record.id,
    )

    assert isinstance(attached.reader, DosOutputReader)
    assert attached.reader.input_parameters is None


def test_persistent_attach_rejects_unknown_effective_task(tmp_path):
    record = _record(tmp_path, "UNKNOWN_TASK", input_definition="DOS")

    with pytest.raises(
        PersistentResultUnavailableError,
        match="unknown persisted effective task 'UNKNOWN_TASK'",
    ):
        attach_persistent_run(_Manager(record), record.id)


def test_persistent_attach_rejects_unknown_input_definition(tmp_path):
    input_file = tmp_path / "run.inp"
    input_file.write_text(_DOS_INPUT)
    record = _record(
        tmp_path,
        "DOS",
        input_definition="UNKNOWN_SCHEMA",
        input_file=input_file,
    )

    with pytest.raises(
        PersistentResultUnavailableError,
        match="unknown persisted input definition 'UNKNOWN_SCHEMA'",
    ):
        attach_persistent_run(_Manager(record), record.id)


def test_persistent_attach_rejects_input_with_different_effective_task(tmp_path):
    input_file = tmp_path / "run.inp"
    input_file.write_text(_DOS_INPUT)
    record = _record(
        tmp_path,
        "SCF",
        input_definition="DOS",
        input_file=input_file,
    )

    with pytest.raises(
        PersistentResultUnavailableError,
        match="persisted effective task does not match historical input",
    ):
        attach_persistent_run(_Manager(record), record.id)


def test_persistent_attach_does_not_hide_task_registry_programming_errors(
    monkeypatch, tmp_path
):
    class BrokenRegistry:
        def __getitem__(self, name):
            raise KeyError(f"broken task definition for {name}")

    monkeypatch.setattr(InputParameters, "task_registry", BrokenRegistry())
    record = _record(tmp_path, "DOS", input_definition="DOS")

    with pytest.raises(KeyError, match="broken task definition"):
        attach_persistent_run(_Manager(record), record.id)


def test_persistent_attach_only_translates_missing_task_entries(monkeypatch, tmp_path):
    class MissingRegistry:
        def __getitem__(self, name):
            raise RegistryEntryNotFound(name)

    monkeypatch.setattr(InputParameters, "task_registry", MissingRegistry())
    record = _record(tmp_path, "DOS", input_definition="DOS")

    with pytest.raises(PersistentResultUnavailableError):
        attach_persistent_run(_Manager(record), record.id)


def test_persistent_attach_respects_explicit_reader(tmp_path):
    record = _record(tmp_path, "BSFEK", input_definition="BSF")
    reader = object()

    attached = attach_persistent_run(_Manager(record), record.id, reader=reader)

    assert attached.reader is reader

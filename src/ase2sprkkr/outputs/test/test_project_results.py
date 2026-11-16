from functools import cache
from types import SimpleNamespace
import pytest
from ase2sprkkr.common.module_registry import RegistryEntryNotFound
from kkrunner.files import FileManifest, FileEntry
from kkrunner.runs import ExecutionMetadata, WorkingDirectory, RunRecord
from kkrunner.client.transports import LocalTransport
from kkrunner import RunState
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.outputs.project_results import ProjectResultUnavailableError, attach_project_run, _reader_for_record
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

    def attach(self, record, *, process_result, client=None, callback=None, output_listener=None):
        return SimpleNamespace(record=record, process_result=process_result, executor=client)

def _record(tmp_path, task, *, input_definition=None, input_file=None):
    files = FileManifest() if input_file is None else FileManifest(inputs=(FileEntry(local_path=str(input_file)),))
    return RunRecord(id='run', state=RunState.RUNNING, created_at='2026-01-01T00:00:00+00:00', updated_at='2026-01-01T00:00:00+00:00', execution=ExecutionMetadata(LocalTransport(), WorkingDirectory(str(tmp_path)), backend_id='1'), files=files, metadata={'ase2sprkkr': {'task': str(task).upper(), 'input_definition': str(input_definition).upper()}}, processing_required=True)

@cache
def _input_text(task):
    parameters = InputParameters.create(task)
    parameters.CONTROL.POTFIL = 'Fe.pot'
    if task == 'BSFEK':
        parameters.TASK.KPATH = 3
    return parameters.to_string(validate=True)

def _input_file(tmp_path, task):
    path = tmp_path / 'run.inp'
    path.write_text(_input_text(task))
    return path

@pytest.mark.parametrize(('task', 'schema', 'reader_type', 'with_input'), (('BSFEK', 'BSF', BsfOutputReader, False), ('DOS', 'DOS', DosOutputReader, False), ('BSFEK', 'BSF', BsfOutputReader, True), ('DOS', 'DOS', DosOutputReader, True)), ids=('bsf-no-input', 'dos-no-input', 'bsf-input', 'dos-input'))
def test_project_attach_uses_persisted_task_and_schema(tmp_path, task, schema, reader_type, with_input):
    input_file = _input_file(tmp_path, task) if with_input else None
    record = _record(tmp_path, task, input_definition=schema.lower(), input_file=input_file)
    reader = _reader_for_record(_Manager(record), record)
    assert isinstance(reader, reader_type)
    assert reader.task_name == task
    parameters = reader.input_parameters
    if not with_input:
        assert parameters is None
        return
    assert parameters._definition is InputParameters.input_parameters_registry[schema]
    assert parameters.definition_name == schema
    assert parameters.requested_task_name == task
    assert parameters.task_name == task

def test_project_attach_does_not_parse_stale_input(tmp_path):
    input_file = tmp_path / 'run.inp'
    input_file.write_text('not valid SPR-KKR input')
    record = _record(tmp_path, 'DOS', input_definition='dos', input_file=input_file)
    reader = _reader_for_record(_Manager(record, input_current=False), record)
    assert isinstance(reader, DosOutputReader)
    assert reader.input_parameters is None

def test_project_attach_rejects_unknown_task(tmp_path):
    record = _record(tmp_path, 'UNKNOWN_TASK', input_definition='dos')
    with pytest.raises(ProjectResultUnavailableError, match="unknown persisted task 'UNKNOWN_TASK'"):
        attach_project_run(_Manager(record), record.id)

def test_project_attach_rejects_unknown_input_definition(tmp_path):
    input_file = _input_file(tmp_path, 'DOS')
    record = _record(tmp_path, 'DOS', input_definition='UNKNOWN_SCHEMA', input_file=input_file)
    with pytest.raises(ProjectResultUnavailableError, match="unknown persisted input definition 'UNKNOWN_SCHEMA'"):
        attach_project_run(_Manager(record), record.id)

def test_project_attach_rejects_input_with_different_task(tmp_path):
    input_file = _input_file(tmp_path, 'DOS')
    record = _record(tmp_path, 'SCF', input_definition='dos', input_file=input_file)
    with pytest.raises(ProjectResultUnavailableError, match='persisted task does not match historical input'):
        attach_project_run(_Manager(record), record.id)

def test_project_attach_does_not_hide_task_registry_programming_errors(monkeypatch, tmp_path):

    class BrokenRegistry:

        def __getitem__(self, name):
            raise KeyError(f'broken task definition for {name}')
    monkeypatch.setattr(InputParameters, 'task_registry', BrokenRegistry())
    record = _record(tmp_path, 'DOS', input_definition='dos')
    with pytest.raises(KeyError, match='broken task definition'):
        attach_project_run(_Manager(record), record.id)

def test_project_attach_only_translates_missing_task_entries(monkeypatch, tmp_path):

    class MissingRegistry:

        def __getitem__(self, name):
            raise RegistryEntryNotFound(name)
    monkeypatch.setattr(InputParameters, 'task_registry', MissingRegistry())
    record = _record(tmp_path, 'DOS', input_definition='dos')
    with pytest.raises(ProjectResultUnavailableError):
        attach_project_run(_Manager(record), record.id)

def test_project_attach_passes_executor_to_manager(tmp_path):
    record = _record(tmp_path, 'BSFEK', input_definition='bsf')
    executor = object()
    attached = attach_project_run(_Manager(record), record.id, executor=executor)
    assert attached.executor is executor
    assert callable(attached.process_result)


def test_load_historical_scf_with_missing_auxiliary_output(tmp_path, monkeypatch):
    """Only the authoritative output is required to reconstruct an SCF result."""
    from types import SimpleNamespace
    from ase2sprkkr.outputs import project_results
    from kkrunner.files import FileVersion
    from kkrunner.files.manifest import FileEntry

    output = tmp_path / "scf.out"
    output.write_text("SCF output")
    missing_dos = tmp_path / "scf_SCFSTART.dos"
    files = FileManifest(
        stdout=FileEntry(str(output), availability_state="available"),
        outputs=(FileEntry(str(missing_dos), optional=True, availability_state="failed",
                           availability_message="Historical artifact was not retrieved"),),
        output_versions=(FileVersion.capture(output, require_digest=True),),
        captured_at="2026-11-13T00:00:00+00:00",
    )
    record = _record(tmp_path, "SCF")
    record = __import__("dataclasses").replace(record, state=RunState.FINISHED, files=files, return_code=0)
    result = SimpleNamespace()
    calls = []

    class Manager(_Manager):
        def retrieve_result_files(self, rec):
            calls.append("retrieve")
            return rec

        def output_file_is_current(self, rec):
            return True

    class Reader:
        input_parameters = None

        def read_from_file(self, output_file, **kwargs):
            calls.append(("parse", output_file))
            return result

    manager = Manager(record)
    monkeypatch.setattr(project_results, "_reader_for_record", lambda _manager, _record: Reader())
    loaded = project_results.load_project_result(manager, record.id)
    assert loaded is result
    assert loaded.run_record is record
    assert calls == ["retrieve", ("parse", str(output))]


@pytest.mark.parametrize("auxiliary_state", ("current", "overwritten", "missing", "unrecorded"))
def test_historical_auxiliary_actions_require_original_file(tmp_path, auxiliary_state):
    """Use the real OUT reader and Grace parser across the persistence boundary."""
    from dataclasses import replace
    from pathlib import Path

    from kkrunner import RunManager, RunStore
    from kkrunner.files import FileVersion
    from kkrunner.runs import HistoricalFileMismatchError

    from ase2sprkkr.outputs import capture_result_file_versions
    from ase2sprkkr.outputs.project_results import load_project_result

    output = tmp_path / "scf.out"
    text = (Path(__file__).parents[1] / "examples" / "Al_SCF.out").read_text()
    output.write_text(text.replace("Al.pot_new\n", "Al.pot_new\n          curve : (9) curve.agr\n", 1))
    auxiliary = tmp_path / "curve.agr"
    auxiliary.write_text("@type xy\n0 1\n&\n")
    output_version = FileVersion.capture(output, require_digest=True)
    auxiliary_version = FileVersion.capture(auxiliary, require_digest=True)
    versions = (output_version,) if auxiliary_state == "unrecorded" else (output_version, auxiliary_version)
    files = FileManifest(
        stdout=FileEntry(str(output)), outputs=(FileEntry(str(auxiliary), optional=True),),
        output_versions=versions, captured_at="2026-11-16T00:00:00+00:00",
    )
    manager = RunManager(RunStore(tmp_path / "runs"))
    record = replace(manager.new_record(directory=tmp_path, files=files,
                     metadata={"ase2sprkkr": {"task": "SCF"}}), state=RunState.FINISHED)
    manager.store.save(record)
    if auxiliary_state == "overwritten":
        auxiliary.write_text("@type xy\n0 999\n&\n")
    elif auxiliary_state == "missing":
        auxiliary.unlink()

    result = load_project_result(manager, record.id)
    assert result.run_record.return_code is None  # Unknown OS status still permits OUT parsing.
    option = result.files["curve"]
    captured = capture_result_file_versions(result)
    if auxiliary_state == "unrecorded":
        assert "curve" not in captured
    else:
        assert captured["curve"] == auxiliary_version
    if auxiliary_state == "current":
        assert "plot" in option.actions()
        assert option.parsed_file() is not None
    else:
        assert option.actions() == ()
        with pytest.raises(HistoricalFileMismatchError):
            option.parsed_file()
        with pytest.raises(HistoricalFileMismatchError):
            option.open()
        with pytest.raises(HistoricalFileMismatchError):
            option.save()

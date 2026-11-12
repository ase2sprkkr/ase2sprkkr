"""Project-scoped SPR-KKR runs using the KKRunner EndpointClient/Run lifecycle."""

from dataclasses import replace
from types import SimpleNamespace
import sys

from kkrunner import RunState
from kkrunner.execution import MpiSpec
from kkrunner.projects import ProjectRunManager, ProjectRepository

from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.input_parameters.task_definition import TaskDefinition
from ase2sprkkr.sprkkr import calculation as calculation_module
from ase2sprkkr.sprkkr.calculation import PreparedCalculation
from ase2sprkkr.sprkkr.calculator import SPRKKR
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


class _Prepared:
    create_run = PreparedCalculation.create_run

    def __init__(self, parameters, input_file, output_file, directory):
        self.input_parameters = parameters
        self.input_file = input_file
        self.potential_file = None
        self.output_file = output_file
        self.directory = str(directory)
        self.closed = False

    def close(self):
        if not self.closed:
            self.closed = True
            self.input_file.close()
            self.output_file.close()


def _parameters(monkeypatch, task="DOS"):
    parameters = InputParameters.create(task)
    parameters.CONTROL.DATASET = "sample"
    parameters.CONTROL.POTFIL = "sample.pot"
    monkeypatch.setattr(
        TaskDefinition, "output_reader",
        lambda self, *args, **kwargs: SimpleNamespace(
            read_from_file=lambda *args, **kwargs: SimpleNamespace(return_code=0)
        ),
    )
    monkeypatch.setattr(InputParameters, "input_files", lambda self, directory: ())
    monkeypatch.setattr(InputParameters, "output_files", lambda self, directory: ())
    monkeypatch.setattr(InputParameters, "cache_files", lambda self, directory: ())
    return parameters


def _manager(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Calculation project")
    return ProjectRunManager(project_repository=projects).for_project(project)


def _use_python_as_test_executable(monkeypatch):
    """Test the real client lifecycle without requiring an SPR-KKR binary."""
    original = SPRKKRProfile.create_execution

    def create_execution(self, *args, **kwargs):
        spec = original(self, *args, **kwargs)
        return replace(
            spec, argv=(sys.executable, "-c", "print('test calculation')"),
            mpi=MpiSpec(), mpi_argv=None, mpi_uses_stdin=True, stdin=None,
        )

    monkeypatch.setattr(SPRKKRProfile, "create_execution", create_execution)


def _runner(parameters, tmp_path, manager):
    input_file = (tmp_path / "test.inp").open("w+")
    input_file.write("input")
    input_file.flush()
    output_file = (tmp_path / "test.out").open("wb")
    prepared = _Prepared(parameters, input_file, output_file, tmp_path)
    return prepared.create_run(
        result_context=SimpleNamespace(directory=str(tmp_path)),
        profile=SPRKKRProfile(mpi=False), project=manager,
    )


def _project_record(parameters, tmp_path, monkeypatch):
    _use_python_as_test_executable(monkeypatch)
    manager = _manager(tmp_path)
    runner = _runner(parameters, tmp_path, manager)
    assert runner.id
    runner.run()
    return manager.store.get(runner.id)


def test_project_persists_task_and_dataset(tmp_path, monkeypatch):
    parameters = _parameters(monkeypatch)
    record = _project_record(parameters, tmp_path, monkeypatch)
    assert record.metadata["ase2sprkkr"]["task"] == parameters.task_name
    assert record.metadata["ase2sprkkr"]["input_definition"] == parameters.definition_name.upper()
    assert record.metadata["ase2sprkkr"]["dataset"] == "sample"
    assert record.state == RunState.FINISHED
    assert record.project_id is not None


def test_project_preserves_native_task_identity(tmp_path, monkeypatch):
    parameters = _parameters(monkeypatch, "ARPES")
    parameters.change_task("SPLEED")
    record = _project_record(parameters, tmp_path, monkeypatch)
    assert record.metadata["ase2sprkkr"]["task"] == "SPLEED"


def test_distinct_task_identity_with_shared_schema(tmp_path, monkeypatch):
    parameters = _parameters(monkeypatch, "BSFEK")
    record = _project_record(parameters, tmp_path, monkeypatch)
    assert parameters.TASK.TASK() == "BSF"
    assert record.metadata["ase2sprkkr"]["task"] == "BSFEK"
    assert record.metadata["ase2sprkkr"]["input_definition"] == "BSF"


def test_bound_manager_is_used_as_persistence_authority(tmp_path, monkeypatch):
    parameters = _parameters(monkeypatch)
    _use_python_as_test_executable(monkeypatch)
    manager = _manager(tmp_path)
    runner = _runner(parameters, tmp_path, manager)
    runner.run()
    assert manager.store.get(runner.id).project_id == manager.project_id


def test_calculator_forwards_persistence_settings(tmp_path, monkeypatch):
    manager = _manager(tmp_path)
    client = object()
    calculator = SPRKKR(directory=tmp_path, executor=client, project=manager)
    captured = {}

    class Prepared:
        def create_run(self, **kwargs):
            captured.update(kwargs)
            return object()

    monkeypatch.setattr(calculator, "prepare", lambda **kwargs: Prepared())
    calculator.create_run()
    assert captured["project"] is manager
    assert captured["executor"] is client


def test_explicit_run_profile_selects_client_over_calculator_client(tmp_path, monkeypatch):
    calculator_client = object()
    selected_client = object()
    profile = SPRKKRProfile(executor="ssh", host="cluster")
    calculator = SPRKKR(directory=tmp_path, executor=calculator_client)
    captured = {}

    class Prepared:
        def create_run(self, **kwargs):
            captured.update(kwargs)
            return object()

    monkeypatch.setattr(calculator, "prepare", lambda **kwargs: Prepared())
    calculator.create_run(profile=profile)
    assert captured["profile"] is profile
    assert captured["executor"] is None

    monkeypatch.setattr(calculation_module, "executor_from_profile", lambda selected: selected_client)
    assert calculation_module.executor_from_profile(profile) is selected_client

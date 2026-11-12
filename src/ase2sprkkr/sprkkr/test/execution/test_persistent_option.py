"""Project-only persistence and per-run execution selection contracts."""
from kkrunner import RunStore
from kkrunner.projects import ProjectRunManager, ProjectRepository

from ase2sprkkr.sprkkr.calculator import SPRKKR
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


def _manager(tmp_path):
    return ProjectRunManager(RunStore(tmp_path / 'runs'), project_repository=ProjectRepository(tmp_path / 'projects'))


def test_calculator_forwards_default_project_and_executor(tmp_path, monkeypatch):
    manager = _manager(tmp_path)
    project = manager.projects.create('Calculation project')
    bound = manager.for_project(project)
    executor = object()
    calculator = SPRKKR(directory=tmp_path, executor=executor, project=bound)
    captured = {}

    class Prepared:

        def create_run(self, **kwargs):
            captured.update(kwargs)
            return object()
    monkeypatch.setattr(calculator, 'prepare', lambda **kwargs: Prepared())
    calculator.create_run()
    assert captured['project'] is bound
    assert captured['executor'] is executor


def test_run_level_profile_overrides_calculator_executor(tmp_path, monkeypatch):
    executor = object()
    calculator = SPRKKR(directory=tmp_path, executor=executor)
    captured = {}

    class Prepared:

        def create_run(self, **kwargs):
            captured.update(kwargs)
            return object()
    monkeypatch.setattr(calculator, 'prepare', lambda **kwargs: Prepared())
    profile = SPRKKRProfile(executor='ssh', host='cluster')
    calculator.create_run(profile=profile)
    assert captured['executor'] is None
    assert captured['profile'] is profile

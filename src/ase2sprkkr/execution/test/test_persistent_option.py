"""Project-only persistence and per-run execution selection contracts."""

import pytest

from ase2sprkkr.execution import ExecutionManager, ProjectRepository, RunStore
from ase2sprkkr.execution.persistence import resolve_project_manager
from ase2sprkkr.sprkkr.calculator import SPRKKR
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


def _manager(tmp_path):
    return ExecutionManager(
        RunStore(tmp_path / "runs"),
        project_repository=ProjectRepository(tmp_path / "projects"),
    )


def test_persistence_requires_an_existing_project(tmp_path):
    manager = _manager(tmp_path)
    assert resolve_project_manager(None) is None
    assert resolve_project_manager(False) is None
    with pytest.raises(ValueError, match="bound to a project"):
        resolve_project_manager(manager)
    with pytest.raises(ValueError, match="require a project"):
        manager.new_record()


def test_project_bound_manager_keeps_store_and_project_identity(tmp_path):
    manager = _manager(tmp_path)
    project = manager.projects.create("Calculation project")
    bound = manager.for_project(project)

    assert bound.project_id == project.id
    assert bound.store is manager.store
    assert bound.for_project(project) is bound
    assert resolve_project_manager(bound) is bound


def test_project_manager_cannot_be_rebound(tmp_path):
    manager = _manager(tmp_path)
    first = manager.projects.create("First")
    second = manager.projects.create("Second")

    with pytest.raises(ValueError, match="Cannot replace"):
        manager.for_project(first).for_project(second)


def test_calculator_forwards_default_project_and_executor(tmp_path, monkeypatch):
    manager = _manager(tmp_path)
    project = manager.projects.create("Calculation project")
    bound = manager.for_project(project)
    executor = object()
    calculator = SPRKKR(directory=tmp_path, executor=executor, project=bound)
    captured = {}

    class Prepared:
        def create_run(self, **kwargs):
            captured.update(kwargs)
            return object()

    monkeypatch.setattr(calculator, "prepare", lambda **kwargs: Prepared())
    calculator.create_run()

    assert captured["project"] is bound
    assert captured["executor"] is executor


def test_run_level_profile_overrides_calculator_executor(tmp_path, monkeypatch):
    executor = object()
    calculator = SPRKKR(directory=tmp_path, executor=executor)
    captured = {}

    class Prepared:
        def create_run(self, **kwargs):
            captured.update(kwargs)
            return object()

    monkeypatch.setattr(calculator, "prepare", lambda **kwargs: Prepared())
    profile = SPRKKRProfile(executor="ssh", host="cluster")
    calculator.create_run(profile=profile)

    assert captured["executor"] is None
    assert captured["profile"] is profile

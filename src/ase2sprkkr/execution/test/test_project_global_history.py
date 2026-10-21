"""Global run history and cross-project namespace regression tests."""

from dataclasses import replace

import pytest

from ..base import RunState
from ..manager import DuplicateSubmissionError, RunStore
from ..project_manager import ExecutionManager
from ..projects import ProjectRepository, ProjectDeleteError


def test_different_projects_share_one_store_and_keep_project_identity(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    root = projects.create("Root")
    first = projects.create("First", parent_id=root.id)
    second = projects.create("Second", parent_id=root.id)
    a = ExecutionManager(project_repository=projects, project=first)
    b = ExecutionManager(project_repository=projects, project=second)
    assert a.store.directory == b.store.directory
    first_run = a.new_record(task="DOS", directory=tmp_path / "work", dataset="same")
    second_run = b.new_record(task="DOS", directory=tmp_path / "work", dataset="same")
    assert first_run.project_id == first.id
    assert second_run.project_id == second.id
    a.store.save(first_run)
    b.store.save(second_run)
    assert {r.project_id for r in a.store.runs()} == {first.id, second.id}


def test_cross_project_active_conflict_and_project_deletion(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    root = projects.create("Root")
    first = projects.create("First", parent_id=root.id)
    second = projects.create("Second", parent_id=root.id)
    a = ExecutionManager(project_repository=projects, project=first)
    b = ExecutionManager(project_repository=projects, project=second)
    first_record = a.new_record(task="SCF", directory=tmp_path / "work")
    second_record = b.new_record(task="DOS", directory=tmp_path / "work")
    a.reserve(first_record)
    with pytest.raises(DuplicateSubmissionError):
        b.reserve(second_record)
    with pytest.raises(ProjectDeleteError, match="active"):
        projects.delete(first.id, expected_generation=first.generation)
    a.complete(first_record, state=RunState.FAILED)
    projects.delete(first.id, expected_generation=first.generation)
    assert all(r.project_id != first.id for r in b.store.runs())


def test_for_project_preserves_explicit_run_store(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Project")
    custom_store = RunStore(tmp_path / "custom-runs")
    original = ExecutionManager(
        store=custom_store, project_repository=projects, submitting_grace=17,
    )

    bound = original.for_project(project)

    assert bound.store is custom_store
    assert bound.submitting_grace == 17
    assert bound.executor_pool is original.executor_pool
    assert bound.for_project(project) is bound

    record = bound.new_record(task="DOS", directory=tmp_path / "work")
    bound.store.save(record)

    assert custom_store.get(record.id) == record
    assert record.project_id == project.id
    assert not (tmp_path / "projects_runs").exists()


def test_for_project_without_explicit_store_uses_shared_project_store(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Project")
    unbound = ExecutionManager(project_repository=projects)

    bound = unbound.for_project(project)

    assert unbound.store is None
    assert bound.project_id == project.id
    assert bound.store.directory == tmp_path / "projects_runs"



def test_unknown_run_prevents_project_deletion_until_terminal(tmp_path):
    """UNKNOWN is recoverable, so deleting its project must not erase history."""
    projects = ProjectRepository(tmp_path / "projects")
    root = projects.create("Root")
    project = projects.create("Child", parent_id=root.id)
    manager = ExecutionManager(project_repository=projects, project=project)
    record = replace(
        manager.new_record(task="DOS", directory=tmp_path / "work"),
        state=RunState.UNKNOWN,
    )
    manager.store.save(record)

    with pytest.raises(ProjectDeleteError, match="active"):
        projects.delete(project.id, expected_generation=project.generation)
    assert manager.store.get(record.id).state == RunState.UNKNOWN

    manager.store.update(
        record, lambda latest: replace(latest, state=RunState.FINISHED),
    )
    projects.delete(project.id, expected_generation=project.generation)
    assert manager.store.runs() == ()

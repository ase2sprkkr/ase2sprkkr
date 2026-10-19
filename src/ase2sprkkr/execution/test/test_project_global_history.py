"""Global run history and cross-project namespace regression tests."""

import pytest

from ..base import RunState
from ..manager import DuplicateSubmissionError
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


def test_deleted_project_cannot_recreate_history_from_stale_runner(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    root = projects.create("Root")
    child = projects.create("Child", parent_id=root.id)
    manager = ExecutionManager(project_repository=projects, project=child)
    record = manager.new_record(task="DOS", directory=tmp_path / "work")
    manager.store.save(record)
    record = manager.complete(record, state=RunState.FAILED)

    projects.delete(child.id, expected_generation=child.generation)
    with pytest.raises(KeyError, match="no longer exists"):
        manager.store.save(record)
    with pytest.raises(KeyError, match="no longer exists"):
        manager.transition(record, state=RunState.FINISHED)
    assert manager.store.runs() == ()

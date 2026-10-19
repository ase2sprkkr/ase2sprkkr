"""Project tracking is opt-in by project identity, never by a persistent flag."""

import pytest

from ase2sprkkr.execution.persistence import resolve_project_manager
from ase2sprkkr.execution.project_manager import ExecutionManager
from ase2sprkkr.execution.projects import ProjectRepository


def test_untracked_run_has_no_store():
    assert resolve_project_manager(None) is None
    assert resolve_project_manager(False) is None
    with pytest.raises(TypeError):
        resolve_project_manager(True)


def test_project_uses_shared_project_runs_store(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("test")
    manager = ExecutionManager(project_repository=projects).for_project(project)
    assert manager.project_id == project.id
    assert manager.store.directory == tmp_path / "projects_runs"
    with pytest.raises(ValueError, match="existing project"):
        ExecutionManager(project_repository=projects).create_run(
            object(), object(), process_result=lambda x: x, task="SCF",
        )


def test_project_schema_has_no_run_reference_list(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("test")
    assert not hasattr(project, "persistent_refs")
    text = (tmp_path / "projects" / f"{project.id}.json").read_text()
    assert "persistent_refs" not in text

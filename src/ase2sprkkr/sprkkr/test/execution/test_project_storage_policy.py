"""Existing project definitions and new runner records have independent storage."""

from kkrunner.projects import ProjectRepository

from ase2sprkkr.sprkkr.projects import default_project_repository, project_execution_manager, resolve_project_manager


def test_default_projects_remain_visible_without_loading_old_development_runs(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_STATE_HOME", str(tmp_path))
    projects = ProjectRepository(tmp_path / "ase2sprkkr/projects")
    project = projects.create("Existing", directory=tmp_path)
    old = tmp_path / "ase2sprkkr/projects_runs/2026-01-01/old.json"
    old.parent.mkdir(parents=True)
    old.write_text('{"schema_version":3}')

    manager = resolve_project_manager(project.id)
    assert manager.project_id == project.id
    assert default_project_repository().project(project.id) == project
    assert manager.store.directory == tmp_path / "kkrunner/projects_runs"
    assert manager.store.runs() == ()
    assert old.read_text() == '{"schema_version":3}'


def test_explicit_project_repository_keeps_its_explicit_store(tmp_path):
    projects = ProjectRepository(tmp_path / "custom/projects")
    manager = project_execution_manager(projects)
    assert manager.store.directory == tmp_path / "custom/projects_runs"

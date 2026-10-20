"""Global run history and cross-project namespace regression tests."""

from concurrent.futures import ThreadPoolExecutor
from threading import Event

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


def test_for_project_preserves_explicit_canonical_store(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Project")
    shared_store = RunStore(projects.runs_directory)
    original = ExecutionManager(
        store=shared_store, project_repository=projects, submitting_grace=17,
    )

    bound = original.for_project(project)
    assert bound.store is shared_store
    assert bound.store.projects_directory == projects.directory
    assert bound.submitting_grace == 17
    assert bound.executor_pool is original.executor_pool
    assert bound.for_project(project) is bound

    record = bound.new_record(task="DOS", directory=tmp_path / "work")
    bound.store.save(record)
    assert shared_store.get(record.id) == record
    with pytest.raises(ProjectDeleteError, match="active"):
        projects.delete(project.id, expected_generation=project.generation)


def test_project_bound_store_rejects_noncanonical_location(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Project")
    custom_store = RunStore(tmp_path / "custom-runs")
    with pytest.raises(ValueError, match="shared projects_runs"):
        ExecutionManager(store=custom_store, project_repository=projects, project=project)

    original = ExecutionManager(store=custom_store, project_repository=projects)
    with pytest.raises(ValueError, match="shared projects_runs"):
        original.for_project(project)
    assert not custom_store.directory.exists()


def test_explicit_store_must_belong_to_same_projects_repository(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Project")
    other = tmp_path / "another-projects"
    store = RunStore(projects.runs_directory, projects_directory=other)
    with pytest.raises(ValueError, match="another project repository"):
        ExecutionManager(store=store, project_repository=projects, project=project)


def test_for_project_without_explicit_store_uses_shared_project_store(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Project")
    unbound = ExecutionManager(project_repository=projects)

    bound = unbound.for_project(project)

    assert unbound.store is None
    assert bound.store.directory == tmp_path / "projects_runs"


def test_delete_waits_for_completion_write_before_checking_runs(tmp_path, monkeypatch):
    """Deletion must observe the state written under the project guard."""
    projects = ProjectRepository(tmp_path / "projects")
    root = projects.create("Root")
    child = projects.create("Child", parent_id=root.id)
    manager = ExecutionManager(project_repository=projects, project=child)
    record = manager.new_record(task="DOS", directory=tmp_path / "work")
    manager.store.save(record)
    writing, release, deleting = Event(), Event(), Event()
    original_write = manager.store._write_unlocked

    def pause_write(value):
        writing.set()
        assert release.wait(10)
        return original_write(value)

    monkeypatch.setattr(manager.store, "_write_unlocked", pause_write)

    def delete():
        deleting.set()
        return projects.delete(child.id, expected_generation=child.generation)

    with ThreadPoolExecutor(max_workers=2) as pool:
        completion = pool.submit(manager.complete, record, state=RunState.FAILED)
        try:
            assert writing.wait(10)
            deletion = pool.submit(delete)
            assert deleting.wait(10)
            assert not deletion.done()
        finally:
            release.set()
        assert completion.result(timeout=10).state is RunState.FAILED
        deletion.result(timeout=10)

    assert all(item.project_id != child.id for item in manager.store.runs())
    with pytest.raises(KeyError, match="no longer exists"):
        manager.store.save(record)


@pytest.mark.parametrize("operation", ("save", "update"))
def test_stale_writer_waiting_for_delete_cannot_restore_run(
    tmp_path, monkeypatch, operation,
):
    """Independent RunStores must agree on the same cross-instance guard."""
    projects = ProjectRepository(tmp_path / "projects")
    root = projects.create("Root")
    child = projects.create("Child", parent_id=root.id)
    manager = ExecutionManager(project_repository=projects, project=child)
    other = ExecutionManager(project_repository=projects, project=child)
    assert other.store is not manager.store
    record = manager.new_record(task="DOS", directory=tmp_path / "work")
    manager.store.save(record)
    terminal = manager.complete(record, state=RunState.FAILED)
    deleting, release, writing = Event(), Event(), Event()
    original_projects = projects._read_projects

    def pause_delete():
        # ProjectRepository.delete calls this only after acquiring its run guard.
        deleting.set()
        assert release.wait(10)
        return original_projects()

    monkeypatch.setattr(projects, "_read_projects", pause_delete)

    def stale_write():
        writing.set()
        if operation == "save":
            return other.store.save(terminal)
        return other.store.update(terminal, lambda latest: latest)

    with ThreadPoolExecutor(max_workers=2) as pool:
        deletion = pool.submit(projects.delete, child.id, expected_generation=child.generation)
        try:
            assert deleting.wait(10)
            writer = pool.submit(stale_write)
            assert writing.wait(10)
            assert not writer.done()
        finally:
            release.set()
        deletion.result(timeout=10)
        with pytest.raises(KeyError, match="no longer exists"):
            writer.result(timeout=10)

    assert manager.store.runs() == ()
    assert not (projects.runs_directory / manager.store._bucket(terminal.created_at)
                / f"{terminal.id}.json").exists()

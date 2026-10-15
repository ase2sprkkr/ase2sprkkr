"""Tests for persistent ASE2SPRKKR projects."""

from __future__ import annotations

import json
import multiprocessing
import time

import pytest

from ..projects import (
    Project,
    ProjectConflictError,
    ProjectDeleteError,
    ProjectRepository,
    ProjectRepositoryError,
)

ROOT_A = "11111111111111111111111111111111"
ROOT_B = "22222222222222222222222222222222"
CHILD = "33333333333333333333333333333333"


def _write_project(directory, project_id, **values):
    data = {
        "schema_version": 1,
        "id": project_id,
        "generation": 1,
        "name": project_id,
        "parent_id": None,
        "persistent_refs": [],
    }
    data.update(values)
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"{project_id}.json").write_text(
        json.dumps(data), encoding="utf-8"
    )


def _save_while_locked(directory, project_id, started, finished):
    repository = ProjectRepository(directory)
    project = repository.project(project_id).rename("Worker")
    started.set()
    repository.save(project, expected_generation=project.generation)
    finished.set()


def test_project_editing_methods_keep_snapshot_generation():
    project = Project(id=ROOT_A, generation=7, name="Original")

    edited = project.rename(" Renamed ").add_run("run-1").move_to(ROOT_B)

    assert edited.generation == 7
    assert edited.name == "Renamed"
    assert edited.parent_id == ROOT_B
    assert edited.persistent_refs == ("run-1",)
    assert edited.remove_run("run-1").persistent_refs == ()
    assert project.name == "Original"


def test_project_rejects_invalid_local_edits():
    project = Project(id=ROOT_A, generation=1, name="Project")

    with pytest.raises(ValueError, match="name"):
        project.rename(" ")
    with pytest.raises(ValueError, match="own parent"):
        project.move_to(project.id)
    with pytest.raises(ValueError, match="reference"):
        project.add_run("")


def test_repository_uses_one_generation_checked_json_per_project(tmp_path):
    repository = ProjectRepository(tmp_path / "projects")
    parent = repository.create("Parent")
    child = repository.create("Child", parent_id=parent.id)

    json_files = sorted(repository.directory.glob("*.json"))
    assert [path.stem for path in json_files] == sorted([parent.id, child.id])
    assert repository.project(child.id).parent_id == parent.id

    edited = parent.rename("Renamed").add_run("run-17")
    saved = repository.save(edited, expected_generation=parent.generation)

    assert saved.generation == 2
    assert saved.name == "Renamed"
    assert saved.persistent_refs == ("run-17",)


def test_stale_save_conflicts_and_explicit_overwrite_replaces_snapshot(tmp_path):
    repository = ProjectRepository(tmp_path / "projects")
    original = repository.create("Original")
    stale = repository.project(original.id).add_run("run-17")
    renamed = repository.save(
        original.rename("Renamed elsewhere"),
        expected_generation=original.generation,
    )

    with pytest.raises(ProjectConflictError) as exc_info:
        repository.save(stale, expected_generation=original.generation)
    assert exc_info.value.current == renamed

    overwritten = repository.save(
        stale,
        expected_generation=original.generation,
        overwrite=True,
    )
    assert overwritten.name == "Original"
    assert overwritten.persistent_refs == ("run-17",)
    assert overwritten.generation == renamed.generation + 1


def test_move_validates_parent_and_rejects_cycles(tmp_path):
    directory = tmp_path / "projects"
    _write_project(directory, ROOT_A, name="A")
    _write_project(directory, ROOT_B, name="B")
    _write_project(directory, CHILD, name="Child", parent_id=ROOT_A)
    repository = ProjectRepository(directory)

    child = repository.project(CHILD)
    moved = repository.save(
        child.move_to(ROOT_B), expected_generation=child.generation
    )
    assert moved.parent_id == ROOT_B

    root = repository.project(ROOT_B)
    with pytest.raises(ProjectRepositoryError, match="cycle"):
        repository.save(
            root.move_to(CHILD), expected_generation=root.generation
        )


def test_delete_rejects_root_and_project_with_children(tmp_path):
    repository = ProjectRepository(tmp_path / "projects")
    root = repository.create("Root")
    parent = repository.create("Parent", parent_id=root.id)
    child = repository.create("Child", parent_id=parent.id)

    with pytest.raises(ProjectDeleteError, match="Root"):
        repository.delete(root.id, expected_generation=root.generation)
    with pytest.raises(ProjectDeleteError, match="child"):
        repository.delete(parent.id, expected_generation=parent.generation)

    repository.delete(child.id, expected_generation=child.generation)
    with pytest.raises(KeyError):
        repository.project(child.id)


def test_repository_publishes_local_changes(tmp_path):
    repository = ProjectRepository(tmp_path / "projects")
    changed = []
    repository.subscribe(changed.append)

    root = repository.create("Root")
    child = repository.create("Child", parent_id=root.id)
    child = repository.save(
        child.rename("Renamed"), expected_generation=child.generation
    )
    repository.delete(child.id, expected_generation=child.generation)

    assert changed == [root.id, child.id, child.id, child.id]


def test_refresh_publishes_external_project_change(tmp_path):
    directory = tmp_path / "projects"
    first = ProjectRepository(directory)
    project = first.create("Original")
    first.projects()
    changed = []
    first.subscribe(changed.append)

    second = ProjectRepository(directory)
    second_snapshot = second.project(project.id)
    second.save(
        second_snapshot.rename("Elsewhere"),
        expected_generation=second_snapshot.generation,
    )

    refreshed = first.refresh()

    assert refreshed[0].name == "Elsewhere"
    assert changed == [project.id]


def test_project_lock_serializes_generation_check_and_write(tmp_path):
    if multiprocessing.get_start_method(allow_none=True) == "spawn":
        pytest.skip("lock timing test relies on inherited process primitives")

    repository = ProjectRepository(tmp_path / "projects")
    project = repository.create("Original")
    started = multiprocessing.Event()
    finished = multiprocessing.Event()

    with repository._project_lock(project.id):
        process = multiprocessing.Process(
            target=_save_while_locked,
            args=(repository.directory, project.id, started, finished),
        )
        process.start()
        assert started.wait(timeout=2)
        time.sleep(0.1)
        assert not finished.is_set()

    process.join(timeout=3)
    assert process.exitcode == 0
    assert finished.is_set()
    assert repository.project(project.id).name == "Worker"

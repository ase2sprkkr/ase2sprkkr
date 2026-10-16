from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    ExecutionManager,
    Executor,
    ProjectRepository,
    RunReference,
    RunState,
    RunStore,
)


class _Runner:
    def __init__(self, callback=None):
        self.callback = callback
        self.state = RunState.NEW
        self.reference = None
        self.result = SimpleNamespace(return_code=None)
        self.exception = None
        self.on_submit = None
        self.submit_calls = 0

    def submit(self):
        self.submit_calls += 1
        if self.on_submit is not None:
            self.on_submit()
        self.reference = RunReference(executor="project-test", native_id="native-1")
        self.state = RunState.QUEUED

    def wait(self):
        return self.result

    def cancel(self):
        self.state = RunState.CANCELLED

    def add_output_listener(self, callback):
        return None

    def add_output_error_listener(self, callback):
        return None


class _Executor(Executor):
    name = "project-test"

    def __init__(self):
        self.runners = []

    def create_run(self, reader, spec, callback=None):
        runner = _Runner(callback=callback)
        self.runners.append(runner)
        return runner


def _spec(tmp_path):
    return SimpleNamespace(directory=str(tmp_path))


def _manager(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    manager = ExecutionManager(
        RunStore(tmp_path / "runs"),
        project_repository=projects,
    )
    return manager, projects


def test_project_membership_is_written_after_reserve_before_backend_submit(tmp_path):
    manager, projects = _manager(tmp_path)
    project = projects.create("Project")
    executor = _Executor()
    run = manager.create_run(
        executor,
        object(),
        _spec(tmp_path),
        task="DOS",
        input_definition="dos",
        project=project.id,
    )
    backend = executor.runners[0]
    seen = []
    backend.on_submit = lambda: seen.append(
        (
            manager.store.get(run.id).state,
            projects.project(project.id).persistent_refs,
        )
    )

    assert projects.project(project.id).persistent_refs == ()

    run.submit()

    assert seen == [(RunState.SUBMITTING, (run.id,))]
    assert projects.project(project.id).persistent_refs == (run.id,)
    stored = manager.store.get(run.id)
    assert stored.task == "DOS"
    assert stored.input_definition == "DOS"


def test_project_assignment_replays_on_latest_snapshot_after_conflict(tmp_path, monkeypatch):
    manager, projects = _manager(tmp_path)
    project = projects.create("Original")
    executor = _Executor()
    run = manager.create_run(
        executor,
        object(),
        _spec(tmp_path),
        task="SCF",
        project=project.id,
    )

    original_save = projects.save
    conflicted = False

    def save_with_external_edit(candidate, *, expected_generation, overwrite=False):
        nonlocal conflicted
        if not conflicted:
            conflicted = True
            external = ProjectRepository(projects.directory)
            current = external.project(project.id)
            external.save(
                current.rename("Changed elsewhere"),
                expected_generation=current.generation,
            )
        return original_save(
            candidate,
            expected_generation=expected_generation,
            overwrite=overwrite,
        )

    monkeypatch.setattr(projects, "save", save_with_external_edit)

    run.submit()

    saved = projects.project(project.id)
    assert saved.name == "Changed elsewhere"
    assert saved.persistent_refs == (run.id,)


def test_missing_project_prevents_backend_submit_and_marks_run_failed(tmp_path):
    manager, projects = _manager(tmp_path)
    parent = projects.create("Parent")
    project = projects.create("Child", parent_id=parent.id)
    executor = _Executor()
    run = manager.create_run(
        executor,
        object(),
        _spec(tmp_path),
        task="DOS",
        project=project.id,
    )
    projects.delete(project.id, expected_generation=project.generation)

    with pytest.raises(KeyError):
        run.submit()

    assert executor.runners[0].submit_calls == 0
    stored = manager.store.get(run.id)
    assert stored.state == RunState.FAILED
    assert "Cannot assign persistent run to project" in stored.message


def test_unknown_project_is_rejected_before_runner_creation(tmp_path):
    manager, _projects = _manager(tmp_path)
    executor = _Executor()

    with pytest.raises((KeyError, ValueError)):
        manager.create_run(
            executor,
            object(),
            _spec(tmp_path),
            task="DOS",
            project="missing",
        )

    assert executor.runners == []

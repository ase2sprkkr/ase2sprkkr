"""Project membership is the run-store namespace, not a reference list."""

from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    ExecutionManager, ExecutionSpec, Executor, ProjectRepository, RunState,
)


class _Runner:
    def __init__(self):
        self.state = RunState.NEW
        self.backend_id = None
        self.result = SimpleNamespace(return_code=0)
        self.exception = None
        self.failure_phase = None
        self.failure_message = None
        self.submit_calls = 0
        self.on_submit = None

    def submit(self):
        self.submit_calls += 1
        if self.on_submit:
            self.on_submit()
        self.backend_id = "native-1"
        self.state = RunState.FINISHED

    def wait(self):
        return self.result

    def cancel(self):
        self.state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass


class _Executor(Executor):
    name = "project-test"

    def create_project_run(self, spec, *, run_id, callback=None):
        self.runner = _Runner()
        return self.runner


def _fixture(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Project")
    manager = ExecutionManager(project_repository=projects).for_project(project)
    return manager, projects, project


def test_run_is_durable_in_project_namespace_before_backend_submit(tmp_path):
    manager, projects, project = _fixture(tmp_path)
    executor = _Executor()
    run = manager.create_run(
        executor, ExecutionSpec(argv=(), input_file=None, output_file=None, directory=str(tmp_path)),
        process_result=lambda outcome: outcome, task="DOS",
        input_definition="dos",
    )
    observed = []
    executor.runner.on_submit = lambda: observed.append(
        manager.store.get(run.id).state
    )

    run.submit()

    assert observed == [RunState.SUBMITTING]
    stored = manager.store.get(run.id)
    assert stored.input_definition == "DOS"
    assert stored.execution.backend_id == "native-1"
    assert manager.store.directory == tmp_path / "projects_runs"
    assert not hasattr(projects.project(project.id), "persistent_refs")


def test_missing_project_cannot_create_tracked_run(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    manager = ExecutionManager(project_repository=projects)
    executor = _Executor()
    with pytest.raises(ValueError, match="existing project"):
        manager.create_run(
            executor, ExecutionSpec(argv=(), input_file=None, output_file=None, directory=str(tmp_path)),
            process_result=lambda outcome: outcome, task="DOS",
        )
    assert not hasattr(executor, "runner")


def test_invalid_project_is_rejected_before_backend_is_created(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    manager = ExecutionManager(project_repository=projects)
    executor = _Executor()
    with pytest.raises((KeyError, ValueError)):
        manager.create_run(
            executor, ExecutionSpec(argv=(), input_file=None, output_file=None, directory=str(tmp_path)),
            process_result=lambda outcome: outcome, task="DOS",
            project="missing",
        )
    assert not hasattr(executor, "runner")

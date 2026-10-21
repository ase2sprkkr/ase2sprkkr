"""Project membership is the run-store namespace, not a reference list."""

from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import (
    ExecutionManager, Executor, ProjectRepository, RunState,
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
        self.run_id = run_id
        self.runner = _Runner()
        return self.runner

    def create_run(self, spec, callback=None):
        raise AssertionError("project run must use create_project_run")


def _fixture(tmp_path):
    projects = ProjectRepository(tmp_path / "projects")
    project = projects.create("Project")
    manager = ExecutionManager(project_repository=projects).for_project(project)
    return manager, projects, project


def test_run_is_durable_in_project_namespace_before_backend_submit(tmp_path):
    manager, projects, project = _fixture(tmp_path)
    executor = _Executor()
    run = manager.create_run(
        executor, SimpleNamespace(directory=str(tmp_path)),
        process_result=lambda outcome: outcome, task="DOS",
        input_definition="dos",
    )
    observed = []
    executor.runner.on_submit = lambda: observed.append(
        manager.store.get(run.id).state
    )

    run.submit()

    assert executor.run_id == run.id
    assert observed == [RunState.SUBMITTING]
    stored = manager.store.get(run.id)
    assert stored.input_definition == "DOS"
    assert stored.execution.backend_id == "native-1"
    assert manager.store.directory == tmp_path / "projects_runs"
    assert not hasattr(projects.project(project.id), "persistent_refs")


@pytest.mark.parametrize("explicit", (False, True), ids=("missing", "invalid"))
def test_project_must_exist_before_backend_is_created(tmp_path, explicit):
    projects = ProjectRepository(tmp_path / "projects")
    manager = ExecutionManager(project_repository=projects)
    executor = _Executor()
    options = {"project": "missing"} if explicit else {}
    errors = (KeyError, ValueError) if explicit else ValueError
    message = None if explicit else "existing project"

    with pytest.raises(errors, match=message):
        manager.create_run(
            executor, SimpleNamespace(directory=str(tmp_path)),
            process_result=lambda outcome: outcome, task="DOS",
            **options,
        )
    assert not hasattr(executor, "runner")

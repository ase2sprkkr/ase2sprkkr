"""ExecutionManager scoped to a project and its project-owned run store."""

from __future__ import annotations

from .manager import ExecutionManager as _ExecutionManager, RunStore
from .projects import Project, ProjectRepository


class ExecutionManager(_ExecutionManager):
    """Track a run only when it belongs to an existing project.

    Runs are stored under projects_runs/<project-id>/; projects do not maintain
    duplicate run-reference lists.  A manager without a bound project can
    resolve projects but cannot create tracked runs.
    """

    def __init__(
        self, store=None, *, executor_pool=None, submitting_grace=60.0,
        project_repository=None, project=None,
    ):
        self.projects = project_repository if project_repository is not None else ProjectRepository()
        self._project_id = None
        if project is not None:
            identity = project.id if isinstance(project, Project) else str(project)
            current = self.projects.project(identity)
            self._project_id = current.id
            if store is None:
                root = self.projects.directory.parent / "projects_runs"
                store = RunStore(root / current.id)
        super().__init__(
            store, executor_pool=executor_pool,
            submitting_grace=submitting_grace,
        )

    @property
    def project_id(self):
        return self._project_id

    def for_project(self, project):
        identity = project.id if isinstance(project, Project) else str(project)
        current = self.projects.project(identity)
        if self._project_id is not None and current.id != self._project_id:
            raise ValueError("Cannot replace the project of a project-bound execution manager")
        if self._project_id == current.id:
            return self
        return type(self)(
            executor_pool=self.executor_pool,
            submitting_grace=self.submitting_grace,
            project_repository=self.projects,
            project=current,
        )

    def create_run(
        self, executor, spec, *, process_result, task,
        input_definition=None, dataset=None, callback=None,
        run_id=None, project=None,
    ):
        manager = self if project is None else self.for_project(project)
        if manager.project_id is None:
            raise ValueError("A tracked run must belong to an existing project")
        return _ExecutionManager.create_run(
            manager, executor, spec,
            process_result=process_result, task=task,
            input_definition=input_definition, dataset=dataset,
            callback=callback, run_id=run_id,
        )

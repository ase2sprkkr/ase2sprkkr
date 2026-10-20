"""Project-bound facade over the shared global run store."""

from __future__ import annotations

from .manager import ExecutionManager as _ExecutionManager, RunStore
from .projects import Project, ProjectRepository


class ExecutionManager(_ExecutionManager):
    """Track a run only when it belongs to an existing project.

    All projects share the projects_runs/ RunStore. Each RunRecord contains
    project_id, so reservations are atomic across project boundaries.
    A manager without a bound project cannot create tracked runs.
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
                store = RunStore(self.projects.runs_directory)
        if self._project_id is not None:
            if store.directory.resolve() != self.projects.runs_directory.resolve():
                raise ValueError(
                    "Project-bound runs must use the shared projects_runs store"
                )
            if (store.projects_directory is not None
                    and store.projects_directory.resolve() != self.projects.directory.resolve()):
                raise ValueError("RunStore is bound to another project repository")
            store.projects_directory = self.projects.directory
        super().__init__(
            store, executor_pool=executor_pool,
            submitting_grace=submitting_grace,
        )

    @property
    def project_id(self):
        return self._project_id

    def new_record(self, **kwargs):
        if self._project_id is None:
            raise ValueError("Tracked runs require a project")
        return super().new_record(project_id=self._project_id, **kwargs)

    def for_project(self, project):
        identity = project.id if isinstance(project, Project) else str(project)
        current = self.projects.project(identity)
        if self._project_id is not None and current.id != self._project_id:
            raise ValueError("Cannot replace the project of a project-bound execution manager")
        if self._project_id == current.id:
            return self
        return type(self)(
            store=self.store,
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

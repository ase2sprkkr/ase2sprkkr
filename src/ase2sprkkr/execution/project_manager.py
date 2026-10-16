"""Project-aware persistent execution manager."""

from __future__ import annotations

from .base import RunState
from .manager import ExecutionManager as _ExecutionManager
from .projects import Project, ProjectConflictError, ProjectRepository


_PROJECT_ASSIGNMENT_ATTEMPTS = 3


class ExecutionManager(_ExecutionManager):
    """Persist runs and optionally assign them to an ASE2SPRKKR project.

    Project membership is part of creating a persistent run, not GUI
    post-processing. A project-bound manager shares the run store, executor pool
    and project repository with its parent. The run is first reserved as
    ``SUBMITTING`` and then added to the project before backend submission starts.
    """

    def __init__(
        self,
        store=None,
        *,
        executor_pool=None,
        submitting_grace=60.0,
        project_repository=None,
        project=None,
    ):
        super().__init__(
            store,
            executor_pool=executor_pool,
            submitting_grace=submitting_grace,
        )
        self.projects = (
            project_repository
            if project_repository is not None
            else ProjectRepository()
        )
        self._project_id = None if project is None else self._project_identity(project)

    @staticmethod
    def _project_identity(project):
        return project.id if isinstance(project, Project) else str(project)

    @property
    def project_id(self):
        """Project associated with runs created through this manager context."""
        return self._project_id

    def for_project(self, project):
        """Return a manager context assigning newly persisted runs to *project*."""
        current = self.projects.project(self._project_identity(project))
        if self._project_id is not None and current.id != self._project_id:
            raise ValueError(
                "Cannot replace the project of an already project-bound execution manager"
            )
        if self._project_id == current.id:
            return self
        return type(self)(
            self.store,
            executor_pool=self.executor_pool,
            submitting_grace=self.submitting_grace,
            project_repository=self.projects,
            project=current,
        )

    def create_run(
        self,
        executor,
        reader,
        spec,
        *,
        task,
        input_definition=None,
        dataset=None,
        callback=None,
        run_id=None,
        project=None,
    ):
        """Create a persistent runner, optionally assigning it to *project*."""
        manager = self if project is None else self.for_project(project)
        return _ExecutionManager.create_run(
            manager,
            executor,
            reader,
            spec,
            task=task,
            input_definition=input_definition,
            dataset=dataset,
            callback=callback,
            run_id=run_id,
        )

    def _assign_run_to_project(self, run_id):
        if self._project_id is None:
            return

        project = self.projects.project(self._project_id)
        for attempt in range(_PROJECT_ASSIGNMENT_ATTEMPTS):
            if run_id in project.persistent_refs:
                return
            try:
                self.projects.save(
                    project.add_run(run_id),
                    expected_generation=project.generation,
                )
                return
            except ProjectConflictError as conflict:
                project = conflict.current
                if run_id in project.persistent_refs:
                    return
                if attempt + 1 == _PROJECT_ASSIGNMENT_ATTEMPTS:
                    raise

    def reserve(self, record, executor, spec):
        """Reserve file roles and establish project membership before submit.

        The run record and project live in separate atomic files, so there is no
        filesystem transaction spanning both. File-role materialization and
        conflict detection happen in the base manager's atomic submission
        section; project assignment follows before backend submission starts.
        """
        reserved = super().reserve(record, executor, spec)
        if self._project_id is None:
            return reserved
        try:
            self._assign_run_to_project(reserved.id)
        except Exception as exc:
            try:
                super().transition(
                    reserved,
                    state=RunState.FAILED,
                    message=(
                        f"Cannot assign persistent run to project "
                        f"{self._project_id}: {exc}"
                    ),
                )
            except Exception:
                pass
            raise
        return reserved

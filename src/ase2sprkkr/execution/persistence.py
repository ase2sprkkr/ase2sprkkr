"""Resolve explicit project-backed execution without global persistence."""

from __future__ import annotations

from .project_manager import ExecutionManager
from .projects import Project


def resolve_project_manager(project):
    """Only an existing project enables durable run tracking."""
    if project is None or project is False:
        return None
    if isinstance(project, ExecutionManager):
        if project.project_id is None:
            raise ValueError("ExecutionManager must be bound to a project")
        return project
    if isinstance(project, (Project, str)):
        return ExecutionManager(project=project)
    raise TypeError("project must be a Project, project ID, or project-bound ExecutionManager")

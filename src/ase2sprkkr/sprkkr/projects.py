"""SPR-KKR project storage policy, independent of KKRunner's defaults."""

from pathlib import Path

from platformdirs import user_state_dir

from kkrunner import Project, ProjectRepository, ProjectExecutionManager, RunStore
from kkrunner.persistence import resolve_project_manager as _resolve_project_manager


def default_project_repository():
    """Keep existing user project definitions in the ASE2SPRKKR namespace."""
    return ProjectRepository(Path(user_state_dir("ase2sprkkr")) / "projects")


def project_execution_manager(project_repository=None, *, project=None):
    """Retain project definitions while keeping current run records separate.

    Old development RunRecord formats in the ASE2SPRKKR state directory are
    intentionally neither loaded nor migrated. Explicit repositories retain
    their own colocated run store, including isolated application/test stores.
    """
    projects = project_repository if project_repository is not None else default_project_repository()
    legacy = Path(user_state_dir("ase2sprkkr")) / "projects"
    parent = Path(user_state_dir("kkrunner")) if projects.directory == legacy else projects.directory.parent
    return ProjectExecutionManager(
        project_repository=projects, project=project, store=RunStore(parent / "projects_runs")
    )


def resolve_project_manager(project):
    """Apply application storage policy only when constructing a new manager."""
    if isinstance(project, (Project, str)):
        return project_execution_manager(project=project)
    return _resolve_project_manager(project)

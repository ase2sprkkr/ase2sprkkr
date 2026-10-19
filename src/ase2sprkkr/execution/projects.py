"""Persistent projects grouping ASE2SPRKKR calculation runs."""

from __future__ import annotations

import json
import logging
from contextlib import contextmanager
from dataclasses import dataclass, replace
from pathlib import Path
from threading import RLock
from typing import Callable, Iterator
from uuid import UUID, uuid4

from platformdirs import user_state_dir

from ..common.file_locks import FileLock
from ._durability import durable_mkdir, durable_write_json, fsync_directory

logger = logging.getLogger(__name__)

_PROJECT_SCHEMA_VERSION = 2
_HIERARCHY_LOCK = ".hierarchy.lock"
ProjectChangedListener = Callable[[str], None]


class ProjectRepositoryError(RuntimeError):
    """Base error for persistent project metadata."""


class ProjectConflictError(ProjectRepositoryError):
    """A project changed after the caller loaded it."""

    def __init__(
        self,
        project_id: str,
        *,
        expected_generation: int,
        current: "Project",
    ) -> None:
        self.project_id = project_id
        self.expected_generation = expected_generation
        self.current = current
        super().__init__(
            f"Project {project_id!r} changed from generation "
            f"{expected_generation} to {current.generation}."
        )


class ProjectDeleteError(ProjectRepositoryError):
    """A project cannot be deleted while its hierarchy role requires it."""


@dataclass(frozen=True)
class Project:
    """One immutable project snapshot.

    Projects provide the namespace for tracked ASE2SPRKKR runs independently of any user
    interface. Local editing methods return another snapshot at the same
    generation; :class:`ProjectRepository` advances ``generation`` only after a
    successful persistent save.
    """

    id: str
    generation: int
    name: str
    parent_id: str | None = None

    def __post_init__(self) -> None:
        project_id = _project_id(self.id)
        if not isinstance(self.generation, int) or self.generation < 1:
            raise ValueError("Project generation must be a positive integer")
        name = _project_name(self.name)
        parent_id = None if self.parent_id is None else _project_id(self.parent_id)
        if parent_id == project_id:
            raise ValueError("Project cannot be its own parent")
        object.__setattr__(self, "id", project_id)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "parent_id", parent_id)

    def rename(self, name: str) -> "Project":
        """Return this project with a validated display name."""
        normalized = _project_name(name)
        return self if normalized == self.name else replace(self, name=normalized)

    def move_to(self, parent_id: str | None) -> "Project":
        """Return this project with a different direct parent."""
        parent = None if parent_id is None else _project_id(parent_id)
        if parent == self.id:
            raise ValueError("Project cannot be its own parent")
        return self if parent == self.parent_id else replace(self, parent_id=parent)

def default_projects_directory() -> Path:
    """Return the ASE2SPRKKR state directory containing project JSON files."""
    return Path(user_state_dir("ase2sprkkr")) / "projects"


def _project_id(value: str) -> str:
    text = str(value).strip()
    try:
        return UUID(text).hex
    except (ValueError, AttributeError) as exc:
        raise ValueError(f"Invalid project UUID {value!r}") from exc


def _project_name(value: str) -> str:
    name = str(value).strip()
    if not name:
        raise ValueError("Project name must not be empty")
    return name


def validate_hierarchy(projects: tuple[Project, ...]) -> None:
    """Validate parent references and reject hierarchy cycles."""
    by_id: dict[str, Project] = {}
    for project in projects:
        if project.id in by_id:
            raise ProjectRepositoryError(f"Duplicate project id {project.id!r}")
        by_id[project.id] = project

    for project in projects:
        if project.parent_id is not None and project.parent_id not in by_id:
            raise ProjectRepositoryError(
                f"Project {project.id!r} references missing parent "
                f"{project.parent_id!r}"
            )

    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(project_id: str) -> None:
        if project_id in visited:
            return
        if project_id in visiting:
            raise ProjectRepositoryError("Project hierarchy contains a cycle")
        visiting.add(project_id)
        parent_id = by_id[project_id].parent_id
        if parent_id is not None:
            visit(parent_id)
        visiting.remove(project_id)
        visited.add(project_id)

    for project_id in by_id:
        visit(project_id)


class ProjectRepository:
    """Persist and publish generation-checked project JSON documents.

    Files remain authoritative. The in-memory snapshot exists only to report
    changes discovered by :meth:`refresh`; every save still performs its
    generation check against the current on-disk document while holding the
    corresponding stable sidecar lock.
    """

    def __init__(self, directory: str | Path | None = None) -> None:
        self.directory = (
            Path(directory) if directory is not None else default_projects_directory()
        )
        self._state_lock = RLock()
        self._listeners: set[ProjectChangedListener] = set()
        self._known: dict[str, Project] = {}

    def subscribe(self, listener: ProjectChangedListener) -> Callable[[], None]:
        """Subscribe to local or refreshed project changes."""
        with self._state_lock:
            self._listeners.add(listener)

        def unsubscribe() -> None:
            with self._state_lock:
                self._listeners.discard(listener)

        return unsubscribe

    def _notify(self, project_id: str) -> None:
        with self._state_lock:
            listeners = tuple(self._listeners)
        for listener in listeners:
            try:
                listener(project_id)
            except Exception:  # noqa: BLE001 - subscriber isolation boundary
                logger.exception("Project listener failed for %s", project_id)

    def _remember(self, projects: tuple[Project, ...]) -> None:
        with self._state_lock:
            self._known = {project.id: project for project in projects}

    def _remember_one(self, project: Project) -> None:
        with self._state_lock:
            self._known[project.id] = project

    def projects(self) -> tuple[Project, ...]:
        """Load and remember the complete project forest."""
        projects = self._read_projects()
        self._remember(projects)
        return projects

    def project(self, project_id: str) -> Project:
        """Load and remember one project by UUID."""
        canonical = _project_id(project_id)
        project = self._read_path(self._json_path(canonical))
        self._remember_one(project)
        return project

    def refresh(self) -> tuple[Project, ...]:
        """Reload project files and publish ids whose snapshots changed."""
        projects = self._read_projects()
        current = {project.id: project for project in projects}
        with self._state_lock:
            previous = self._known
            self._known = current
        changed = tuple(
            project_id
            for project_id in set(previous) | set(current)
            if previous.get(project_id) != current.get(project_id)
        )
        for project_id in changed:
            self._notify(project_id)
        return projects

    def create(self, name: str, *, parent_id: str | None = None) -> Project:
        """Create a root or child project at generation one."""
        parent = None if parent_id is None else _project_id(parent_id)
        project = Project(
            id=uuid4().hex,
            generation=1,
            name=name,
            parent_id=parent,
        )
        with self._hierarchy_lock():
            projects = self._read_projects()
            if parent is not None and all(item.id != parent for item in projects):
                raise KeyError(parent)
            with self._project_lock(project.id):
                path = self._json_path(project.id)
                if path.exists():
                    raise ProjectRepositoryError(
                        f"Generated project id already exists: {project.id}"
                    )
                self._write(project)
        self._remember_one(project)
        self._notify(project.id)
        return project

    def save(
        self,
        project: Project,
        *,
        expected_generation: int,
        overwrite: bool = False,
    ) -> Project:
        """Persist an edited project snapshot using optimistic concurrency.

        ``overwrite`` means exactly that: the supplied snapshot replaces newer
        persisted values after validating it against the current project forest.
        It is intentionally not an implicit merge operation.
        """
        if not isinstance(project, Project):
            raise TypeError("project must be a Project")

        canonical = _project_id(project.id)
        with self._project_lock(canonical):
            current = self._read_path(self._json_path(canonical))
            self._check_generation(
                current,
                expected_generation=expected_generation,
                overwrite=overwrite,
            )
            if project.parent_id == current.parent_id:
                saved = self._save_locked(project, current)
            else:
                saved = None

        if saved is None:
            # Parent changes affect a cross-project invariant. Re-enter with the
            # hierarchy lock first so create/delete/move all use one lock order.
            with self._hierarchy_lock(), self._project_lock(canonical):
                current = self._read_path(self._json_path(canonical))
                self._check_generation(
                    current,
                    expected_generation=expected_generation,
                    overwrite=overwrite,
                )
                projects = self._read_projects()
                candidate = replace(project, generation=current.generation)
                snapshot = tuple(
                    candidate if item.id == canonical else item for item in projects
                )
                validate_hierarchy(snapshot)
                saved = self._save_locked(project, current)

        self._remember_one(saved)
        self._notify(saved.id)
        return saved

    def delete(self, project_id: str, *, expected_generation: int) -> None:
        """Delete an inactive leaf project and its globally stored run history."""
        canonical = _project_id(project_id)
        from .manager import RunStore, _FINAL_STATES
        run_store = RunStore(
            self.directory.parent / "projects_runs", projects_directory=self.directory,
        )
        # Reserve submissions first; all run writes then use the per-project
        # guard. A callback arriving after deletion cannot recreate history.
        with self._hierarchy_lock(), run_store.submission_lock(), self._file_lock(
            self.directory / f".runs-{canonical}.lock"
        ), self._project_lock(canonical):
            current = self._read_path(self._json_path(canonical))
            self._check_generation(
                current,
                expected_generation=expected_generation,
                overwrite=False,
            )
            if current.parent_id is None:
                raise ProjectDeleteError("Root projects cannot be deleted")
            projects = self._read_projects()
            if any(item.parent_id == canonical for item in projects):
                raise ProjectDeleteError("Projects with child projects cannot be deleted")
            histories = tuple(
                record for record in run_store.runs()
                if record.project_id == canonical
            )
            if any(record.state not in _FINAL_STATES for record in histories):
                raise ProjectDeleteError("Cannot delete a project with active runs")
            try:
                self._json_path(canonical).unlink()
                fsync_directory(self.directory)
                for record in histories:
                    record_path = run_store._path(record)
                    record_path.unlink(missing_ok=True)
                    fsync_directory(record_path.parent)
            except OSError as exc:
                raise ProjectRepositoryError(
                    f"Cannot delete project {canonical!r}: {exc}"
                ) from exc

        with self._state_lock:
            self._known.pop(canonical, None)
        self._notify(canonical)

    def _read_projects(self) -> tuple[Project, ...]:
        try:
            paths = sorted(self.directory.glob("*.json"))
        except OSError as exc:
            raise ProjectRepositoryError(
                f"Cannot list project directory {self.directory}: {exc}"
            ) from exc
        projects = tuple(self._read_path(path) for path in paths)
        validate_hierarchy(projects)
        return projects

    def _save_locked(self, project: Project, current: Project) -> Project:
        updated = replace(project, generation=current.generation + 1)
        self._write(updated)
        return updated

    @staticmethod
    def _check_generation(
        current: Project,
        *,
        expected_generation: int,
        overwrite: bool,
    ) -> None:
        if current.generation != expected_generation and not overwrite:
            raise ProjectConflictError(
                current.id,
                expected_generation=expected_generation,
                current=current,
            )

    def _read_path(self, path: Path) -> Project:
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except FileNotFoundError as exc:
            raise KeyError(path.stem) from exc
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise ProjectRepositoryError(f"Cannot read project {path}: {exc}") from exc
        if not isinstance(data, dict):
            raise ProjectRepositoryError(
                f"Project file {path} must contain a JSON object"
            )
        schema_version = data.pop("schema_version", None)
        if schema_version != _PROJECT_SCHEMA_VERSION:
            raise ProjectRepositoryError(
                f"Unsupported project schema version {schema_version!r} in {path}"
            )
        try:
            project = Project(**data)
        except (TypeError, ValueError) as exc:
            raise ProjectRepositoryError(f"Invalid project file {path}: {exc}") from exc
        if path.stem != project.id:
            raise ProjectRepositoryError(
                f"Project file name {path.stem!r} does not match id {project.id!r}"
            )
        return project

    def _write(self, project: Project) -> None:
        try:
            durable_mkdir(self.directory, parents=True, exist_ok=True, mode=0o700)
            durable_write_json(
                self._json_path(project.id),
                {
                    "schema_version": _PROJECT_SCHEMA_VERSION,
                    "id": project.id,
                    "generation": project.generation,
                    "name": project.name,
                    "parent_id": project.parent_id,
                    },
                mode=0o600,
                indent=2,
            )
        except OSError as exc:
            raise ProjectRepositoryError(
                f"Cannot write project {project.id!r}: {exc}"
            ) from exc

    def _json_path(self, project_id: str) -> Path:
        return self.directory / f"{project_id}.json"

    def _lock_path(self, project_id: str) -> Path:
        return self.directory / f"{project_id}.lock"

    @contextmanager
    def _file_lock(self, path: Path) -> Iterator[None]:
        try:
            with FileLock(path):
                yield
        except OSError as exc:
            raise ProjectRepositoryError(
                f"Cannot lock project repository at {path}: {exc}"
            ) from exc

    def _project_lock(self, project_id: str):
        return self._file_lock(self._lock_path(project_id))

    def _hierarchy_lock(self):
        return self._file_lock(self.directory / _HIERARCHY_LOCK)

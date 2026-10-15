"""Resolution of optional persistent execution bookkeeping."""

from __future__ import annotations

from threading import RLock

from .manager import ExecutionManager as BaseExecutionManager
from .project_manager import ExecutionManager


_DEFAULT_MANAGER = None
_DEFAULT_MANAGER_LOCK = RLock()


def default_execution_manager():
    """Return the process-wide manager used by ``persistent=True``."""
    global _DEFAULT_MANAGER
    with _DEFAULT_MANAGER_LOCK:
        if _DEFAULT_MANAGER is None:
            _DEFAULT_MANAGER = ExecutionManager()
        return _DEFAULT_MANAGER


def resolve_execution_manager(persistent):
    """Resolve the public ``persistent`` option to an execution manager.

    ``False`` and ``None`` disable persistence, ``True`` selects the process-wide
    default manager, and an explicit execution manager is used as-is.
    """
    if persistent is None or persistent is False:
        return None
    if persistent is True:
        return default_execution_manager()
    if isinstance(persistent, BaseExecutionManager):
        return persistent
    raise TypeError(
        "persistent must be False, True, None, or an ExecutionManager instance"
    )

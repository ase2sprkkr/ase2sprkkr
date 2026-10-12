"""Calculation execution subsystem."""

from .base import (
    ExecutionSpec,
    Executor,
    LaunchSpec,
    MpiNotFoundWarning,
    Runner,
    RunReference,
    RunState,
    RunStatus,
    RunTarget,
)
from .local import LocalExecutor, LocalRunner
from .manager import (
    CorruptRunRecord,
    ExecutionManager,
    PersistentRunner,
    RunRecord,
    RunStore,
)
from .pool import ExecutorPool, default_executor_pool
from .psij import PsijExecutor, PsijRunner
from .remote import FileSyncPolicy, RemoteRunner
from .remote_endpoint import RemoteEndpointExecutor
from .tracking import (
    RunConflictError,
    RunRegistry,
    TrackedRunner,
    TrackingExecutor,
    execution_target,
)
from .transport import LocalTransport, RunningCommand, SshTransport, Transport

# Public remote executors represent reusable endpoint identity only. Per-run
# filesystem mapping and synchronization are carried by LaunchSpec.
RemoteExecutor = RemoteEndpointExecutor

__all__ = [
    "ExecutionSpec",
    "Executor",
    "LaunchSpec",
    "MpiNotFoundWarning",
    "Runner",
    "RunReference",
    "RunState",
    "RunStatus",
    "RunTarget",
    "RunRecord",
    "RunStore",
    "CorruptRunRecord",
    "ExecutionManager",
    "PersistentRunner",
    "ExecutorPool",
    "default_executor_pool",
    "RunConflictError",
    "RunRegistry",
    "TrackedRunner",
    "TrackingExecutor",
    "execution_target",
    "LocalExecutor",
    "LocalRunner",
    "PsijExecutor",
    "PsijRunner",
    "FileSyncPolicy",
    "RemoteExecutor",
    "RemoteEndpointExecutor",
    "RemoteRunner",
    "Transport",
    "RunningCommand",
    "LocalTransport",
    "SshTransport",
]

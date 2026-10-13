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
)
from .endpoint import ExecutorEndpoint
from .local import LocalExecutor, LocalRunner
from .manager import (
    CorruptRunRecord,
    DuplicateSubmissionError,
    ExecutionManager,
    FileFingerprint,
    FileRetrievalState,
    PersistentRunner,
    ResultFilesRecord,
    ResultFilesState,
    RunRecord,
    RunStore,
    SubmissionKey,
)
from .pool import ExecutorPool, default_executor_pool
from .psij import PsijExecutor, PsijRunner
from .remote import FileSyncPolicy, RemoteRunner
from .remote_endpoint import RemoteEndpointExecutor
from .transport import LocalTransport, RunningCommand, SshTransport, Transport

# Public remote executors represent reusable endpoint identity only. Per-run
# filesystem mapping and synchronization are carried by LaunchSpec.
RemoteExecutor = RemoteEndpointExecutor

__all__ = [
    "ExecutionSpec",
    "Executor",
    "ExecutorEndpoint",
    "LaunchSpec",
    "MpiNotFoundWarning",
    "Runner",
    "RunReference",
    "RunState",
    "RunStatus",
    "SubmissionKey",
    "DuplicateSubmissionError",
    "FileFingerprint",
    "FileRetrievalState",
    "ResultFilesState",
    "ResultFilesRecord",
    "RunRecord",
    "RunStore",
    "CorruptRunRecord",
    "ExecutionManager",
    "PersistentRunner",
    "ExecutorPool",
    "default_executor_pool",
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
"""Calculation execution subsystem."""

from ..common.file_fingerprints import FileFingerprint
from ..common.file_versions import FileState, FileVersion, FileVersionCaptureError
from .base import (
    ExecutionResources,
    ExecutionSpec,
    Executor,
    LaunchSpec,
    MpiNotFoundWarning,
    MpiSpec,
    ResourceOwnedRunner,
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
    FileRetrievalState,
    PersistentRunner,
    ResultFilesRecord,
    RunRecord,
    RunStore,
    SubmissionKey,
)
from .pool import ExecutorPool, default_executor_pool
from .projects import (
    Project,
    ProjectConflictError,
    ProjectDeleteError,
    ProjectRepository,
    ProjectRepositoryError,
    default_projects_directory,
)
from .project_manager import ExecutionManager
from .psij import PsijExecutor, PsijRunner
from .remote import FileSyncPolicy, RemoteRunner
from .remote_endpoint import RemoteEndpointExecutor
from .transport import LocalTransport, RunningCommand, SshTransport, Transport

# Public remote executors represent reusable endpoint identity only. Per-run
# filesystem mapping and synchronization are carried by LaunchSpec.
RemoteExecutor = RemoteEndpointExecutor

__all__ = [
    "ExecutionResources",
    "ExecutionSpec",
    "Executor",
    "ExecutorEndpoint",
    "LaunchSpec",
    "MpiSpec",
    "MpiNotFoundWarning",
    "ResourceOwnedRunner",
    "Runner",
    "RunReference",
    "RunState",
    "RunStatus",
    "SubmissionKey",
    "DuplicateSubmissionError",
    "FileFingerprint",
    "FileState",
    "FileVersion",
    "FileVersionCaptureError",
    "FileRetrievalState",
    "ResultFilesRecord",
    "RunRecord",
    "RunStore",
    "CorruptRunRecord",
    "ExecutionManager",
    "PersistentRunner",
    "Project",
    "ProjectRepository",
    "ProjectRepositoryError",
    "ProjectConflictError",
    "ProjectDeleteError",
    "default_projects_directory",
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
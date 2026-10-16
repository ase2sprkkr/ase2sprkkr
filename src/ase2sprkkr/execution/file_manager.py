"""File-manifest policy for persistent execution.

The lifecycle implementation in :mod:`ase2sprkkr.execution.manager` is kept
independent of SPR-KKR task identities. This layer gives persistent runs a
concrete input/output/cache manifest and uses only that manifest for concurrent
submission protection.
"""

from __future__ import annotations

import fnmatch
import glob
import json
import os
import posixpath
from dataclasses import dataclass, replace
from uuid import uuid4

from ..common.file_versions import FileState, FileVersion
from .base import RunState
from .endpoint import ExecutorEndpoint, endpoint_for_executor
from .manager import (
    CorruptRunRecord,
    DuplicateSubmissionError,
    ExecutionManager as _LifecycleExecutionManager,
    FileRetrievalState,
    PersistentRunner,
    ResultFilesRecord as _LegacyResultFilesRecord,
    RunRecord,
    RunStore as _LegacyRunStore,
    _now,
    _unique,
)


@dataclass(frozen=True)
class ResultFilesRecord(_LegacyResultFilesRecord):
    """Concrete run file manifest plus immutable historical snapshots.

    ``input_files`` and ``output_files`` are the access roles chosen for this
    concrete run. ``cache_files`` is orthogonal provenance: every cache belongs
    to exactly one of the two access roles after reservation, but remains tagged
    so later cache reuse does not make an old result look overwritten.
    """

    cache_files: tuple[str, ...] = ()

    def __post_init__(self):
        super().__post_init__()
        object.__setattr__(
            self,
            "cache_files",
            tuple(os.fspath(value) for value in self.cache_files),
        )

    @property
    def materialized(self):
        """Whether all cache declarations have a concrete input/output role."""
        roles = set(self.input_files) | set(self.output_files)
        return all(cache in roles for cache in self.cache_files)

    @property
    def historical_input_files(self):
        caches = set(self.cache_files)
        return tuple(path for path in self.input_files if path not in caches)

    @property
    def historical_output_files(self):
        caches = set(self.cache_files)
        return tuple(path for path in self.output_files if path not in caches)

    @property
    def manifests(self):
        """Files belonging to the immutable historical result."""
        return _unique((*self.historical_input_files, *self.historical_output_files))


@dataclass(frozen=True)
class SubmissionKey:
    """Backend location needed for persistent submission recovery.

    Collision identity deliberately does not live here. Concurrent-run
    conflicts are derived from :class:`ResultFilesRecord` input/output roles.
    """

    endpoint: ExecutorEndpoint
    directory: str

    def __post_init__(self):
        if isinstance(self.endpoint, dict):
            object.__setattr__(self, "endpoint", ExecutorEndpoint(**self.endpoint))
        # Backend paths may intentionally contain ``~`` (notably SSH before the
        # target HOME has been resolved), so never reinterpret them locally.
        object.__setattr__(self, "directory", os.fspath(self.directory))


class RunStore(_LegacyRunStore):
    """Run store decoding both legacy and file-manifest version-1 records."""

    @staticmethod
    def _decode(path):
        try:
            with path.open() as file:
                data = json.load(file)
            version = data.get("schema_version")
            if version != 1:
                raise ValueError(
                    f"unsupported schema version {version!r}; expected 1"
                )
            data.pop("requested_task", None)

            submission = data.get("submission")
            if isinstance(submission, dict):
                # Version-1 records may contain the retired task/dataset fields.
                data["submission"] = SubmissionKey(
                    endpoint=submission["endpoint"],
                    directory=submission["directory"],
                )
            result_files = data.get("result_files")
            if isinstance(result_files, dict):
                data["result_files"] = ResultFilesRecord(**result_files)
            return RunRecord(**data)
        except CorruptRunRecord:
            raise
        except (json.JSONDecodeError, TypeError, ValueError, KeyError) as exc:
            raise CorruptRunRecord(path, str(exc)) from exc


def _path_conflicts(first, second):
    """Conservatively decide whether two normalized file paths may overlap."""
    first = os.fspath(first)
    second = os.fspath(second)
    first_magic = glob.has_magic(first)
    second_magic = glob.has_magic(second)
    if not first_magic and not second_magic:
        return first == second
    if first_magic and not second_magic:
        return fnmatch.fnmatchcase(second, first)
    if second_magic and not first_magic:
        return fnmatch.fnmatchcase(first, second)
    if first == second:
        return True

    # Exact glob intersection is not generally cheap. Different fixed prefixes
    # cannot overlap; otherwise choose the safe answer and serialize the runs.
    def fixed_prefix(value):
        positions = [
            position
            for token in ("*", "?", "[")
            if (position := value.find(token)) >= 0
        ]
        return value[: min(positions)] if positions else value

    a = fixed_prefix(first)
    b = fixed_prefix(second)
    return a.startswith(b) or b.startswith(a)


def _sets_conflict(writes, reads_or_writes):
    return any(
        _path_conflicts(write, other)
        for write in writes
        for other in reads_or_writes
    )


class ExecutionManager(_LifecycleExecutionManager):
    """Persistent lifecycle with file-derived concurrent-run protection."""

    def __init__(self, store=None, *, executor_pool=None, submitting_grace=60.0):
        super().__init__(
            RunStore() if store is None else store,
            executor_pool=executor_pool,
            submitting_grace=submitting_grace,
        )

    @staticmethod
    def _result_files_for_spec(spec):
        directory = os.path.abspath(os.fspath(spec.directory))

        def path(value):
            if value is None:
                return None
            value = getattr(value, "name", value)
            try:
                value = os.fspath(value)
            except TypeError:
                return None
            if isinstance(value, int):
                return None
            if not os.path.isabs(value):
                value = os.path.join(directory, value)
            return os.path.abspath(value)

        input_file = path(getattr(spec, "input_file", None))
        output_file = path(getattr(spec, "output_file", None))
        input_files = _unique(
            (input_file, *(path(value) for value in getattr(spec, "input_files", ())))
        )
        output_files = _unique(
            (output_file, *(path(value) for value in getattr(spec, "output_files", ())))
        )
        cache_files = _unique(
            path(value) for value in getattr(spec, "cache_files", ())
        )
        return ResultFilesRecord(
            input_file=input_file,
            output_file=output_file,
            input_files=input_files,
            output_files=output_files,
            cache_files=cache_files,
        )

    def new_record(
        self,
        *,
        task,
        directory,
        dataset=None,
        endpoint=None,
        submission=None,
        result_files=None,
        input_definition=None,
        run_id=None,
    ):
        now = _now()
        return RunRecord(
            id=run_id or uuid4().hex,
            state=RunState.NEW,
            created_at=now,
            updated_at=now,
            finished_at=None,
            task=str(task).upper(),
            dataset=None if dataset is None else str(dataset),
            directory=os.path.abspath(os.fspath(directory)),
            endpoint=endpoint,
            submission=submission,
            result_files=result_files,
            input_definition=input_definition,
        )

    @staticmethod
    def _merged_result_files(latest, requested):
        from .manager import _UNCHANGED

        if requested is _UNCHANGED or requested is None:
            return latest.result_files
        if isinstance(requested, dict):
            requested = ResultFilesRecord(**requested)
        elif isinstance(requested, _LegacyResultFilesRecord) and not isinstance(
            requested, ResultFilesRecord
        ):
            requested = ResultFilesRecord(
                input_file=requested.input_file,
                output_file=requested.output_file,
                input_files=requested.input_files,
                output_files=requested.output_files,
                input_versions=requested.input_versions,
                output_versions=requested.output_versions,
                captured_at=requested.captured_at,
                retrieval_state=requested.retrieval_state,
                retrieval_message=requested.retrieval_message,
            )
        current = latest.result_files
        if current is None:
            return requested
        if current.captured_at is not None and requested.captured_at is not None:
            return current
        return requested

    @staticmethod
    def _backend_path(record, path):
        """Map a local run path into the persisted submission namespace."""
        path = os.path.abspath(os.fspath(path))
        directory = os.path.abspath(os.fspath(record.directory))
        if record.submission is None:
            return os.path.normcase(path)
        try:
            relative = os.path.relpath(path, directory)
        except ValueError:
            return os.path.normcase(path)

        backend_directory = os.fspath(record.submission.directory)
        if record.endpoint is not None and record.endpoint.host is not None:
            parts = () if relative == "." else relative.split(os.sep)
            return posixpath.normpath(posixpath.join(backend_directory, *parts))
        return os.path.normcase(os.path.normpath(os.path.join(backend_directory, relative)))

    def _cache_exists(self, record, cache):
        """Test whether a reusable exact cache is available to this run."""
        if glob.has_magic(cache):
            # A pattern can be extended by the run even when it already has
            # matches, so it is conservatively a writer.
            return False
        if record.endpoint is None or record.endpoint.host is None:
            return os.path.isfile(cache)

        # For a remote run, the backend filesystem is authoritative. A local
        # cache that has not been staged yet is itself a future backend write.
        executor = self.executor_for(record)
        transport = getattr(executor, "transport", None)
        resolve = getattr(executor, "_resolve_submission_directory", None)
        if transport is None or resolve is None:
            return False
        from .remote import _run_sync

        async def probe():
            path = await resolve(self._backend_path(record, cache))
            return await transport.file_signature(path) is not None

        return _run_sync(probe())

    def _materialize_cache_roles(self, record):
        result_files = record.result_files
        if result_files is None or result_files.materialized:
            return result_files
        inputs = list(result_files.input_files)
        outputs = list(result_files.output_files)
        for cache in result_files.cache_files:
            if cache in inputs or cache in outputs:
                continue
            if self._cache_exists(record, cache):
                inputs.append(cache)
            else:
                outputs.append(cache)
        return replace(
            result_files,
            input_files=_unique(inputs),
            output_files=_unique(outputs),
        )

    def _records_conflict(self, first, second):
        if first.endpoint != second.endpoint:
            return False
        a = first.result_files
        b = second.result_files
        if a is None or b is None:
            return False
        a_input = tuple(self._backend_path(first, path) for path in a.input_files)
        a_output = tuple(self._backend_path(first, path) for path in a.output_files)
        b_input = tuple(self._backend_path(second, path) for path in b.input_files)
        b_output = tuple(self._backend_path(second, path) for path in b.output_files)
        return _sets_conflict(a_output, (*b_input, *b_output)) or _sets_conflict(
            b_output, a_input
        )

    def reserve(self, record):
        """Atomically classify caches, check file locks and persist SUBMITTING."""
        if record.submission is None:
            raise ValueError("Persistent run lacks submission metadata")

        attempted_recovery = set()
        while True:
            with self.store.submission_lock():
                # Existence is intentionally sampled only inside the same
                # critical section that reserves the chosen access roles.
                record = replace(
                    record,
                    result_files=self._materialize_cache_roles(record),
                )
                conflicts = tuple(
                    active
                    for active in self.store.runs(active=True)
                    if active.id != record.id and self._records_conflict(record, active)
                )
                if not conflicts:
                    return self.transition(
                        record,
                        state=RunState.SUBMITTING,
                        reference=None,
                        result_files=record.result_files,
                    )
                stale = tuple(
                    conflict
                    for conflict in conflicts
                    if conflict.id not in attempted_recovery
                    and self._is_stale_submitting(conflict)
                )

            if not stale:
                raise DuplicateSubmissionError(record, conflicts)
            progressed = False
            for conflict in stale:
                attempted_recovery.add(conflict.id)
                recovered = self.recover_submitting(conflict)
                if recovered.state != RunState.SUBMITTING or recovered.reference is not None:
                    progressed = True
            if not progressed:
                raise DuplicateSubmissionError(record, conflicts)

    def _capture_result_files(self, result_files):
        if result_files is None or result_files.captured_at is not None:
            return result_files
        return replace(
            result_files,
            input_versions=self._result_files_snapshot(
                result_files.historical_input_files
            ),
            output_versions=self._result_files_snapshot(
                result_files.historical_output_files
            ),
            captured_at=_now(),
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
    ):
        endpoint = endpoint_for_executor(executor)
        submission = SubmissionKey(
            endpoint=endpoint,
            directory=executor.submission_directory(spec),
        )
        record = self.new_record(
            task=task,
            input_definition=input_definition,
            dataset=dataset,
            directory=spec.directory,
            endpoint=endpoint,
            submission=submission,
            result_files=self._result_files_for_spec(spec),
            run_id=run_id,
        )
        self.register_executor(executor, endpoint)
        return self._wrap_runner(
            executor,
            record,
            lambda completed: executor.create_persistent_run(
                reader,
                spec,
                run_id=record.id,
                callback=completed,
            ),
            callback,
        )


__all__ = [
    "ExecutionManager",
    "ResultFilesRecord",
    "RunRecord",
    "RunStore",
    "SubmissionKey",
    "PersistentRunner",
    "DuplicateSubmissionError",
    "CorruptRunRecord",
    "FileRetrievalState",
    "FileState",
    "FileVersion",
]

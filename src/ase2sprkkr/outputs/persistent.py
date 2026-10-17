"""Rehydrate and check SPR-KKR results from persistent execution records."""

from __future__ import annotations

import pyparsing as pp

from ..common.module_registry import RegistryEntryNotFound
from ..common.process_output_reader import OutputParseError
from ..common.warnings import DataValidityError
from ..execution import FileState, RunState
from ..input_parameters.input_parameters import InputParameters
from ..input_parameters.task_definition import TaskMismatchError
from .reader_registry import reader_registry


class PersistentResultUnavailableError(RuntimeError):
    """A stored run cannot currently be reconstructed as its historical TaskResult."""

    def __init__(self, record, message):
        self.record = record
        super().__init__(f"Run {record.id}: {message}")


class PersistentResultMismatchError(PersistentResultUnavailableError):
    """Persisted result data cannot be reconstructed consistently by current code."""


class UnrecoverablePersistentResultError(PersistentResultMismatchError):
    """Persisted metadata cannot be reconstructed by the current code version."""


def _record(manager, run):
    return manager.store.get(run) if isinstance(run, str) else manager.store.get(run.id)


def _task_definition_for_record(record):
    """Return the task definition persisted for this run."""
    try:
        return InputParameters.task_registry[record.task]
    except RegistryEntryNotFound:
        raise UnrecoverablePersistentResultError(
            record, f"unknown persisted task {record.task!r}"
        ) from None


def _input_parameters_for_record(manager, record, task_definition):
    """Reconstruct parameters through the common strict task-read contract."""
    result_files = record.result_files
    input_file = None if result_files is None else result_files.input_file
    if input_file is None or not manager.input_file_is_current(record):
        return None
    if record.input_definition is None:
        return None

    input_registry = InputParameters.input_parameters_registry
    try:
        persisted_definition = input_registry[record.input_definition]
    except RegistryEntryNotFound:
        raise UnrecoverablePersistentResultError(
            record,
            f"unknown persisted input definition {record.input_definition!r}",
        ) from None

    try:
        return input_registry.read_for_task(
            task_definition,
            input_file,
            input_definition=persisted_definition,
        )
    except TaskMismatchError as exc:
        raise UnrecoverablePersistentResultError(
            record,
            "persisted task does not match historical input "
            f"({record.task!r} != {exc.task_name!r})",
        ) from exc
    except (pp.ParseBaseException, DataValidityError) as exc:
        raise PersistentResultMismatchError(
            record,
            f"persisted input cannot be reconstructed: {exc}",
        ) from exc


def _reader_for_record(manager, record):
    """Create the reader declared by the persisted task."""
    task_definition = _task_definition_for_record(record)
    parameters = _input_parameters_for_record(manager, record, task_definition)
    reader_class = reader_registry.class_for_task(task_definition)

    return reader_class(
        parameters,
        None,
        record.directory,
        task_name=record.task,
    )


def attach_persistent_run(
    manager,
    run,
    *,
    reader=None,
    executor=None,
    callback=None,
    output_listener=None,
):
    """Attach to one active stored run without reconstructing its result.

    Reattachment belongs to execution/file finalization. When no live task
    reader is supplied the backend uses a parser-free capture adapter; the
    historical task/input reader is reconstructed only after the run reaches
    AVAILABLE.
    """
    record = _record(manager, run)
    if record.state == RunState.FINALIZATION_FAILED:
        record = manager.retry_finalization(record)
    if record.state in (
        RunState.AVAILABLE,
        RunState.CHECKED,
        RunState.FINISHED,
        RunState.EXECUTION_FAILED,
        RunState.FAILED,
        RunState.CANCELLED,
    ):
        raise RuntimeError(
            f"Run {record.id} is {record.state.value}; use result loading/checking, not attach"
        )

    return manager.attach(
        record,
        reader=reader,
        executor=executor,
        callback=callback,
        output_listener=output_listener,
    )


def _load_result(manager, record, reader=None):
    """Parse one locally available persistent result without changing run state."""
    record = manager.retrieve_result_files(record)
    record = manager.check_result_files(record)
    result_files = record.result_files
    if result_files is None or result_files.output_file is None:
        raise PersistentResultUnavailableError(
            record, "no persisted output result file"
        )

    file_state = result_files.state
    if file_state is not FileState.CURRENT:
        if file_state is FileState.OVERWRITTEN:
            message = "result files have been overwritten"
        elif file_state is FileState.MISSING:
            message = "result files are missing"
        else:
            message = "result files are unavailable"
        raise PersistentResultUnavailableError(record, message)

    output_file = result_files.output_file
    if not manager.output_file_is_current(record):
        raise PersistentResultUnavailableError(
            record, "output result file is unavailable"
        )

    reader_supplied = reader is not None
    if reader is None:
        reader = _reader_for_record(manager, record)

    return_code = record.return_code
    if return_code is None:
        return_code = 0

    input_file = None
    if manager.input_file_is_current(record) and (
        reader_supplied or getattr(reader, "input_parameters", None) is not None
    ):
        input_file = result_files.input_file

    try:
        result = reader.read_from_file(
            output_file,
            return_code=return_code,
            input_file=input_file,
            directory=record.directory,
        )
    except OutputParseError as exc:
        raise PersistentResultMismatchError(
            record, f"persisted output cannot be parsed: {exc}"
        ) from exc
    return record, result


def check_persistent_result(manager, run, *, reader=None):
    """Check/parse an AVAILABLE result and persist the outcome.

    Known reconstruction and availability failures are retryable: the run stays
    AVAILABLE and ``failure_reason`` records the last failed check. Unexpected
    exceptions make the run FAILED. A successful check advances it to CHECKED
    and clears the previous failure reason.
    """
    record = _record(manager, run)
    if record.state == RunState.FINISHED:
        record = manager.mark_available(record)
    if record.state not in (RunState.AVAILABLE, RunState.CHECKED):
        raise PersistentResultUnavailableError(
            record, f"result cannot be checked while run is {record.state.value}"
        )

    try:
        loaded_record, result = _load_result(manager, record, reader=reader)
    except PersistentResultUnavailableError as exc:
        manager.record_result_check_failure(record, exc, retryable=True)
        raise
    except Exception as exc:
        manager.record_result_check_failure(record, exc, retryable=False)
        raise

    checked = manager.mark_checked(loaded_record)
    result.run_record = checked
    return result


def load_persistent_result(manager, run, *, reader=None):
    """Load a persistent TaskResult, checking it on demand when necessary.

    The first successful result retrieval is checked automatically by the live
    persistent runner. Calling this function is an explicit/on-demand check and
    therefore also retries an AVAILABLE result whose previous check failed.
    """
    record = _record(manager, run)
    if record.state == RunState.FINALIZATION_FAILED:
        raise PersistentResultUnavailableError(
            record, "result finalization failed; retry finalization first"
        )
    if record.state in (
        RunState.EXECUTION_FAILED,
        RunState.FAILED,
        RunState.CANCELLED,
    ):
        reason = f": {record.failure_reason}" if record.failure_reason else ""
        raise PersistentResultUnavailableError(
            record,
            f"run ended as {record.state.value}; result checking is disabled{reason}",
        )

    return check_persistent_result(manager, record, reader=reader)

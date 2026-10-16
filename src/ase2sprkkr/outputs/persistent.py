"""Rehydrate SPR-KKR results from persistent execution records."""

from __future__ import annotations

import pyparsing as pp

from ..common.module_registry import RegistryEntryNotFound
from ..common.warnings import DataValidityError
from ..execution import FileState, RunState
from ..input_parameters.input_parameters import InputParameters
from ..input_parameters.task_definition import TaskMismatchError
from .reader_registry import reader_registry


class PersistentResultUnavailableError(RuntimeError):
    """A stored run cannot be reconstructed as its historical TaskResult."""

    def __init__(self, record, message):
        self.record = record
        super().__init__(f"Run {record.id}: {message}")


class PersistentResultMismatchError(PersistentResultUnavailableError):
    """Persisted result data cannot be reconstructed consistently by current code."""


class UnrecoverablePersistentResultError(PersistentResultMismatchError):
    """Retrying reconstruction with unchanged persisted data cannot succeed."""


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
    """Attach to or explicitly retry result processing of one stored run."""
    record = _record(manager, run)
    if record.state in (RunState.FINALIZATION_FAILED, RunState.RESULT_MISMATCH):
        record = manager.retry_finalization(record)

    if reader is None:
        try:
            reader = _reader_for_record(manager, record)
        except PersistentResultMismatchError as exc:
            if record.state == RunState.FINALIZING:
                manager.mark_result_mismatch(record, exc)
            raise

    # Attach/transport failures are deliberately not classified as result
    # mismatches. They may be transient and leave FINALIZING retryable.
    return manager.attach(
        record,
        reader=reader,
        executor=executor,
        callback=callback,
        output_listener=output_listener,
    )


def load_persistent_result(manager, run, *, reader=None):
    """Reconstruct the historical TaskResult of one terminal persistent run.

    Reader selection is determined by the persisted task. When the historical
    input file is still current and its exact input definition was persisted,
    parameters are reconstructed from that definition.
    """
    record = _record(manager, run)
    if record.state == RunState.FINALIZATION_FAILED:
        raise PersistentResultUnavailableError(
            record, "result finalization failed; retry result processing first"
        )
    if record.state == RunState.RESULT_MISMATCH:
        raise PersistentResultUnavailableError(
            record, "persisted result mismatch; retry result processing first"
        )

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
    output_available = manager.output_file_is_current(record)
    if not output_available and record.state == RunState.FINISHED:
        raise PersistentResultUnavailableError(
            record, "output result file is unavailable"
        )

    reader_supplied = reader is not None
    if reader is None:
        reader = _reader_for_record(manager, record)

    return_code = record.return_code
    if return_code is None:
        return_code = (
            0
            if record.state in (RunState.FINALIZING, RunState.FINISHED)
            else 1
        )

    input_file = None
    if manager.input_file_is_current(record) and (
        reader_supplied or getattr(reader, "input_parameters", None) is not None
    ):
        input_file = result_files.input_file

    if output_available:
        result = reader.read_from_file(
            output_file,
            return_code=return_code,
            input_file=input_file,
            directory=record.directory,
        )
    else:
        result = reader.create_result(output_file, input_file)
        result.complete(None, return_code)

    result.run_record = record
    return result
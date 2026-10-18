"""Rehydrate SPR-KKR results from persistent execution records."""

from __future__ import annotations

from ..common.module_registry import RegistryEntryNotFound
from ..execution import FailurePhase, FileState, RunState
from ..input_parameters.input_parameters import InputParameters
from ..input_parameters.task_definition import TaskMismatchError
from .reader_registry import reader_registry


class ProjectResultUnavailableError(RuntimeError):
    """A stored run cannot be reconstructed as its historical TaskResult."""

    def __init__(self, record, message):
        self.record = record
        super().__init__(f"Run {record.id}: {message}")


class UnrecoverableProjectResultError(ProjectResultUnavailableError):
    """The persisted metadata is insufficient to reconstruct this result."""


def _record(manager, run):
    return manager.store.get(run) if isinstance(run, str) else manager.store.get(run.id)


def _task_definition_for_record(record):
    try:
        return InputParameters.task_registry[record.task]
    except RegistryEntryNotFound:
        raise UnrecoverableProjectResultError(
            record, f"unknown persisted task {record.task!r}"
        ) from None


def _input_parameters_for_record(manager, record, task_definition):
    result_files = record.files
    input_file = None if result_files is None else result_files.input_file
    if input_file is None or not manager.input_file_is_current(record):
        return None
    if record.input_definition is None:
        return None

    input_registry = InputParameters.input_parameters_registry
    try:
        definition = input_registry[record.input_definition]
    except RegistryEntryNotFound:
        raise UnrecoverableProjectResultError(
            record,
            f"unknown persisted input definition {record.input_definition!r}",
        ) from None

    try:
        return input_registry.read_for_task(
            task_definition,
            input_file,
            input_definition=definition,
        )
    except TaskMismatchError as exc:
        raise UnrecoverableProjectResultError(
            record,
            "persisted task does not match historical input "
            f"({record.task!r} != {exc.task_name!r})",
        ) from exc


def _reader_for_record(manager, record):
    """Create the task-owned output reader declared by the persisted task."""
    task_definition = _task_definition_for_record(record)
    parameters = _input_parameters_for_record(manager, record, task_definition)
    reader_class = reader_registry.class_for_task(task_definition)
    return reader_class(
        parameters,
        None,
        record.execution.directory.local,
        task_name=record.task,
    )


def _result_processor_for_record(manager, record):
    """Build the generic execution→TaskResult boundary in the outputs layer."""
    reader = _reader_for_record(manager, record)

    def process_result(outcome):
        return reader.read_from_file(
            outcome.output_file,
            outcome.error_file,
            return_code=(0 if outcome.return_code is None else outcome.return_code),
            input_file=outcome.input_file,
            directory=record.execution.directory.local,
        )

    return process_result


def attach_project_run(
    manager,
    run,
    *,
    executor=None,
    callback=None,
    output_listener=None,
):
    """Reconnect an active run and resume its client-side lifecycle."""
    record = _record(manager, run)
    if record.state in (RunState.FINISHED, RunState.FAILED, RunState.CANCELLED):
        raise ProjectResultUnavailableError(
            record, f"run is already terminal ({record.state.value})"
        )

    try:
        process_result = _result_processor_for_record(manager, record)
    except UnrecoverableProjectResultError as exc:
        # If execution is already waiting for result reconstruction, failure to
        # construct the task-owned processor is itself a postprocessing failure.
        if record.state == RunState.POSTPROCESSING:
            manager.complete(
                record,
                state=RunState.FAILED,
                message=str(exc),
                failure_phase=FailurePhase.POSTPROCESSING,
                artifacts_available=True,
            )
        raise

    return manager.attach(
        record,
        process_result=process_result,
        executor=executor,
        callback=callback,
        output_listener=output_listener,
    )


def load_project_result(manager, run):
    """Reconstruct the historical TaskResult of one terminal project run."""
    record = _record(manager, run)
    if record.state not in (RunState.FINISHED, RunState.FAILED, RunState.CANCELLED):
        raise ProjectResultUnavailableError(
            record, f"run is still active ({record.state.value})"
        )

    record = manager.retrieve_result_files(record)
    record = manager.check_result_files(record)
    result_files = record.files
    if result_files is None or result_files.output_file is None:
        raise ProjectResultUnavailableError(
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
        raise ProjectResultUnavailableError(record, message)

    output_file = result_files.output_file
    output_available = manager.output_file_is_current(record)
    if not output_available:
        raise ProjectResultUnavailableError(
            record, "output result file is unavailable"
        )

    reader = _reader_for_record(manager, record)
    return_code = record.return_code
    if return_code is None:
        return_code = 0 if record.state == RunState.FINISHED else 1

    input_file = None
    if manager.input_file_is_current(record) and (
        getattr(reader, "input_parameters", None) is not None
    ):
        input_file = result_files.input_file

    result = reader.read_from_file(
        output_file,
        return_code=return_code,
        input_file=input_file,
        directory=record.execution.directory.local,
    )
    result.run_record = record
    return result

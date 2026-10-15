"""Rehydrate SPR-KKR results from persistent execution records."""

from __future__ import annotations

import os

from ..common.module_registry import RegistryEntryNotFound
from ..execution import ResultFilesState, RunState
from ..input_parameters.input_parameters import InputParameters
from .reader_registry import reader_registry
from .readers.default import DefaultOutputReader


class PersistentResultUnavailableError(RuntimeError):
    """A stored run cannot be reconstructed as its historical TaskResult."""

    def __init__(self, record, message):
        self.record = record
        super().__init__(f"Run {record.id}: {message}")


def _record(manager, run):
    return manager.store.get(run) if isinstance(run, str) else manager.store.get(run.id)


def _task_definition_for_record(record):
    """Return the best reconstructable bound task definition.

    Schema-2 records persist the bound task explicitly. Older records only have
    ``task``; that value is still worth trying because historically it was often
    a normal task name or alias. A concrete native submode simply misses the
    task registry and falls through to native reader discovery.
    """
    name = getattr(record, "task_definition", None) or record.task
    try:
        return InputParameters.task_registry[name]
    except RegistryEntryNotFound:
        return None


def _input_parameters_for_record(manager, record, task_definition):
    """Reconstruct input parameters only from verified historical input."""
    result_files = record.result_files
    input_file = None if result_files is None else result_files.input_file
    if input_file is None or not manager.input_file_is_current(record):
        return None
    if task_definition is not None:
        return task_definition.read_input_parameters(input_file)
    return InputParameters.from_file(input_file)


def _reader_for_record(manager, record):
    """Create the best reader reconstructable from persisted task context."""
    task_definition = _task_definition_for_record(record)
    parameters = _input_parameters_for_record(manager, record, task_definition)
    if parameters is not None and parameters.task_definition is not None:
        task_definition = parameters.task_definition

    if task_definition is not None:
        reader_class = task_definition.output_reader_class(
            parameters,
            execution_task=record.task,
        )
    else:
        try:
            reader_class = reader_registry[record.task]
        except RegistryEntryNotFound:
            reader_class = DefaultOutputReader

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
    """Attach to a stored run using the best reconstructable output reader."""
    record = _record(manager, run)
    if reader is None:
        reader = _reader_for_record(manager, record)
    return manager.attach(
        record,
        reader=reader,
        executor=executor,
        callback=callback,
        output_listener=output_listener,
    )


def load_persistent_result(manager, run, *, reader=None):
    """Reconstruct the historical TaskResult of one terminal persistent run.

    Result-file retrieval and freshness remain responsibilities of the execution
    layer. When the historical input file is still current, reader creation
    receives reconstructed task-bound input parameters. If that context is
    unavailable or stale, reconstruction degrades to the persisted
    task-definition/native-task identities instead of making the output result
    unavailable.
    """
    record = manager.retrieve_result_files(run)
    record = manager.check_result_files(record)
    result_files = record.result_files
    if result_files is None or result_files.output_file is None:
        raise PersistentResultUnavailableError(
            record, "no persisted output result file"
        )

    if result_files.state == ResultFilesState.OVERWRITTEN:
        raise PersistentResultUnavailableError(
            record, "result files have been overwritten"
        )

    output_file = result_files.output_file
    output_exists = os.path.isfile(output_file)
    if not output_exists and record.state == RunState.FINISHED:
        raise PersistentResultUnavailableError(
            record, "output result file is unavailable"
        )

    if reader is None:
        reader = _reader_for_record(manager, record)

    return_code = record.return_code
    if return_code is None:
        return_code = 0 if record.state == RunState.FINISHED else 1

    if output_exists:
        result = reader.read_from_file(
            output_file,
            return_code=return_code,
            input_file=(
                result_files.input_file
                if manager.input_file_is_current(record)
                else None
            ),
            directory=record.directory,
        )
    else:
        # Failed/cancelled runs may have no parseable output at all. A shell
        # result still conveys their terminal outcome without pretending that a
        # successful calculation was reconstructed.
        result = reader.create_result(
            output_file,
            result_files.input_file
            if manager.input_file_is_current(record)
            else None,
        )
        result.complete(None, return_code)

    result.run_record = record
    return result

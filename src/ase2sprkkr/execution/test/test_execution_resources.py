import io
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution.base import (
    ExecutionResources,
    ResourceOwnedRunner,
    RunState,
    SubmissionNotStarted,
)
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


class _Parameters:
    _definition = SimpleNamespace(mpi=False, executable="sprkkr")
    CONTROL = SimpleNamespace(
        POTFIL=SimpleNamespace(result=None),
        DATASET=SimpleNamespace(result=None),
    )
    _execution_task_name = "SCF"

    def output_reader(self, calculator, print_output, directory):
        return object()

    def input_files(self, directory):
        return ()

    def output_files(self, directory):
        return ()

    def cache_files(self, directory):
        return ()


class _Calculator:
    persistent = False

    @staticmethod
    def value_or_default(name, value):
        return False if value is None else value


class _Runner:
    result = SimpleNamespace(return_code=None)
    reference = None

    def __init__(self, callback=None, *, submit_error=None):
        self.callback = callback
        self.submit_error = submit_error
        self._state = RunState.NEW

    @property
    def state(self):
        return self._state

    def submit(self):
        if self.submit_error is not None:
            raise self.submit_error
        self._state = RunState.FINISHED
        if self.callback:
            self.callback(self.result)

    def wait(self):
        return self.result

    def cancel(self):
        self._state = RunState.CANCELLED

    def add_output_listener(self, callback):
        pass

    def add_output_error_listener(self, callback):
        pass


class _Executor:
    def __init__(self, *, factory_error=None, submit_error=None):
        self.factory_error = factory_error
        self.submit_error = submit_error

    def create_run(self, reader, spec, callback=None):
        if self.factory_error is not None:
            raise self.factory_error
        return _Runner(callback, submit_error=self.submit_error)


def _run_process(executor, input_file, output_file):
    return InputParameters.run_process(
        _Parameters(),
        _Calculator(),
        input_file,
        output_file,
        profile=SPRKKRProfile(mpi=False),
        executor=executor,
        run_async=True,
    )


def test_files_are_closed_when_runner_creation_fails():
    input_file = io.BytesIO(b"input")
    output_file = io.BytesIO()

    with pytest.raises(RuntimeError, match="factory failed"):
        _run_process(
            _Executor(factory_error=RuntimeError("factory failed")),
            input_file,
            output_file,
        )

    assert input_file.closed
    assert output_file.closed


def test_files_are_closed_when_submission_fails_before_start():
    input_file = io.BytesIO(b"input")
    output_file = io.BytesIO()
    runner = _run_process(
        _Executor(submit_error=SubmissionNotStarted("staging failed")),
        input_file,
        output_file,
    )

    assert isinstance(runner, ResourceOwnedRunner)
    assert not input_file.closed
    assert not output_file.closed

    with pytest.raises(SubmissionNotStarted, match="staging failed"):
        runner.submit()

    assert input_file.closed
    assert output_file.closed


def test_completion_releases_files_without_wait():
    input_file = io.BytesIO(b"input")
    output_file = io.BytesIO()
    runner = _run_process(_Executor(), input_file, output_file)

    runner.submit()

    assert runner.state == RunState.FINISHED
    assert input_file.closed
    assert output_file.closed


def test_unstarted_runner_can_release_owned_files_explicitly():
    input_file = io.BytesIO(b"input")
    output_file = io.BytesIO()
    runner = _run_process(_Executor(), input_file, output_file)

    runner.close()
    runner.close()

    assert input_file.closed
    assert output_file.closed


def test_cancel_started_run_does_not_close_resources_before_monitor_finishes():
    input_file = io.BytesIO(b"input")
    underlying = _Runner()
    underlying._state = RunState.RUNNING
    resources = ExecutionResources(input_file)
    runner = ResourceOwnedRunner(underlying, resources)

    runner.cancel()

    assert runner.state == RunState.CANCELLED
    assert not resources.closed
    assert not input_file.closed


def test_cancel_unstarted_run_releases_resources():
    input_file = io.BytesIO(b"input")
    resources = ExecutionResources(input_file)
    runner = ResourceOwnedRunner(_Runner(), resources)

    runner.cancel()

    assert runner.state == RunState.CANCELLED
    assert resources.closed
    assert input_file.closed

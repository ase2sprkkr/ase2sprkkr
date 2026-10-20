"""Shared monitor waiting contract for local-project and remote runners."""

import asyncio

import pytest

from ..base import ExecutionSpec, RunState
from ..local_project import LocalProjectRunner
from ..remote import RemoteRunner


@pytest.mark.parametrize("runner_class", (LocalProjectRunner, RemoteRunner))
@pytest.mark.parametrize("state, message", (
    (RunState.NEW, "Run has not been submitted"),
    (RunState.RUNNING, "Attached run has not been started"),
    (RunState.UNKNOWN, "Attached run has not been started"),
    (RunState.FINISHED, "Attached run has not been started"),
))
def test_wait_without_started_monitor_fails_immediately(
    tmp_path, runner_class, state, message,
):
    spec = ExecutionSpec(
        argv=(), input_file=None, output_file=None, directory=str(tmp_path),
    )
    runner = runner_class(None, spec, state=state)
    with pytest.raises(RuntimeError, match=message):
        runner.wait()


@pytest.mark.parametrize("runner_class", (LocalProjectRunner, RemoteRunner))
def test_wait_and_run_async_return_completed_monitor_result(tmp_path, runner_class):
    spec = ExecutionSpec(
        argv=(), input_file=None, output_file=None, directory=str(tmp_path),
    )
    runner = runner_class(None, spec)
    runner._thread = object()  # A monitor is already running/attached.
    runner._done.set()
    assert runner.wait() is runner.result
    assert asyncio.run(runner.run_async()) is runner.result

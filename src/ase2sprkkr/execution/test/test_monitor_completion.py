"""Shared monitor completion mechanics, independently of backend state machines."""

from threading import Event
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution.base import finish_monitor, wait_for_monitor


def _runner(callback=None, error=None):
    return SimpleNamespace(
        callback=callback,
        result=object(),
        _exception=error,
        _callback_exception=None,
        _done=Event(),
    )


def test_monitor_callback_receives_result_before_waking_waiters():
    observed = []
    runner = _runner(callback=lambda result: observed.append(result))
    finish_monitor(runner)

    assert observed == [runner.result]
    assert runner._done.is_set()
    assert wait_for_monitor(runner) is runner.result


def test_callback_failure_does_not_block_monitor_completion():
    def failed(_result):
        raise ValueError("callback failed")

    runner = _runner(callback=failed)
    finish_monitor(runner)

    assert runner._done.is_set()
    with pytest.raises(ValueError, match="callback failed"):
        wait_for_monitor(runner)


def test_execution_error_has_priority_over_callback_error():
    def failed(_result):
        raise ValueError("callback failed")

    original_error = RuntimeError("execution failed")
    runner = _runner(callback=failed, error=original_error)
    finish_monitor(runner)

    with pytest.raises(RuntimeError, match="execution failed") as error:
        wait_for_monitor(runner)
    assert error.value is original_error
    assert isinstance(runner._callback_exception, ValueError)

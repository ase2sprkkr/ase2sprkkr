"""Shared monitor completion mechanics, independently of backend state machines."""

from threading import Event, Lock
from types import SimpleNamespace

import pytest

from ase2sprkkr.execution.base import RunState, finish_monitor, start_monitor, wait_for_monitor


def _runner(callback=None, error=None):
    return SimpleNamespace(
        callback=callback,
        result=object(),
        _exception=error,
        _callback_exception=None,
        _done=Event(),
        _thread=object(),  # Model a monitor already started and completed.
        _state=RunState.FINISHED,
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


def test_start_monitor_is_idempotent():
    entered, release = Event(), Event()
    calls = []

    def monitor():
        calls.append(1)
        entered.set()
        assert release.wait(2)

    runner = SimpleNamespace(
        _lock=Lock(), _thread=None, _backend_id="job-1", _thread_main=monitor,
    )
    try:
        start_monitor(runner, missing_backend="No job ID")
        assert entered.wait(2)
        thread = runner._thread
        start_monitor(runner, missing_backend="No job ID")
        assert runner._thread is thread
    finally:
        release.set()
        if runner._thread is not None:
            runner._thread.join(timeout=2)

    assert calls == [1]
    assert not runner._thread.is_alive()


def test_start_monitor_requires_backend_identity():
    runner = SimpleNamespace(
        _lock=Lock(), _thread=None, _backend_id=None, _thread_main=lambda: None,
    )
    with pytest.raises(RuntimeError, match="No backend identity"):
        start_monitor(runner, missing_backend="No backend identity")
    assert runner._thread is None

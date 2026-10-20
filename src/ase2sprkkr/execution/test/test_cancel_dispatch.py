"""Backend cancellation must not block a calling GUI thread."""

from threading import Event, Lock
from types import SimpleNamespace

import pytest

from ..base import dispatch_cancellation, request_cancellation


def test_cancel_dispatch_does_not_wait_for_backend():
    entered, release = Event(), Event()

    def cancel():
        entered.set()
        assert release.wait(2)
        return "accepted"

    future = dispatch_cancellation(cancel)
    assert entered.wait(2)
    assert not future.done()
    release.set()
    assert future.result(timeout=2) == "accepted"


def test_cancel_dispatch_preserves_backend_error():
    def cancel():
        raise OSError("scheduler unavailable")

    with pytest.raises(OSError, match="scheduler unavailable"):
        dispatch_cancellation(cancel).result(timeout=2)


def test_request_cancellation_dispatches_once_and_reuses_future():
    entered, release = Event(), Event()
    calls = []
    runner = SimpleNamespace(_lock=Lock(), _cancel_future=None)

    def cancel():
        calls.append(1)
        entered.set()
        assert release.wait(2)
        return "cancelled"

    try:
        first = request_cancellation(runner, cancel)
        assert entered.wait(2)
        second = request_cancellation(
            runner, lambda: pytest.fail("Cancellation was dispatched twice"),
        )
        assert first is second
        assert not second.done()
    finally:
        release.set()
    assert first.result(timeout=2) == "cancelled"
    assert calls == [1]

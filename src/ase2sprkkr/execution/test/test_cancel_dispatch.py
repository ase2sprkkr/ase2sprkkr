"""Backend cancellation must not block a calling GUI thread."""

from threading import Event

import pytest

from ..base import dispatch_cancellation


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

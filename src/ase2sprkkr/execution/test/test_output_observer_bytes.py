"""One byte-output decoding and observer error policy across execution backends."""

import pytest

from ase2sprkkr.execution.base import OutputObservers
from ase2sprkkr.execution.local import LocalRunner
from ase2sprkkr.execution.local_project import _ObserverSink
from ase2sprkkr.execution.psij import PsijRunner
from ase2sprkkr.execution.remote import RemoteRunner


@pytest.mark.parametrize(
    "data, expected",
    [
        (b"", []),
        (b"ASCII", [("ASCII", "err")]),
        ("Příliš".encode("utf8"), [("Příliš", "err")]),
        (b"broken:\xff", [("broken:\ufffd", "err")]),
    ],
)
def test_emit_bytes_decodes_nonempty_chunks(data, expected):
    received = []
    observers = OutputObservers(lambda text, kind: received.append((text, kind)))
    observers.emit_bytes(data, "err")
    assert received == expected


@pytest.mark.parametrize("runner_cls", (LocalRunner, PsijRunner, RemoteRunner))
def test_runner_output_methods_delegate_to_shared_bytes_policy(runner_cls):
    received = []
    runner = object.__new__(runner_cls)
    runner._output_observers = OutputObservers(
        lambda text, kind: received.append((text, kind))
    )
    runner._emit_output(b"warning:\xff", "err")
    runner._emit_output(b"", "out")
    assert received == [("warning:\ufffd", "err")]


def test_local_project_sink_uses_shared_bytes_policy():
    received = []
    observers = OutputObservers(
        lambda text, kind: received.append((text, kind))
    )
    sink = _ObserverSink(observers, "out")
    sink.write("výstup".encode("utf8"))
    sink.write(b"")
    sink.flush()
    assert received == [("výstup", "out")]


def test_observer_error_does_not_prevent_other_byte_listeners():
    observed, failures = [], []
    observers = OutputObservers()
    def broken(text, kind):
        raise RuntimeError("listener failed")
    observers.add(broken)
    observers.add(lambda text, kind: observed.append((text, kind)))
    observers.add_error_listener(lambda callback, exc: failures.append((callback, str(exc))))
    observers.emit_bytes(b"data", "out")
    observers.emit_bytes(b"next", "out")
    assert observed == [("data", "out"), ("next", "out")]
    assert failures == [(broken, "listener failed")]

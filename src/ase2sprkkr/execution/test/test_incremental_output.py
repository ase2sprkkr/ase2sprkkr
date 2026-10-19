"""Observer decoding across arbitrary binary stream boundaries."""

from ..base import OutputObservers


def test_utf8_chunks_and_independent_channels():
    seen = []
    observers = OutputObservers(lambda text, kind: seen.append((text, kind)))
    observers.emit_bytes(b"\xc5", "out")
    observers.emit_bytes(b"error", "err")
    observers.emit_bytes(b"\xbe", "out")
    observers.finish_bytes()
    assert seen == [("error", "err"), ("ž", "out")]


def test_incomplete_utf8_replaced_at_eof_once():
    seen = []
    observers = OutputObservers(lambda text, kind: seen.append((text, kind)))
    observers.emit_bytes(b"prefix:\xc5", "out")
    assert seen == [("prefix:", "out")]
    observers.finish_bytes()
    observers.finish_bytes()
    assert seen == [("prefix:", "out"), ("\ufffd", "out")]

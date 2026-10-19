import asyncio
import threading

import pytest

from ..async_utils import run_coro_sync
from ..process_output_reader import ProcessOutputParser


class _OutputParser(ProcessOutputParser):
    async def read_output(self, stdout, *args):
        out = b""
        while line := await stdout.readline():
            out += line
        return out

    async def read_error(self, stderr, *args):
        out = b""
        while line := await stderr.readline():
            out += line
        return out


def test_read_output_file(tmp_path):
    output = tmp_path / "stdout"
    error = tmp_path / "stderr"
    output.write_bytes(b"one\ntwo\n")
    error.write_bytes(b"warning\n")

    out, err, return_code = _OutputParser().read_from_file(output, error)

    assert out == b"one\ntwo\n"
    assert err == b"warning\n"
    assert return_code == 0


def test_run_coro_sync_outside_event_loop():
    async def operation():
        return 42

    assert run_coro_sync(operation()) == 42


def test_run_coro_sync_inside_event_loop():
    async def check():
        event_loop_thread = threading.get_ident()

        async def operation():
            return threading.get_ident()

        assert run_coro_sync(operation()) != event_loop_thread

    asyncio.run(check())


@pytest.mark.parametrize("error", (ValueError("failed"), KeyboardInterrupt("interrupted")))
def test_run_coro_sync_propagates_worker_errors(error):
    async def check():
        async def operation():
            raise error

        with pytest.raises(type(error)) as raised:
            run_coro_sync(operation())
        assert raised.value is error

    asyncio.run(check())


def test_process_parser_can_read_from_an_existing_event_loop(tmp_path):
    output = tmp_path / "stdout"
    error = tmp_path / "stderr"
    output.write_bytes(b"out\n")
    error.write_bytes(b"err\n")

    async def check():
        assert _OutputParser().read_from_file(output, error) == (
            b"out\n", b"err\n", 0,
        )

    asyncio.run(check())

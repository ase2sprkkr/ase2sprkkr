import asyncio
import sys
import tempfile
import threading
import time
from types import SimpleNamespace

import pytest

from ...execution import ExecutionSpec, LocalExecutor, LocalRunner, RunState
from ..profile import SPRKKRProfile


class _Result:
    def complete(self, error, return_code):
        self.error = error
        self.return_code = return_code


class _Parser:
    print_output = False

    async def read_output(self, stdout, result):
        out = b""
        while data := await stdout.readline():
            out += data
        result.output = out
        return out

    async def read_error(self, stderr, result):
        out = b""
        while data := await stderr.readline():
            out += data
        result.error_output = out
        return out

    def result(self, output, error, return_code):
        if return_code:
            exc = ValueError(f"Process returned {return_code}")
            exc.output = output
            exc.error = error
            raise exc
        return output, error, return_code


class _Reader:
    print_output = False

    def __init__(self):
        self.parser = _Parser()

    def create_result(self, output_file, input_file):
        return _Result()

    async def parse_streams(self, stdout, stderr, result):
        return await asyncio.gather(
            self.parser.read_output(stdout, result),
            self.parser.read_error(stderr, result),
        )

    def complete_result(self, result, output, error, return_code):
        _, error, return_code = self.parser.result(output, error, return_code)
        result.complete(error, return_code)
        return result


def test_profile_builds_serial_command():
    input_file = SimpleNamespace(name="/tmp/Fe.inp")
    output_file = object()
    potential_file = "/tmp/Fe.pot"
    task_output = "/tmp/Fe_DOS.dos"
    profile = SPRKKRProfile(
        executable_dir=False,
        executable_suffix="",
        mpi=False,
    )

    spec = profile.create_execution(
        "kkrscf",
        input_file,
        output_file,
        "/tmp",
        input_files=(input_file, potential_file),
        output_files=(output_file, task_output),
        mpi_supported=True,
    )

    assert spec.argv == ("kkrscf",)
    assert spec.stdin is input_file
    assert spec.input_file is input_file
    assert spec.output_file is output_file
    assert spec.input_files == (input_file, potential_file)
    assert spec.output_files == (output_file, task_output)


def test_profile_builds_mpi_command():
    input_file = SimpleNamespace(name="/tmp/Fe.inp")
    profile = SPRKKRProfile(
        executable_dir="/opt/sprkkr",
        executable_suffix="_v97",
        mpi=["mpirun", "-np", "4"],
    )

    spec = profile.create_execution(
        "kkrscf",
        input_file,
        object(),
        "/tmp",
        mpi_supported=True,
    )

    assert spec.argv == (
        "mpirun",
        "-np",
        "4",
        "/opt/sprkkr/kkrscf_v97MPI",
        "/tmp/Fe.inp",
    )
    assert spec.stdin is None


def test_non_mpi_task_does_not_use_profile_mpi():
    input_file = SimpleNamespace(name="/tmp/Fe.inp")
    profile = SPRKKRProfile(
        executable_dir=False,
        executable_suffix="",
        mpi=["mpirun", "-np", "4"],
    )

    spec = profile.create_execution(
        "kkrgen",
        input_file,
        object(),
        "/tmp",
        mpi_supported=False,
    )

    assert spec.argv == ("kkrgen",)
    assert spec.stdin is input_file


def _local_run(code):
    directory = tempfile.TemporaryDirectory()
    output = open(f"{directory.name}/run.out", "wb")
    spec = ExecutionSpec(
        argv=(sys.executable, "-S", "-c", code),
        input_file=SimpleNamespace(name=f"{directory.name}/run.inp"),
        output_file=output,
        directory=directory.name,
    )
    runner = LocalExecutor().create_run(_Reader(), spec)
    runner._test_directory = directory
    return runner


def test_local_executor_returns_unstarted_runner():
    runner = _local_run("pass")
    assert isinstance(runner, LocalRunner)
    assert runner.state == RunState.NEW
    assert runner.reference is None


def test_local_runner_streams_output_to_parser_and_listeners():
    runner = _local_run(
        "import sys; print('hello', flush=True); "
        "print('problem', file=sys.stderr, flush=True)"
    )
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((kind, text)))

    result = runner.run()

    assert runner.state == RunState.FINISHED
    assert result.output == b"hello\n"
    assert result.error_output == b"problem\n"
    assert result.return_code == 0
    assert ("out", "hello\n") in seen
    assert ("err", "problem\n") in seen


def test_local_runner_submit_then_wait():
    runner = _local_run("print('ok')")
    runner.submit()
    result = runner.wait()
    assert result.output == b"ok\n"
    assert runner.state == RunState.FINISHED


def test_local_runner_cancel_is_thread_safe():
    runner = _local_run("import time; time.sleep(30)")
    caller_thread = threading.get_ident()
    stopped_from = []
    original_stop = runner._stop_process

    async def stop():
        stopped_from.append(threading.get_ident())
        await original_stop()

    runner._stop_process = stop
    runner.submit()

    deadline = time.monotonic() + 5
    while runner.state != RunState.RUNNING and time.monotonic() < deadline:
        time.sleep(0.01)
    assert runner.state == RunState.RUNNING

    runner.cancel()
    with pytest.raises(ValueError):
        runner.wait()

    assert runner.state == RunState.CANCELLED
    assert stopped_from
    assert stopped_from[0] != caller_thread

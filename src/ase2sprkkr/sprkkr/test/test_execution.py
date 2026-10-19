import asyncio
import sys
import tempfile
import threading
import time
from types import SimpleNamespace

import pytest

from ...execution import (
    ExecutionSpec,
    LaunchSpec,
    LocalExecutor,
    LocalRunner,
    MpiSpec,
    RunState,
)
from ..profile import SPRKKRProfile


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


def test_profile_builds_portable_mpi_launch_intent():
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

    assert spec.argv == ("/opt/sprkkr/kkrscf_v97",)
    assert spec.mpi_argv == ("/opt/sprkkr/kkrscf_v97MPI", "/tmp/Fe.inp")
    assert spec.launch == LaunchSpec(
        mpi=MpiSpec(required=True, launcher=("mpirun", "-np", "4"))
    )
    assert spec.stdin is input_file


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
        error_file=open(f"{directory.name}/run.err", "wb"),
        directory=directory.name,
    )
    runner = LocalExecutor().create_run(spec)
    runner._test_directory = directory
    return runner


def test_local_executor_returns_unstarted_runner():
    runner = _local_run("pass")
    assert isinstance(runner, LocalRunner)
    assert runner.state == RunState.NEW
    assert runner.backend_id is None


def test_local_runner_streams_output_to_parser_and_listeners():
    runner = _local_run(
        "import sys; print('hello', flush=True); "
        "print('problem', file=sys.stderr, flush=True)"
    )
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((kind, text)))

    result = runner.run()

    assert runner.state == RunState.FINISHED
    from pathlib import Path
    assert Path(result.output_file).read_bytes() == b"hello\n"
    assert Path(result.error_file).read_bytes() == b"problem\n"
    assert result.return_code == 0
    assert ("out", "hello\n") in seen
    assert ("err", "problem\n") in seen


def test_local_runner_submit_then_wait():
    runner = _local_run("print('ok')")
    runner.submit()
    result = runner.wait()
    from pathlib import Path
    assert Path(result.output_file).read_bytes() == b"ok\n"
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
    runner.wait()

    assert runner.state == RunState.CANCELLED
    assert stopped_from
    assert stopped_from[0] != caller_thread


def test_local_binary_output_does_not_abort_observers():
    runner = _local_run("import os; os.write(1, b'\\xff\\n')")
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((text, kind)))
    assert runner.run().return_code == 0
    assert seen == [("\\ufffd\\n", "out")]

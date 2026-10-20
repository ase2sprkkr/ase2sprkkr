"""End-to-end supervisor lifecycle tests, independent of the SSH transport."""

import asyncio
import os
import sys

from .. import remote_process_worker
from ..remote_supervisor_stream import follow_socket


def _spec(tmp_path, socket_fd, *, detached, argv=None):
    return {
        "argv": argv or [sys.executable, "-c", "import time; time.sleep(20)"],
        "directory": str(tmp_path),
        "output_file": str(tmp_path / "stdout"),
        "error_file": str(tmp_path / "stderr"),
        "receipt_id": "a" * 32,
        "socket_fd": socket_fd,
        "detached": detached,
    }


def test_attached_run_stops_when_follow_connection_drops(tmp_path, monkeypatch):
    async def check():
        captured = []
        monkeypatch.setattr(remote_process_worker, "publish_outcome",
                            lambda *args: captured.append(args))
        rd, wr = os.pipe()
        task = asyncio.create_task(remote_process_worker._run(_spec(tmp_path, wr, detached=False)))
        with os.fdopen(rd, "rb") as pipe:
            path = (await asyncio.to_thread(pipe.readline)).decode().strip()
        reader, writer = await asyncio.open_unix_connection(path)
        writer.write(b'{"command":"follow","offsets":{"out":0,"err":0}}\n')
        await writer.drain()
        assert await reader.readexactly(8) == b"A2SFOLL1"
        writer.close()
        await writer.wait_closed()
        assert await asyncio.wait_for(task, timeout=5) != 0
        assert captured[0][1].value == "cancelled"
    asyncio.run(check())


def test_detached_run_survives_follow_disconnect(tmp_path, monkeypatch):
    async def check():
        captured = []
        monkeypatch.setattr(remote_process_worker, "publish_outcome",
                            lambda *args: captured.append(args))
        rd, wr = os.pipe()
        spec = _spec(tmp_path, wr, detached=True, argv=[
            sys.executable, "-c",
            "import time; print('first', flush=True); time.sleep(.2); print('last', flush=True)",
        ])
        task = asyncio.create_task(remote_process_worker._run(spec))
        with os.fdopen(rd, "rb") as pipe:
            path = (await asyncio.to_thread(pipe.readline)).decode().strip()
        reader, writer = await asyncio.open_unix_connection(path)
        writer.write(b'{"command":"follow","offsets":{"out":0,"err":0}}\n')
        await writer.drain()
        assert await reader.readexactly(8) == b"A2SFOLL1"
        writer.close()
        await writer.wait_closed()
        assert await asyncio.wait_for(task, timeout=5) == 0
        assert (tmp_path / "stdout").read_bytes() == b"first\nlast\n"
        assert captured[0][1].value == "finished"
    asyncio.run(check())


def test_cancel_is_delivered_to_supervisor_not_just_supervisor_pid(tmp_path, monkeypatch):
    async def check():
        captured = []
        monkeypatch.setattr(remote_process_worker, "publish_outcome",
                            lambda *args, **kwargs: captured.append(args))
        rd, wr = os.pipe()
        task = asyncio.create_task(remote_process_worker._run(
            _spec(tmp_path, wr, detached=True)
        ))
        with os.fdopen(rd, "rb") as pipe:
            path = (await asyncio.to_thread(pipe.readline)).decode().strip()
        reader, writer = await asyncio.open_unix_connection(path)
        writer.write(b'{"command":"cancel"}\n')
        await writer.drain()
        assert b'"accepted":true' in await reader.readline()
        writer.close()
        await writer.wait_closed()
        await asyncio.wait_for(task, 5)
        assert captured[0][1].value == "cancelled"
    asyncio.run(check())


def test_failed_detached_launch_clears_previous_output(tmp_path, monkeypatch):
    async def check():
        captured = []
        monkeypatch.setattr(
            remote_process_worker, "publish_outcome",
            lambda *args, **kwargs: captured.append(args),
        )
        (tmp_path / "stdout").write_bytes(b"previous")
        rd, wr = os.pipe()
        spec = _spec(tmp_path, wr, detached=True, argv=[
            str(tmp_path / "missing-command"),
        ])
        result = await remote_process_worker._run(spec)
        with os.fdopen(rd, "rb") as pipe:
            assert pipe.readline()
        assert result != 0
        assert (tmp_path / "stdout").read_bytes() == b""
        assert captured[0][1].value == "failed"
        assert not list(tmp_path.glob(".*ase2sprkkr-*"))
    asyncio.run(check())


def test_output_spool_write_failure_is_reported_without_losing_diagnosis(tmp_path, monkeypatch):
    from ..base import FailurePhase, RunState
    from ..remote_supervisor_stream import SupervisorStreams

    reported = []
    monkeypatch.setattr(
        remote_process_worker, "publish_outcome",
        lambda *args, **kwargs: reported.append((args, kwargs)),
    )

    def fail_write(self, kind, data):
        raise OSError("disk full during output")

    monkeypatch.setattr(SupervisorStreams, "publish", fail_write)

    async def check():
        rd, wr = os.pipe()
        spec = _spec(tmp_path, wr, detached=True, argv=[
            sys.executable, "-c", "print('payload', flush=True)",
        ])
        await asyncio.wait_for(remote_process_worker._run(spec), 5)
        with os.fdopen(rd, "rb") as read_pipe:
            assert read_pipe.readline()

    asyncio.run(check())
    assert len(reported) == 1
    args, kwargs = reported[0]
    assert args[1] == RunState.FAILED
    assert args[3] == FailurePhase.EXECUTION
    assert "Failed to write process output: disk full during output" in kwargs["error"]


def test_failed_final_output_sync_is_reported_as_artifact_failure(tmp_path, monkeypatch):
    from ..base import FailurePhase, RunState
    from ..remote_supervisor_stream import SupervisorStreams

    reported = []
    monkeypatch.setattr(
        remote_process_worker, "publish_outcome",
        lambda *args, **kwargs: reported.append((args, kwargs)),
    )

    def fail_sync(self):
        raise OSError("cannot fsync output")

    monkeypatch.setattr(SupervisorStreams, "sync", fail_sync)

    async def check():
        rd, wr = os.pipe()
        spec = _spec(tmp_path, wr, detached=True, argv=[
            sys.executable, "-c", "print('payload', flush=True)",
        ])
        await asyncio.wait_for(remote_process_worker._run(spec), 5)
        with os.fdopen(rd, "rb") as read_pipe:
            assert read_pipe.readline()

    asyncio.run(check())
    assert len(reported) == 1
    args, kwargs = reported[0]
    assert args[1] == RunState.FAILED
    assert args[2] == 0
    assert args[3] == FailurePhase.ARTIFACTS
    assert kwargs["error"] == "Failed to sync process output: cannot fsync output"


def test_receipt_write_failure_never_emits_durable_terminal_claim(tmp_path, monkeypatch):
    from ..remote_supervisor_stream import SupervisorStreams

    observed = []
    original_terminal = SupervisorStreams.terminal

    def record_terminal(self, outcome):
        observed.append(outcome)
        original_terminal(self, outcome)

    monkeypatch.setattr(SupervisorStreams, "terminal", record_terminal)

    def fail_receipt(*args, **kwargs):
        raise OSError("mailbox fsync failed")

    monkeypatch.setattr(remote_process_worker, "publish_outcome", fail_receipt)
    monkeypatch.setattr(remote_process_worker, "_write_outcome", fail_receipt)

    async def check():
        rd, wr = os.pipe()
        spec = _spec(tmp_path, wr, detached=True, argv=[
            sys.executable, "-c", "print('payload', flush=True)",
        ])
        assert await asyncio.wait_for(remote_process_worker._run(spec), 5) == 0
        with os.fdopen(rd, "rb") as read_pipe:
            assert read_pipe.readline()

    asyncio.run(check())
    assert observed == [{
        "state": "unknown",
        "return_code": None,
        "failure_phase": None,
        "error": "Unable to persist terminal receipt: mailbox fsync failed",
    }]
    assert (tmp_path / "stdout").read_bytes() == b"payload\n"



def test_supervisor_pumps_both_binary_channels_to_separate_files(tmp_path, monkeypatch):
    reported = []
    monkeypatch.setattr(
        remote_process_worker, "publish_outcome",
        lambda *args, **kwargs: reported.append((args, kwargs)),
    )

    async def run():
        rd, wr = os.pipe()
        spec = _spec(
            tmp_path, wr, detached=True,
            argv=[
                sys.executable, "-c",
                "import sys; "
                "sys.stdout.buffer.write(b'A' * 65536 + bytes([0, 255])); "
                "sys.stdout.buffer.flush(); "
                "sys.stderr.buffer.write(b'B' * 65536 + bytes([1, 254])); "
                "sys.stderr.buffer.flush()",
            ],
        )
        return_code = await asyncio.wait_for(remote_process_worker._run(spec), 5)
        with os.fdopen(rd, "rb") as pipe:
            assert pipe.readline()
        return return_code

    assert asyncio.run(run()) == 0
    assert (tmp_path / "stdout").read_bytes() == b"A" * 65536 + bytes([0, 255])
    assert (tmp_path / "stderr").read_bytes() == b"B" * 65536 + bytes([1, 254])
    assert reported[0][0][1].value == "finished"


def test_managed_calculation_has_separate_group_but_same_session(tmp_path):
    import subprocess

    async def check():
        spec = _spec(tmp_path, -1, detached=True, argv=[
            sys.executable, "-c",
            "import os, time; print(os.getsid(0), os.getpgrp(), flush=True); time.sleep(1)",
        ])
        process, guard, write_fd = await remote_process_worker._launch_managed_process(
            spec, subprocess.DEVNULL, asyncio.subprocess.PIPE,
        )
        try:
            assert os.getsid(process.pid) == os.getsid(0)
            assert os.getpgid(process.pid) == guard.pid
            assert os.getpgrp() != guard.pid
            await asyncio.wait_for(process.wait(), 5)
        finally:
            os.close(write_fd)
            await asyncio.wait_for(guard.wait(), 5)

    asyncio.run(check())


def test_watchdog_kills_calculation_group_when_supervisor_dies(tmp_path):
    import signal
    import subprocess

    async def check():
        ready = tmp_path / "child-pid"
        script = (
            "import os, signal, time, pathlib\n"
            "signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
            f"pathlib.Path({str(ready)!r}).write_text(str(os.getpid()))\n"
            "while True: time.sleep(10)\n"
        )
        spec = _spec(tmp_path, -1, detached=True, argv=[
            sys.executable, "-c", script,
        ])
        process, guard, write_fd = await remote_process_worker._launch_managed_process(
            spec, subprocess.DEVNULL, asyncio.subprocess.PIPE,
        )
        try:
            for _ in range(100):
                if ready.exists():
                    break
                await asyncio.sleep(0.02)
            assert ready.exists(), "Calculation never started"
            os.close(write_fd)  # The supervisor has disappeared.
            write_fd = None
            assert await asyncio.wait_for(process.wait(), 5) == -signal.SIGKILL
            await asyncio.wait_for(guard.wait(), 5)
        finally:
            if write_fd is not None:
                os.close(write_fd)
            if process.returncode is None:
                os.killpg(guard.pid, signal.SIGKILL)
                await process.wait()
            if guard.returncode is None:
                await guard.wait()

    asyncio.run(check())


def test_supervisor_sigkill_stops_calculation_and_grandchild(tmp_path):
    """Real supervisor death must not leave a calculation process tree running."""
    import json
    import signal
    import subprocess
    import time

    pid_file = tmp_path / "calculation-pids"
    calculation = (
        "import os, subprocess, sys, time, pathlib\n"
        "sub = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])\n"
        f"pathlib.Path({str(pid_file)!r}).write_text(f'{{os.getpid()}} {{sub.pid}}')\n"
        "time.sleep(60)\n"
    )
    rd, wr = os.pipe()
    spec = _spec(tmp_path, wr, detached=True, argv=[
        sys.executable, "-c", calculation,
    ])
    socket_dir = tmp_path / "socket-dir"
    socket_dir.mkdir()  # The submitter creates this directory before starting the worker.
    spec["socket_dir"] = str(socket_dir)
    env = os.environ.copy()
    env[remote_process_worker._SPEC_ENV] = json.dumps(spec)
    supervisor = subprocess.Popen(
        [sys.executable, "-m", "ase2sprkkr.execution.remote_process_worker"],
        cwd=str(tmp_path), env=env, start_new_session=True,
        pass_fds=(wr,), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
    )
    os.close(wr)
    pids = []
    try:
        with os.fdopen(rd, "rb") as reader:
            assert reader.readline(), "Supervisor did not publish follow socket"
        for _ in range(150):
            if pid_file.is_file():
                pids = list(map(int, pid_file.read_text().split()))
                break
            if supervisor.poll() is not None:
                raise AssertionError("Supervisor died before launching calculation")
            time.sleep(0.02)
        assert len(pids) == 2
        os.kill(supervisor.pid, signal.SIGKILL)
        supervisor.wait(timeout=5)

        def is_alive(pid):
            try:
                state = open(f"/proc/{pid}/stat").read().rsplit(")", 1)[1].split()[0]
            except FileNotFoundError:
                return False
            return state != "Z"

        for _ in range(150):
            if not any(is_alive(pid) for pid in pids):
                break
            time.sleep(0.03)
        assert not any(is_alive(pid) for pid in pids)
    finally:
        if supervisor.poll() is None:
            supervisor.kill()
            supervisor.wait(timeout=5)
        for pid in pids:
            if os.path.exists(f"/proc/{pid}"):
                try:
                    os.kill(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass


def test_launch_failure_after_child_creation_waits_for_calculation(tmp_path, monkeypatch):
    """The watchdog stopping a group does not reap the calculation child."""
    import subprocess

    calls = []

    class Watchdog:
        pid = 1234

        class stdout:
            @staticmethod
            async def readexactly(_n):
                return b"R"

        async def wait(self):
            calls.append("watchdog")

    class Calculation:
        async def wait(self):
            calls.append("calculation")

    async def spawn(*argv, **kwargs):
        return Calculation() if "--exec" in argv else Watchdog()

    async def cancelled_handshake(*args):
        raise asyncio.CancelledError

    monkeypatch.setattr(remote_process_worker.asyncio, "create_subprocess_exec", spawn)
    monkeypatch.setattr(remote_process_worker.asyncio, "to_thread", cancelled_handshake)

    async def check():
        spec = _spec(tmp_path, -1, detached=True)
        try:
            await remote_process_worker._launch_managed_process(
                spec, subprocess.DEVNULL, asyncio.subprocess.PIPE,
            )
        except asyncio.CancelledError:
            pass
        else:
            assert False, "The simulated cancellation must propagate"

    asyncio.run(check())
    assert calls == ["watchdog", "calculation"]


def test_watchdog_eof_during_ready_handshake_cleans_up(tmp_path, monkeypatch):
    """A watchdog dying before readiness must not launch a calculation."""
    import pytest
    import subprocess

    events = []

    class FailedWatchdog:
        pid = 1234

        class stdout:
            @staticmethod
            async def readexactly(_size):
                raise asyncio.IncompleteReadError(b"", 1)

        async def wait(self):
            events.append("reaped")

    async def spawn(*argv, **kwargs):
        events.append("started")
        return FailedWatchdog()

    monkeypatch.setattr(remote_process_worker.asyncio, "create_subprocess_exec", spawn)

    async def check():
        with pytest.raises(asyncio.IncompleteReadError):
            await remote_process_worker._launch_managed_process(
                _spec(tmp_path, -1, detached=True),
                subprocess.DEVNULL, asyncio.subprocess.PIPE,
            )

    asyncio.run(check())
    assert events == ["started", "reaped"]

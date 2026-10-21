"""SIGTERM/SIGKILL escalation shared by direct and supervisor runners."""

import asyncio
import os
import select
import subprocess
import signal
import sys

import pytest

from ..local_process import stop_local_process


@pytest.mark.skipif(os.name != "posix", reason="Process-group signals require POSIX")
@pytest.mark.parametrize(
    "process_group,reuse_wait",
    ((False, False), (False, True), (True, True)),
)
def test_stop_ignoring_sigterm_escalates_without_cancelling_wait(
    process_group, reuse_wait,
):
    async def exercise():
        process = await asyncio.create_subprocess_exec(
            sys.executable, "-u", "-c",
            "import signal, time; "
            "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
            "print('ready', flush=True); time.sleep(30)",
            stdout=asyncio.subprocess.PIPE,
            start_new_session=True,
        )
        waiting = asyncio.create_task(process.wait()) if reuse_wait else None
        try:
            assert await asyncio.wait_for(process.stdout.readline(), 5) == b"ready\n"
            await asyncio.wait_for(
                stop_local_process(
                    process, process_group=process_group,
                    timeout=0.05, waiting=waiting,
                ),
                5,
            )
            assert process.returncode == -signal.SIGKILL
            if waiting is not None:
                assert waiting.done()
                assert not waiting.cancelled()
                assert waiting.result() == -signal.SIGKILL
        finally:
            if process.returncode is None:
                process.kill()
            await process.wait()

    asyncio.run(exercise())


@pytest.mark.skipif(os.name != "posix", reason="Process-group signals require POSIX")
def test_stop_process_group_signals_surviving_descendants():
    # The leader exits normally but leaves a SIGTERM-resistant descendant.
    # The descendant inherits stdout, so an EOF proves that it was killed.
    descendant = (
        "import signal, time; "
        "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
        "print('ready', flush=True); time.sleep(30)"
    )
    leader_code = (
        "import subprocess, sys; "
        f"subprocess.Popen([sys.executable, '-u', '-c', {descendant!r}])"
    )
    leader = subprocess.Popen(
        [sys.executable, "-u", "-c", leader_code],
        stdout=subprocess.PIPE,
        start_new_session=True,
        bufsize=0,
    )
    try:
        assert select.select([leader.stdout], [], [], 5)[0], "descendant did not start"
        assert leader.stdout.readline() == b"ready\n"
        assert leader.wait(timeout=5) == 0

        class FinishedLeader:
            pid = leader.pid
            returncode = 0

            async def wait(self):
                return 0

        asyncio.run(stop_local_process(FinishedLeader(), process_group=True, timeout=0.05))
        assert select.select([leader.stdout], [], [], 5)[0], "descendant survived SIGTERM"
        assert leader.stdout.read(1) == b""
    finally:
        try:
            os.killpg(leader.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        leader.stdout.close()
        leader.wait(timeout=5)

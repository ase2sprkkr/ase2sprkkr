"""Regression coverage for process-group termination."""

import asyncio
import signal

import pytest

from .. import local_process


@pytest.mark.parametrize("leader_already_exited", [False, True])
def test_group_escalation_does_not_stop_at_leader_exit(monkeypatch, leader_already_exited):
    signals = []

    def kill_group(pgid, signum):
        signals.append((pgid, signum))

    monkeypatch.setattr(local_process.os, "killpg", kill_group)

    class Leader:
        pid = 12345
        returncode = -signal.SIGTERM if leader_already_exited else None

        async def wait(self):
            self.returncode = -signal.SIGTERM
            return self.returncode

    asyncio.run(local_process.stop_local_process(
        Leader(), process_group=True, timeout=0.01,
    ))

    assert signals == [
        (Leader.pid, signal.SIGTERM),
        (Leader.pid, signal.SIGKILL),
    ]

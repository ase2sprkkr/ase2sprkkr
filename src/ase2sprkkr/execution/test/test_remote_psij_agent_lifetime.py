import time

from ase2sprkkr.execution.remote_psij import _ActiveRun, _RemotePsijAgent


def test_psij_agent_does_not_idle_out_while_job_is_active(tmp_path, monkeypatch):
    agent = _RemotePsijAgent(tmp_path)
    agent._active["run"] = _ActiveRun(object(), object())

    polls = 0

    def poll_active():
        nonlocal polls
        polls += 1
        if polls == 3:
            agent._active.clear()
        return False

    monkeypatch.setattr(agent, "_discover", lambda: False)
    monkeypatch.setattr(agent, "_poll_active", poll_active)
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_psij._AGENT_IDLE_TIMEOUT", 0.0
    )
    monkeypatch.setattr(
        "ase2sprkkr.execution.remote_psij._AGENT_POLL_INTERVAL", 0.0
    )
    monkeypatch.setattr(time, "sleep", lambda _: None)

    agent.run()

    # The agent may need one extra idle poll after the active run disappears.
    # The invariant is that it must not exit while the run is still active.
    assert polls >= 3

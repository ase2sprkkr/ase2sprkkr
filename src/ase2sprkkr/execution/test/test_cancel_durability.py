from contextlib import nullcontext
import json

import ase2sprkkr.execution.detached as detached


def _prepare_request(tmp_path, run_id="run"):
    run_directory = detached._run_directory(tmp_path, run_id)
    run_directory.mkdir(parents=True)
    detached._request_path(tmp_path, run_id).write_text(
        json.dumps({"executor": "local"})
    )
    return run_id


def _fake_state(monkeypatch, state, writes):
    monkeypatch.setattr(detached, "_state_lock", lambda *_args: nullcontext())
    monkeypatch.setattr(
        detached, "_read_state_unlocked", lambda *_args: dict(state)
    )
    monkeypatch.setattr(detached, "_read_state", lambda *_args: dict(state))

    def write(_directory, _run_id, **changes):
        state.update(changes)
        writes.append(dict(changes))
        return dict(state)

    monkeypatch.setattr(detached, "_write_state_unlocked", write)


def test_running_cancel_is_not_terminal_before_process_group_stops(
    tmp_path, monkeypatch
):
    run_id = _prepare_request(tmp_path)
    state = {"state": "running", "pid": 123, "finished_at": None}
    writes = []
    _fake_state(monkeypatch, state, writes)

    def terminate(*_args):
        assert state["state"] == "running"
        return dict(state), True, True

    monkeypatch.setattr(detached, "_terminate_process_group", terminate)

    result = detached.cancel_detached(tmp_path, run_id)

    assert result["state"] == "cancelled"
    assert [change.get("state") for change in writes] == ["cancelled"]


def test_running_cancel_stays_active_when_process_group_survives(
    tmp_path, monkeypatch
):
    run_id = _prepare_request(tmp_path)
    state = {"state": "running", "pid": 123, "finished_at": None}
    writes = []
    _fake_state(monkeypatch, state, writes)
    monkeypatch.setattr(
        detached,
        "_terminate_process_group",
        lambda *_args: (dict(state), False, True),
    )

    result = detached.cancel_detached(tmp_path, run_id)

    assert result["state"] == "running"
    assert writes == []


def test_prelaunch_cancel_is_safe_while_holding_state_lock(tmp_path, monkeypatch):
    run_id = _prepare_request(tmp_path)
    state = {"state": "submitting", "pid": None, "finished_at": None}
    writes = []
    _fake_state(monkeypatch, state, writes)

    def unexpected_terminate(*_args):
        raise AssertionError("no process exists to terminate")

    monkeypatch.setattr(detached, "_terminate_process_group", unexpected_terminate)

    result = detached.cancel_detached(tmp_path, run_id)

    assert result["state"] == "cancelled"
    assert [change.get("state") for change in writes] == ["cancelled"]

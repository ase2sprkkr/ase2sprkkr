from uuid import uuid4

import pytest

from ase2sprkkr.execution import FailurePhase, RunState
from ase2sprkkr.execution import remote_mailbox


@pytest.fixture(autouse=True)
def isolated_mailbox(tmp_path, monkeypatch):
    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path))


def test_terminal_outcome_survives_reads_until_ack():
    key = uuid4().hex
    assert remote_mailbox.read_outcome(key) is None
    assert remote_mailbox.publish_outcome(key, RunState.FAILED, 7, FailurePhase.EXECUTION)
    expected = {
        "state": "failed",
        "return_code": 7,
        "failure_phase": "execution",
        "error": None,
    }
    assert remote_mailbox.read_outcome(key) == expected
    assert remote_mailbox.read_outcome(key) == expected
    assert remote_mailbox.acknowledge_outcome(key)
    assert not remote_mailbox.acknowledge_outcome(key)
    assert remote_mailbox.read_outcome(key) is None


def test_first_writer_wins_even_when_cancel_races_with_completion():
    key = uuid4().hex
    assert remote_mailbox.publish_outcome(key, RunState.FINISHED, 0)
    assert not remote_mailbox.publish_outcome(key, RunState.CANCELLED, None)
    assert remote_mailbox.read_outcome(key)["state"] == "finished"


def test_rejects_invalid_id_and_nonterminal_state(tmp_path):
    with pytest.raises(ValueError, match="receipt ID"):
        remote_mailbox.read_outcome("../outside")
    with pytest.raises(ValueError, match="terminal"):
        remote_mailbox.publish_outcome(uuid4().hex, RunState.RUNNING, None)
    assert list(tmp_path.rglob("*.json")) == []


def test_receipt_uses_configured_state_root(tmp_path, monkeypatch):
    key = uuid4().hex
    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path / "explicit"))
    assert remote_mailbox._path(key) == (
        tmp_path / "explicit" / "remote-outcomes" / f"{key}.json"
    )


def test_receipt_uses_xdg_state_home_when_no_override(tmp_path, monkeypatch):
    key = uuid4().hex
    monkeypatch.delenv("ASE2SPRKKR_REMOTE_STATE_DIR", raising=False)
    monkeypatch.setenv("XDG_STATE_HOME", str(tmp_path / "xdg"))
    assert remote_mailbox._path(key) == (
        tmp_path / "xdg" / "ase2sprkkr" / "remote-outcomes" / f"{key}.json"
    )

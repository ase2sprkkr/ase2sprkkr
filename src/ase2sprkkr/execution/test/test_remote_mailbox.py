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


@pytest.mark.parametrize("override", (True, False), ids=("explicit", "xdg"))
def test_receipt_state_root_respects_override_and_xdg(tmp_path, monkeypatch, override):
    key = uuid4().hex
    if override:
        monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path / "explicit"))
        root = tmp_path / "explicit"
    else:
        monkeypatch.delenv("ASE2SPRKKR_REMOTE_STATE_DIR", raising=False)
        monkeypatch.setenv("XDG_STATE_HOME", str(tmp_path / "xdg"))
        root = tmp_path / "xdg" / "ase2sprkkr"
    assert remote_mailbox._path(key) == root / "remote-outcomes" / f"{key}.json"

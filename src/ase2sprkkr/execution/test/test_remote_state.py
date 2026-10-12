from ..remote import _read_state, _run_directory, _write_state


def test_remote_state_accepts_directory_field(tmp_path):
    run_id = "run"
    _run_directory(tmp_path, run_id).mkdir(parents=True)

    state = _write_state(
        tmp_path,
        run_id,
        state="queued",
        directory="/remote/workspace",
    )

    assert state["directory"] == "/remote/workspace"
    assert _read_state(tmp_path, run_id)["directory"] == "/remote/workspace"

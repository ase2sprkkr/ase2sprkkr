import io
import sys

from .. import detached


def test_detached_worker_and_follow_handle_merged_output_once(tmp_path, monkeypatch):
    output = tmp_path / "combined.log"
    request = {
        "id": "merged",
        "executor": "local",
        "directory": str(tmp_path),
        "argv": [
            sys.executable,
            "-c",
            "import sys; print('stdout'); print('stderr', file=sys.stderr)",
        ],
        "input_file": None,
        "stdin": None,
        "output_file": str(output),
        "error_file": str(output),
        "launch": {},
        "mpi_argv": None,
    }

    monkeypatch.setattr(detached, "_spawn_detached_worker", lambda *args: None)
    detached.submit_detached(request)

    detached.run_detached_worker(tmp_path, "merged")

    combined = output.read_bytes()
    assert combined.count(b"stdout\n") == 1
    assert combined.count(b"stderr\n") == 1

    stdout = io.BytesIO()
    stderr = io.BytesIO()
    state = detached.follow_detached(
        tmp_path,
        "merged",
        stdout,
        stderr,
        poll_interval=0,
    )

    assert state["state"] == "finished"
    assert stdout.getvalue() == combined
    assert stderr.getvalue() == b""

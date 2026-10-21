"""Local projects and SSH use the same detached supervisor and submission ticket."""
import os
import sys
import time
from dataclasses import replace

import pytest

from ase2sprkkr.execution import (
    Execution, ExecutionSpec, ExecutorEndpoint, FileEntry, FileManifest,
    LocalExecutor, RunRecord, RunState, WorkingDirectory,
)
from ase2sprkkr.execution.local_project import local_manifest
from ase2sprkkr.execution.remote_service import (
    _decode_handle, _encode_handle, _process_running, _recover_submission,
    acknowledge_submission, status_remote,
)
from ase2sprkkr.execution.base import RunStatus, FailurePhase

RUN_ID = "b" * 32

pytestmark = pytest.mark.skipif(
    not sys.platform.startswith("linux"),
    reason="The shared supervisor uses Linux /proc process identity",
)


@pytest.fixture(autouse=True)
def private_backend_state(tmp_path, monkeypatch):
    monkeypatch.setenv("ASE2SPRKKR_REMOTE_STATE_DIR", str(tmp_path / "backend-state"))


def _spec(tmp_path, command=None, *, out="out.txt", err="err.txt"):
    return ExecutionSpec(
        argv=command or (sys.executable, "-c", "print('output', flush=True)"),
        input_file=None, output_file=None if out is None else str(tmp_path / out),
        error_file=None if err is None else str(tmp_path / err),
        directory=str(tmp_path),
    )


def _record(tmp_path, backend_id, *, stdout=True, stderr=True):
    return RunRecord(
        id=RUN_ID, state=RunState.RUNNING, task="DOS",
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            ExecutorEndpoint(), WorkingDirectory(str(tmp_path)), backend_id,
        ),
        files=FileManifest(
            stdout=FileEntry(local_path=str(tmp_path / "out.txt")) if stdout else None,
            stderr=FileEntry(local_path=str(tmp_path / "err.txt")) if stderr else None,
        ),
    )


def test_project_launch_reconnect_and_shared_ticket(tmp_path):
    executor = LocalExecutor()
    runner = executor.create_project_run(_spec(tmp_path), run_id=RUN_ID)
    runner.submit()
    handle = runner.backend_id
    assert _decode_handle(handle)["kind"] == "process"
    assert _recover_submission(RUN_ID) == {"remote_id": handle}
    assert runner.wait().return_code == 0
    assert runner.state == RunState.FINISHED
    assert not _process_running({**_decode_handle(handle), "start_time": "0"})
    assert (tmp_path / "out.txt").read_text() == "output\n"
    assert executor.status(_record(tmp_path, handle))[0].state == RunState.FINISHED
    assert executor.recover_submission(RUN_ID, _record(tmp_path, None)).backend_id == handle

    executor.acknowledge_submission(_record(tmp_path, handle))
    assert _recover_submission(RUN_ID) is None
    assert executor.status(_record(tmp_path, handle))[0].state == RunState.FINISHED

    seen = []
    attached = executor.attach(
        _record(tmp_path, handle),
        output_listener=lambda text, kind: seen.append((text, kind)),
    )
    assert attached.wait().return_code == 0
    assert seen == [("output\n", "out")]
    executor.acknowledge(_record(tmp_path, handle))
    assert runner.state == RunState.FINISHED
    assert runner.return_code == 0


def test_project_file_objects_are_not_serialized(tmp_path):
    spec = _spec(tmp_path)
    with (tmp_path / "extra.out").open("w+b") as extra, (
        tmp_path / "cache.pot"
    ).open("w+b") as cache:
        spec = replace(spec, output_files=(extra,), cache_files=(cache,))
        runner = LocalExecutor().create_project_run(spec, run_id=RUN_ID)
        runner.submit()
        assert runner.wait().return_code == 0


def test_process_failure_preserves_execution_phase(tmp_path):
    runner = LocalExecutor().create_project_run(
        _spec(tmp_path, (sys.executable, "-c", "import sys; sys.exit(7)")),
        run_id=RUN_ID,
    )
    runner.submit()
    assert runner.wait().return_code == 7
    assert runner.state == RunState.FAILED
    assert runner.failure_phase == FailurePhase.EXECUTION


def test_cancel_is_supervisor_owned_and_nonblocking(tmp_path):
    runner = LocalExecutor().create_project_run(
        _spec(tmp_path, (
            sys.executable, "-c",
            "import time; print('started', flush=True); time.sleep(30)",
        )),
        run_id=RUN_ID,
    )
    runner.submit()
    started = time.monotonic()
    runner.cancel()
    assert time.monotonic() - started < 2
    runner.wait()
    assert runner.state == RunState.CANCELLED
    assert status_remote(runner.backend_id)["state"] == "cancelled"


def test_missing_ticket_keeps_submission_uncertain(tmp_path):
    with pytest.raises(NotImplementedError, match="not recoverable"):
        LocalExecutor().recover_submission(RUN_ID, _record(tmp_path, None))


def test_missing_supervisor_outcome_is_unknown(tmp_path):
    handle = _encode_handle({
        "kind": "process", "pid": 9999999, "start_time": "0",
        "receipt_id": "a" * 32, "socket_path": "/missing/socket",
    })
    assert LocalExecutor().status(_record(tmp_path, handle))[0].state == RunState.UNKNOWN


def test_wildcards_track_concrete_file_availability(tmp_path):
    (tmp_path / "one.dat").write_text("x")
    files = FileManifest(outputs=(
        FileEntry(local_path=str(tmp_path / "*.dat"), availability_state="pending"),
    ))
    status = RunStatus(RunState.FINISHED, backend_id="test")
    updated = local_manifest(files, status)
    assert len(updated.outputs) == 2
    assert updated.outputs[0].availability_state == "available"
    assert updated.outputs[1].source_pattern == str(tmp_path / "*.dat")
    assert updated.outputs[1].availability_state == "available"


def test_failed_launch_clears_stale_stdout(tmp_path):
    out = tmp_path / "out.txt"
    out.write_bytes(b"previous")
    runner = LocalExecutor().create_project_run(
        _spec(tmp_path, (str(tmp_path / "missing-program"),)), run_id=RUN_ID,
    )
    runner.submit()
    runner.wait()
    assert runner.state == RunState.FAILED
    assert runner.failure_phase == FailurePhase.LAUNCH
    assert out.read_bytes() == b""


def test_none_stream_paths_are_live_only(tmp_path):
    runner = LocalExecutor().create_project_run(
        _spec(tmp_path, out=None, err=None), run_id=RUN_ID,
    )
    seen = []
    runner.add_output_listener(lambda text, kind: seen.append((text, kind)))
    runner.submit()
    assert runner.wait().return_code == 0
    assert not (tmp_path / "out.txt").exists()
    assert not (tmp_path / "err.txt").exists()
    attached = LocalExecutor().attach(
        _record(tmp_path, runner.backend_id, stdout=False, stderr=False),
        start=False,
    )
    assert attached.backend_id == runner.backend_id
    assert attached.result.output_file is None
    assert attached.result.error_file is None
    assert attached.state == RunState.FINISHED


def test_terminal_without_ticket_is_still_accessible_by_handle(tmp_path):
    runner = LocalExecutor().create_project_run(_spec(tmp_path), run_id=RUN_ID)
    runner.submit()
    runner.wait()
    assert acknowledge_submission(RUN_ID)["acknowledged"]
    assert LocalExecutor().status(_record(tmp_path, runner.backend_id))[0].state == RunState.FINISHED


def test_local_prelaunch_error_is_known_submission_failure(tmp_path, monkeypatch):
    from ase2sprkkr.execution import remote_service
    from ase2sprkkr.execution.base import SubmissionNotStarted
    from ase2sprkkr.execution import FailurePhase

    def failed_launch(*args, **kwargs):
        raise OSError("cannot spawn supervisor")

    monkeypatch.setattr(remote_service.subprocess, "Popen", failed_launch)
    runner = LocalExecutor().create_project_run(_spec(tmp_path), run_id=RUN_ID)
    with pytest.raises(SubmissionNotStarted, match="cannot spawn supervisor"):
        runner.submit()
    assert runner.backend_id is None
    assert runner.state == RunState.FAILED
    assert runner.failure_phase == FailurePhase.SUBMISSION
    assert not _recover_submission(RUN_ID)


def test_cancelled_run_preserves_partial_stream_availability(tmp_path):
    out = tmp_path / "out.txt"
    out.write_bytes(b"partial output\n")
    files = FileManifest(
        stdout=FileEntry(local_path=str(out), availability_state="pending"),
    )
    manifest = local_manifest(files, RunStatus(RunState.CANCELLED))
    assert manifest.stdout.availability_state == "available"
    assert out.read_bytes() == b"partial output\n"



def test_terminal_local_runner_never_regresses_to_running(tmp_path, monkeypatch):
    from ase2sprkkr.execution import local_project
    runner = LocalExecutor().create_project_run(_spec(tmp_path), run_id=RUN_ID)
    runner._backend_id = "opaque-handle"
    runner._state = RunState.FINISHED
    runner._return_code = 0
    runner.result.return_code = 0
    monkeypatch.setattr(local_project, "local_status", lambda *args: RunStatus(RunState.RUNNING))
    assert runner.state == RunState.FINISHED
    assert runner.return_code == 0


def test_project_request_flushes_each_declared_file_once(tmp_path):
    class NamedFile:
        def __init__(self, name):
            self.name = str(tmp_path / name)
            self.flush_count = 0

        def flush(self):
            self.flush_count += 1

    output = NamedFile("supplement.out")
    inp = NamedFile("supplement.in")
    cache = NamedFile("restart.pot")
    spec = replace(
        _spec(tmp_path),
        output_files=(None, output, "a.dat", object()),
        input_files=(None, inp, "b.inp", 1),
        cache_files=(None, cache, "restart.pot", object()),
    )
    runner = LocalExecutor().create_project_run(spec, run_id=RUN_ID)
    request = runner._request()

    assert request["output_files"] == [output.name, str(tmp_path / "a.dat")]
    assert request["input_files"] == [inp.name, str(tmp_path / "b.inp")]
    assert request["cache_files"] == [cache.name, str(tmp_path / "restart.pot")]
    assert (output.flush_count, inp.flush_count, cache.flush_count) == (1, 1, 1)


def test_submission_cleanup_uses_parsed_manifest_and_protects_caches(tmp_path):
    from ase2sprkkr.execution.remote_service import (
        _clean_request_outputs, _request_spec,
    )

    for name in ("stale.dat", "input.dat", "stdin.dat", "cache.dat", "sync-cache.dat"):
        (tmp_path / name).write_text(name)

    request = {
        "directory": str(tmp_path),
        "argv": [sys.executable, "-c", "pass"],
        "input_file": "input.dat",
        "stdin": "stdin.dat",
        "output_files": ["*.dat"],
        "cache_files": ["cache.dat"],
        "client_sync": {
            "cache_files": [{"path": "sync-cache.dat", "before": None}],
        },
    }
    spec = _request_spec(request)
    _clean_request_outputs(spec, tmp_path)

    assert not (tmp_path / "stale.dat").exists()
    for name in ("input.dat", "stdin.dat", "cache.dat", "sync-cache.dat"):
        assert (tmp_path / name).read_text() == name


def test_local_project_request_uses_common_command_fields(tmp_path):
    from ase2sprkkr.execution import LaunchSpec
    from ase2sprkkr.execution.remote_service import _request_spec

    spec = replace(
        _spec(tmp_path),
        mpi_argv=("mpi-program", "input.inp"),
        mpi_uses_stdin=False,
        launch=LaunchSpec.from_mpi(3),
    )
    request = LocalExecutor().create_project_run(spec, run_id=RUN_ID)._request()
    for key, value in spec.command_request().items():
        assert request[key] == value
    reconstructed = _request_spec(request)
    assert reconstructed.mpi_uses_stdin is False
    assert reconstructed.mpi_argv == spec.mpi_argv


def test_terminal_project_runner_does_not_query_backend(monkeypatch, tmp_path):
    import pytest
    from ase2sprkkr.execution.local_project import LocalProjectRunner
    from ase2sprkkr.execution.base import RunState

    spec = _spec(tmp_path)
    runner = LocalProjectRunner(
        LocalExecutor(), spec, backend_id="old-supervisor", state=RunState.FINISHED,
    )
    monkeypatch.setattr(
        "ase2sprkkr.execution.local_project.local_status",
        lambda *args: pytest.fail("Terminal runner must not query backend"),
    )
    assert runner.state is RunState.FINISHED


def test_local_project_invalid_output_is_known_non_submission(tmp_path):
    """Output/cache collisions fail preflight without destroying existing files."""
    from ase2sprkkr.execution.base import SubmissionNotStarted

    existing = tmp_path / "out.txt"
    extra = tmp_path / "other.out"
    existing.write_text("restart cache")
    extra.write_text("old output")
    spec = replace(
        _spec(tmp_path),
        cache_files=(str(existing),),
        output_files=(str(extra),),
    )
    runner = LocalExecutor().create_project_run(spec, run_id=RUN_ID)
    with pytest.raises(SubmissionNotStarted, match="overlaps protected"):
        runner.submit()
    assert runner.state is RunState.FAILED
    assert runner.failure_phase is FailurePhase.SUBMISSION
    assert runner.backend_id is None
    assert existing.read_text() == "restart cache"
    assert extra.read_text() == "old output"
    assert _recover_submission(RUN_ID) is None

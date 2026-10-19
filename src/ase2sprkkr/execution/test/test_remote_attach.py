import sys
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from ..base import ExecutionSpec, LaunchSpec, RunState
from ..manager import RunRecord
from ..run_model import Execution, WorkingDirectory, FileManifest
from ..endpoint import endpoint_for_executor
from ..postprocessing import PostprocessingRunner
from ..remote_endpoint import RemoteEndpointExecutor
from ..transport import LocalTransport


def _executor(*, backend="local"):
    return RemoteEndpointExecutor(
        LocalTransport(),
        backend=backend,
        command=(sys.executable, "-m", "ase2sprkkr.tools.ase2sprkkr", "remote"),
    )


def _launch(local_root, remote_root):
    return LaunchSpec(local_root=local_root, remote_root=remote_root)


def _record(runner):
    return RunRecord(
        id=runner.backend_id, task="SCF", state=RunState.RUNNING,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=runner.execution,
        files=runner.files,
    )


def _processed(raw):
    def process(outcome):
        return SimpleNamespace(
            return_code=outcome.return_code,
            output=(
                b""
                if outcome.output_file is None
                else Path(outcome.output_file).read_bytes()
            ),
            error=(
                b""
                if outcome.error_file is None
                else Path(outcome.error_file).read_bytes()
            ),
        )

    runner = PostprocessingRunner(raw, process)
    runner.start()
    return runner


def _wait_for_state(runner, expected, timeout=5):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if runner.state == expected:
            return
        time.sleep(0.01)
    pytest.fail(f"remote run did not reach {expected.value}; current state={runner.state.value}")


def test_remote_executor_attach_reconnects_running_run_and_restores_stage_out(tmp_path):
    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)

    input_file = (workdir / "run.inp").open("w+b")
    input_file.write(b"input")
    input_file.seek(0)
    output_file = (workdir / "run.out").open("wb")
    generated = workdir / "generated.dat"
    cache = workdir / "cache.dat"
    cache.write_text("before")

    first_executor = _executor()
    spec = ExecutionSpec(
        argv=(
            sys.executable,
            "-c",
            (
                "import pathlib,time\n"
                "release=pathlib.Path('release')\n"
                "while not release.exists():\n"
                "    time.sleep(0.01)\n"
                "pathlib.Path('generated.dat').write_text('generated')\n"
                "pathlib.Path('cache.dat').write_text('after')\n"
                "print('done')\n"
            ),
        ),
        input_file=input_file,
        output_file=output_file,
        directory=str(workdir),
        stdin=input_file,
        input_files=(input_file,),
        output_files=(output_file, generated),
        cache_files=(cache,),
        launch=_launch(local_root, remote_root),
    )
    original_raw = first_executor.create_run(spec)
    original = PostprocessingRunner(original_raw, lambda outcome: outcome)
    original.submit()
    _wait_for_state(original, RunState.RUNNING)

    run = _record(original)
    assert run.execution.backend_id is not None

    # A fresh endpoint object reconstructs per-run staging policy solely from
    # persisted run manifest; no task reader or live runner state is used.
    attached = _processed(_executor().attach(run, start=False))
    (remote_root / "calculation" / "release").write_text("go")

    assert attached.backend_id == run.execution.backend_id
    result = attached.wait()

    assert attached.state == RunState.FINISHED
    assert result.return_code == 0
    assert result.output == b"done\n"
    assert (workdir / "run.out").read_bytes() == b"done\n"
    assert (workdir / "generated.dat").read_text() == "generated"
    assert (workdir / "cache.dat").read_text() == "after"

    assert original.wait().return_code == 0


def test_remote_executor_attach_replays_finished_run(tmp_path):
    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    workdir = local_root / "calculation"
    workdir.mkdir(parents=True)

    input_file = (workdir / "run.inp").open("w+b")
    input_file.write(b"input")
    input_file.seek(0)
    output_file = (workdir / "run.out").open("wb")

    first_executor = _executor()
    spec = ExecutionSpec(
        argv=(sys.executable, "-c", "print('finished')"),
        input_file=input_file,
        output_file=output_file,
        directory=str(workdir),
        stdin=input_file,
        input_files=(input_file,),
        output_files=(output_file,),
        launch=_launch(local_root, remote_root),
    )
    original = first_executor.create_run(spec)
    assert original.run().return_code == 0
    run = _record(original)
    assert run.execution.backend_id is not None

    attached = _processed(_executor().attach(run, start=False))
    result = attached.wait()

    assert attached.state == RunState.FINISHED
    assert result.return_code == 0
    assert result.output == b"finished\n"
    assert (workdir / "run.out").read_bytes() == b"finished\n"


def test_remote_executor_attach_rejects_reference_for_other_backend(tmp_path):
    executor = _executor()
    run = RunRecord(
        id="wrong-backend", task="SCF", state=RunState.RUNNING,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            endpoint_for_executor(_executor(backend="psij")),
            WorkingDirectory(str(tmp_path), "/remote/calculation"),
            "job",
        ),
        files=FileManifest(),
    )
    with pytest.raises(ValueError, match="different remote endpoint"):
        executor.attach(run)


def test_remote_attach_preserves_existing_partial_stream(tmp_path, monkeypatch):
    import asyncio
    from ..run_model import FileEntry
    local_root = tmp_path / "local"
    remote_root = tmp_path / "remote"
    local_root.mkdir()
    remote_root.mkdir()
    local_output = local_root / "run.out"
    local_output.write_bytes(b"partial")
    endpoint = _executor()
    context = endpoint._context(endpoint, _launch(local_root, remote_root))
    run = RunRecord(
        id="partial", task="SCF", state=RunState.RUNNING,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(endpoint_for_executor(endpoint),
                            WorkingDirectory(str(local_root), str(remote_root)),
                            "opaque-handle"),
        files=FileManifest(stdout=FileEntry(
            local_path=str(local_output),
            remote_path=str(remote_root / "run.out"),
        )),
    )

    async def status(*args, **kwargs):
        return {"state": "running"}
    monkeypatch.setattr(context, "_status", status)
    runner = asyncio.run(context.attach(run, start=False))
    assert runner.spec.output_file.tell() == len(b"partial")
    runner.spec.output_file.close()
    assert local_output.read_bytes() == b"partial"

"""Crash boundary: remote staging metadata must be persisted before submit."""

import pytest

from ase2sprkkr.execution.remote_endpoint import _RemoteRun
from ase2sprkkr.execution import (
    ExecutionSpec, FileEntry, FileManifest, FileMapping,
    LaunchSpec, RemoteEndpointExecutor, RunState, RunStore, LocalTransport,
)

from ase2sprkkr.execution.manager import ExecutionManager


def test_prepared_manifest_is_durable_before_remote_backend_submit(tmp_path, monkeypatch):
    executor = RemoteEndpointExecutor(LocalTransport())
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    directory = tmp_path / "local"
    directory.mkdir()
    output = directory / "run.out"
    spec = ExecutionSpec(
        argv=("ignored",),
        directory=str(directory),
        input_file=None,
        output_file=str(output),
        launch=LaunchSpec(
            local_root=str(tmp_path),
            remote_root=str(tmp_path / "target"),
        ),
    )
    run = manager.create_run(
        executor, spec, task="DOS", process_result=lambda outcome: outcome,
    )
    expected = FileManifest(
        stdout=FileEntry(str(output), "/target/prepared/run.out", availability_state="pending"),
        mapping=FileMapping(str(tmp_path), "/target/prepared"),
    )

    async def prepare(_runner):
        return {"directory": "/target/prepared", "id": run.id}, expected

    submitted = []

    async def submit(_request):
        current = manager.store.get(run.id)
        submitted.append(current)
        return {"state": "running", "remote_id": "opaque-job"}

    monkeypatch.setattr(_RemoteRun, "_prepare_async", prepare)
    monkeypatch.setattr("ase2sprkkr.execution.remote_endpoint.start_monitor", lambda *args, **kwargs: None)
    executor._submit = submit

    async def status(_directory, _backend_id, _files=None):
        return {"state": "running"}

    # The test substitutes an opaque backend id and disables the monitor;
    # supply the corresponding fake status instead of invoking the CLI.
    executor._status = status
    # This test is exclusively about the pre-submit durability boundary.
    # Avoid starting a background completion worker for a synthetic backend
    # which intentionally has no process to wait for.
    monkeypatch.setattr(run._runner, "_start_completion", lambda: None)
    run.submit()

    assert len(submitted) == 1
    snapshot = submitted[0]
    assert snapshot.state == RunState.SUBMITTING
    assert snapshot.execution.backend_id is None
    assert snapshot.execution.directory.remote == "/target/prepared"
    assert snapshot.files == expected
    assert manager.store.get(run.id).execution.backend_id == "opaque-job"


def test_remote_manifest_uses_same_declarations_as_local(tmp_path):
    import asyncio
    from ase2sprkkr.execution.run_model import FileManifest

    directory = tmp_path / "local"
    directory.mkdir(exist_ok=True)
    remote = tmp_path / "target"
    spec = ExecutionSpec(
        argv=("ignored",), directory=str(directory),
        input_file="main.inp", input_files=("main.inp", "main.inp", "extra.inp", "extra.inp"),
        output_file="main.out", error_file="main.err",
        output_files=("main.out", "main.err", "one.dat", "one.dat", "two.dat"),
        cache_files=("restart.pot", "restart.pot"),
        launch=LaunchSpec(local_root=str(tmp_path), remote_root=str(remote)),
    )
    local = FileManifest.from_spec(spec)
    executor = RemoteEndpointExecutor(LocalTransport())
    runner = executor.create_run(spec)
    manifest = asyncio.run(runner._submitted_manifest({
        "client_sync": {"output_file": str(remote / "local" / "main.out"),
                        "error_file": str(remote / "local" / "main.err")},
    }))
    assert tuple(e.local_path for e in manifest.inputs) == (
        str(directory / "main.inp"), str(directory / "extra.inp"),
    )
    assert tuple(e.local_path for e in manifest.inputs) == local.input_files
    assert tuple(e.remote_path for e in manifest.inputs) == (
        str(remote / "local" / "main.inp"),
        str(remote / "local" / "extra.inp"),
    )
    assert tuple(e.local_path for e in manifest.outputs) == tuple(e.local_path for e in local.outputs)
    assert tuple(e.local_path for e in manifest.outputs) == (
        str(directory / "one.dat"), str(directory / "two.dat"),
    )
    assert manifest.output_files == local.output_files
    assert manifest.stderr.remote_path == str(remote / "local" / "main.err")
    assert manifest.stdout.remote_path == str(remote / "local" / "main.out")
    assert manifest.outputs[0].remote_path == str(remote / "local" / "one.dat")
    assert manifest.caches[0].file.remote_path == str(remote / "local" / "restart.pot")
    assert len(manifest.caches) == 1


def test_remote_file_preparation_deduplicates_paths_before_mapping(tmp_path):
    """Exact and wildcard declarations share the same normalization rules."""
    import asyncio

    directory = tmp_path / "work"
    directory.mkdir()
    remote = tmp_path / "remote"
    executor = RemoteEndpointExecutor(LocalTransport())
    spec = ExecutionSpec(
        argv=("ignored",), directory=str(directory),
        input_file=None, output_file="stdout.out",
        output_files=(
            None, "stdout.out", "one.dat", "./one.dat",
            "cut-*.dat", "./cut-*.dat",
        ),
        cache_files=(
            None, "restart.dat", "./restart.dat",
            "checkpoint-*.dat", "./checkpoint-*.dat",
        ),
        launch=LaunchSpec(local_root=str(tmp_path), remote_root=str(remote)),
    )
    runner = executor.create_run(spec)

    async def prepare():
        used = {}
        await runner._prepare_cache_files(used)
        await runner._prepare_output_files()
        return used

    used = asyncio.run(prepare())
    assert runner._output_targets == {
        str(directory / "one.dat"): str(remote / "work" / "one.dat"),
    }
    assert runner._output_patterns == (str(remote / "work" / "cut-*.dat"),)
    assert runner._cache_targets == {
        str(directory / "restart.dat"): (
            str(remote / "work" / "restart.dat"), None,
        ),
    }
    assert runner._cache_patterns == (
        str(remote / "work" / "checkpoint-*.dat"),
    )
    assert used == {str(remote / "work" / "restart.dat"): str(directory / "restart.dat")}


@pytest.mark.parametrize("method", ("submit", "start"))
def test_postprocessing_rejected_submit_preserves_error_and_does_not_monitor(method):
    from ase2sprkkr.execution.base import SubmissionNotStarted
    from ase2sprkkr.execution.postprocessing import PostprocessingRunner

    failure = SubmissionNotStarted("preflight refused")

    class RejectedRunner:
        state = RunState.NEW

        def submit(self):
            raise failure

    runner = PostprocessingRunner(RejectedRunner(), lambda result: result)
    with pytest.raises(SubmissionNotStarted, match="preflight refused") as caught:
        getattr(runner, method)()
    assert caught.value is failure
    assert runner._started is False
    assert runner._thread is None

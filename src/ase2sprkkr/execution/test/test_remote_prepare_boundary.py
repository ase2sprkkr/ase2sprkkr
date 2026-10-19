"""Crash boundary: remote staging metadata must be persisted before submit."""

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
    monkeypatch.setattr(_RemoteRun, "_start_monitor", lambda _runner: None)
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


def test_remote_manifest_deduplicates_declared_inputs_and_outputs(tmp_path):
    import asyncio

    directory = tmp_path / "local"
    directory.mkdir()
    executor = RemoteEndpointExecutor(LocalTransport())
    remote = tmp_path / "remote"
    inp = directory / "calc.inp"
    out = directory / "calc.out"
    auxiliary = directory / "dos.dat"
    spec = ExecutionSpec(
        argv=("ignored",), directory=str(directory),
        input_file=str(inp), input_files=(str(inp), str(inp)),
        output_file=str(out), output_files=(str(out), str(auxiliary), str(auxiliary)),
        launch=LaunchSpec(local_root=str(tmp_path), remote_root=str(remote)),
    )
    runner = executor.create_run(spec)

    async def check():
        request = {
            "client_sync": {
                "output_file": str(remote / "local" / "calc.out"),
            }
        }
        return await runner._submitted_manifest(request)

    manifest = asyncio.run(check())
    assert [entry.local_path for entry in manifest.inputs] == [str(inp)]
    assert [entry.local_path for entry in manifest.outputs] == [str(auxiliary)]
    assert manifest.stdout.local_path == str(out)

"""Streamless runs remain attachable for status and artifact retrieval."""

import asyncio
from types import SimpleNamespace

from ..base import LaunchSpec, RunState
from ..endpoint import ExecutorEndpoint
from ..local_project import attach_local
from ..run_model import Execution, FileManifest, WorkingDirectory
from ..remote_endpoint import RemoteEndpointExecutor
from ..transport import LocalTransport


def test_local_attach_without_replay_streams(monkeypatch, tmp_path):
    from .. import local_project
    monkeypatch.setattr(local_project, "_decode_handle", lambda _: {"kind": "process"})
    monkeypatch.setattr(
        local_project, "local_status",
        lambda *_: SimpleNamespace(state=RunState.RUNNING, failure_phase=None,
                                   message=None, return_code=None, files=None),
    )
    run = SimpleNamespace(
        id="run", files=FileManifest(),
        execution=Execution(ExecutorEndpoint(),
                            WorkingDirectory(str(tmp_path)), "handle"),
    )
    runner = attach_local(SimpleNamespace(), run, start=False)
    assert runner.backend_id == "handle"
    assert runner.files == run.files


def test_remote_attach_without_streams(monkeypatch, tmp_path):
    root = tmp_path / "local"
    root.mkdir()
    endpoint = RemoteEndpointExecutor(LocalTransport(), backend="local")
    context = endpoint._context(
        endpoint, LaunchSpec(local_root=str(root),
                             remote_root=str(tmp_path / "remote")),
    )

    async def status(*_):
        return {"state": "running"}
    monkeypatch.setattr(context, "_status", status)
    run = SimpleNamespace(
        files=FileManifest(),
        execution=SimpleNamespace(
            backend_id="handle",
            directory=WorkingDirectory(str(root), str(tmp_path / "remote")),
        ),
    )
    raw = asyncio.run(context.attach(run, start=False))
    assert raw.backend_id == "handle"
    assert raw.spec.output_file is None
    assert raw.spec.error_file is None

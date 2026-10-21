import asyncio

import pytest

from ase2sprkkr.execution import Executor, LaunchSpec
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor


class _Transport:
    host = "cluster"



def test_run_context_is_staging_policy_not_executor(tmp_path):
    endpoint = RemoteEndpointExecutor(_Transport())
    context = endpoint._context(
        endpoint,
        LaunchSpec(
            local_root=tmp_path,
            remote_root="/scratch/a",
            input_sync="always",
            output_sync="missing",
        ),
    )

    assert not isinstance(context, Executor)
    assert context.endpoint is endpoint
    assert context.transport is endpoint.transport
    assert context.local_root == tmp_path.resolve()
    assert context.remote_root == "/scratch/a"



def test_different_run_staging_policies_share_one_endpoint(tmp_path):
    endpoint = RemoteEndpointExecutor(_Transport())
    first = endpoint._context(
        endpoint,
        LaunchSpec(local_root=tmp_path, remote_root="/scratch/one"),
    )
    second = endpoint._context(
        endpoint,
        LaunchSpec(local_root=tmp_path, remote_root="/scratch/two"),
    )

    assert first is not second
    assert first.endpoint is endpoint
    assert second.endpoint is endpoint
    assert first.remote_root == "/scratch/one"
    assert second.remote_root == "/scratch/two"


def test_attach_restores_typed_manifest_without_sync_dictionary(tmp_path, monkeypatch):
    from ase2sprkkr.execution import RunState
    from ase2sprkkr.execution.endpoint import endpoint_for_executor
    from ase2sprkkr.execution.manager import RunRecord
    from ase2sprkkr.execution.run_model import (
        FileEntry, FileManifest, CacheEntry, Execution, WorkingDirectory,
    )
    from ase2sprkkr.execution.transport import LocalTransport

    local = tmp_path / "local"
    remote = tmp_path / "remote"
    local.mkdir()
    remote.mkdir()

    def entry(name, **kwargs):
        return FileEntry(str(local / name), str(remote / name), **kwargs)

    files = FileManifest(
        stdout=entry("stdout.txt"),
        outputs=(
            entry("exact.out"),
            entry("data-*.out"),
            entry("data-a.out", source_pattern=str(local / "data-*.out")),
            FileEntry(str(local / "local-only.out")),
        ),
        caches=(
            CacheEntry(entry("cache.dat"), (2, 3)),
            CacheEntry(entry("restart-*.pot")),
            CacheEntry(entry("restart-a.pot",
                             source_pattern=str(local / "restart-*.pot"))),
            CacheEntry(FileEntry(remote_path=str(remote / "remote-only.dat")), (5, 7)),
            CacheEntry(FileEntry(str(local / "local-only-cache.dat"))),
        ),
    )
    endpoint = RemoteEndpointExecutor(LocalTransport())
    context = endpoint._context(endpoint, LaunchSpec(local_root=local, remote_root=remote))
    run = RunRecord(
        id="test-manifest", task="SCF", state=RunState.RUNNING,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            endpoint_for_executor(endpoint), WorkingDirectory(str(local), str(remote)),
            "opaque-id",
        ),
        files=files,
    )

    async def status(*args, **kwargs):
        return {"state": "running"}

    monkeypatch.setattr(context, "_status", status)
    runner = asyncio.run(context.attach(run, start=False))
    try:
        assert runner._output_targets == {
            str(local / "exact.out"): str(remote / "exact.out"),
        }
        assert runner._output_patterns == (str(remote / "data-*.out"),)
        assert runner._cache_targets == {
            str(local / "cache.dat"): (str(remote / "cache.dat"), (2, 3)),
        }
        assert runner._cache_patterns == (str(remote / "restart-*.pot"),)
        assert runner._cache_before == {
            str(remote / "remote-only.dat"): (5, 7),
        }
    finally:
        runner.spec.output_file.close()


def test_absolute_remote_root_needs_no_remote_home_query(tmp_path, monkeypatch):
    endpoint = RemoteEndpointExecutor(_Transport())
    context = endpoint._context(
        endpoint,
        LaunchSpec(local_root=tmp_path, remote_root="/scratch/../work"),
    )
    calls = []

    async def remote_home():
        calls.append("home")
        raise AssertionError("Absolute remote roots must not query remote HOME")

    monkeypatch.setattr(endpoint, "_remote_home", remote_home)
    directory = tmp_path / "calculation"
    directory.mkdir()
    assert asyncio.run(context.map_path(directory)) == "/work/calculation"
    assert asyncio.run(context.map_directory(directory)) == "/work/calculation"
    assert asyncio.run(context.unmap_path("/work/calculation")) == directory
    assert asyncio.run(context._resolved_remote_root()) == "/work"
    assert calls == []


@pytest.mark.parametrize(
    "root, expected",
    [
        (None, "/home/user"),
        ("~", "/home/user"),
        ("~/projects", "/home/user/projects"),
        ("projects", "/home/user/projects"),
    ],
)
def test_home_relative_remote_root_still_resolves_home_once(
    tmp_path, monkeypatch, root, expected,
):
    endpoint = RemoteEndpointExecutor(_Transport())
    context = endpoint._context(
        endpoint, LaunchSpec(local_root=tmp_path, remote_root=root),
    )
    calls = []

    async def remote_home():
        calls.append("home")
        return "/home/user"

    monkeypatch.setattr(endpoint, "_remote_home", remote_home)
    assert asyncio.run(context._resolved_remote_root()) == expected
    assert calls == ["home"]


def test_attach_closes_stdout_if_opening_stderr_fails(tmp_path, monkeypatch):
    """A failed second open must release the first stream immediately."""
    import builtins
    from ase2sprkkr.execution import RunState
    from ase2sprkkr.execution.endpoint import endpoint_for_executor
    from ase2sprkkr.execution.manager import RunRecord
    from ase2sprkkr.execution.run_model import (
        Execution, FileEntry, FileManifest, WorkingDirectory,
    )
    from ase2sprkkr.execution.transport import LocalTransport
    from ase2sprkkr.execution import remote_endpoint

    local, remote = tmp_path / "local", tmp_path / "remote"
    local.mkdir()
    remote.mkdir()
    stdout = local / "stdout.txt"
    stderr = local / "stderr.txt"
    endpoint = RemoteEndpointExecutor(LocalTransport())
    context = endpoint._context(endpoint, LaunchSpec(
        local_root=local, remote_root=remote,
    ))
    run = RunRecord(
        id="stream-open-failure", task="SCF", state=RunState.RUNNING,
        created_at="2026-10-13T08:00:00+00:00",
        updated_at="2026-10-13T08:00:00+00:00",
        execution=Execution(
            endpoint_for_executor(endpoint),
            WorkingDirectory(str(local), str(remote)), "remote-id",
        ),
        files=FileManifest(
            stdout=FileEntry(str(stdout), str(remote / "stdout.txt")),
            stderr=FileEntry(str(stderr), str(remote / "stderr.txt")),
        ),
    )

    async def status(*args, **kwargs):
        return {"state": "running"}

    monkeypatch.setattr(context, "_status", status)
    opened = []
    actual_open = builtins.open

    def failing_open(path, *args, **kwargs):
        if str(path) == str(stderr):
            raise PermissionError("stderr not writable")
        stream = actual_open(path, *args, **kwargs)
        if str(path) == str(stdout):
            opened.append(stream)
        return stream

    monkeypatch.setattr(remote_endpoint, "open", failing_open, raising=False)
    with pytest.raises(PermissionError, match="stderr not writable"):
        asyncio.run(context.attach(run, start=False))
    assert len(opened) == 1
    assert opened[0].closed


def test_unmap_rejects_symlink_parent_outside_local_root(tmp_path):
    root = tmp_path / "local"
    outside = tmp_path / "outside"
    root.mkdir()
    outside.mkdir()
    (root / "redirect").symlink_to(outside, target_is_directory=True)
    endpoint = RemoteEndpointExecutor(_Transport())
    context = endpoint._context(endpoint, LaunchSpec(local_root=root, remote_root="/remote"))
    with pytest.raises(ValueError, match="outside mapped root"):
        asyncio.run(context.unmap_path("/remote/redirect/output.dat"))


def test_unmap_preserves_internal_symlink_spelling(tmp_path):
    root = tmp_path / "local"
    actual = root / "actual"
    actual.mkdir(parents=True)
    (root / "alias").symlink_to(actual, target_is_directory=True)
    endpoint = RemoteEndpointExecutor(_Transport())
    context = endpoint._context(endpoint, LaunchSpec(local_root=root, remote_root="/remote"))
    assert asyncio.run(context.unmap_path("/remote/alias/output.dat")) == (
        root / "alias" / "output.dat"
    )

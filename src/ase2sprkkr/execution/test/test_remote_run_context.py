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


def test_manifest_sync_groups_declarations_without_replaying_pattern_children():
    from ase2sprkkr.execution.run_model import FileEntry, FileManifest, CacheEntry
    from ase2sprkkr.execution.remote_endpoint import _RemoteRunContext

    manifest = FileManifest(
        outputs=(
            FileEntry("/local/exact.out", "/remote/exact.out"),
            FileEntry("/local/data-*.out", "/remote/data-*.out"),
            FileEntry("/local/data-a.out", "/remote/data-a.out",
                      source_pattern="/local/data-*.out"),
            FileEntry("/local/local-only.out"),
        ),
        caches=(
            CacheEntry(FileEntry("/local/cache.dat", "/remote/cache.dat"), (2, 3)),
            CacheEntry(FileEntry("/local/restart-*.pot", "/remote/restart-*.pot")),
            CacheEntry(FileEntry("/local/restart-a.pot", "/remote/restart-a.pot",
                                 source_pattern="/local/restart-*.pot")),
            CacheEntry(FileEntry(remote_path="/remote/remote-only.dat"), (5, 7)),
            CacheEntry(FileEntry(local_path="/local/local-only-cache.dat")),
        ),
    )
    assert _RemoteRunContext._manifest_sync(manifest) == {
        "output_files": ["/remote/exact.out"],
        "output_patterns": ["/remote/data-*.out"],
        "cache_files": [{"path": "/remote/cache.dat", "before": (2, 3)}],
        "cache_patterns": ["/remote/restart-*.pot"],
        "cache_before": [{"path": "/remote/remote-only.dat", "signature": (5, 7)}],
    }


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

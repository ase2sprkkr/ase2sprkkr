import asyncio

import pytest

from ase2sprkkr.execution import ExecutionSpec, LaunchSpec
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor


class _Transport:
    host = "cluster"


def test_remote_submission_directory_uses_target_namespace_without_io(tmp_path):
    executor = RemoteEndpointExecutor(_Transport())
    directory = tmp_path / "project" / "case"
    directory.mkdir(parents=True)

    first = ExecutionSpec(
        argv=("true",),
        input_file=None,
        output_file=None,
        directory=directory,
        launch=LaunchSpec(local_root=tmp_path, remote_root="/scratch/a"),
    )
    second = ExecutionSpec(
        argv=("true",),
        input_file=None,
        output_file=None,
        directory=directory,
        launch=LaunchSpec(local_root=tmp_path, remote_root="/scratch/b"),
    )

    assert executor.submission_directory(first) == "/scratch/a/project/case"
    assert executor.submission_directory(second) == "/scratch/b/project/case"


def test_remote_home_variants_have_same_symbolic_namespace(tmp_path):
    executor = RemoteEndpointExecutor(_Transport())
    directory = tmp_path / "case"
    directory.mkdir()

    implicit = ExecutionSpec(
        argv=("true",),
        input_file=None,
        output_file=None,
        directory=directory,
        launch=LaunchSpec(local_root=tmp_path),
    )
    explicit = ExecutionSpec(
        argv=("true",),
        input_file=None,
        output_file=None,
        directory=directory,
        launch=LaunchSpec(local_root=tmp_path, remote_root="~"),
    )

    assert executor.submission_directory(implicit) == "~/case"
    assert executor.submission_directory(explicit) == "~/case"


@pytest.mark.parametrize("name", ("case", "case-alias"))
def test_context_and_submission_directory_share_local_path_resolution(
    tmp_path, name,
):
    local_root = tmp_path / "workspace"
    directory = local_root / "case"
    directory.mkdir(parents=True)
    (local_root / "case-alias").symlink_to(directory, target_is_directory=True)

    executor = RemoteEndpointExecutor(_Transport())
    launch = LaunchSpec(local_root=local_root, remote_root="/scratch/project")
    spec = ExecutionSpec(
        argv=("true",), directory=local_root / name,
        input_file=None, output_file=None, launch=launch,
    )
    context = executor._context(executor, launch)

    async def home():
        return "/home/cluster"

    executor._remote_home = home
    expected = "/scratch/project/case"
    assert executor.submission_directory(spec) == expected
    assert asyncio.run(context.map_path(local_root / name)) == expected


def test_context_and_submission_directory_reject_symlink_escaping_root(
    tmp_path,
):
    local_root = tmp_path / "workspace"
    local_root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (local_root / "escape").symlink_to(
        outside, target_is_directory=True,
    )
    executor = RemoteEndpointExecutor(_Transport())
    launch = LaunchSpec(local_root=local_root, remote_root="/scratch/project")
    spec = ExecutionSpec(
        argv=("true",), directory=local_root / "escape",
        input_file=None, output_file=None, launch=launch,
    )
    context = executor._context(executor, launch)
    with pytest.raises(ValueError, match="outside mapped root"):
        executor.submission_directory(spec)
    with pytest.raises(ValueError, match="outside mapped root"):
        asyncio.run(context.map_path(local_root / "escape"))


def test_context_and_submission_directory_map_root_itself(tmp_path):
    executor = RemoteEndpointExecutor(_Transport())
    launch = LaunchSpec(
        local_root=tmp_path, remote_root="/scratch/root",
    )
    spec = ExecutionSpec(
        argv=("true",), directory=tmp_path,
        input_file=None, output_file=None, launch=launch,
    )
    context = executor._context(executor, launch)

    async def home():
        return "/home/cluster"

    executor._remote_home = home
    assert executor.submission_directory(spec) == "/scratch/root"
    assert asyncio.run(context.map_path(tmp_path)) == "/scratch/root"

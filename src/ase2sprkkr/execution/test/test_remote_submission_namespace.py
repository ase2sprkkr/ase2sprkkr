import asyncio

import pytest

from ase2sprkkr.execution import ExecutionSpec, LaunchSpec
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor


class _Transport:
    host = "cluster"


@pytest.mark.parametrize(
    ("relative", "remote_root", "expected"),
    (
        ("project/case", "/scratch/a", "/scratch/a/project/case"),
        ("project/case", "/scratch/b", "/scratch/b/project/case"),
        ("case", None, "~/case"),
        ("case", "~", "~/case"),
    ),
    ids=("target-a", "target-b", "implicit-home", "explicit-home"),
)
def test_remote_submission_directory_uses_target_namespace_without_io(
    tmp_path, relative, remote_root, expected,
):
    executor = RemoteEndpointExecutor(_Transport())
    directory = tmp_path / relative
    directory.mkdir(parents=True)

    options = {} if remote_root is None else {"remote_root": remote_root}
    spec = ExecutionSpec(
        argv=("true",),
        input_file=None,
        output_file=None,
        directory=directory,
        launch=LaunchSpec(local_root=tmp_path, **options),
    )
    assert executor.submission_directory(spec) == expected


@pytest.mark.parametrize(
    ("name", "suffix"),
    (("case", "/case"), ("case-alias", "/case"), (".", "")),
    ids=("regular", "symlink-alias", "root"),
)
def test_context_and_submission_directory_share_local_path_resolution(
    tmp_path, name, suffix,
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
    expected = "/scratch/project" + suffix
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

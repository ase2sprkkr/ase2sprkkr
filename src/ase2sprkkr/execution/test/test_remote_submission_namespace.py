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

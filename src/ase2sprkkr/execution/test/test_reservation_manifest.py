"""Persistent reservation must use the same remote paths as staging."""

from types import SimpleNamespace

from ase2sprkkr.execution import (
    ExecutionSpec, FileManifest, LaunchSpec, RemoteEndpointExecutor,
)


class _Transport:
    def __init__(self):
        self.commands = []

    async def run(self, argv, **kwargs):
        self.commands.append(tuple(argv))
        return SimpleNamespace(stdout=b"/home/test")


def test_remote_manifest_resolves_files_without_staging(tmp_path):
    transport = _Transport()
    executor = RemoteEndpointExecutor(transport)
    spec = ExecutionSpec(
        argv=("true",), input_file="main.in", output_file="main.out",
        directory=str(tmp_path / "work"), cache_files=("TAU",),
        launch=LaunchSpec(local_root=str(tmp_path), remote_root="/remote"),
    )
    files = executor.reservation_manifest(spec, FileManifest.from_spec(spec))
    assert files.inputs[0].remote_path == "/remote/work/main.in"
    assert files.stdout.remote_path == "/remote/work/main.out"
    assert files.caches[0].file.remote_path == "/remote/work/TAU"
    assert files.mapping.remote_root == "/remote"
    assert all(command[:1] != ("mkdir",) for command in transport.commands)


def test_different_local_roots_map_to_identical_remote_target(tmp_path):
    executor = RemoteEndpointExecutor(_Transport())

    def manifest(root):
        spec = ExecutionSpec(
            argv=("true",), input_file=None, output_file="result.out",
            directory=str(root / "work"),
            launch=LaunchSpec(local_root=str(root), remote_root="/shared"),
        )
        return executor.reservation_manifest(spec, FileManifest.from_spec(spec))

    first = manifest(tmp_path / "first")
    second = manifest(tmp_path / "second")
    assert first.stdout.local_path != second.stdout.local_path
    assert first.stdout.remote_path == second.stdout.remote_path == "/shared/work/result.out"

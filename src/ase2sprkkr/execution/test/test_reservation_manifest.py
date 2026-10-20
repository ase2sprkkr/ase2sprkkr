"""Persistent reservation must use the same remote paths as staging."""

from types import SimpleNamespace

import pytest

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
    files = executor.reservation_context(spec, FileManifest.from_spec(spec))[0]
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
        return executor.reservation_context(spec, FileManifest.from_spec(spec))[0]

    first = manifest(tmp_path / "first")
    second = manifest(tmp_path / "second")
    assert first.stdout.local_path != second.stdout.local_path
    assert first.stdout.remote_path == second.stdout.remote_path == "/shared/work/result.out"



def test_reservation_maps_stdin_for_target_side_conflicts(tmp_path):
    executor = RemoteEndpointExecutor(_Transport())
    spec = ExecutionSpec(
        argv=("cat",), input_file=None, output_file="result.out",
        directory=str(tmp_path / "work"), stdin="shared.in",
        launch=LaunchSpec(local_root=str(tmp_path), remote_root="/remote"),
    )
    files = executor.reservation_context(spec, FileManifest.from_spec(spec))[0]
    assert files.input_file is None
    assert files.stdin.local_path == str(tmp_path / "work" / "shared.in")
    assert files.stdin.remote_path == "/remote/work/shared.in"


def test_cache_mapping_yields_real_cache_entries(tmp_path):
    executor = RemoteEndpointExecutor(_Transport())
    spec = ExecutionSpec(
        argv=("true",), input_file=None, output_file=None,
        directory=str(tmp_path / "work"), cache_files=("TAU",),
        launch=LaunchSpec(local_root=str(tmp_path), remote_root="/remote"),
    )
    manifest = executor.reservation_context(spec, FileManifest.from_spec(spec))[0]
    assert len(manifest.caches) == 1
    assert manifest.caches[0].file.remote_path == "/remote/work/TAU"


def test_create_run_is_offline_and_resolves_home_only_when_submitted(tmp_path, monkeypatch):
    from ase2sprkkr.execution import ExecutionManager, RunStore

    transport = _Transport()
    executor = RemoteEndpointExecutor(transport)
    manager = ExecutionManager(RunStore(tmp_path / "runs"))
    spec = ExecutionSpec(
        argv=("true",), input_file="primary.in", output_file="result.out",
        directory=str(tmp_path / "case"),
        launch=LaunchSpec(local_root=str(tmp_path), remote_root="~"),
    )
    run = manager.create_run(
        executor, spec, task="DOS", process_result=lambda result: result,
    )
    assert transport.commands == []
    assert run.record.execution.directory.remote is None
    assert run.record.files.mapping is None

    captured = []

    def intercept_reservation(record):
        captured.append(record)
        raise RuntimeError("stop before staging")

    monkeypatch.setattr(manager, "reserve", intercept_reservation)
    with pytest.raises(RuntimeError, match="stop before staging"):
        run.submit()
    assert captured[0].execution.directory.remote == "/home/test/case"
    assert captured[0].files.inputs[0].remote_path == "/home/test/case/primary.in"
    assert transport.commands


def test_tilde_and_absolute_roots_protect_same_scf_directory(tmp_path):
    from ase2sprkkr.execution import (
        DuplicateSubmissionError, ExecutionManager, RunStore,
    )
    from ase2sprkkr.execution.endpoint import endpoint_for_executor

    executor = RemoteEndpointExecutor(_Transport())
    manager = ExecutionManager(RunStore(tmp_path / "runs"))

    def record(task, root):
        spec = ExecutionSpec(
            argv=("true",), input_file=None, output_file=task + ".out",
            directory=str(tmp_path / "case"),
            launch=LaunchSpec(local_root=str(tmp_path), remote_root=root),
        )
        files, directory = executor.reservation_context(
            spec, FileManifest.from_spec(spec),
        )
        return manager.new_record(
            task=task, directory=spec.directory, remote_directory=directory,
            endpoint=endpoint_for_executor(executor), files=files,
        )

    first, second = record("SCF", "~"), record("DOS", "/home/test")
    assert first.execution.directory.remote == second.execution.directory.remote
    manager.reserve(first)
    with pytest.raises(DuplicateSubmissionError):
        manager.reserve(second)


@pytest.mark.parametrize("policy,conflicts", [
    ("changed", True), ("always", True), ("missing", True), ("never", False),
])
@pytest.mark.parametrize("role", ("input_file", "input_files", "stdin"))
def test_remote_uploads_require_write_reservation(tmp_path, policy, conflicts, role):
    from ase2sprkkr.execution import (
        DuplicateSubmissionError, ExecutionManager, RunStore,
    )
    from ase2sprkkr.execution.endpoint import endpoint_for_executor

    executor = RemoteEndpointExecutor(_Transport())
    manager = ExecutionManager(RunStore(tmp_path / "runs"))

    def record(task, subdir):
        local_root = tmp_path / subdir
        sources = {"input_file": None, "input_files": (), "stdin": None}
        sources[role] = ("shared.in",) if role == "input_files" else "shared.in"
        spec = ExecutionSpec(
            argv=("true",), output_file=task + ".out",
            directory=str(local_root / "case"), **sources,
            launch=LaunchSpec(
                local_root=str(local_root), remote_root="/remote",
                input_sync=policy,
            ),
        )
        files, directory = executor.reservation_context(
            spec, FileManifest.from_spec(spec),
        )
        return manager.new_record(
            task=task, directory=spec.directory, remote_directory=directory,
            endpoint=endpoint_for_executor(executor), files=files,
        )

    manager.reserve(record("DOS", "left"))
    second = record("BSF", "right")
    if conflicts:
        with pytest.raises(DuplicateSubmissionError):
            manager.reserve(second)
    else:
        manager.reserve(second)

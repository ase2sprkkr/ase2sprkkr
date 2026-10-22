import asyncio

from ase2sprkkr.execution import (
    ExecutionSpec,
    LaunchSpec,
    LocalTransport,
    MpiSpec,
    RemoteExecutor,
    FileManifest,
)


def test_remote_request_carries_target_executor_options(tmp_path):
    local_root = tmp_path / "local"
    directory = local_root / "calculation"
    directory.mkdir(parents=True)
    input_file = (directory / "run.inp").open("w+b")
    input_file.write(b"input")
    input_file.flush()
    potential = directory / "start.pot"
    potential.write_text("potential")
    output_file = (directory / "run.out").open("wb")

    executor = RemoteExecutor(
        LocalTransport(),
        backend="psij",
        backend_options={"backend": "slurm", "url": None},
    )
    launch = LaunchSpec(
        mpi=MpiSpec(required=True, autodetect=True, process_count=16),
        local_root=local_root,
        remote_root=tmp_path / "remote",
    )
    spec = ExecutionSpec(
        argv=("program",),
        mpi_argv=("programMPI", input_file.name),
        launch=launch,
        input_file=input_file,
        output_file=output_file,
        directory=str(directory),
        stdin=input_file,
        input_files=(input_file, potential),
    )
    runner = executor.create_run(spec)
    context = runner.executor

    remote_directory = asyncio.run(context.map_directory(directory))
    request = asyncio.run(runner._prepare_request("test-run", remote_directory, FileManifest.from_spec(spec)))

    assert request["executor"] == "psij"
    assert request["executor_options"] == {"backend": "slurm", "url": None}
    assert request["argv"] == ["program"]
    assert request["mpi_uses_stdin"] is True
    assert request["mpi_argv"][0] == "programMPI"
    assert request["mpi_argv"][1].endswith("/calculation/run.inp")
    assert request["launch"] == {
        "mpi": {
            "required": True,
            "autodetect": True,
            "launcher": (),
            "process_count": 16,
            "warning": True,
        },
        "local_root": str(local_root),
        "remote_root": str(tmp_path / "remote"),
        "input_sync": "changed",
        "output_sync": "changed",
    }
    remote_input = tmp_path / "remote" / "calculation" / "run.inp"
    remote_potential = tmp_path / "remote" / "calculation" / "start.pot"
    assert "output_files" not in request
    assert "cache_files" not in request
    assert request["client_sync"]["output_files"] == []
    assert request["client_sync"]["cache_files"] == []
    assert request["input_file"] == str(remote_input)
    assert remote_input.read_bytes() == b"input"
    assert remote_potential.read_text() == "potential"

    input_file.close()
    output_file.close()

def test_nonproject_remote_request_does_not_publish_project_identity(tmp_path):
    root = tmp_path / "local"
    directory = root / "work"
    directory.mkdir(parents=True)
    executor = RemoteExecutor(LocalTransport(), backend="psij")
    runner = executor.create_run(ExecutionSpec(
        argv=("program",), directory=str(directory), input_file=None,
        output_file=None, launch=LaunchSpec(
            local_root=root, remote_root=tmp_path / "remote",
        ),
    ))
    request, _manifest = asyncio.run(runner._prepare_async())
    assert request["id"] is None


def test_remote_sync_declarations_have_one_wire_authority(tmp_path):
    from ase2sprkkr.execution.remote_service import _request_spec

    root = tmp_path / "local"
    directory = root / "work"
    directory.mkdir(parents=True)
    output = directory / "run.out"
    auxiliary = directory / "plot-*.dat"
    cache = directory / "restart.dat"
    cache.write_text("restart")
    executor = RemoteExecutor(LocalTransport())
    spec = ExecutionSpec(
        argv=("program",), directory=str(directory), input_file=None,
        output_file=str(output),
        output_files=(auxiliary,),
        cache_files=(cache,),
        launch=LaunchSpec(local_root=root, remote_root=tmp_path / "remote"),
    )
    run = executor.create_run(spec)
    remote_dir = tmp_path / "remote" / "work"
    request = asyncio.run(
        run._prepare_request("test-run", str(remote_dir), FileManifest.from_spec(spec))
    )
    assert "output_files" not in request
    assert "cache_files" not in request
    restored = _request_spec(request)
    assert str(remote_dir / "plot-*.dat") in restored.output_files
    assert str(remote_dir / "restart.dat") in restored.cache_files

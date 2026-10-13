import asyncio
from types import SimpleNamespace

from ase2sprkkr.execution import (
    ExecutionSpec,
    LaunchSpec,
    LocalTransport,
    RemoteExecutor,
)


class _Reader:
    print_output = False

    def create_result(self, output_file, input_file):
        return SimpleNamespace()


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
        mpi=True,
        process_count=16,
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
    runner = executor.create_run(_Reader(), spec)
    context = runner.executor

    remote_directory = asyncio.run(context.map_directory(directory))
    request = asyncio.run(runner._prepare_request("test-run", remote_directory))

    assert request["executor"] == "psij"
    assert request["executor_options"] == {"backend": "slurm", "url": None}
    assert request["argv"] == ["program"]
    assert request["mpi_argv"][0] == "programMPI"
    assert request["mpi_argv"][1].endswith("/calculation/run.inp")
    assert request["launch"] == {
        "mpi": True,
        "process_count": 16,
        "mpi_warning": True,
        "local_root": str(local_root),
        "remote_root": str(tmp_path / "remote"),
        "input_sync": "changed",
        "output_sync": "changed",
    }
    assert request["client_sync"]["input_files"] == [
        str(tmp_path / "remote" / "calculation" / "run.inp"),
        str(tmp_path / "remote" / "calculation" / "start.pot"),
    ]

    input_file.close()
    output_file.close()
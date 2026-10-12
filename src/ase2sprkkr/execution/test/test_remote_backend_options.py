import asyncio
from types import SimpleNamespace

from ase2sprkkr.execution import ExecutionSpec, LocalTransport, RemoteExecutor


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
    output_file = (directory / "run.out").open("wb")

    executor = RemoteExecutor(
        LocalTransport(),
        backend="psij",
        backend_options={"backend": "slurm", "url": None},
        local_root=local_root,
        remote_root=tmp_path / "remote",
    )
    spec = ExecutionSpec(
        argv=("program",),
        input_file=input_file,
        output_file=output_file,
        directory=str(directory),
        stdin=input_file,
    )
    runner = executor.create_run(_Reader(), spec)

    remote_directory = asyncio.run(executor.map_directory(directory))
    request = asyncio.run(runner._prepare_request("test-run", remote_directory))

    assert request["executor"] == "psij"
    assert request["executor_options"] == {"backend": "slurm", "url": None}

    input_file.close()
    output_file.close()

from types import SimpleNamespace

from ...execution import ExecutionSpec, LocalExecutor
from ..profile import SPRKKRProfile


class _Reader:
    def __init__(self):
        self.call = None

    def create_process(self, command, output_file, **kwargs):
        self.call = command, output_file, kwargs
        return "run"


def test_profile_builds_serial_command():
    input_file = SimpleNamespace(name="/tmp/Fe.inp")
    output_file = object()
    profile = SPRKKRProfile(
        executable_dir=False,
        executable_suffix="",
        mpi=False,
    )

    spec = profile.create_execution(
        "kkrscf",
        input_file,
        output_file,
        "/tmp",
        mpi_supported=True,
    )

    assert spec.argv == ("kkrscf",)
    assert spec.stdin is input_file
    assert spec.input_file is input_file
    assert spec.output_file is output_file


def test_profile_builds_mpi_command():
    input_file = SimpleNamespace(name="/tmp/Fe.inp")
    profile = SPRKKRProfile(
        executable_dir="/opt/sprkkr",
        executable_suffix="_v97",
        mpi=["mpirun", "-np", "4"],
    )

    spec = profile.create_execution(
        "kkrscf",
        input_file,
        object(),
        "/tmp",
        mpi_supported=True,
    )

    assert spec.argv == (
        "mpirun",
        "-np",
        "4",
        "/opt/sprkkr/kkrscf_v97MPI",
        "/tmp/Fe.inp",
    )
    assert spec.stdin is None


def test_non_mpi_task_does_not_use_profile_mpi():
    input_file = SimpleNamespace(name="/tmp/Fe.inp")
    profile = SPRKKRProfile(
        executable_dir=False,
        executable_suffix="",
        mpi=["mpirun", "-np", "4"],
    )

    spec = profile.create_execution(
        "kkrgen",
        input_file,
        object(),
        "/tmp",
        mpi_supported=False,
    )

    assert spec.argv == ("kkrgen",)
    assert spec.stdin is input_file


def test_local_executor_does_not_modify_command():
    reader = _Reader()
    input_file = SimpleNamespace(name="/tmp/Fe.inp")
    output_file = object()
    spec = ExecutionSpec(
        argv=("custom-launcher", "--flag", "program"),
        input_file=input_file,
        output_file=output_file,
        directory="/tmp",
        stdin=None,
    )

    run = LocalExecutor().create_run(reader, spec)

    assert run == "run"
    command, created_output, kwargs = reader.call
    assert command == ["custom-launcher", "--flag", "program"]
    assert created_output is output_file
    assert kwargs["stdin"] is None
    assert kwargs["input_file"] == input_file.name

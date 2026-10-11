from types import SimpleNamespace

from ...input_parameters.input_parameters import InputParameters
from ..base import ExecutionSpec


def test_input_parameters_resolve_mpi_on_executor(tmp_path, monkeypatch):
    reader = object()
    parameters = object.__new__(InputParameters)
    object.__setattr__(
        parameters,
        "_definition",
        SimpleNamespace(mpi=True, executable="kkrscf"),
    )
    monkeypatch.setattr(
        InputParameters,
        "output_reader",
        lambda self, calculator, print_output, directory: reader,
    )

    seen = {}
    runner = object()

    class Profile:
        mpi = "auto"

        def create_execution(self, executable, input_file, output_file, directory, **kwargs):
            seen["profile_mpi"] = kwargs["mpi"]
            return ExecutionSpec((executable,), input_file, output_file, directory)

    class Executor:
        def resolve_mpi(self, mpi):
            seen["executor_mpi"] = mpi
            return ["remote-mpirun"]

        def create_run(self, output_reader, spec, callback=None):
            assert output_reader is reader
            seen["spec"] = spec
            return runner

    calculator = SimpleNamespace(
        profile=Profile(),
        executor=Executor(),
        value_or_default=lambda name, value: value,
    )
    input_file = (tmp_path / "run.inp").open("w+")
    output_file = (tmp_path / "run.out").open("wb")

    created = parameters.run_process(
        calculator,
        input_file,
        output_file,
        directory=str(tmp_path),
        run_async=True,
    )

    assert created is runner
    assert seen["executor_mpi"] == "auto"
    assert seen["profile_mpi"] == ["remote-mpirun"]

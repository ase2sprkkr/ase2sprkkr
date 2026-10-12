from types import SimpleNamespace

from ...input_parameters.input_parameters import InputParameters
from ...sprkkr.profile import SPRKKRProfile


def test_input_parameters_resolve_mpi_on_executor(tmp_path, monkeypatch):
    reader = object()
    parameters = InputParameters.create("SCF")
    monkeypatch.setattr(
        InputParameters,
        "output_reader",
        lambda self, calculator, print_output, directory: reader,
    )

    seen = {}
    runner = object()

    class Executor:
        def resolve_mpi(self, mpi):
            seen["executor_mpi"] = mpi
            return ("remote-mpirun",)

        def create_run(self, output_reader, spec, callback=None):
            assert output_reader is reader
            seen["spec"] = spec
            return runner

    profile = SPRKKRProfile(mpi="auto")
    executor = Executor()
    calculator = SimpleNamespace(
        value_or_default=lambda name, value: value,
    )

    with (tmp_path / "run.inp").open("w+") as input_file, (
        tmp_path / "run.out"
    ).open("wb") as output_file:
        created = parameters.run_process(
            calculator,
            input_file,
            output_file,
            directory=str(tmp_path),
            run_async=True,
            profile=profile,
            executor=executor,
        )

        assert created is runner
        assert seen["executor_mpi"] == "auto"
        assert seen["spec"].argv == (
            "remote-mpirun",
            "kkrscfMPI",
            input_file.name,
        )
        assert seen["spec"].stdin is None

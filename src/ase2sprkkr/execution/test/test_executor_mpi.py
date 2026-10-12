from types import SimpleNamespace
import warnings

import pytest

from ...execution import (
    ExecutionSpec,
    LaunchSpec,
    LocalExecutor,
    MPIRunnerNotFoundWarning,
)
from ...execution.local import materialize_local_spec
from ...input_parameters.input_parameters import InputParameters
from ...sprkkr.profile import SPRKKRProfile


def test_input_parameters_preserve_mpi_launch_intent(tmp_path, monkeypatch):
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
            return LaunchSpec.from_mpi(mpi)

        def create_run(self, output_reader, spec, callback=None):
            assert output_reader is reader
            seen["spec"] = spec
            return runner

    profile = SPRKKRProfile(mpi="auto", mpi_warning=False)
    calculator = SimpleNamespace(value_or_default=lambda name, value: value)

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
            executor=Executor(),
        )

        assert created is runner
        assert seen["executor_mpi"] == "auto"
        assert seen["spec"].argv == ("kkrscf",)
        assert seen["spec"].mpi_argv == ("kkrscfMPI", input_file.name)
        assert seen["spec"].launch == LaunchSpec(
            mpi="auto", mpi_warning=False
        )
        assert seen["spec"].stdin is input_file


def test_mpi_process_count_is_portable_launch_intent():
    assert LocalExecutor().resolve_mpi(16) == LaunchSpec(
        mpi=True, process_count=16
    )


def test_local_executor_materializes_mpi_on_execution_host(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "ase2sprkkr.configuration.mpi_runner",
        lambda mpi: ["srun"] if mpi is True else None,
    )
    input_file = (tmp_path / "run.inp").open("w+")
    spec = ExecutionSpec(
        argv=("kkrscf",),
        mpi_argv=("kkrscfMPI", input_file.name),
        launch=LaunchSpec(mpi=True, process_count=16),
        input_file=input_file,
        output_file=None,
        directory=str(tmp_path),
        stdin=input_file,
    )

    resolved = materialize_local_spec(spec)

    assert resolved.argv == ("srun", "-n", "16", "kkrscfMPI", input_file.name)
    assert resolved.stdin is None
    input_file.close()


def test_missing_explicit_mpi_emits_specific_warning(monkeypatch, tmp_path):
    monkeypatch.setattr("ase2sprkkr.configuration.mpi_runner", lambda mpi: None)
    input_file = (tmp_path / "run.inp").open("w+")
    spec = ExecutionSpec(
        argv=("kkrscf",),
        mpi_argv=("kkrscfMPI", input_file.name),
        launch=LaunchSpec(mpi=True),
        input_file=input_file,
        output_file=None,
        directory=str(tmp_path),
        stdin=input_file,
    )

    with pytest.warns(MPIRunnerNotFoundWarning):
        with pytest.raises(FileNotFoundError, match="no MPI launcher"):
            materialize_local_spec(spec)
    input_file.close()


def test_mpi_warning_can_be_suppressed(monkeypatch, tmp_path):
    monkeypatch.setattr("ase2sprkkr.configuration.mpi_runner", lambda mpi: None)
    input_file = (tmp_path / "run.inp").open("w+")
    spec = ExecutionSpec(
        argv=("kkrscf",),
        mpi_argv=("kkrscfMPI", input_file.name),
        launch=LaunchSpec(mpi=True, mpi_warning=False),
        input_file=input_file,
        output_file=None,
        directory=str(tmp_path),
        stdin=input_file,
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with pytest.raises(FileNotFoundError, match="no MPI launcher"):
            materialize_local_spec(spec)
    assert not [w for w in caught if issubclass(w.category, MPIRunnerNotFoundWarning)]
    input_file.close()


def test_local_auto_falls_back_to_serial_without_warning(monkeypatch, tmp_path):
    monkeypatch.setattr("ase2sprkkr.configuration.mpi_runner", lambda mpi: None)
    input_file = (tmp_path / "run.inp").open("w+")
    spec = ExecutionSpec(
        argv=("kkrscf",),
        mpi_argv=("kkrscfMPI", input_file.name),
        launch=LaunchSpec(mpi="auto"),
        input_file=input_file,
        output_file=None,
        directory=str(tmp_path),
        stdin=input_file,
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        resolved = materialize_local_spec(spec)

    assert resolved.argv == ("kkrscf",)
    assert resolved.stdin is input_file
    assert not [w for w in caught if issubclass(w.category, MPIRunnerNotFoundWarning)]
    input_file.close()


def test_gdb_forces_serial_execution(tmp_path):
    profile = SPRKKRProfile(mpi="auto", mpi_warning=False)
    input_file = (tmp_path / "run.inp").open("w+")

    spec = profile.create_execution(
        "kkrscf",
        input_file,
        None,
        str(tmp_path),
        mpi_supported=True,
        gdb=True,
    )

    assert spec.argv == ("gdb", "kkrscf")
    assert spec.mpi_argv is None
    assert spec.launch.mpi is False
    assert spec.launch.mpi_warning is False
    input_file.close()


def test_explicit_local_executor_overrides_remote_profile(tmp_path, monkeypatch):
    reader = object()
    parameters = InputParameters.create("SCF")
    monkeypatch.setattr(
        InputParameters,
        "output_reader",
        lambda self, calculator, print_output, directory: reader,
    )

    runner = object()
    seen = {}

    class ExplicitLocalExecutor(LocalExecutor):
        def resolve_mpi(self, mpi):
            seen["mpi"] = mpi
            return LaunchSpec()

        def create_run(self, output_reader, spec, callback=None):
            assert output_reader is reader
            seen["spec"] = spec
            return runner

    def must_not_resolve_profile(profile):
        raise AssertionError("explicit executor was replaced by the profile backend")

    monkeypatch.setattr(
        "ase2sprkkr.input_parameters.input_parameters.executor_from_profile",
        must_not_resolve_profile,
    )

    profile = SPRKKRProfile(executor="ssh", host="cluster", mpi="auto")
    executor = ExplicitLocalExecutor()
    calculator = SimpleNamespace(
        value_or_default=lambda name, value: value,
        executor=None,
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
        assert seen["mpi"] == "auto"
        assert seen["spec"].argv == ("kkrscf",)
        assert seen["spec"].mpi_argv is None
        assert seen["spec"].launch == LaunchSpec()

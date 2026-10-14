import asyncio
import warnings
from types import SimpleNamespace

import pytest

from ...execution import (
    ExecutionSpec,
    LaunchSpec,
    LocalExecutor,
    MpiNotFoundWarning,
    RemoteEndpointExecutor,
    RunState,
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
    runner = SimpleNamespace(state=RunState.NEW, result=object())

    class Executor:
        def resolve_mpi(self, mpi):
            seen["executor_mpi"] = mpi
            return LaunchSpec.from_mpi(mpi)

        def create_run(self, output_reader, spec, callback=None):
            assert output_reader is reader
            seen["spec"] = spec
            return runner

    profile = SPRKKRProfile(mpi="auto")
    calculator = SimpleNamespace(
        value_or_default=lambda name, value: value,
        persistent=False,
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
            executor=Executor(),
        )

        assert created.result is runner.result
        assert seen["executor_mpi"] == "auto"
        assert seen["spec"].argv == ("kkrscf",)
        assert seen["spec"].mpi_argv == ("kkrscfMPI", input_file.name)
        assert seen["spec"].launch == LaunchSpec(mpi="auto")
        assert seen["spec"].stdin is input_file


def test_mpi_process_count_is_portable_launch_intent():
    assert LocalExecutor().resolve_mpi(16) == LaunchSpec(
        mpi=True, process_count=16
    )


def test_local_executor_materializes_mpi_on_execution_host(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "ase2sprkkr.execution.local.resolve_mpi_command",
        lambda mpi: ("srun",) if mpi is True else None,
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


def test_local_auto_warns_and_falls_back_to_serial(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "ase2sprkkr.execution.local.resolve_mpi_command", lambda mpi: None
    )
    input_file = (tmp_path / "run.inp").open("w+")
    spec = ExecutionSpec(
        argv=("kkrscf",),
        mpi_argv=("kkrscfMPI", input_file.name),
        launch=LaunchSpec(mpi="auto", mpi_warning=True),
        input_file=input_file,
        output_file=None,
        directory=str(tmp_path),
        stdin=input_file,
    )

    with pytest.warns(MpiNotFoundWarning, match="without MPI"):
        resolved = materialize_local_spec(spec)

    assert resolved.argv == ("kkrscf",)
    assert resolved.stdin is input_file
    input_file.close()


def test_local_auto_can_disable_missing_mpi_warning(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "ase2sprkkr.execution.local.resolve_mpi_command", lambda mpi: None
    )
    input_file = (tmp_path / "run.inp").open("w+")
    spec = ExecutionSpec(
        argv=("kkrscf",),
        mpi_argv=("kkrscfMPI", input_file.name),
        launch=LaunchSpec(mpi="auto", mpi_warning=False),
        input_file=input_file,
        output_file=None,
        directory=str(tmp_path),
        stdin=input_file,
    )

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        resolved = materialize_local_spec(spec)

    assert not caught
    assert resolved.argv == ("kkrscf",)
    assert resolved.stdin is input_file
    input_file.close()


def test_remote_endpoint_warns_from_remote_autodetection(monkeypatch):
    endpoint = object.__new__(RemoteEndpointExecutor)
    endpoint.backend = "local"
    calls = []

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        if action == "resolve-mpi":
            return {"mpi": None}
        return {"state": "queued", "directory": "/remote/run"}

    monkeypatch.setattr(endpoint, "_json_request", request)
    payload = {
        "directory": "/remote/run",
        "launch": {"mpi": "auto", "mpi_warning": True},
    }

    with pytest.warns(MpiNotFoundWarning, match="remote host"):
        asyncio.run(endpoint._submit(payload))

    assert [call[0] for call in calls] == ["resolve-mpi", "submit"]


def test_remote_endpoint_skips_warning_probe_when_disabled(monkeypatch):
    endpoint = object.__new__(RemoteEndpointExecutor)
    endpoint.backend = "local"
    calls = []

    async def request(action, *args, payload=None):
        calls.append((action, args, payload))
        return {"state": "queued", "directory": "/remote/run"}

    monkeypatch.setattr(endpoint, "_json_request", request)
    payload = {
        "directory": "/remote/run",
        "launch": {"mpi": "auto", "mpi_warning": False},
    }

    asyncio.run(endpoint._submit(payload))

    assert [call[0] for call in calls] == ["submit"]


def test_explicit_local_executor_overrides_remote_profile(tmp_path, monkeypatch):
    reader = object()
    parameters = InputParameters.create("SCF")
    monkeypatch.setattr(
        InputParameters,
        "output_reader",
        lambda self, calculator, print_output, directory: reader,
    )

    runner = SimpleNamespace(state=RunState.NEW, result=object())
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
        persistent=False,
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

        assert created.result is runner.result
        assert seen["mpi"] == "auto"
        assert seen["spec"].argv == ("kkrscf",)
        assert seen["spec"].mpi_argv is None
        assert seen["spec"].launch == LaunchSpec()

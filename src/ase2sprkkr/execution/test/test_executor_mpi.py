import asyncio
import warnings
from types import SimpleNamespace

import pytest

from ...execution import (
    ExecutionSpec,
    LaunchSpec,
    LocalExecutor,
    MpiNotFoundWarning,
    MpiSpec,
    RemoteEndpointExecutor,
    RunState,
)
from ...execution.local import materialize_local_spec
from ...input_parameters.input_parameters import InputParameters
from ...input_parameters.task_definition import TaskDefinition
from ...sprkkr import calculation as calculation_module
from ...sprkkr.calculation import PreparedCalculation
from ...sprkkr.profile import SPRKKRProfile


class _Prepared:
    create_run = PreparedCalculation.create_run

    def __init__(self, parameters, input_file, output_file, directory):
        self.input_parameters = parameters
        self.input_file = input_file
        self.potential_file = None
        self.output_file = output_file
        self.directory = str(directory)
        self.closed = False

    def close(self):
        if self.closed:
            return
        self.closed = True
        self.input_file.close()
        self.output_file.close()


def test_input_parameters_preserve_mpi_launch_intent(tmp_path, monkeypatch):
    reader = object()
    parameters = InputParameters.create("SCF")
    monkeypatch.setattr(
        TaskDefinition,
        "output_reader",
        lambda self, **kwargs: reader,
    )

    seen = {}
    runner = SimpleNamespace(
        state=RunState.NEW,
        result=object(),
        reference=None,
        exception=None,
    )

    class Executor:
        def resolve_mpi(self, mpi):
            seen["executor_mpi"] = mpi
            return LaunchSpec.from_mpi(mpi)

        def create_run(self, output_reader, spec, callback=None):
            assert output_reader is reader
            seen["spec"] = spec
            return runner

    profile = SPRKKRProfile(mpi="auto")
    calculator = SimpleNamespace(directory=str(tmp_path))

    with (tmp_path / "run.inp").open("w+") as input_file, (
        tmp_path / "run.out"
    ).open("wb") as output_file:
        created = _Prepared(parameters, input_file, output_file, tmp_path).create_run(
            result_context=calculator,
            profile=profile,
            executor=Executor(),
        )

        assert created.result is runner.result
        assert seen["executor_mpi"] == "auto"
        assert seen["spec"].argv == ("kkrscf",)
        assert seen["spec"].mpi_argv == ("kkrscfMPI", input_file.name)
        assert seen["spec"].launch == LaunchSpec.from_mpi("auto")
        assert seen["spec"].stdin is input_file
        created.close()


def test_mpi_process_count_is_portable_launch_intent():
    assert LocalExecutor().resolve_mpi(16) == LaunchSpec(
        mpi=MpiSpec(required=True, autodetect=True, process_count=16)
    )


def test_local_executor_materializes_mpi_on_execution_host(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "ase2sprkkr.execution.mpi.find_default_mpi_launcher",
        lambda: ("srun",),
    )
    input_file = (tmp_path / "run.inp").open("w+")
    spec = ExecutionSpec(
        argv=("kkrscf",),
        mpi_argv=("kkrscfMPI", input_file.name),
        launch=LaunchSpec.from_mpi(16),
        input_file=input_file,
        output_file=None,
        directory=str(tmp_path),
        stdin=input_file,
    )

    resolved = materialize_local_spec(spec)

    assert resolved.argv == ("srun", "-n", "16", "kkrscfMPI", input_file.name)
    assert resolved.stdin is None
    input_file.close()


def test_required_local_mpi_without_launcher_fails(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "ase2sprkkr.execution.mpi.find_default_mpi_launcher", lambda: None
    )
    input_file = (tmp_path / "run.inp").open("w+")
    spec = ExecutionSpec(
        argv=("kkrscf",),
        mpi_argv=("kkrscfMPI", input_file.name),
        launch=LaunchSpec.from_mpi(True),
        input_file=input_file,
        output_file=None,
        directory=str(tmp_path),
        stdin=input_file,
    )

    with pytest.raises(FileNotFoundError, match="no MPI launcher"):
        materialize_local_spec(spec)
    input_file.close()


def test_local_auto_warns_and_falls_back_to_serial(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "ase2sprkkr.execution.local.resolve_mpi_command", lambda mpi: None
    )
    input_file = (tmp_path / "run.inp").open("w+")
    spec = ExecutionSpec(
        argv=("kkrscf",),
        mpi_argv=("kkrscfMPI", input_file.name),
        launch=LaunchSpec.from_mpi("auto"),
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
        launch=LaunchSpec(
            mpi=MpiSpec(autodetect=True, warning=False)
        ),
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
        "launch": {
            "mpi": {
                "required": False,
                "autodetect": True,
                "launcher": [],
                "process_count": None,
                "warning": True,
            }
        },
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
        "launch": {
            "mpi": {
                "required": False,
                "autodetect": True,
                "launcher": [],
                "process_count": None,
                "warning": False,
            }
        },
    }

    asyncio.run(endpoint._submit(payload))

    assert [call[0] for call in calls] == ["submit"]


def test_explicit_local_executor_overrides_remote_profile(tmp_path, monkeypatch):
    reader = object()
    parameters = InputParameters.create("SCF")
    monkeypatch.setattr(
        TaskDefinition,
        "output_reader",
        lambda self, **kwargs: reader,
    )

    runner = SimpleNamespace(
        state=RunState.NEW,
        result=object(),
        reference=None,
        exception=None,
    )
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

    monkeypatch.setattr(calculation_module, "executor_from_profile", must_not_resolve_profile)

    profile = SPRKKRProfile(executor="ssh", host="cluster", mpi="auto")
    executor = ExplicitLocalExecutor()
    calculator = SimpleNamespace(directory=str(tmp_path))

    with (tmp_path / "run.inp").open("w+") as input_file, (
        tmp_path / "run.out"
    ).open("wb") as output_file:
        created = _Prepared(parameters, input_file, output_file, tmp_path).create_run(
            result_context=calculator,
            profile=profile,
            executor=executor,
        )

        assert created.result is runner.result
        assert seen["mpi"] == "auto"
        assert seen["spec"].argv == ("kkrscf",)
        assert seen["spec"].mpi_argv is None
        assert seen["spec"].launch == LaunchSpec()
        created.close()

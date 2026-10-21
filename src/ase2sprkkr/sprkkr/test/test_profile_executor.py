from types import SimpleNamespace

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.execution import ExecutorPool, LocalExecutor, RunState, SshTransport
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.input_parameters.task_definition import TaskDefinition
from ase2sprkkr.profiles import ProfilesConfigurationContainer
from ase2sprkkr.sprkkr import calculation as calculation_module
from ase2sprkkr.sprkkr.calculation import PreparedCalculation
from ase2sprkkr.sprkkr.profile import (
    SPRKKRProfile,
    executor_from_profile,
    resolve_sprkkr_profile,
)


def fresh_config(monkeypatch):
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: None),
    )
    return definition.create_object()


def test_named_ssh_profile_resolves_executor_settings(monkeypatch):
    configuration = fresh_config(monkeypatch)
    configuration.profiles.set(
        {
            "cluster": {
                "executor": "ssh",
                "host": "user@cluster",
                "local_root": "/home/user/work",
                "remote_root": "/scratch/user/work",
                "ssh_options": ["-p", "2222"],
                "ssh_zstd": True,
                "remote_command": ["/opt/ase2sprkkr/bin/ase2sprkkr", "remote"],
                "input_sync": "missing",
                "output_sync": "always",
                "mpi": 16,
            }
        }
    )

    profile = resolve_sprkkr_profile("cluster", configuration=configuration)

    assert profile.executor == "ssh"
    assert profile.host == "user@cluster"
    assert profile.local_root == "/home/user/work"
    assert profile.remote_root == "/scratch/user/work"
    assert profile.ssh_options == ("-p", "2222")
    assert profile.ssh_zstd is True
    assert profile.remote_command == (
        "/opt/ase2sprkkr/bin/ase2sprkkr",
        "remote",
    )
    assert profile.input_sync == "missing"
    assert profile.output_sync == "always"
    assert profile.mpi == 16


def test_local_profile_creates_local_executor():
    executor = executor_from_profile(SPRKKRProfile())

    assert isinstance(executor, LocalExecutor)


def test_ssh_profile_creates_endpoint_only_executor():
    profile = SPRKKRProfile(
        executor="ssh",
        host="user@cluster",
        local_root="/home/user/work",
        remote_root="/scratch/user/work",
        ssh_options=("-p", "2222"),
        ssh_zstd=True,
        remote_command=("custom-ase2sprkkr", "remote"),
        input_sync="missing",
        output_sync="always",
    )

    executor = executor_from_profile(profile)

    assert isinstance(executor, RemoteEndpointExecutor)
    assert isinstance(executor.transport, SshTransport)
    assert executor.transport.host == "user@cluster"
    assert executor.transport.options == ("-p", "2222")
    assert executor.transport.zstd is True
    assert executor.command == ("custom-ase2sprkkr", "remote")
    input_file = SimpleNamespace(name="/home/user/work/test.inp")
    spec = profile.create_execution(
        "kkrscf",
        input_file,
        object(),
        "/home/user/work",
        mpi_supported=True,
    )
    assert spec.launch.local_root == "/home/user/work"
    assert spec.launch.remote_root == "/scratch/user/work"
    assert spec.launch.input_sync == "missing"
    assert spec.launch.output_sync == "always"


def test_ssh_profile_requires_host():
    with pytest.raises(ValueError, match="requires 'host'"):
        executor_from_profile(SPRKKRProfile(executor="ssh"))


def test_psij_profile_creates_configured_executor(monkeypatch):
    created = []

    class FakePsijExecutor:
        def __init__(self, backend, *, url=None):
            self.backend = backend
            self.url = url
            created.append((backend, url))

    monkeypatch.setattr("ase2sprkkr.execution.psij.PsijExecutor", FakePsijExecutor)
    monkeypatch.setattr(
        "ase2sprkkr.execution.pool.default_executor_pool",
        lambda: ExecutorPool(),
    )

    executor = executor_from_profile(
        SPRKKRProfile(
            executor="psij",
            psij_backend="slurm",
            psij_url="local://scheduler",
        )
    )

    assert isinstance(executor, FakePsijExecutor)
    assert created == [("slurm", "local://scheduler")]


def test_unknown_runtime_executor_is_rejected():
    with pytest.raises(ValueError, match="unsupported execution endpoint"):
        executor_from_profile(SPRKKRProfile(executor="custom"))


def test_profile_executor_is_selected_before_mpi_resolution(monkeypatch, tmp_path):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "test"
    parameters.CONTROL.POTFIL = "test.pot"
    backend_runner = SimpleNamespace(
        state=RunState.NEW,
        result=object(),
        reference=None,
        exception=None,
        add_output_listener=lambda callback: None,
    )
    calls = []

    class FakeExecutor:
        def resolve_mpi(self, value):
            calls.append(("resolve_mpi", value))
            return False

        def create_run(self, spec, callback=None):
            calls.append(("create_run", spec.argv))
            return backend_runner

    fake_executor = FakeExecutor()
    monkeypatch.setattr(
        calculation_module,
        "executor_from_profile",
        lambda profile: calls.append(("factory", profile.executor)) or fake_executor,
    )
    monkeypatch.setattr(
        TaskDefinition,
        "output_reader",
        lambda self, **kwargs: SimpleNamespace(print_output=False),
    )

    profile = SPRKKRProfile(executor="ssh", host="cluster", mpi=8)
    # Neither stream is read or written: the mocked executor only inspects
    # the execution specification. Avoid filesystem I/O for this ordering test.
    closed = []
    input_file = SimpleNamespace(
        name=str(tmp_path / "test.inp"), close=lambda: closed.append("input"),
    )
    output_file = SimpleNamespace(
        name=str(tmp_path / "test.out"), close=lambda: closed.append("output"),
    )
    prepared = SimpleNamespace(
        input_parameters=parameters,
        input_file=input_file,
        potential_file=None,
        output_file=output_file,
        directory=str(tmp_path),
        close=lambda: (input_file.close(), output_file.close()),
    )

    result = PreparedCalculation.create_run(
        prepared,
        result_context=SimpleNamespace(directory=str(tmp_path)),
        profile=profile,
        error_file=True,  # Stderr policy is tested separately; avoid opening it.
    )
    try:
        assert result.result is None  # Postprocessing starts only after submit
        assert calls[0] == ("factory", "ssh")
        assert calls[1] == ("resolve_mpi", 8)
        assert calls[2][0] == "create_run"
    finally:
        result.close()

    assert closed == ["input", "output"]
    assert not (tmp_path / "test.inp").exists()
    assert not (tmp_path / "test.out").exists()

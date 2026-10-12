from types import SimpleNamespace

import pytest

from ase2sprkkr.configuration import definition
from ase2sprkkr.execution import LocalExecutor, RemoteExecutor, SshTransport
from ase2sprkkr.input_parameters import input_parameters as input_parameters_module
from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.profiles import ProfilesConfigurationContainer
from ase2sprkkr.sprkkr import profile as profile_module
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


def test_ssh_profile_creates_remote_executor():
    profile = SPRKKRProfile(
        executor="ssh",
        host="user@cluster",
        local_root="/home/user/work",
        remote_root="/scratch/user/work",
        ssh_options=("-p", "2222"),
        remote_command=("custom-ase2sprkkr", "remote"),
        input_sync="missing",
        output_sync="always",
    )

    executor = executor_from_profile(profile)

    assert isinstance(executor, RemoteExecutor)
    assert isinstance(executor.transport, SshTransport)
    assert executor.transport.host == "user@cluster"
    assert executor.transport.options == ("-p", "2222")
    assert executor.local_root.as_posix() == "/home/user/work"
    assert executor.remote_root == "/scratch/user/work"
    assert executor.command == ("custom-ase2sprkkr", "remote")
    assert executor.input_files.value == "missing"
    assert executor.output_files.value == "always"


def test_ssh_profile_requires_host():
    with pytest.raises(ValueError, match="requires 'host'"):
        executor_from_profile(SPRKKRProfile(executor="ssh"))


def test_psij_profile_creates_configured_executor(monkeypatch):
    created = []

    class FakePsijExecutor:
        def __init__(self, backend, *, url=None):
            created.append((backend, url))

    monkeypatch.setattr(profile_module, "PsijExecutor", FakePsijExecutor)

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
    with pytest.raises(ValueError, match="Unknown execution backend"):
        executor_from_profile(SPRKKRProfile(executor="custom"))


def test_profile_executor_is_selected_before_mpi_resolution(monkeypatch, tmp_path):
    parameters = InputParameters.create("DOS")
    parameters.CONTROL.DATASET = "test"
    parameters.CONTROL.POTFIL = "test.pot"
    runner = object()
    calls = []

    class FakeExecutor:
        def resolve_mpi(self, value):
            calls.append(("resolve_mpi", value))
            return False

        def create_run(self, reader, spec, callback=None):
            calls.append(("create_run", spec.argv))
            return runner

    fake_executor = FakeExecutor()
    monkeypatch.setattr(
        input_parameters_module,
        "executor_from_profile",
        lambda profile: calls.append(("factory", profile.executor)) or fake_executor,
    )
    monkeypatch.setattr(
        InputParameters,
        "output_reader",
        lambda self, *args, **kwargs: SimpleNamespace(print_output=False),
    )

    profile = SPRKKRProfile(executor="ssh", host="cluster", mpi=8)
    calculator = SimpleNamespace(
        profile=profile,
        executor=LocalExecutor(),
        value_or_default=lambda name, value: value,
    )

    input_file = (tmp_path / "test.inp").open("w+")
    output_file = (tmp_path / "test.out").open("wb")
    try:
        result = parameters.run_process(
            calculator,
            input_file,
            output_file,
            directory=str(tmp_path),
            profile=profile,
            run_async=True,
        )
    finally:
        input_file.close()
        output_file.close()

    assert result is runner
    assert calls[0] == ("factory", "ssh")
    assert calls[1] == ("resolve_mpi", 8)
    assert calls[2][0] == "create_run"

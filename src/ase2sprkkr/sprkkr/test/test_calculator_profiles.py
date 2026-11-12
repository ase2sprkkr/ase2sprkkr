from types import SimpleNamespace

from ase2sprkkr.configuration import definition
from kkrunner.client import EndpointClient
from kkrunner.client.transports import SshTransport
from kkrunner.runs import RunState
from ase2sprkkr.profiles import ProfilesConfigurationContainer
from ase2sprkkr.sprkkr import profile as profile_module
from ase2sprkkr.sprkkr.calculator import SPRKKR
from ase2sprkkr.sprkkr.profile import SPRKKRProfile, executor_from_profile


def fresh_profile_config(monkeypatch):
    configuration = definition.create_object()
    monkeypatch.setattr(profile_module, "config", configuration)
    monkeypatch.setattr(
        ProfilesConfigurationContainer,
        "_ase_config",
        staticmethod(lambda: None),
    )
    return configuration


def test_calculator_uses_configured_default_profile(monkeypatch):
    configuration = fresh_profile_config(monkeypatch)
    configuration.profiles.set(
        {
            "cluster": {
                "mpi": 8,
                "mpi_warning": False,
                "executable_dir": "/cluster/sprkkr",
                "executable_suffix": ".x",
            }
        }
    )
    configuration.default_profile.set("cluster")

    calculator = SPRKKR()

    assert calculator.profile == SPRKKRProfile(
        executable_dir="/cluster/sprkkr",
        executable_suffix=".x",
        mpi=8,
        mpi_warning=False,
    )


def test_calculator_accepts_named_profile(monkeypatch):
    configuration = fresh_profile_config(monkeypatch)
    configuration.profiles.set(
        {
            "local": {"executable_dir": "/local"},
            "cluster": {"executable_dir": "/cluster", "mpi": 4},
        }
    )
    configuration.default_profile.set("local")

    calculator = SPRKKR(profile="cluster")

    assert calculator.profile.executable_dir == "/cluster"
    assert calculator.profile.mpi == 4
    assert calculator.profile.mpi_warning is True


def test_named_profile_without_mpi_uses_auto(monkeypatch):
    configuration = fresh_profile_config(monkeypatch)
    configuration.profiles.set({"local": {"executable_dir": "/local"}})

    calculator = SPRKKR(profile="local")

    assert calculator.profile.mpi == "auto"
    assert calculator.profile.mpi_warning is True


def test_calculator_init_overrides_copy_explicit_profile(monkeypatch):
    fresh_profile_config(monkeypatch)
    profile = SPRKKRProfile(
        executable_dir="/old",
        executable_suffix=".old",
        mpi=2,
    )

    calculator = SPRKKR(
        profile=profile,
        executable_dir="/new",
        executable_suffix=".new",
        mpi=6,
    )

    assert calculator.profile == SPRKKRProfile(
        executable_dir="/new",
        executable_suffix=".new",
        mpi=6,
    )
    assert calculator.profile is not profile
    assert profile == SPRKKRProfile("/old", ".old", 2)


def test_calculate_temporary_override_preserves_calculator_profile(monkeypatch):
    fresh_profile_config(monkeypatch)
    calculator = SPRKKR(profile=SPRKKRProfile("/calculator", ".calc", 1))
    profile = SPRKKRProfile("/cluster", ".x", 4)
    captured = []

    class Runner:
        state = RunState.FINISHED

        def run(self):
            return SimpleNamespace(task_name="DOS")

    def create_run(self, **kwargs):
        captured.append(kwargs["profile"])
        return Runner()

    monkeypatch.setattr(SPRKKR, "create_run", create_run)

    calculator.calculate(
        profile=profile,
        mpi=2,
        executable_dir="/override",
        executable_suffix=False,
    )

    assert captured[-1] == SPRKKRProfile("/override", False, 2)
    assert captured[-1] is not profile
    assert profile == SPRKKRProfile("/cluster", ".x", 4)

    # The one-off override must not affect the next calculation.
    calculator.calculate()
    assert captured[-1] is calculator.profile
    assert calculator.profile == SPRKKRProfile("/calculator", ".calc", 1)


def test_ssh_profile_uses_managed_remote_execution_by_default():
    client = executor_from_profile(
        SPRKKRProfile(executor="ssh", host="user@cluster")
    )
    assert isinstance(client, EndpointClient)
    assert isinstance(client.transport, SshTransport)
    assert client.backend.name == "process"


def test_ssh_profile_can_use_psij_on_target():
    client = executor_from_profile(
        SPRKKRProfile(
            executor="ssh", host="user@cluster", backend="psij",
            psij_backend="slurm", psij_url="local://scheduler",
        )
    )
    assert isinstance(client, EndpointClient)
    assert client.backend.name == "psij"
    assert client.backend.options == {
        "executor_name": "slurm", "url": "local://scheduler",
    }


def test_direct_psij_profile_defaults_to_slurm():
    client = executor_from_profile(SPRKKRProfile(executor="psij"))
    assert isinstance(client, EndpointClient)
    assert client.backend.name == "psij"
    assert client.backend.options == {"executor_name": "slurm", "url": None}

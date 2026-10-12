from types import SimpleNamespace

from ase2sprkkr.configuration import definition
from ase2sprkkr.profiles import ProfilesConfigurationContainer
from ase2sprkkr.sprkkr import profile as profile_module
from ase2sprkkr.sprkkr.calculator import SPRKKR
from ase2sprkkr.sprkkr.profile import SPRKKRProfile


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


def test_calculate_overrides_copy_explicit_profile(monkeypatch):
    fresh_profile_config(monkeypatch)
    calculator = SPRKKR(profile=SPRKKRProfile("/calculator", ".calc", 1))
    profile = SPRKKRProfile("/cluster", ".x", 4)
    captured = {}

    class Runner:
        def run(self):
            return SimpleNamespace(task_name="DOS")

    def create_run(self, **kwargs):
        captured.update(kwargs)
        return Runner()

    monkeypatch.setattr(SPRKKR, "create_run", create_run)

    calculator.calculate(
        profile=profile,
        mpi=2,
        executable_dir="/override",
        executable_suffix=False,
    )

    assert captured["profile"] == SPRKKRProfile("/override", False, 2)
    assert captured["profile"] is not profile
    assert profile == SPRKKRProfile("/cluster", ".x", 4)


def test_calculate_reuses_profile_when_no_override_is_given(monkeypatch):
    fresh_profile_config(monkeypatch)
    profile = SPRKKRProfile("/cluster", ".x", 4)
    calculator = SPRKKR(profile=profile)
    captured = {}

    class Runner:
        def run(self):
            return SimpleNamespace(task_name="DOS")

    def create_run(self, **kwargs):
        captured.update(kwargs)
        return Runner()

    monkeypatch.setattr(SPRKKR, "create_run", create_run)

    calculator.calculate()

    assert captured["profile"] is calculator.profile

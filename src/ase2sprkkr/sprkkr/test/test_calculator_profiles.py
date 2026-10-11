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


def test_calculator_profile_fields_are_single_authority(monkeypatch):
    fresh_profile_config(monkeypatch)
    calculator = SPRKKR(
        profile=SPRKKRProfile(
            executable_dir="/old",
            executable_suffix=".old",
            mpi=2,
        )
    )

    calculator.executable_dir = "/new"
    calculator.executable_suffix = ".new"
    calculator.mpi = 6

    assert calculator.profile == SPRKKRProfile(
        executable_dir="/new",
        executable_suffix=".new",
        mpi=6,
    )

import json
from dataclasses import asdict

import pytest

from ...execution import mpi as mpi_runtime
from ...execution.base import LaunchSpec


@pytest.mark.parametrize(
    "slurm, available, expected, probed",
    [
        (True, "srun", ("srun",), ("srun",)),
        (True, "mpirun", ("mpirun",), ("srun", "mpirun")),
        (False, "mpiexec", ("mpiexec",),
         ("mpirun", "mpirun.opmpirun", "mpirun.mpich", "mpiexec")),
    ],
)
def test_mpi_launcher_selection_and_probe_order(
    monkeypatch, slurm, available, expected, probed,
):
    if slurm:
        monkeypatch.setenv("SLURM_NODELIST", "node[01-02]")
    else:
        monkeypatch.delenv("SLURM_NODELIST", raising=False)
    seen = []

    def which(name):
        seen.append(name)
        return f"/usr/bin/{name}" if name == available else None

    monkeypatch.setattr(mpi_runtime.shutil, "which", which)
    assert mpi_runtime.find_default_mpi_launcher() == expected
    assert seen == list(probed)


def test_mpi_spec_distinguishes_required_auto_and_explicit_launchers():
    required = LaunchSpec.from_mpi(True).mpi
    automatic = LaunchSpec.from_mpi("auto").mpi
    explicit = LaunchSpec.from_mpi(("custom-mpi", "--bind")).mpi

    assert required.required is True
    assert required.autodetect is True
    assert required.allows_serial_fallback is False
    assert automatic.required is False
    assert automatic.autodetect is True
    assert automatic.allows_serial_fallback is True
    assert explicit.required is True
    assert explicit.autodetect is False
    assert explicit.launcher == ("custom-mpi", "--bind")


def test_launch_spec_json_round_trip_uses_nested_mpi_spec():
    launch = LaunchSpec.from_mpi(4)
    encoded = json.loads(json.dumps(asdict(launch)))

    assert encoded["mpi"] == {
        "required": True,
        "autodetect": True,
        "launcher": [],
        "process_count": 4,
        "warning": True,
    }
    assert "process_count" not in {key for key in encoded if key != "mpi"}
    assert "mpi_warning" not in encoded
    assert LaunchSpec(**encoded) == launch


def test_auto_disables_mpi_when_only_one_cpu_is_available(monkeypatch):
    monkeypatch.setattr(mpi_runtime.os, "sched_getaffinity", lambda _pid: {3})

    def must_not_probe():
        raise AssertionError("single-CPU auto mode must not probe MPI launchers")

    monkeypatch.setattr(mpi_runtime, "find_default_mpi_launcher", must_not_probe)
    spec = LaunchSpec.from_mpi("auto").mpi

    assert mpi_runtime.resolve_mpi_command(spec) is None


@pytest.mark.parametrize(
    "affinity",
    [None, {2, 7}],
    ids=("no-sched-getaffinity", "multiple-cpus"),
)
def test_auto_mpi_probes_host_when_cpu_count_is_not_one(monkeypatch, affinity):
    if affinity is None:
        monkeypatch.delattr(mpi_runtime.os, "sched_getaffinity", raising=False)
    else:
        monkeypatch.setattr(
            mpi_runtime.os, "sched_getaffinity", lambda _pid: affinity,
        )
    monkeypatch.setattr(
        mpi_runtime, "find_default_mpi_launcher", lambda: ("mpirun",),
    )
    spec = LaunchSpec.from_mpi("auto").mpi

    assert mpi_runtime.resolve_mpi_command(spec) == ("mpirun",)


def test_explicit_mpi_command_does_not_probe_host(monkeypatch):
    def must_not_probe():
        raise AssertionError("explicit MPI command must not probe the host")

    monkeypatch.setattr(mpi_runtime, "find_default_mpi_launcher", must_not_probe)
    spec = LaunchSpec.from_mpi(("custom-mpi", "--bind")).mpi

    assert mpi_runtime.resolve_mpi_command(spec) == (
        "custom-mpi",
        "--bind",
    )


@pytest.mark.parametrize(
    ("launcher", "count_option"),
    (("mpirun", "-np"), ("srun", "-n")),
)
def test_process_count_is_materialized_for_launcher(monkeypatch, launcher, count_option):
    monkeypatch.setattr(
        mpi_runtime, "find_default_mpi_launcher", lambda: (launcher,),
    )
    spec = LaunchSpec.from_mpi(4).mpi

    assert mpi_runtime.resolve_mpi_command(spec) == (launcher, count_option, "4")


def test_resolver_rejects_unparsed_mpi_values():
    try:
        mpi_runtime.resolve_mpi_command("auto")
    except TypeError as exc:
        assert "MpiSpec" in str(exc)
    else:
        raise AssertionError("raw MPI values must be normalized before resolution")

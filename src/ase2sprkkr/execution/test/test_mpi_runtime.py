from ...execution import mpi as mpi_runtime


def test_slurm_allocation_prefers_srun(monkeypatch):
    monkeypatch.setenv("SLURM_NODELIST", "node[01-02]")
    seen = []

    def which(name):
        seen.append(name)
        return f"/usr/bin/{name}"

    monkeypatch.setattr(mpi_runtime.shutil, "which", which)

    assert mpi_runtime.find_default_mpi_runner() == ("srun",)
    assert seen == ["srun"]


def test_default_launcher_uses_legacy_probe_order(monkeypatch):
    monkeypatch.delenv("SLURM_NODELIST", raising=False)
    seen = []

    def which(name):
        seen.append(name)
        return "/usr/bin/mpiexec" if name == "mpiexec" else None

    monkeypatch.setattr(mpi_runtime.shutil, "which", which)

    assert mpi_runtime.find_default_mpi_runner() == ("mpiexec",)
    assert seen == ["mpirun", "mpirun.opmpirun", "mpirun.mpich", "mpiexec"]


def test_auto_disables_mpi_when_only_one_cpu_is_available(monkeypatch):
    monkeypatch.setattr(mpi_runtime.os, "sched_getaffinity", lambda _pid: {3})

    def must_not_probe():
        raise AssertionError("single-CPU auto mode must not probe MPI launchers")

    monkeypatch.setattr(mpi_runtime, "find_default_mpi_runner", must_not_probe)

    assert mpi_runtime.resolve_mpi_command("auto") is None


def test_auto_probes_host_when_multiple_cpus_are_available(monkeypatch):
    monkeypatch.setattr(mpi_runtime.os, "sched_getaffinity", lambda _pid: {2, 7})
    monkeypatch.setattr(
        mpi_runtime, "find_default_mpi_runner", lambda: ("mpirun",)
    )

    assert mpi_runtime.resolve_mpi_command("auto") == ("mpirun",)


def test_explicit_mpi_command_does_not_probe_host(monkeypatch):
    def must_not_probe():
        raise AssertionError("explicit MPI command must not probe the host")

    monkeypatch.setattr(mpi_runtime, "find_default_mpi_runner", must_not_probe)

    assert mpi_runtime.resolve_mpi_command(("custom-mpi", "--bind")) == (
        "custom-mpi",
        "--bind",
    )

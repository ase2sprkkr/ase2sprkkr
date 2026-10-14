"""Host-specific MPI launcher detection and materialization."""

from __future__ import annotations

import os
import shutil


_MPI_LAUNCHERS = ("mpirun", "mpirun.opmpirun", "mpirun.mpich", "mpiexec")


def find_default_mpi_runner():
    """Return the preferred MPI launcher available on the current host.

    ``srun`` is preferred inside a SLURM allocation. Outside SLURM, launchers
    are considered in the same order as the legacy configuration-layer probe.
    Host detection is intentionally performed at execution time rather than
    cached in user configuration state.
    """
    if os.environ.get("SLURM_NODELIST") and shutil.which("srun"):
        return ("srun",)
    for launcher in _MPI_LAUNCHERS:
        if shutil.which(launcher):
            return (launcher,)
    return None


def _single_cpu_available():
    affinity = getattr(os, "sched_getaffinity", None)
    return affinity is not None and len(affinity(0)) == 1


def resolve_mpi_command(mpi):
    """Materialize portable MPI intent for the current execution host.

    ``mpi`` is the portable value carried by :class:`LaunchSpec`: ``False`` or
    ``None`` for serial execution, ``True`` for required MPI, ``"auto"`` for
    host-dependent MPI with serial fallback, or an explicit launcher argv.
    This function deliberately does not read user configuration.
    """
    if mpi is None or mpi is False:
        return None
    if mpi is True:
        return find_default_mpi_runner()
    if mpi == "auto":
        if _single_cpu_available():
            return None
        return find_default_mpi_runner()
    if isinstance(mpi, str):
        return (mpi,)
    return tuple(os.fspath(part) for part in mpi)

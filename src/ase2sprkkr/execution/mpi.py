"""Host-specific MPI launcher materialization."""

from __future__ import annotations

import os
import shutil
from typing import Optional, Tuple

from .base import MpiSpec


_MPI_LAUNCHERS = ("mpirun", "mpirun.opmpirun", "mpirun.mpich", "mpiexec")


def find_default_mpi_launcher() -> Optional[Tuple[str, ...]]:
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


def _single_cpu_available() -> bool:
    affinity = getattr(os, "sched_getaffinity", None)
    return affinity is not None and len(affinity(0)) == 1


def _with_process_count(
    launcher: Tuple[str, ...], process_count: Optional[int]
) -> Tuple[str, ...]:
    if process_count is None:
        return launcher
    executable = os.path.basename(os.fspath(launcher[0]))
    count_option = "-n" if executable == "srun" else "-np"
    return (*launcher, count_option, str(process_count))


def resolve_mpi_command(spec: MpiSpec) -> Optional[Tuple[str, ...]]:
    """Materialize normalized MPI intent for the current execution host.

    The returned tuple is the complete MPI launcher prefix, including a
    launcher-specific process-count option when requested. This function does
    not read user configuration and does not decide warning/error policy.
    """
    if not isinstance(spec, MpiSpec):
        raise TypeError("resolve_mpi_command() requires an MpiSpec")
    if not spec.requested:
        return None

    launcher = spec.launcher
    if not launcher:
        if spec.allows_serial_fallback and _single_cpu_available():
            return None
        launcher = find_default_mpi_launcher()
        if launcher is None:
            return None

    return _with_process_count(launcher, spec.process_count)

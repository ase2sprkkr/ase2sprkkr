"""Machine-specific SPR-KKR execution profile.

The split follows the same idea as ASE's calculator profiles: the profile knows
how an installed program is invoked, while the executor only knows how to run
an already prepared command.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any

from ..configuration import config
from ..execution import ExecutionSpec
from ..common.misc import first_non_none


@dataclass(frozen=True)
class SPRKKRProfile:
    """Describe how SPR-KKR executables are invoked on a computer."""

    executable_dir: Any = None
    executable_suffix: Any = None
    mpi: Any = None

    def _executable_suffix(self, override):
        suffix = first_non_none(override, self.executable_suffix)
        if suffix is False:
            return ""
        if suffix is True or suffix is None:
            return config.executables.suffix() or ""
        return suffix

    def _executable_dir(self, override):
        directory = first_non_none(override, self.executable_dir)
        if directory is None:
            directory = config.executables.dir()
        return directory

    def create_execution(
        self,
        executable,
        input_file,
        output_file,
        directory,
        *,
        mpi_supported=False,
        executable_suffix=None,
        executable_dir=None,
        mpi=None,
        gdb=False,
    ):
        """Create a concrete execution specification for one SPR-KKR task.

        ``mpi`` is an already resolved MPI command.  Resolving ``True``,
        ``"auto"`` or a process count belongs to the executor, because only
        the executor knows the machine on which the command will run.
        """

        executable += self._executable_suffix(executable_suffix)
        executable_dir = self._executable_dir(executable_dir)
        if executable_dir is not False and executable_dir:
            executable = os.path.join(executable_dir, executable)

        mpi_command = mpi if mpi_supported else None

        if mpi_command:
            input_name = input_file.name
            # SPR-KKR limits the input-file command-line argument to 80 chars.
            if len(input_name) > 80:
                input_name = os.path.relpath(input_name, directory)
            argv = [*mpi_command, executable + "MPI", input_name]
            stdin = None
        else:
            argv = [executable]
            stdin = input_file

        if gdb:
            stdin = None
            if mpi_command:
                print(f"run {input_file.name}")
                argv = argv[:-1]
            else:
                print(f"run <{input_file.name}")
            argv = ["gdb", *argv]

        return ExecutionSpec(
            argv=tuple(argv),
            input_file=input_file,
            output_file=output_file,
            directory=directory,
            stdin=stdin,
        )

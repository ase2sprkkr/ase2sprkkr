"""ASE calculator facade for SPR-KKR.

The calculator stores persistent user/ASE defaults.  Materialization of one
concrete calculation lives in :mod:`ase2sprkkr.sprkkr.calculation`; process
execution is delegated to the generic execution subsystem.
"""

from __future__ import annotations

import copy
import os
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Union

from ase.calculators.calculator import BaseCalculator, all_changes

from ..common.file_utils import Directory
from ..common.misc import config_property, first_non_none
from ..configuration import config
from ..potentials.potentials import Potential
from .profile import SPRKKRProfile, resolve_sprkkr_profile
from .sprkkr_atoms import SPRKKRAtoms


class SPRKKR(BaseCalculator):
    """ASE-facing calculator for SPR-KKR.

    ``SPRKKR`` is a mutable facade holding defaults shared by calculations:
    structure/potential, input parameters, filename templates, working directory
    and execution policy.  Each call snapshots those defaults into a separate
    :class:`PreparedCalculation` before execution.

    A calculator instance is not intended to be prepared concurrently from
    multiple threads because filename-template counters are mutable facade
    state.  Prepared runners themselves own independent resources and may run
    concurrently.
    """

    implemented_properties = ["energy"]

    def __init__(
        self,
        label=None,
        directory=".",
        atoms=None,
        potential=None,
        input_file=True,
        output_file=True,
        potential_file=True,
        input_parameters=None,
        options=None,
        print_output=None,
        mpi: Optional[bool] = None,
        empty_spheres: Optional[Union[str, bool, Dict]] = None,
        executable_suffix=None,
        executable_dir: Optional[Union[str, bool]] = None,
        profile=None,
        executor=None,
        persistent=False,
        **kwargs,
    ):
        if kwargs:
            print(f"Warning - unknown arguments in the SPRKKR calculator: {kwargs}")

        if potential and not isinstance(potential, (bool, str)):
            if atoms:
                raise ValueError(
                    "You can not specify both atoms and potential. "
                    "If you specify one of them, the other will be generated from it."
                )
            atoms = potential.atoms
        elif atoms:
            atoms = SPRKKRAtoms.promote_ase_atoms(atoms)

        self._potential = None
        self.atoms = atoms
        self.potential = potential

        super().__init__(use_cache=False)
        if directory is False:
            self._directory = directory
        else:
            self.directory = directory
        self.atoms = atoms
        self.potential = potential
        self.input_file = input_file
        self.output_file = output_file
        self.potential_file = potential_file
        self.input_parameters = input_parameters
        if options:
            self.input_parameters.set(options, unknown="find")
        self.print_output = print_output
        self.empty_spheres = empty_spheres
        self.profile = resolve_sprkkr_profile(profile).with_overrides(
            mpi=mpi,
            executable_dir=executable_dir,
            executable_suffix=executable_suffix,
        )
        self.executor = executor
        self.persistent = persistent
        self.label = label
        self._counter = 0

    print_output = config_property("print_output", config.running.print_output)
    empty_spheres = config_property("empty_spheres", config.running.empty_spheres)

    @property
    def input_parameters(self):
        """Default :class:`InputParameters` used for new calculations."""
        if self._input_parameters is None:
            self.input_parameters = "SCF"
        return self._input_parameters

    @input_parameters.setter
    def input_parameters(self, value):
        if value is not None:
            value = InputParameters.create_input_parameters(value)
        self._input_parameters = value

    def set(
        self,
        options: Union[Dict[str, Any], str, None] = None,
        value=None,
        *,
        unknown="find",
        **kwargs,
    ):
        """Set values on the calculator's persistent input-parameter defaults."""
        if isinstance(options, str):
            options = {options: value}
        elif value is not None:
            raise ValueError(
                "If value argument is given to SPRKKR.set, the options argument have to be string name of the value"
            )
        if kwargs or options:
            self.input_parameters.set(options, **kwargs, unknown=unknown)

    def get(self, name):
        """Return an input-parameter value by name or ``SECTION.NAME`` path."""
        return self.input_parameters.get(name)

    @property
    def potential(self) -> Potential:
        """Potential associated with the calculator, generated lazily from atoms."""
        if self._potential is None:
            if not self._atoms:
                return None
            self._potential = Potential.from_atoms(self._atoms)
        return self._potential

    @potential.setter
    def potential(self, pot):
        self._potential = pot
        if pot and not isinstance(pot, (bool, str)):
            atoms = pot.atoms
            if atoms:
                self._atoms = atoms
            else:
                pot._atoms = self._atoms

    @property
    def atoms(self) -> SPRKKRAtoms:
        """Structure associated with the calculator."""
        if self._atoms is None:
            if self._potential is None or isinstance(self._potential, bool):
                return None
            if isinstance(self._potential, str):
                self._potential = Potential.from_file(self._potential)
            self._atoms = self._potential.atoms
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        atoms = SPRKKRAtoms.promote_ase_atoms(atoms)
        self._atoms = atoms
        if self._potential:
            if isinstance(self._potential, str):
                self._potential = None
            elif self._potential is not True:
                self._potential.atoms = atoms

    @property
    def directory(self):
        """Default working directory, or ``False`` for a temporary run directory."""
        return self._directory

    @directory.setter
    def directory(self, directory: Union[bool, str, os.PathLike]):
        if not isinstance(directory, bool):
            directory = str(Path(directory))
        self._directory = directory

    def _advance_counter(self):
        """Advance the filename-template counter used by the legacy templator."""
        self._counter += 1
        return str(self._counter)

    def _prepare_calculation(self, **kwargs):
        """Snapshot defaults and materialize one concrete calculation."""
        prepared = PreparedCalculation.prepare(
            default_atoms=self._atoms,
            default_potential=self._potential,
            default_input_parameters=self.input_parameters,
            default_input_file=self.input_file,
            default_potential_file=self.potential_file,
            default_output_file=self.output_file,
            default_directory=self._directory,
            default_empty_spheres=self.empty_spheres,
            label=self.label,
            counter=self._counter,
            **kwargs,
        )
        self._counter = prepared.counter
        return prepared

    @staticmethod
    def _is_temporary_directory(directory):
        if directory is False:
            return True
        return isinstance(directory, Directory) and directory.dir is False

    def save_input(
        self,
        atoms=None,
        input_parameters=None,
        potential=None,
        input_file=None,
        potential_file=None,
        output_file=False,
        directory: Union[str, bool, Directory, None] = None,
        create_subdirs: bool = False,
        options=None,
        task=None,
        empty_spheres=None,
        return_files=False,
        restart_scf: bool = False,
    ):
        """Materialize SPR-KKR input files without starting an executable.

        Preparation always works on an independent snapshot of
        ``InputParameters``; generated POTFIL/DATASET values therefore never
        mutate calculator defaults or caller-owned parameter objects.

        ``False`` as ``input_file`` or ``potential_file`` requests a managed
        named temporary artifact. Temporary artifacts are useful for an
        executing runner, but ``save_input`` cannot return ownership of a
        temporary POTFIL because its public tuple exposes only the pathname.
        Supply a persistent POTFIL name when ``return_files=True``.

        ``directory=False`` is intentionally rejected here: a temporary
        directory is owned naturally by a runner, whereas returning filenames
        from ``save_input`` after deleting their directory would be invalid.
        Use :meth:`create_run`/``calculate`` for temporary run directories.
        """
        selected_directory = self._directory if directory is None else directory
        if self._is_temporary_directory(selected_directory):
            raise ValueError(
                "save_input() requires a persistent working directory; "
                "directory=False is supported by create_run()/calculate() instead"
            )

        prepared = self._prepare_calculation(
            atoms=atoms,
            input_parameters=input_parameters,
            potential=potential,
            input_file=input_file,
            potential_file=potential_file,
            output_file=output_file,
            directory=directory,
            create_subdirs=create_subdirs,
            options=options,
            task=task,
            empty_spheres=empty_spheres,
            allow_temporary_files=return_files,
            restart_scf=restart_scf,
        )

        if return_files:
            if prepared.temporary_potential:
                prepared.close()
                raise ValueError(
                    "save_input(return_files=True) cannot transfer ownership of "
                    "a temporary POTFIL pathname; supply potential_file as a path"
                )
            prepared.release_file_ownership()
            prepared.close_directory_context()
            return (
                prepared.input_parameters,
                prepared.input_file,
                prepared.potential_file,
                prepared.output_file,
            )

        input_name = prepared.input_file.name
        output_name = prepared.output_file.name if prepared.output_file else None
        potential_name = prepared.potential_file
        input_parameters = prepared.input_parameters
        prepared.close()
        return input_parameters, input_name, potential_name, output_name

    def create_run(
        self,
        atoms=None,
        input_parameters=None,
        potential=None,
        input_file=None,
        potential_file=None,
        output_file=None,
        directory: Union[str, bool, Directory, None] = None,
        create_subdirs: bool = False,
        options=None,
        task=None,
        empty_spheres: Optional[str | bool] = None,
        print_output=None,
        read_callback: Optional[callable] = None,
        callback=None,
        gdb=False,
        profile: Optional[SPRKKRProfile] = None,
        executor=None,
        persistent=None,
        provenance=None,
        restart_scf: bool = False,
    ):
        """Prepare a calculation and return an unstarted owning runner.

        An explicit per-run executor has highest precedence. Otherwise an
        explicit per-run profile selects its backend; if neither is supplied,
        the calculator-level executor (when present) overrides the calculator
        profile backend.

        ``provenance`` is opaque caller-owned metadata. ASE2SPRKKR neither
        interprets nor persists it; a live result receives the exact same object
        before the terminal callback is invoked.

        ``callback`` is a terminal runner callback. It can therefore be invoked
        with a partially populated result when backend execution or parsing
        fails. ``read_callback`` observes streamed stdout/stderr.
        """
        explicit_profile = profile is not None
        run_profile = resolve_sprkkr_profile(
            self.profile if profile is None else profile
        )
        print_output = first_non_none(print_output, self.print_output)
        if persistent is None:
            persistent = self.persistent
        if executor is None and not explicit_profile:
            executor = self.executor

        prepared = self._prepare_calculation(
            atoms=atoms,
            input_parameters=input_parameters,
            potential=potential,
            input_file=input_file,
            potential_file=potential_file,
            output_file=output_file,
            directory=directory,
            create_subdirs=create_subdirs,
            empty_spheres=empty_spheres,
            options=options,
            task=task,
            allow_temporary_files=True,
            restart_scf=restart_scf,
        )

        def finalized(result):
            result.provenance = provenance
            if callback:
                callback(result)

        return _create_sprkkr_run(
            prepared,
            result_context=self,
            profile=run_profile,
            executor=executor,
            persistent=persistent,
            print_output=print_output,
            read_callback=read_callback,
            callback=finalized,
            gdb=gdb,
        )

    def value_or_default(self, name, value):
        """Return calculator default *name* when *value* is ``None``."""
        if value is None:
            value = getattr(self, name)
        return value

    def _update_from_result(self, result, *, adopt_potential=False):
        """Update calculator state from a completed synchronous result."""
        if hasattr(result, "energy"):
            self.results["energy"] = result.energy

        if adopt_potential:
            try:
                potential = result.potential
            except (ValueError, FileNotFoundError):
                return
            if potential:
                self.potential = potential

    def calculate(
        self,
        atoms=None,
        properties=("energy",),
        system_changes=all_changes,
        input_parameters=None,
        potential=None,
        input_file=None,
        potential_file=None,
        output_file=None,
        directory: Union[str, bool, Directory, None] = None,
        create_subdirs: bool = False,
        empty_spheres: Optional[str | bool] = None,
        mpi: bool = None,
        options=None,
        task=None,
        print_output=None,
        read_callback: Optional[callable] = None,
        executable_suffix=None,
        callback=None,
        executable_dir=None,
        gdb=False,
        profile=None,
        executor=None,
        persistent=None,
        provenance=None,
        restart_scf: bool = False,
    ):
        """Run one SPR-KKR task synchronously and return its parsed result.

        ``properties`` and ``system_changes`` are accepted for ASE compatibility;
        task selection is currently controlled by ``input_parameters``/``task``.
        Per-call ``mpi``, ``executable_dir`` and ``executable_suffix`` values are
        immutable overrides of the selected execution profile.

        ``provenance`` is opaque runtime metadata forwarded unchanged to the
        returned result. It is intentionally not part of persistent run state.

        The callback of this synchronous convenience method is intentionally a
        success callback: it runs only after the result has updated calculator
        state. This differs from the lower-level terminal callback accepted by
        :meth:`create_run`.
        """
        if executor is None and profile is None:
            executor = self.executor
        run_profile = resolve_sprkkr_profile(
            self.profile if profile is None else profile
        ).with_overrides(
            mpi=mpi,
            executable_dir=executable_dir,
            executable_suffix=executable_suffix,
        )
        runner = self.create_run(
            atoms=atoms,
            input_parameters=input_parameters,
            potential=potential,
            input_file=input_file,
            potential_file=potential_file,
            output_file=output_file,
            directory=directory,
            create_subdirs=create_subdirs,
            options=options,
            task=task,
            empty_spheres=empty_spheres,
            print_output=print_output,
            read_callback=read_callback,
            gdb=gdb,
            profile=run_profile,
            executor=executor,
            persistent=persistent,
            provenance=provenance,
            restart_scf=restart_scf,
        )
        result = runner.run()
        self._update_from_result(
            result,
            adopt_potential=(
                result.task_name.upper() == "SCF" and not potential and not atoms
            ),
        )
        if callback:
            callback(result)
        return result

    def scf(self, *args, **kwargs):
        """Run an SCF calculation."""
        self.calculate(input_parameters="SCF", *args, **kwargs)

    def phagen(self, *args, **kwargs):
        """Run a PHAGEN calculation."""
        self.calculate(input_parameters="PHAGEN", *args, **kwargs)

    def kkrgen(self, *args, **kwargs):
        """Run a KKRGEN calculation."""
        self.calculate(input_parameters="KKRGEN", *args, **kwargs)

    def kkrspec(self, *args, **kwargs):
        """Run a KKRSPEC calculation."""
        self.calculate(input_parameters="KKRSPEC", *args, **kwargs)

    def kkrch(self, *args, **kwargs):
        """Run a KKRCH calculation."""
        self.calculate(input_parameters="KKRCH", *args, **kwargs)

    def run_xband(self, executable="xband"):
        """Start the external ``xband`` executable detached from calculator IO."""
        subprocess.Popen(executable, stdout=subprocess.DEVNULL)

    def copy_with_potential(self, potential):
        """Return a shallow calculator copy using *potential*."""
        out = copy.copy(self)
        out.potential = potential
        out.results = out.results.copy()
        return out

    @property
    def potential_object(self):
        """Return the calculator potential as a :class:`Potential` object."""
        if isinstance(self.potential, str):
            self.potential = Potential.from_file(self.potential)
        return self.potential

    def change_task(self, task, retain_values=False):
        """Change the persistent default input-parameter task."""
        self.input_parameters.change_task(task, retain_values)


# Keep these imports below the calculator class. Importing InputParameters while
# sprkkr.configuration is still initializing creates a cycle through potentials;
# calculation also imports InputParameters, so it must be delayed for the same
# reason. Method bodies resolve these globals only when called, after module init.
from ..input_parameters.input_parameters import InputParameters  # NOQA: E402
from .calculation import (  # NOQA: E402
    FilenameTemplator as _FilenameTemplator,
    PreparedCalculation,
    create_sprkkr_run as _create_sprkkr_run,
)


class FilenameTemplator(_FilenameTemplator):
    """Backward-compatible facade for the historical calculator-based API."""

    def __init__(self, calculator, task, symbols):
        self.calculator = calculator
        super().__init__(calculator.label, calculator._counter, task, symbols)

    def __call__(self, template):
        rendered = super().__call__(template)
        self.calculator._counter = self.counter
        return rendered


if os.path.basename(sys.argv[0]) != "sphinx-build":
    SPRKKR.InputParameters = InputParameters
    SPRKKR.Potential = Potential
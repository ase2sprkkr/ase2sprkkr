"""Containers for configuration parameters of SPR-KKR tasks."""

from __future__ import annotations

import io
import os
import pyparsing as pp
from typing import Union

from . import definitions
from .task_definition import TaskRegistry
from ..sprkkr.configuration import ConfigurationFile, ConfigurationSection
from ..common.decorators import cached_class_property
from ..common.warnings import DataValidityError
from ..configuration import config
from ..potentials.potentials import Potential
from ..sprkkr.profile import SPRKKRProfile, executor_from_profile
from ..execution.base import ExecutionResources, ResourceOwnedRunner
from ..execution.persistence import resolve_execution_manager


class InputSection(ConfigurationSection):
    """Input parameters sections have nothing special yet."""


class InputParameters(ConfigurationFile):
    """Concrete values of an SPR-KKR input file.

    Canonical task identity is resolved through :meth:`get_task_definition`.
    The default implementation follows the input schema back to its owning
    TaskDefinition. ``task_name`` is only the concrete task view used to create
    these parameters, so aliases such as BSFEK/BSFKK need no second task object.
    """

    task_registry = TaskRegistry(definitions)

    def __init__(
        self,
        definition,
        inputfile=None,
        outputfile=False,
        *,
        task_name=None,
    ):
        super().__init__(definition)
        self._task_name = str(task_name or definition.name).upper()
        self._inputfile = inputfile
        self._outputfile = inputfile if outputfile is False else outputfile

    def get_task_definition(self):
        """Return the canonical task governing these parameters, if any.

        The lookup is intentionally overridable instead of being duplicated as
        instance state.
        """
        return getattr(self._definition, "task_definition", None)

    @property
    def task_definition(self):
        return self.get_task_definition()

    def _require_task_definition(self):
        task = self.get_task_definition()
        if task is None:
            raise RuntimeError(
                "These input parameters are not associated with an SPR-KKR task definition"
            )
        return task

    @property
    def potential(self):
        pot = getattr(self, "_potential", None)
        potfil = self.CONTROL.POTFIL.result
        if pot:
            if potfil == pot[0]:
                return pot[1]
        out = Potential.from_file(potfil) if potfil else None
        self._potential = potfil, out
        return out

    def resolve_executable_suffix(self, suffix: Union[str, bool]):
        if suffix is False:
            return ""
        if suffix is True:
            return config.executables.suffix()
        return suffix

    def copy(self, copy_values: bool = False):
        """Copy parameters while preserving their concrete task view."""
        values = self.as_dict(copy=copy_values, only_changed=True)
        out = self.__class__(definition=self._definition, task_name=self.task_name)
        out.set(values, unknown="add", retain_invalid="all", report_invalid="ignore")
        return out

    @property
    def task_name(self):
        """Concrete requested/view task name used for persistence and UX.

        Canonical identity is ``get_task_definition().name``. A low-level native
        execution mode such as ARPES->SPLEED is interpreted by TaskDefinition and
        is deliberately not persisted as task identity.
        """
        return self._task_name

    def uses_separate_single_site_contour(self):
        switches = (
            ("ENERGY", "SPLITSS"),
            ("CONTROL", "SPLITSS"),
            ("CONTROL", "FSOHFF"),
        )
        for section, name in switches:
            try:
                if self[section][name]():
                    return True
            except (KeyError, AttributeError):
                pass
        return self.task_name.upper() == "COMPTON"

    def is_mpi(self, mpi=True):
        return mpi and self._require_task_definition().mpi

    def input_files(self, directory):
        return self._require_task_definition().input_files(self, directory)

    def output_files(self, directory):
        return self._require_task_definition().output_files(self, directory)

    def cache_files(self, directory):
        return self._require_task_definition().cache_files(self, directory)

    def run_process(
        self,
        calculator,
        input_file,
        output_file,
        directory=".",
        print_output=None,
        read_callback=None,
        run_async=False,
        callback=None,
        gdb=False,
        profile=None,
        executor=None,
        persistent=None,
    ):
        """Run the SPR-KKR task described by these input parameters."""
        task = self._require_task_definition()
        resources = ExecutionResources(input_file, output_file)
        handed_off = False

        try:
            print_output = calculator.value_or_default("print_output", print_output)
            reader = self.output_reader(calculator, print_output, directory)
            explicit_profile = profile is not None
            profile = profile or getattr(calculator, "profile", None) or SPRKKRProfile()
            if executor is None and not explicit_profile:
                executor = getattr(calculator, "executor", None)
            if executor is None:
                executor = executor_from_profile(profile)
            if persistent is None:
                persistent = calculator.persistent
            manager = resolve_execution_manager(persistent)

            execution_profile = profile
            if task.mpi:
                execution_profile = profile.with_overrides(
                    mpi=executor.resolve_mpi(profile.mpi)
                )
        except BaseException:
            resources.close(suppress_errors=True)
            raise

        try:
            def completed(out):
                resources.close(suppress_errors=True)
                if callback:
                    callback(out)

            potential_file = self.CONTROL.POTFIL.result
            if potential_file and not os.path.isabs(potential_file):
                potential_file = os.path.join(directory, potential_file)

            input_files = tuple(
                file
                for file in (
                    input_file,
                    potential_file,
                    *self.input_files(directory),
                )
                if file is not None
            )
            output_files = tuple(
                file
                for file in (
                    output_file,
                    *self.output_files(directory),
                )
                if file is not None
            )
            cache_files = tuple(
                file for file in self.cache_files(directory) if file is not None
            )

            spec = execution_profile.create_execution(
                task.executable,
                input_file,
                output_file,
                directory,
                input_files=input_files,
                output_files=output_files,
                cache_files=cache_files,
                mpi_supported=task.mpi,
                gdb=gdb,
            )
            if manager is None:
                runner = executor.create_run(reader, spec, callback=completed)
            else:
                runner = manager.create_run(
                    executor,
                    reader,
                    spec,
                    task=self.task_name,
                    dataset=self.CONTROL.DATASET.result,
                    callback=completed,
                )
            if read_callback:
                runner.add_output_listener(read_callback)

            runner = ResourceOwnedRunner(runner, resources)
            handed_off = True
            if run_async:
                return runner
            return runner.run()
        except FileNotFoundError as e:
            add = (
                "Cannot find SPRKKR executable. Check executable_dir and executable_suffix "
                "in the selected SPR-KKR profile.\n"
            )
            if execution_profile.mpi:
                add += (
                    "Also, please check the MPI setting in the selected profile, or explicitly "
                    "disable MPI with mpi=False in SPRKKR.calculate().\n\n"
                )
            e.strerror = add + "SPRKKR cannot be run due to the following error: \n" + e.strerror
            raise
        finally:
            if not handed_off:
                resources.close(suppress_errors=True)

    def output_reader(self, calculator=None, print_output=False, directory=None):
        if not directory and calculator and calculator.directory:
            directory = calculator.directory

        task_definition = self.get_task_definition()
        if task_definition is not None:
            cls = task_definition.output_reader_class(self)
        else:
            from ..outputs.task_result import KkrOutputReader
            task = self.TASK.TASK().lower()
            cls = KkrOutputReader.class_for_task(task)
        return cls(self, calculator, directory, print_output)

    def read_output_from_file(self, filename, directory=None):
        directory = directory or os.path.dirname(filename)
        return self.output_reader(directory=directory).read_from_file(filename)

    def executable_params(self, directory=None):
        """Return ``([executable], stdin, stdout)`` parameters for ``Popen``."""

    def __str__(self):
        return f"<Configuration container {self._get_path()}>"

    def set_option(self, name, value):
        try:
            option = self.get_member(name, is_option=True)
        except KeyError:
            raise KeyError(f"No option with name {name} in any of the members") from None
        option.set(value)

    @classmethod
    def definition(cls, name):
        return cls.task_registry[name].input_parameters_definition

    @cached_class_property
    def definitions(cls):
        return {
            name: cls.task_registry[name].input_parameters_definition
            for name in cls.task_registry.names()
        }

    @classmethod
    def is_it_an_input_parameters_name(cls, name):
        name = str(name).upper()
        return name if name in cls.task_registry else False

    @classmethod
    def create_input_parameters(cls, arg):
        if isinstance(arg, str):
            name = cls.is_it_an_input_parameters_name(arg)
            if name:
                return cls.create(name)
            return cls.from_file(arg)
        if isinstance(arg, io.IOBase):
            return cls.from_file(arg)
        return arg

    @classmethod
    def create(cls, name):
        requested = cls.task_registry.normalize_name(name)
        return cls.task_registry[requested].create_input_parameters(task_name=requested)

    @classmethod
    def default_parameters(cls):
        return cls.create("SCF")

    @classmethod
    def from_file(cls, filename, allow_dangerous: bool = False):
        """Parse using the overridable schema collection, deduplicating aliases."""
        last = None
        is_stream = hasattr(filename, "seek")
        seen = set()
        for definition in cls.definitions.values():
            identity = id(definition)
            if identity in seen:
                continue
            seen.add(identity)
            if is_stream:
                filename.seek(0)
            try:
                return definition.read_from_file(
                    filename, allow_dangerous=allow_dangerous
                )
            except (pp.ParseBaseException, DataValidityError) as e:
                last = e
        if last is None:
            raise ValueError("No SPR-KKR input-parameter definitions are available")
        raise last

    def __repr__(self):
        return self._definition.configuration_type_name + " for task " + self.task_name.upper()

    def change_task(self, task, retain_values=False):
        if retain_values:
            vals = self.to_dict(only_changed=True)

        requested = self.task_registry.normalize_name(task)
        task_definition = self.task_registry[requested]
        self._definition = task_definition.input_parameters_definition
        self._task_name = requested
        self._init_members_from_the_definition()

        if retain_values:
            self.set(
                vals,
                unknown="ignore",
                retain_invalid="none",
                report_invalid="ignore",
            )

    def save_to_file(self, file, atoms=None, *, validate="save"):
        if self._definition.save_hook:
            self._definition.save_hook(getattr(file, "name", None), atoms)
        super().save_to_file(file, atoms, validate=validate)

"""Containers for configuration parameters of SPR-KKR tasks."""

from __future__ import annotations

import io
import os
from typing import TYPE_CHECKING, Union
import pyparsing as pp

from ..outputs.task_result import KkrOutputReader
from ..sprkkr.configuration import ConfigurationFile, ConfigurationSection
from ..common.decorators import cached_class_property
from ..common.warnings import DataValidityError
from ..configuration import config
from ..potentials.potentials import Potential
from ..sprkkr.profile import SPRKKRProfile, executor_from_profile
from ..execution.base import ExecutionResources, ResourceOwnedRunner
from ..execution.persistence import resolve_execution_manager

if TYPE_CHECKING:
    from .task_definition import TaskDefinition


class InputSection(ConfigurationSection):
    """Input-parameter sections currently need no additional behaviour."""


class InputParameters(ConfigurationFile):
    """Concrete configuration values for an SPR-KKR input-file schema."""

    task_registry = None
    """Externally supplied registry of effective :class:`TaskDefinition` objects."""

    definition_registry = None
    """Externally supplied registry of canonical input-parameter definitions."""

    @classmethod
    def _require_task_registry(cls):
        registry = cls.task_registry
        if registry is None:
            raise RuntimeError("InputParameters.task_registry has not been configured")
        return registry

    @property
    def potential(self):
        pot = getattr(self, "_potential", None)
        potfil = self.CONTROL.POTFIL.result
        if pot and potfil == pot[0]:
            return pot[1]
        out = Potential.from_file(potfil) if potfil else None
        self._potential = potfil, out
        return out

    def resolve_executable_suffix(self, suffix: Union[str, bool]):
        """Return the suffix appended to the SPR-KKR executable name."""
        if suffix is False:
            return ""
        if suffix is True:
            return config.executables.suffix()
        return suffix

    def __init__(
        self,
        definition,
        container=None,
        *,
        inputfile=None,
        outputfile=False,
        requested_task_definition: TaskDefinition | None = None,
    ):
        self._inputfile = inputfile
        self._outputfile = inputfile if outputfile is False else outputfile
        self._requested_task_definition = requested_task_definition
        super().__init__(definition, container)

    def get_task_definition(self) -> TaskDefinition:
        """Resolve the effective task from this object's schema and current values."""
        return self._definition.task_definition(self)

    @property
    def task_definition(self) -> TaskDefinition:
        """The effective task represented by the current input values."""
        return self.get_task_definition()

    @property
    def requested_task_definition(self) -> TaskDefinition | None:
        """Task explicitly requested when this object was created, if any.

        This is provenance/intention only. Execution uses ``task_definition``,
        which is resolved from the current values and may differ after editing.
        """
        return self._requested_task_definition

    def copy(self, copy_values: bool = False):
        """Copy values while preserving explicit requested-task provenance."""
        out = self._definition.create_object(
            requested_task_definition=self._requested_task_definition,
        )
        values = self.as_dict(copy=copy_values, only_changed=True)
        out.set(values, unknown="add", retain_invalid="all", report_invalid="ignore")
        return out

    @property
    def task_name(self):
        """Name of the effective task represented by the current values."""
        return self.task_definition.name

    @property
    def requested_task_name(self):
        """Name of the explicitly requested task, or ``None`` for parsed inputs."""
        task = self.requested_task_definition
        return task.name if task is not None else None

    @property
    def _execution_task_name(self):
        """Return the effective task name used for execution bookkeeping."""
        return self.task_name

    def uses_separate_single_site_contour(self):
        """Whether SPR-KKR evaluates single-site terms on a separate contour."""
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
        """Return whether this effective task should be run using MPI."""
        return mpi and self.task_definition.mpi

    def input_files(self, directory):
        """Return task-specific input files for this concrete calculation."""
        return self.task_definition.input_files(self, directory)

    def output_files(self, directory):
        """Return task-specific output files for this concrete calculation."""
        return self.task_definition.output_files(self, directory)

    def cache_files(self, directory):
        """Return task-specific reusable cache files for this calculation."""
        return self.task_definition.cache_files(self, directory)

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
        """Run the effective SPR-KKR task described by these input parameters."""
        resources = ExecutionResources(input_file, output_file)
        handed_off = False
        task = self.task_definition

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
                    task=self._execution_task_name,
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
                "Cannot find SPR-KKR executable. Check executable_dir and executable_suffix "
                "in the selected SPR-KKR profile.\n"
            )
            if execution_profile.mpi:
                add += (
                    "Also, please check the MPI setting in the selected profile, or explicitly "
                    "disable MPI with mpi=False in SPRKKR.calculate().\n\n"
                )
            e.strerror = (
                add
                + "SPRKKR cannot be run due to the following error: \n"
                + e.strerror
            )
            raise
        finally:
            if not handed_off:
                resources.close(suppress_errors=True)

    def output_reader(self, calculator=None, print_output=False, directory=None):
        """Return the effective task's output reader/parser."""
        cls = self.task_definition.result_reader

        if not directory and calculator and calculator.directory:
            directory = calculator.directory

        if cls is None:
            task_section = getattr(self, "TASK", None)
            task_option = (
                getattr(task_section, "TASK", None) if task_section is not None else None
            )
            task = task_option() if task_option is not None else self.task_name
            cls = KkrOutputReader.class_for_task(task.lower())
        return cls(self, calculator, directory, print_output)

    def read_output_from_file(self, filename, directory=None):
        """Read and parse output of a previously run task."""
        directory = directory or os.path.dirname(filename)
        return self.output_reader(directory=directory).read_from_file(filename)

    def executable_params(self, directory=None):
        """Return ``([executable], stdin, stdout)`` parameters for Popen."""

    def __str__(self):
        return f"<Configuration container {self._get_path()}>"

    def set_option(self, name, value):
        """Set the first option matching *name* anywhere in this input file."""
        try:
            option = self.get_member(name, is_option=True)
        except KeyError:
            raise KeyError(
                f"No option with name {name} in any of the members"
            ) from None
        option.set(value)

    @classmethod
    def definition(cls, name):
        """Return the schema used by a named effective task."""
        return cls._require_task_registry()[name].input_parameters_definition

    @cached_class_property
    def definitions():
        """Return canonical input-file schemas when their registry is available."""
        registry = InputParameters.definition_registry
        if registry is not None:
            return {name: registry[name] for name in registry.names(include_aliases=False)}
        tasks = InputParameters._require_task_registry()
        return {
            name: task.input_parameters_definition
            for name, task in tasks.definitions.items()
        }

    @classmethod
    def is_it_an_input_parameters_name(cls, name):
        name = str(name).upper()
        registry = cls.task_registry
        return name if registry is not None and name in registry else False

    @classmethod
    def create_input_parameters(cls, arg):
        """Create or load an :class:`InputParameters` object."""
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
        """Create default input parameters for an explicitly requested task."""
        return cls._require_task_registry()[name].create_input_parameters()

    @classmethod
    def default_parameters(cls):
        """Create default SCF input parameters."""
        return cls.create("SCF")

    @classmethod
    def from_file(cls, filename, allow_dangerous: bool = False):
        """Read an input file without inventing requested-task provenance."""
        last = None
        is_stream = hasattr(filename, "seek")

        registry = cls.definition_registry
        if registry is not None:
            definitions = registry.values(include_aliases=False)
        else:
            seen = set()
            definitions = []
            for task in cls._require_task_registry().definitions.values():
                definition = task.input_parameters_definition
                marker = id(definition)
                if marker not in seen:
                    seen.add(marker)
                    definitions.append(definition)

        for definition in definitions:
            if is_stream:
                filename.seek(0)
            try:
                return definition.read_from_file(
                    filename,
                    allow_dangerous=allow_dangerous,
                    inputfile=filename,
                )
            except (pp.ParseBaseException, DataValidityError) as exc:
                last = exc
        if last is None:
            raise ValueError("No SPR-KKR input-parameter definitions are registered")
        raise last

    def calculate(self, *args, **kwargs):
        """Create a calculator and run these input parameters."""
        calc = calculator.SPRKKR()
        calc.calculate(input_parameters=self, *args, **kwargs)

    def __repr__(self):
        requested = self.requested_task_name
        suffix = f" (requested {requested})" if requested and requested != self.task_name else ""
        return f"{self._definition.configuration_type_name} for task {self.task_name}{suffix}"

    def change_task(self, task, retain_values=False):
        """Explicitly request a new task and replace its input schema."""
        if retain_values:
            vals = self.to_dict(only_changed=True)
        new_task = self._require_task_registry()[task]
        self._definition = new_task.input_parameters_definition
        self._requested_task_definition = new_task
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


# Kept at the bottom to avoid the long-standing calculator/input-parameter cycle.
from ..sprkkr import calculator  # NOQA: E402

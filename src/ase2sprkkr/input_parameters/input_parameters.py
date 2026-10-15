"""Containers for configuration parameters of SPR-KKR tasks."""

from __future__ import annotations

import io
import os
import pyparsing as pp

from . import definitions, tasks
from .task_definition import InputParametersRegistry, TaskRegistry
from ..sprkkr.configuration import ConfigurationFile, ConfigurationSection
from ..common.decorators import cached_class_property
from ..common.warnings import DataValidityError


class InputSection(ConfigurationSection):
    """Input parameters sections have nothing special yet."""


class InputParameters(ConfigurationFile):
    """Concrete values of an SPR-KKR input file.

    Input-file schemas and task identities are separate registries. Standard
    schemas carry their canonical task directly; shared schemas resolve the task
    from their concrete values.
    """

    input_parameters_registry = InputParametersRegistry(definitions)
    task_registry = TaskRegistry(input_parameters_registry, tasks)

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

    @property
    def task_definition(self):
        """Return the task represented by the current parameter values."""
        return self._definition.task_definition(self)

    def _require_task_definition(self):
        task = self.task_definition
        if task is None:
            raise RuntimeError(
                "These input parameters are not associated with an SPR-KKR task definition"
            )
        return task

    @property
    def potential(self):
        from ..potentials.potentials import Potential

        pot = getattr(self, "_potential", None)
        potfil = self.CONTROL.POTFIL.result
        if pot:
            if potfil == pot[0]:
                return pot[1]
        out = Potential.from_file(potfil) if potfil else None
        self._potential = potfil, out
        return out

    def copy(self, copy_values: bool = False):
        """Copy parameters while preserving their concrete task view."""
        values = self.as_dict(copy=copy_values, only_changed=True)
        out = self.__class__(definition=self._definition, task_name=self.task_name)
        out.set(values, unknown="add", retain_invalid="all", report_invalid="ignore")
        return out

    @property
    def task_name(self):
        """Task name originally requested when these parameters were created."""
        return self._task_name

    @property
    def executable(self):
        """Executable basename declared by the effective task."""
        return self._require_task_definition().executable

    @property
    def mpi_supported(self):
        """Whether the effective task has an MPI executable."""
        return self._require_task_definition().mpi

    @property
    def execution_task_name(self):
        """Task identity represented by the concrete parameter values."""
        return self._require_task_definition().execution_task_name(self)

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

    def output_reader(self, calculator=None, print_output=False, directory=None):
        if not directory and calculator and calculator.directory:
            directory = calculator.directory

        cls = self._require_task_definition().output_reader_class(self)
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
        name = cls.input_parameters_registry.normalize_name(name)
        if name in cls.input_parameters_registry:
            return cls.input_parameters_registry[name]
        return cls.task_registry.input_parameters_definition(name)

    @cached_class_property
    def definitions(cls):
        return {
            name: cls.input_parameters_registry[name]
            for name in cls.input_parameters_registry.names()
        }

    @classmethod
    def is_it_an_input_parameters_name(cls, name):
        name = cls.task_registry.normalize_name(name)
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
        return cls.task_registry.create_input_parameters(requested)

    @classmethod
    def default_parameters(cls):
        return cls.create("SCF")

    @classmethod
    def from_file(cls, filename, allow_dangerous: bool = False):
        """Try all registered input-file schemas."""
        last = None
        is_stream = hasattr(filename, "seek")
        for definition in cls.definitions.values():
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
        self._definition = self.task_registry.input_parameters_definition(requested)
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

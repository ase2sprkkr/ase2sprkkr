"""Containers for configuration parameters of SPR-KKR tasks."""

from __future__ import annotations

import io
import pyparsing as pp

from . import definitions, tasks
from .task_definition import InputParametersRegistry, TaskRegistry
from ..sprkkr.configuration import ConfigurationFile, ConfigurationSection
from ..common.decorators import cached_class_property
from ..common.warnings import DataValidityError
from ..potentials.potentials import Potential


class InputSection(ConfigurationSection):
    """Input parameters sections have nothing special yet."""


class InputParameters(ConfigurationFile):
    """Concrete values of an SPR-KKR input file."""

    input_parameters_registry = InputParametersRegistry(definitions)
    task_registry = TaskRegistry(input_parameters_registry, tasks)

    def __init__(
        self,
        definition,
        inputfile=None,
        outputfile=False,
        *,
        requested_task_name=None,
        from_result=None,
    ):
        super().__init__(definition)
        self._requested_task_name = (
            None
            if requested_task_name is None
            else str(requested_task_name).upper()
        )
        self._from_result = from_result
        self._inputfile = inputfile
        self._outputfile = inputfile if outputfile is False else outputfile

    @property
    def task_definition(self):
        """Return the task represented by the current values, if any."""
        return self.task_registry.task_for_parameters(self)

    def task_definition_required(self):
        """Return the task represented by the current values, failing if detached."""
        task = self.task_definition
        if task is None:
            raise RuntimeError(
                "These input parameters are not associated with an SPR-KKR task definition"
            )
        return task

    @property
    def task_name(self):
        """Return the ASE2SPRKKR task represented by the current values."""
        task = self.task_definition
        return None if task is None else task.name

    @property
    def requested_task_name(self):
        """Task requested when these parameters were created or last retargeted.

        Shared input schemas may use this only as context for callable defaults.
        Concrete parameter values remain authoritative for :attr:`task_name`.
        """
        return self._requested_task_name

    @property
    def from_result(self):
        """Result from which these parameters were derived, if any.

        Like :attr:`requested_task_name`, this is provenance and default-value
        context rather than part of the task identity.
        """
        return self._from_result

    @property
    def definition_name(self):
        """Return the registered input-parameter schema name."""
        return self._definition.name

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

    def copy(self, copy_values: bool = False):
        """Copy parameters while preserving provenance and default context."""
        values = self.as_dict(copy=copy_values, only_changed=True)
        out = self.__class__(
            definition=self._definition,
            requested_task_name=self.requested_task_name,
            from_result=self.from_result,
        )
        out.set(values, unknown="add", retain_invalid="all", report_invalid="ignore")
        return out

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
        task = self.task_definition
        return task is not None and task.name == "COMPTON"

    def _task_files(self, directory, kind):
        task = self.task_definition_required()
        return self.input_parameters_registry.files_for_task(
            task, self, directory, kind
        )

    def input_files(self, directory):
        return self._task_files(directory, "input")

    def output_files(self, directory):
        return self._task_files(directory, "output")

    def cache_files(self, directory):
        return self._task_files(directory, "cache")

    def output_reader(self, calculator=None, print_output=False, directory=None):
        """Create the reader for the effective task represented by these values."""
        return self.task_definition_required().output_reader(
            input_parameters=self,
            calculator=calculator,
            print_output=print_output,
            directory=directory,
        )

    def read_output_from_file(self, filename, directory=None):
        """Read a completed output through the effective task's registered reader."""
        return self.task_definition_required().read_output_from_file(
            filename,
            input_parameters=self,
            directory=directory,
        )

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
        normalized = cls.input_parameters_registry.normalize_name(name)
        if normalized in cls.input_parameters_registry:
            return cls.input_parameters_registry[normalized]
        task = cls.task_registry[normalized]
        return cls.input_parameters_registry.definition_for_task(task)

    @cached_class_property
    def definitions(cls):
        return {
            name: cls.input_parameters_registry[name]
            for name in cls.input_parameters_registry.names()
        }

    @classmethod
    def is_task_name(cls, name):
        """Return the normalized task name when *name* identifies a task."""
        name = cls.task_registry.normalize_name(name)
        return name if name in cls.task_registry else False

    @classmethod
    def create_input_parameters(cls, arg):
        if isinstance(arg, str):
            name = cls.is_task_name(arg)
            if name:
                return cls.create(name)
            return cls.from_file(arg)
        if isinstance(arg, io.IOBase):
            return cls.from_file(arg)
        return arg

    @classmethod
    def create(cls, name):
        task = cls.task_registry[name]
        return cls.input_parameters_registry.create_for_task(task)

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
        task = self.task_definition
        if task is None:
            return self._definition.configuration_type_name + " (detached schema)"
        return self._definition.configuration_type_name + " for task " + task.name

    def change_task(self, task, retain_values=False):
        if retain_values:
            vals = self.to_dict(only_changed=True)

        requested = self.task_registry.normalize_name(task)
        task_definition = self.task_registry[requested]
        self._definition = self.input_parameters_registry.definition_for_task(
            task_definition
        )
        self._requested_task_name = requested
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

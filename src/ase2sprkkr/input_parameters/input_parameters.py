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
        task_name=None,
    ):
        super().__init__(definition)
        self._task_name = None if task_name is None else str(task_name).upper()
        self._inputfile = inputfile
        self._outputfile = inputfile if outputfile is False else outputfile

    @property
    def task_definition(self):
        """Return the effective task represented by the current values."""
        return self.task_registry.effective_task(self)

    def _require_task_definition(self):
        task = self.task_definition
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

    def copy(self, copy_values: bool = False):
        """Copy parameters while preserving requested-task provenance."""
        values = self.as_dict(copy=copy_values, only_changed=True)
        out = self.__class__(definition=self._definition, task_name=self.task_name)
        out.set(values, unknown="add", retain_invalid="all", report_invalid="ignore")
        return out

    @property
    def task_name(self):
        """Originally requested task, or ``None`` for a detached schema object."""
        return self._task_name

    @property
    def execution_task_name(self):
        """Return the concrete SPR-KKR TASK value used for execution."""
        try:
            value = self["TASK"]["TASK"]()
        except KeyError:
            value = None
        if value:
            return str(value).upper()
        return self._require_task_definition().name

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
        task = self._require_task_definition()
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
        from ..outputs.reader_registry import reader_registry

        if not directory and calculator and calculator.directory:
            directory = calculator.directory
        cls = reader_registry.class_for_task(self._require_task_definition())
        return cls(self, calculator, directory, print_output)

    def read_output_from_file(self, filename, directory=None):
        directory = directory or os.path.dirname(filename)
        return self.output_reader(directory=directory).read_from_file(filename)

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
        self._task_name = requested
        self._init_members_from_the_definition()

        if retain_values:
            self.set(
                vals,
                unknown="ignore",
                retain_invalid="none",
                report_invalid="ignore",
            )

        self._definition.initialize_for_task(self, requested)

    def save_to_file(self, file, atoms=None, *, validate="save"):
        if self._definition.save_hook:
            self._definition.save_hook(getattr(file, "name", None), atoms)
        super().save_to_file(file, atoms, validate=validate)

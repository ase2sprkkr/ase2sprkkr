"""Containers for configuration parameters of SPR-KKR tasks."""

from __future__ import annotations

import io
import pyparsing as pp

from . import definitions, tasks
from .task_definition import InputParametersRegistry, TaskRegistry
from ..sprkkr.configuration import (
    ConfigurationFile,
    ConfigurationSection,
    ConfigurationValue,
)
from ..common.configuration_containers import BaseConfigurationContainer
from ..common.repeated_configuration_containers import RepeatedConfigurationContainer
from ..common.decorators import cached_class_property
from ..common.warnings import DataValidityError
from ..potentials.potentials import Potential


class InputSection(ConfigurationSection):
    """Input parameters sections have nothing special yet."""


class InputValue(ConfigurationValue):
    """Input-file option supporting a default specific to one parameter instance."""

    @property
    def default_value(self):
        if hasattr(self, "_default_override"):
            return self._default_override
        return super().default_value

    def set_default(self, value):
        if self._definition.is_generated:
            raise TypeError("Generated options cannot have instance defaults")
        if self._definition.default_value is None:
            raise ValueError(
                "An instance default can only override an existing definition default"
            )
        self._default_override = None if value is None else self._pack_value(value)

    def clear_default(self):
        self.__dict__.pop("_default_override", None)


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
    ):
        super().__init__(definition)
        self._requested_task_name = (
            None
            if requested_task_name is None
            else str(requested_task_name).upper()
        )
        self._inputfile = inputfile
        self._outputfile = inputfile if outputfile is False else outputfile

    @property
    def task_name(self):
        """Return the ASE2SPRKKR task represented by the current values."""
        name = self._definition.task_name(self)
        if name is None:
            return None
        return str(name).strip().upper()

    @property
    def task_definition(self):
        """Return the task represented by the current values, if any."""
        name = self.task_name
        return None if name is None else self.task_registry[name]

    def task_definition_required(self):
        task = self.task_definition
        if task is None:
            raise RuntimeError(
                "These input parameters are not associated with an SPR-KKR task definition"
            )
        return task

    @property
    def requested_task_name(self):
        """Task requested when these parameters were created or last retargeted."""
        return self._requested_task_name

    @property
    def definition_name(self):
        return self.input_parameters_registry.normalize_name(self._definition.name)

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

    @staticmethod
    def _has_instance_defaults(source):
        if isinstance(source, InputValue):
            return hasattr(source, "_default_override")
        if isinstance(source, RepeatedConfigurationContainer):
            return any(
                InputParameters._has_instance_defaults(member)
                for member in source._raw_values()
            )
        if isinstance(source, BaseConfigurationContainer):
            return any(
                InputParameters._has_instance_defaults(member)
                for member in source._members.values()
            )
        return False

    @staticmethod
    def _copy_instance_defaults(source, target, copy_values):
        if isinstance(source, RepeatedConfigurationContainer):
            for key, source_member in source._raw_items():
                try:
                    target_member = target[key]
                except (IndexError, KeyError, TypeError):
                    if not InputParameters._has_instance_defaults(source_member):
                        continue
                    if isinstance(target._values, list):
                        while len(target._values) <= key:
                            target._add()
                    else:
                        target._add(key)
                    target_member = target[key]
                InputParameters._copy_instance_defaults(
                    source_member, target_member, copy_values
                )
            return

        if isinstance(source, BaseConfigurationContainer):
            for name, source_member in source._members.items():
                target_member = target._members.get(name)
                if target_member is None:
                    continue
                InputParameters._copy_instance_defaults(
                    source_member, target_member, copy_values
                )
            return

        if not isinstance(source, InputValue):
            return
        if not hasattr(source, "_default_override"):
            return
        value = source._default_override
        if copy_values and value is not None:
            value = source._definition.copy_value(value, all_values=True)
        target._default_override = value

    def copy(self, copy_values: bool = False):
        values = self.as_dict(copy=copy_values, only_changed="explicit")
        out = self.__class__(
            definition=self._definition,
            requested_task_name=self.requested_task_name,
        )
        out.set(values, unknown="add", retain_invalid="all", report_invalid="ignore")
        self._copy_instance_defaults(self, out, copy_values)
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
        return self.task_definition_required().output_reader(
            input_parameters=self,
            calculator=calculator,
            print_output=print_output,
            directory=directory,
        )

    def read_output_from_file(self, filename, directory=None):
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
        """Try registered schemas, preferring a task matching ``CONTROL.ADSI``."""
        last = None
        first_match = None
        is_stream = hasattr(filename, "seek")
        for definition in cls.definitions.values():
            if is_stream:
                filename.seek(0)
            try:
                parsed = definition.read_from_file(
                    filename, allow_dangerous=allow_dangerous
                )
            except (pp.ParseBaseException, DataValidityError) as e:
                last = e
                continue

            try:
                adsi = parsed.CONTROL.ADSI()
            except (AttributeError, KeyError):
                return parsed
            if not isinstance(adsi, str):
                return parsed

            adsi = cls.task_registry.normalize_name(adsi)
            if adsi not in cls.task_registry:
                return parsed
            if parsed.task_name == adsi:
                return parsed
            if first_match is None:
                first_match = parsed

        if first_match is not None:
            return first_match
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
            vals = self.to_dict(only_changed="explicit")

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

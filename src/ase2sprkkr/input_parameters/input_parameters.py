"""Containers for configuration parameters of SPR-KKR tasks."""

from __future__ import annotations

import io
from collections.abc import Mapping
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
from ..common.module_registry import RegistryEntryNotFound
from ..potentials.potentials import Potential

class _LazyInputParametersDefinitions(Mapping):
    """Read-only compatibility mapping that materializes only requested schemas.

    User-installed definitions can be malformed or require optional packages.
    Merely requesting a different definition must not import those modules.
    """

    def __init__(self, registry):
        self._registry = registry

    def __getitem__(self, name):
        return self._registry[name]

    def __iter__(self):
        return iter(self._registry.names())

    def __len__(self):
        return len(self._registry.names())

    def __contains__(self, name):
        return name in self._registry


class InputSection(ConfigurationSection):
    """Input parameters sections have nothing special yet."""


class InputValue(ConfigurationValue):
    """Input-file option supporting a default specific to one parameter instance."""

    @property
    def default_value(self):
        """Return an instance-specific override of the definition default."""
        if hasattr(self, "_default_override"):
            return self._default_override
        return super().default_value

    def set_default(self, value):
        """Override this instance's definition default without setting a value."""
        if self._definition.is_generated:
            raise TypeError("Generated options cannot have instance defaults")
        if self._definition.default_value is None:
            raise ValueError(
                "An instance default can only override an existing definition default"
            )
        self._default_override = None if value is None else self._pack_value(value)

    def clear_default(self):
        """Restore the default supplied by the input-parameter definition."""
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
        name = self._definition.task_definition_name(self)
        if name is None:
            return None
        return str(name).strip().upper()

    @property
    def task_definition(self):
        """Return the task represented by the current values, if any."""
        name = self.task_name
        return None if name is None else self.task_registry[name]

    def task_definition_required(self):
        """Return the task represented by the current values, failing if detached."""
        task = self.task_definition
        if task is None:
            raise RuntimeError(
                "These input parameters are not associated with an SPR-KKR task definition"
            )
        expected = self.input_parameters_registry.definition_for_task(task)
        if (expected is not self._definition
                and expected.name.upper() != self._definition.derived_from_name.upper()):
            raise RuntimeError(
                f"Task {task.name} uses a different input schema; "
                "use change_task() to switch schemas explicitly."
            )
        return task

    @property
    def requested_task_name(self):
        """Task requested when these parameters were created or last retargeted.

        Shared input schemas may use this only as context for callable defaults.
        Concrete parameter values remain authoritative for :attr:`task_name`.
        """
        return self._requested_task_name

    @property
    def definition_name(self):
        """Return the normalized registered input-parameter schema name."""
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
        """Copy parameters, including instance-specific defaults."""
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
        return _LazyInputParametersDefinitions(cls.input_parameters_registry)

    @classmethod
    def is_task_name(cls, name):
        """Return the normalized task name when *name* identifies a task."""
        name = cls.task_registry.normalize_name(name)
        return name if name in cls.task_registry else False

    @classmethod
    def for_definition(cls, name):
        """Create input parameters using a named input-parameter definition.

        Unlike :meth:`for_task`, this also accepts shared schemas such
        as BSF and user-defined presets.
        """
        return cls.input_parameters_registry[name].create_object()

    @classmethod
    def for_task(cls, name):
        """Create input parameters for a registered SPR-KKR task."""
        task = cls.task_registry[name]
        return cls.input_parameters_registry.create_for_task(task)

    @classmethod
    def create(cls, arg):
        """Create or obtain input parameters from a task, schema, or file.

        Strings naming both a task and a schema are resolved as tasks.
        Unknown strings are treated as filenames. Existing parameter
        instances are returned unchanged.
        """
        if isinstance(arg, str):
            name = cls.task_registry.normalize_name(arg)
            if name in cls.task_registry:
                return cls.for_task(name)
            if name in cls.input_parameters_registry:
                return cls.for_definition(name)
            return cls.from_file(arg)
        if isinstance(arg, io.IOBase):
            return cls.from_file(arg)
        return arg

    @classmethod
    def default_parameters(cls):
        return cls.for_task("SCF")

    @classmethod
    def from_file(cls, filename, allow_dangerous: bool = False, *,
                  allow_incomplete: bool = False, report_invalid=None):
        """Read with the IPD class-level selector autodetection."""
        from .input_parameters_definitions import InputParametersDefinition

        return InputParametersDefinition.from_file(
            filename, allow_dangerous=allow_dangerous,
            allow_incomplete=allow_incomplete, report_invalid=report_invalid,
        )

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

            # Explicit TASK.TASK may select a different canonical schema.
            # Preserve requested provenance while keeping the effective task
            # and its input definition consistent.
            effective = self.task_definition
            if effective is not None:
                canonical = self.input_parameters_registry.definition_for_task(effective)
                if canonical is not self._definition:
                    values = self.to_dict(only_changed="explicit")
                    self._definition = canonical
                    self._init_members_from_the_definition()
                    self.set(
                        values,
                        unknown="ignore",
                        retain_invalid="none",
                        report_invalid="ignore",
                    )

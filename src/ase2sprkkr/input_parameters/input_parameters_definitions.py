"""Definitions of SPR-KKR input-parameter file formats.

An :class:`InputParametersDefinition` describes an SPR-KKR input-file schema.
Execution, result-reader and auxiliary-file policy belong to
:class:`TaskDefinition`.  The common one-task case keeps a direct reference to
its task.  Shared schemas resolve an effective task name from concrete values
and delegate object lookup to the task registry supplied by ``InputParameters``.
"""

import pyparsing as pp
from ..common.container_definitions import dict_from_parsed
from ..sprkkr.configuration import (
    ConfigurationValueDefinition,
    ConfigurationSectionDefinition,
    ConfigurationFileDefinition,
    CustomConfigurationValue,
    CustomConfigurationSection,
)
from ..common.grammar_types import mixed, flag
from ..common.grammar import generate_grammar, delimitedList
from ..common.decorators import cached_class_property, cache
from .input_parameters import InputParameters, InputSection

with generate_grammar():
    section_line_ends = pp.Regex("(\n[ \t]*)*[ \t]")


class InputValueDefinition(ConfigurationValueDefinition):
    """Describe the format of one input-parameter value."""

    @cached_class_property
    def delimiter():
        with generate_grammar():
            return pp.Suppress("=").set_name("=")

    prefix = "\t"
    type_from_type_map = {bool: flag}
    type_of_dangerous = mixed

    def copy(self, **kwargs):
        """Copy input-specific definition metadata not represented by constructor args."""
        copied = super().copy(**kwargs)
        if hasattr(self, "_energy_default_value"):
            copied._energy_default_value = self._energy_default_value
        return copied


class InputSectionDefinition(ConfigurationSectionDefinition):
    """Describe the format of one input-parameter section."""

    child_class = InputValueDefinition
    result_class = InputSection
    custom_class = staticmethod(CustomConfigurationValue.factory(InputValueDefinition, mixed))

    @cached_class_property
    def delimiter():
        with generate_grammar():
            out = (pp.Optional(section_line_ends) + pp.WordStart()).suppress()
            return out.set_name("\n")

    do_not_skip_whitespaces_before_name = True


class InputParametersDefinition(ConfigurationFileDefinition):
    """Describe the syntax and semantics of an SPR-KKR input file.

    A normal schema corresponds to one task and stores a direct reference to
    that :class:`TaskDefinition`.  A shared schema sets ``dynamic_task = True``
    and overrides :meth:`task_definition_name` to derive the effective task from
    concrete values.  Registry construction is intentionally outside this class;
    dynamic resolution uses the registry supplied by the runtime
    :class:`InputParameters` object.
    """

    save_hook = None
    child_class = InputSectionDefinition
    result_class = InputParameters
    custom_class = staticmethod(CustomConfigurationSection.factory(InputSectionDefinition))
    configuration_type_name = "INPUT PARAMETERS"
    dynamic_task = False

    def __init__(self, *args, task_definition=None, **kwargs):
        self._task_definition = task_definition
        super().__init__(*args, **kwargs)

    def _associate_task_definition(self, task_definition):
        """Associate the ordinary one-task schema with its task.

        Shared schemas deliberately keep no local collection of tasks; their
        effective task is obtained through the externally supplied task registry.
        """
        if self.dynamic_task:
            return task_definition
        if self._task_definition is None:
            self._task_definition = task_definition
        elif self._task_definition is not task_definition:
            raise ValueError(
                f"Input-parameter definition {self.name!r} is already associated "
                f"with task {self._task_definition.name}"
            )
        return task_definition

    def task_definition_name(self, parameters):
        """Return the effective task name represented by *parameters*.

        The base implementation is only for the common one-task case. Shared
        definitions override this method and leave ``_task_definition`` unset.
        """
        if self._task_definition is None:
            raise RuntimeError(
                f"Input-parameter definition {self.name!r} must override "
                "task_definition_name(parameters)"
            )
        return self._task_definition.name

    def task_definition(self, parameters):
        """Resolve the effective :class:`TaskDefinition` for *parameters*."""
        if self._task_definition is not None:
            return self._task_definition

        registry = getattr(parameters, "task_registry", None)
        if registry is None:
            raise RuntimeError(
                f"Input-parameter definition {self.name!r} requires an external "
                "task registry to resolve its effective task"
            )
        return registry[self.task_definition_name(parameters)]

    def copy(self, *args, **kwargs):
        """Copy the schema without copying a direct task association."""
        copied = super().copy(*args, **kwargs)
        copied._task_definition = None
        return copied

    def create_object(self, container=None, **kwargs):
        """Create concrete input parameters from this schema."""
        self._finalize()
        return self.result_class(self, container, **kwargs)

    @cached_class_property
    def delimiter():
        with generate_grammar():
            def ws(x):
                return x.set_whitespace_chars("")

            after_consumed_line_end = pp.Empty().add_condition(
                lambda s, loc, _t: loc > 0 and s[:loc].rstrip(" \t\r").endswith("\n")
            )

            out = (
                pp.Optional(section_line_ends)
                + pp.OneOrMore(ws(pp.LineEnd()))
                + pp.FollowedBy(ws(pp.Regex(r"[^\s]")))
                | after_consumed_line_end
            ).suppress()
            return out.set_name("\n")

    @classmethod
    @cache
    def custom_value_grammar(cls):
        value = cls.child_class.custom_member_grammar()
        delim = cls.child_class.delimiter
        return delimitedList(value, delim).set_parse_action(
            lambda x: dict_from_parsed(x.asList())
        )

    def _generic_info(self):
        return f"Input parameters file format {self.name}"

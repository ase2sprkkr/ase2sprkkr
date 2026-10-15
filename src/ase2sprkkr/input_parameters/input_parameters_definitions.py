"""Definitions of SPR-KKR input-parameter file formats."""

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
    @cached_class_property
    def delimiter():
        with generate_grammar():
            return pp.Suppress("=").set_name("=")

    prefix = "\t"
    type_from_type_map = {bool: flag}
    type_of_dangerous = mixed

    def copy(self, **kwargs):
        copied = super().copy(**kwargs)
        if hasattr(self, "_energy_default_value"):
            copied._energy_default_value = self._energy_default_value
        return copied


class InputSectionDefinition(ConfigurationSectionDefinition):
    child_class = InputValueDefinition
    result_class = InputSection
    custom_class = staticmethod(
        CustomConfigurationValue.factory(InputValueDefinition, mixed)
    )

    @cached_class_property
    def delimiter():
        with generate_grammar():
            out = (pp.Optional(section_line_ends) + pp.WordStart()).suppress()
            return out.set_name("\n")

    do_not_skip_whitespaces_before_name = True


class InputParametersDefinition(ConfigurationFileDefinition):
    """Syntax and parameter semantics of one SPR-KKR input-file schema."""

    _task_definition = None
    """Canonical 1:1 task, or None for schemas shared by multiple tasks."""

    save_hook = None
    child_class = InputSectionDefinition
    result_class = InputParameters
    custom_class = staticmethod(
        CustomConfigurationSection.factory(InputSectionDefinition)
    )
    configuration_type_name = "INPUT PARAMETERS"

    def task_definition_name(self, parameters):
        """Return the task name represented by *parameters*.

        Shared schemas override this method. For standard schemas it is normally
        unused because ``_task_definition`` is set directly.
        """
        return str(self.name).upper()

    def task_definition(self, parameters):
        """Resolve the task represented by concrete parameter values."""
        if self._task_definition is not None:
            return self._task_definition
        return parameters.task_registry[self.task_definition_name(parameters)]

    def initialize_for_task(self, parameters, task_name):
        """Initialize an explicit task discriminator when the schema has one.

        Shared schemas commonly encode their concrete task in ``TASK.TASK``.
        When that option accepts the requested task, make the requested task an
        explicit value. Schemas using another discriminator (for example BSF)
        keep their own task-dependent defaults and may override this hook.
        """
        if self._task_definition is not None:
            return
        try:
            option = parameters["TASK"]["TASK"]
            grammar_type = option._definition.type
            keywords = grammar_type.keywords
        except (KeyError, AttributeError):
            return
        value = grammar_type.convert(task_name)
        if value in keywords:
            option.set(value)

    @classmethod
    @cache
    def custom_value_grammar(cls):
        value = cls.child_class.custom_member_grammar()
        delim = cls.child_class.delimiter
        return delimitedList(value, delim).set_parse_action(
            lambda x: dict_from_parsed(x.asList())
        )

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

    def _generic_info(self):
        return f"Input parameters format {self.name}"

"""
This file contains the classes for definitions of InputParameterss:
the list of sections and their allowed (or standard) options
and their value formats. InputParameterss and potentials have a simlilar
structure, so they share common functionalities from
sprkkr.common.configuration_definitions
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
    """This class describes the format of one value of
    a task configuration"""

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
        if hasattr(self, '_energy_default_value'):
            copied._energy_default_value = self._energy_default_value
        return copied


class InputSectionDefinition(ConfigurationSectionDefinition):
    """This class describes the format of one
    value of a task section"""

    child_class = InputValueDefinition
    """ standard child class """

    result_class = InputSection
    """ The standard class for InputParameters section """

    custom_class = staticmethod(CustomConfigurationValue.factory(InputValueDefinition, mixed))
    """ Factory for custom values in the input sections. """

    @cached_class_property
    def delimiter():
        with generate_grammar():
            out = (pp.Optional(section_line_ends) + pp.WordStart()).suppress()
            return out.set_name("\n")

    do_not_skip_whitespaces_before_name = True


class InputParametersDefinition(ConfigurationFileDefinition):
    """Describe only the syntax and semantics of an SPR-KKR input file.

    Task-level concerns such as executable selection, MPI support, result
    readers and auxiliary file manifests belong to :class:`TaskDefinition`.
    """

    save_hook = None
    """Input parameters can have a save hook run before serialization."""

    child_class = InputSectionDefinition
    """ Sections of the :class:`InputParameters` are defined by :class:`InputSectionDefinition` """

    result_class = InputParameters
    """ Parsing a file directly through a schema yields detached InputParameters. """

    custom_class = staticmethod(CustomConfigurationSection.factory(InputSectionDefinition))
    """ The class factory for custom sections in the container """

    configuration_type_name = "INPUT PARAMETERS"
    """ Name of the container type in the runtime documentation """

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
        return delimitedList(value, delim).set_parse_action(lambda x: dict_from_parsed(x.asList()))

    def _generic_info(self):
        return f"Input parameters format {self.name}"

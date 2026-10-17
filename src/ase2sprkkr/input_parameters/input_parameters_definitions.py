"""Definitions of SPR-KKR input-parameter file formats."""

import os
from functools import cache

import pyparsing as pp
from ..common.grammar_types import Keyword
from ..common.decorators import maybeclassmethod
from ..common.warnings import DataValidityError
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
from .input_parameters import InputParameters, InputSection, InputValue

with generate_grammar():
    section_line_ends = pp.Regex("(\n[ \t]*)*[ \t]")


MAX_INPUT_FILE_SIZE = 1024 * 1024


def _read_input_text(file):
    """Return a UTF-8 snapshot no larger than MAX_INPUT_FILE_SIZE bytes."""
    if isinstance(file, (str, bytes, os.PathLike)):
        with open(file, "rb") as stream:
            content = stream.read(MAX_INPUT_FILE_SIZE + 1)
    elif hasattr(file, "read"):
        content = file.read(MAX_INPUT_FILE_SIZE + 1)
    else:
        raise TypeError("Input file must be a path or a readable stream")

    if isinstance(content, str):
        data = content.encode("utf-8")
        if len(data) > MAX_INPUT_FILE_SIZE:
            raise ValueError(
                f"SPR-KKR input exceeds {MAX_INPUT_FILE_SIZE} bytes"
            )
        return content
    if not isinstance(content, bytes):
        raise TypeError("Input stream must return str or bytes")
    if len(content) > MAX_INPUT_FILE_SIZE:
        raise ValueError(
            f"SPR-KKR input exceeds {MAX_INPUT_FILE_SIZE} bytes"
        )
    return content.decode("utf-8")


@cache
def _task_selector_map():
    """A single map for the singleton input-definition registry."""
    selectors = {}
    for definition in InputParameters.definitions.values():
        for token in definition.accepted_task_selectors():
            if not token:
                raise ValueError(f"Empty TASK selector in IPD {definition.name}")
            selectors.setdefault(token, [])
            if all(definition is not other for other in selectors[token]):
                selectors[token].append(definition)
    return {key: tuple(value) for key, value in selectors.items()}


@cache
def _task_preparse_definition():
    """Build the minimal IPD and its complete grammar exactly once."""
    values = tuple(sorted(_task_selector_map()))
    if not values:
        raise ValueError("No registered input-parameter TASK selectors")
    ipd = _CachedPreparseIPD(
        "_TASK_PREPARSE",
        [InputSectionDefinition(
            "TASK",
            [InputValueDefinition(
                "TASK",
                Keyword(*values),
                name_in_grammar=False,
                is_optional=True,
                is_required=False,
                is_repeated=True,
                repeated_delimiter=pp.Empty(),
            )],
            is_optional=True,
        )],
    )
    ipd.grammar()  # Compile and cache the completed minimal grammar now.
    return ipd


def _read_selected_definitions(text, selectors, *, allow_dangerous=False):
    """Try only IPDs nominated by detected native TASK keywords."""
    mapping = _task_selector_map()
    candidates = []
    for selector in selectors:
        selector = str(selector).upper()
        if selector not in mapping:
            raise ValueError(f"Unknown TASK selector: {selector}")
        for ipd in mapping[selector]:
            if all(ipd is not candidate for candidate in candidates):
                candidates.append(ipd)

    def parse(definition):
        return definition.read_from_string(
            text, allow_dangerous=allow_dangerous
        )

    if len(candidates) == 1:
        return parse(candidates[0])

    successful = []
    failures = []
    for definition in candidates:
        try:
            result = parse(definition)
        except (pp.ParseBaseException, DataValidityError) as exc:
            failures.append((definition, exc))
        else:
            successful.append((definition, result))
    if len(successful) == 1:
        return successful[0][1]
    if len(successful) > 1:
        names = ", ".join(definition.name for definition, _ in successful)
        raise ValueError(f"Ambiguous TASK selectors, multiple IPDs accept input: {names}")
    if failures:
        names = ", ".join(definition.name for definition, _ in failures)
        raise ValueError(f"No selected IPD can parse the input: {names}") from failures[-1][1]
    raise ValueError("No input-parameter definition selected")


class InputValueDefinition(ConfigurationValueDefinition):
    result_class = InputValue

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
    """Canonical 1:1 task declared together with this schema.

    This is declaration-time metadata for the simple case where one input
    schema represents exactly one task. Schemas shared by multiple tasks leave
    it as ``None`` and override :meth:`task_definition_name` using their values.

    The association is internal and must not be changed after registration.
    """

    child_class = InputSectionDefinition
    result_class = InputParameters
    custom_class = staticmethod(
        CustomConfigurationSection.factory(InputSectionDefinition)
    )
    configuration_type_name = "INPUT PARAMETERS"

    task_selectors = None
    """Native TASK keywords selecting this IPD. Defaults to the IPD's name."""

    def accepted_task_selectors(self):
        return tuple(
            str(name).upper()
            for name in (self.task_selectors if self.task_selectors is not None
                         else (self.name,))
        )

    @maybeclassmethod
    def from_file(self, cls, file, allow_dangerous=False):
        """Read one input; on the class autodetect an IPD, on an instance use it.

        Read a bounded snapshot of the input once, independently of path/stream.
        A missing TASK selector selects SCF, not an ADSI-derived definition.
        """
        text = _read_input_text(file)
        if self is not None:
            return self.read_from_string(text, allow_dangerous=allow_dangerous)

        preview = _task_preparse_definition().read_from_string(
            text, allow_dangerous=allow_dangerous
        )
        selectors = preview.TASK.TASK()
        if not selectors:
            selectors = ("SCF",)
        return _read_selected_definitions(
            text, selectors, allow_dangerous=allow_dangerous
        )

    def task_definition_name(self, parameters):
        """Return the ASE2SPRKKR task name represented by *parameters*.

        The standard one-schema/one-task case is answered directly by the
        schema declaration. Shared schemas override this method and derive the
        identity from concrete parameter values. Resolving the returned name to
        a :class:`TaskDefinition` remains the responsibility of ``TaskRegistry``.
        """
        task = self._task_definition
        return None if task is None else task.name

    def initialize_from_result(self, parameters, result):
        """Initialize instance-specific defaults derived from a previous result."""

    def input_files(self, parameters, directory):
        """Return schema-specific auxiliary input files."""
        return ()

    def output_files(self, parameters, directory):
        """Return schema-specific output files."""
        return ()

    def cache_files(self, parameters, directory):
        """Return schema-specific optional cache files."""
        return ()

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


class _CachedPreparseIPD(InputParametersDefinition):
    """Preparser only: no mutations after grammar compilation."""

    @cache
    def grammar(self, allow_dangerous=False):
        return super().grammar(allow_dangerous)

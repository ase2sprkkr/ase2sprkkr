"""Definitions of SPR-KKR input-parameter file formats."""

import os
from io import StringIO

import pyparsing as pp
from ..common.grammar_types import Keyword
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
from ..common.decorators import cached_class_property, cache, maybeclassmethod
from ..common.file_utils import filename_from_file
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


def _snapshot_stream(text, filename):
    """Reparse the same bytes without losing the origin of a named input."""
    stream = StringIO(text)
    if filename is not None:
        stream.name = filename
    return stream


class _TaskSelectorMap(dict):
    """Map TASK selectors to IPDs, allowing filtering by their origin."""

    def definitions(self, selector, *, user_defined=None):
        """Return all, built-in, or user-defined candidates for a selector."""
        candidates = self[selector]
        if user_defined is None:
            return candidates
        if type(user_defined) is not bool:
            raise TypeError("user_defined must be True, False, or None")
        return tuple(
            definition for definition in candidates
            if definition.is_user_defined is user_defined
        )


@cache
def _task_selector_map():
    """A single map for the singleton input-definition registry."""
    selectors = {}
    for definition in InputParameters.input_parameters_registry.iter_canonical_values():
        for token in definition.accepted_task_selectors():
            if not token:
                raise ValueError(f"Empty TASK selector in IPD {definition.name}")
            selectors.setdefault(token, [])
            if all(definition is not other for other in selectors[token]):
                selectors[token].append(definition)
    return _TaskSelectorMap({key: tuple(value) for key, value in selectors.items()})


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


def _read_selected_definitions(text, selectors, *, filename=None,
                               allow_dangerous=False, allow_incomplete=False,
                               report_invalid=None):
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
        return definition.read_from_file(
            _snapshot_stream(text, filename),
            allow_dangerous=allow_dangerous,
            **({"allow_incomplete": True} if allow_incomplete else {}),
            **({"report_invalid": report_invalid} if report_invalid is not None else {}),
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
        details = "\n".join(
            f"  {definition.name}: {type(error).__name__}: {error}"
            for definition, error in failures
        )
        raise ValueError(
            f"No selected IPD can parse the input:\n{details}"
        ) from failures[-1][1]
    raise ValueError("No input-parameter definition selected")


class InputValueDefinition(ConfigurationValueDefinition):
    result_class = InputValue

    def _resolve_is_optional(self, is_optional):
        return True if is_optional is None else is_optional

    @cached_class_property
    def delimiter():
        with generate_grammar():
            return pp.Suppress("=").set_name("=")

    prefix = "\t"
    type_from_type_map = {bool: flag}
    type_of_dangerous = mixed


class InputSectionDefinition(ConfigurationSectionDefinition):
    child_class = InputValueDefinition
    result_class = InputSection

    def _resolve_is_optional(self, is_optional):
        return True if is_optional is None else is_optional
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

    is_user_defined = True
    """Whether the definition comes from user code.

    Standalone definitions default to True; registration overrides this
    according to the defining module's origin.
    """

    _derived_from_name = None
    """ If one copies the definition, it should be listed under its new name, but
    retain the original name as well, to retain the information, to which "task"
    it originally belongs """


    child_class = InputSectionDefinition
    result_class = InputParameters
    custom_class = staticmethod(
        CustomConfigurationSection.factory(InputSectionDefinition)
    )
    configuration_type_name = "INPUT PARAMETERS"

    task_selectors = None
    """Native TASK selectors, defaulting to the IPD name.

    A string denotes one selector; an iterable denotes multiple selectors.
    Like other registration metadata, this must not change once the singleton
    registry has constructed the cached selector map and preparse grammar.
    """

    def accepted_task_selectors(self):
        names = self.task_selectors
        if names is None:
            names = (self.name,)
        elif isinstance(names, str):
            names = (names,)
        names = tuple(names)
        if not names:
            raise ValueError(f"IPD {self.name}: task_selectors may not be empty")
        normalized = []
        for name in names:
            if not isinstance(name, str) or not name.strip():
                raise ValueError(
                    f"IPD {self.name}: TASK selectors must be nonempty strings"
                )
            normalized.append(name.strip().upper())
        if len(set(normalized)) != len(normalized):
            raise ValueError(f"IPD {self.name}: duplicate TASK selectors")
        return tuple(normalized)

    @maybeclassmethod
    def from_file(self, cls, file, allow_dangerous=False, *,
                  allow_incomplete=False, report_invalid=None):
        """Read one input; on the class autodetect an IPD, on an instance use it.

        Read a bounded snapshot of the input once, independently of path/stream.
        A missing TASK selector selects SCF, not an ADSI-derived definition.
        """
        filename = filename_from_file(file, None)
        text = _read_input_text(file)
        if self is not None:
            return self.read_from_file(
                _snapshot_stream(text, filename),
                allow_dangerous=allow_dangerous,
                allow_incomplete=allow_incomplete,
                report_invalid=report_invalid,
            )

        preview = _task_preparse_definition().read_from_string(
            text, allow_dangerous=False
        )
        selected = preview.TASK.TASK()
        selectors = tuple(selected) if selected is not None else ()
        if not selectors:
            if allow_incomplete:
                raise ValueError("Cannot autodetect incomplete input without TASK selector; supply an explicit definition or task")
            selectors = ("SCF",)
        return _read_selected_definitions(
            text, selectors, filename=filename,
            allow_dangerous=allow_dangerous,
            allow_incomplete=allow_incomplete,
            report_invalid=report_invalid,
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

    def copy(self, *args, new_schema=False, **kwargs):
        copied = super().copy(**kwargs)
        if not new_schema:
            copied._derived_from_name = self.derived_from_name
        copied.is_user_defined = self.is_user_defined
        if '_task_definition' in self.__dict__:
            copied._task_definition = self._task_definition
        return copied

    @property
    def derived_from_name(self):
        return self._derived_from_name or self.name

    def is_derived_definition(self):
        return self._derived_from_name and self._derived_from_name != self.name

    @property
    def display_name(self):
        name = self.name.upper()
        if self._derived_from_name:
            derived_from = self._derived_from_name.upper()
        else:
            derived_from = name

        return (
            derived_from if derived_from == name
            else f"{name}({derived_from})"
        )

    def info(self, generic=True):
        description = super().info(generic)
        return (
            f"{self.display_name} - {description}"
            if description else self.display_name
        )



class _CachedPreparseIPD(InputParametersDefinition):
    """Preparser only: no mutations after grammar compilation."""

    @cache
    def grammar(self, allow_dangerous=False, *, allow_incomplete=False):
        return super().grammar(allow_dangerous, allow_incomplete=allow_incomplete)

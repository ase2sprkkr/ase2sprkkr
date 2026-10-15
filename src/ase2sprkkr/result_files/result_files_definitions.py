from ..sprkkr.configuration import (
    ConfigurationValueDefinition,
    ConfigurationFileDefinition,
    ConfigurationSectionDefinition,
    ConfigurationSection,
)

from ..common.configuration_definitions import SeparatorDefinition
from ..common.decorators import cached_class_property, cache
from ..common.grammar_types import (
    unsigned,
    Array,
    Table,
    Keyword,
    GrammarType,
    pot_mixed,
    line_string,
    line_end,
    Separator as SeparatorType,
    BlankSeparator as BlankSeparatorType,
)
import pyparsing as pp
from .result_files import ResultFile
from ..common.grammar import generate_grammar


class ResultFileValueDefinition(ConfigurationValueDefinition):
    """This class describes the format of one value of
    a header of a result file"""

    @cached_class_property
    def delimiter():
        with generate_grammar():
            return pp.Empty().set_name("\t")

    prefix = " "

    type_from_type_map = {str: line_string}
    type_of_dangerous = pot_mixed


class ResultFileSection(ConfigurationSection):
    pass


class ResultFileSectionDefinition(ConfigurationSectionDefinition):
    @cached_class_property
    def delimiter():
        with generate_grammar():
            return pp.Suppress("\n").set_whitespace_chars(" \r\t").set_name("\n")
    custom_class = None
    force_order = True


class BlankSeparator(SeparatorDefinition):
    """
    A special (hidden) value, that appears in a result file header

    The separator is a blank line
    """

    separator_type = BlankSeparatorType()


class Separator(SeparatorDefinition):
    """
    A special (hidden) value, that appears in a result file header
    """

    separator_type = SeparatorType(char="#", length=80)


class ResultFileDefinition(ConfigurationFileDefinition):
    """This class describes the format of one
    structure of a structured result file"""

    force_order = True
    """ The order of items in a result file is fixed """

    value_name_format = "{:<12}"

    child_class = ResultFileValueDefinition
    """ standard child class """

    custom_class = None
    """ No custom members in the result files """

    delimiter = line_end
    """ Options are delimited by newlines in result files. """

    result_class = ResultFile

    configuration_type_name = "RESULT FILE"

    def add_ignored(self, out):
        """In the result files, do not ignore doublecrossed lines"""
        return out

    @property
    def extension(self):
        mod = self.__class__.__module__
        return mod.rsplit(".", 1)[1]


@cache
def result_file_header():
    """Return the members of the common result file header, up to
    the first ``KEYWORD`` value (which differs for each result file)"""
    V = ResultFileValueDefinition
    return [
        V("TITLE", str, info="Name of the calculation"),
        V("SYSTEM", str, info="Description of the calculated system"),
        V("NQ_EFF", int, info="Number of atomic sites"),
        V("NT_EFF", int, info="Number of types of atoms"),
        BlankSeparator(),
        V("NE", int),
        V("IREL", int),
        V("EFERMI", float, info="Fermi energy"),
        V("INFO", str),
        BlankSeparator(),
        V("ORBITALS", Table({"IQ": int, "NLQ": unsigned}, named_result=True), name_in_grammar=False),
        V(
            "TYPES",
            Table({"IT": int, "TXT_T": str, "CONC": float, "NAT": int, "IQAT": Array(int)}),
            name_in_grammar=False,
        ),
    ]


def create_result_file_definition(keyword, add, name=None, cls=ResultFileDefinition, **kwargs):
    """Create a definition of a result file. Supply the additional values
    containing the actual data.
    """
    V = ResultFileValueDefinition
    if isinstance(keyword, str):
        if name is None:
            name = keyword
        keyword = V("KEYWORD", Keyword(keyword))
    elif isinstance(keyword, GrammarType):
        keyword = V("KEYWORD", keyword)
    fields = [keyword]
    fields.extend(result_file_header())
    if add:
        fields.extend(add)
    return cls(f"{name}", fields, **kwargs)

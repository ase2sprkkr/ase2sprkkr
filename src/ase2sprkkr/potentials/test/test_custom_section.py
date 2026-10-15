import re

from ...common.grammar import delimitedList, generate_grammar
from ..custom_potential_section import SectionString
from ..potential_definitions import PotentialDefinition


def test_section_delimiter_pattern():
    pattern = re.compile(SectionString.delimiter_pattern)

    for text in (
        "\n****************\n",
        "    \n****************\n",
        " \n \n   \n****************\n",
        "\n****************\n   \n    \n",
    ):
        assert pattern.fullmatch(text) is not None

    for text in (
        "a\n****************\n",
        "\n****************\na",
        "\n******* *********\na",
    ):
        assert pattern.fullmatch(text) is None


def test_custom_member_grammar():
    with generate_grammar():
        grammar = PotentialDefinition.custom_member_grammar()

    section = """HOST MADELUNG POTENTIAL
        IQ      VLMMAD
NLMTOP-POT         4
         1   1 -2.34816691089662E-01
         1   2  0.00000000000000E+00to je konec"""

    grammar.parse_string(section, True)
    assert grammar.parse_string(section + "\n")[0][1].endswith("to je konec")

    with generate_grammar():
        sections = delimitedList(grammar, SectionString.delimiter())

    parsed = sections.parse_string(section + "\n************************\n" + section, True)
    assert len(parsed) == 2

    parsed = sections.parse_string(
        section
        + "\n************************\n"
        + section
        + "     \n****************************       \n   \n"
        + section
        + section,
        True,
    )
    assert len(parsed) == 3

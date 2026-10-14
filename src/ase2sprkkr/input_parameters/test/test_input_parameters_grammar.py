import io
from functools import partial

import numpy as np
import pyparsing as pp
import pytest

from ...common import grammar_types as gt
from ...common.configuration_containers import CustomSection, Section
from ...common.grammar import generate_grammar
from ...common.options import CustomOption, Option
from ...common.warnings import DataValidityError
from .. import input_parameters
from .. import input_parameters_definitions as cd
from ._pytest_helpers import as_array, assert_equal, assert_not_valid, assert_parse

V = cd.InputValueDefinition


@pytest.mark.slow
def test_input_parameters_definition():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("GRID", gt.SetOf(int, length=1), fixed_value=3),
                V("NE", gt.SetOf(int, min_length=1)),
                V("Ime", float, 0.0),
            ],
            "SITES": [V("NL", int)],
        }
    )
    grammar = definition.sections["SITES"].values["NL"].grammar()
    check_parse = partial(assert_parse, grammar=lambda: grammar)
    check_invalid = partial(assert_not_valid, grammar=lambda: grammar)

    check_parse("NL=3", ("NL", 3))
    grammar = definition.sections["ENERGY"].values["NE"].grammar()
    check_parse("NE={3}", ("NE", as_array(3)))
    grammar = definition.sections["ENERGY"].values["Ime"].grammar()
    check_parse("Ime= 0.5", ("Ime", 0.5))
    grammar = definition.sections["ENERGY"].values["GRID"].grammar()
    check_parse("GRID={3}", ("GRID", as_array(3)))
    grammar = definition.sections["ENERGY"].grammar()
    check_parse("ENERGY Ime= 0.5", ("ENERGY", {"Ime": 0.5}))
    check_parse(
        "ENERGY Ime= 0.5 NE={5}",
        ("ENERGY", {"Ime": 0.5, "NE": as_array(5)}),
    )
    check_parse(
        """ENERGY Ime= 0.5
                                   NE={5}""",
        ("ENERGY", {"Ime": 0.5, "NE": as_array(5)}),
    )
    check_parse(
        """ENERGY Ime= 0.5

                                   NE={5}""",
        ("ENERGY", {"Ime": 0.5, "NE": as_array(5)}),
    )
    check_invalid("""ENERGY Ime= 0.5

NE={5}""")

    check_invalid("ENERGY GRID={1}")
    check_invalid(" ENERGY GRID={3}")
    check_parse("ENERGY GRID={3}", ("ENERGY", {"GRID": as_array(3)}))

    grammar = definition.grammar()
    check_parse("ENERGY GRID={3}", {"ENERGY": {"GRID": as_array(3)}})
    check_parse(
        """ENERGY GRID={3}
                     NE={300}


          """,
        {"ENERGY": {"GRID": 3, "NE": 300}},
    )
    check_parse(
        """ENERGY
                     NE={300}
                     GRID={3}


SITES NL=2""",
        {"ENERGY": {"NE": 300, "GRID": 3}, "SITES": {"NL": 2}},
    )

    with generate_grammar():
        grammar = definition["ENERGY"]._grammar_of_values()
    check_parse(
        """GRID={3}
                     NE={300}
                     """,
        {"GRID": as_array(3), "NE": as_array(300)},
    )
    check_parse(
        """GRID={3}
                   NE={300}
                   NXXX=5""",
        {"GRID": as_array(3), "NE": as_array(300), "NXXX": 5},
    )

    grammar = definition.grammar()
    check_parse(
        """ENERGY GRID={3}
                     NE={300}
                     NXXX=5


SITES NL=2""",
        {
            "ENERGY": {"GRID": as_array(3), "NE": as_array(300), "NXXX": 5},
            "SITES": {"NL": 2},
        },
    )
    check_parse(
        """ENERGY GRID={3}
                     NE={300}
                     NXXX

SITES NL=2

              """,
        {
            "ENERGY": {"GRID": as_array(3), "NE": as_array(300), "NXXX": True},
            "SITES": {"NL": 2},
        },
    )

    for indent in ("     ", "  "):
        check_parse(
            f"""ENERGY GRID={{3}}
                     NE={{300}}
                     NXXX

{indent}SITES NL=2

              """,
            {
                "ENERGY": {
                    "GRID": as_array(3),
                    "NE": as_array(300),
                    "NXXX": True,
                    "SITES": True,
                    "NL": 2,
                }
            },
        )

    check_invalid(
        """ENERGY GRID={3}
                     NE={300}
                     NXXX


SITES NL=2

SITES NL=3
              """
    )

    check_parse(
        """ENERGY GRID={3}
                     NE={300}
                     NXXX

    #numbered_arrays

SITES NL=2

XSITES NR=3
              """,
        {
            "ENERGY": {"GRID": as_array(3), "NE": as_array(300), "NXXX": True},
            "SITES": {"NL": 2},
            "XSITES": {"NR": 3},
        },
    )
    check_parse(
        """ENERGY GRID={3}
                     NE={300}
                     NXXX


SITES NL=2

XSITES NR=3 NF=1
                     NZ=5.5
              """,
        {
            "ENERGY": {"GRID": as_array(3), "NE": as_array(300), "NXXX": True},
            "SITES": {"NL": 2},
            "XSITES": {"NR": 3, "NF": 1, "NZ": 5.5},
        },
    )
    check_parse(
        """ENERGY GRID={3}
                     NXXX
                     NE={300}
                     NZZZ=4


SITES NL=2

XSITES NR=3 FLAG
                     FLOAT=3.5
              """,
        {
            "ENERGY": {
                "GRID": as_array(3),
                "NXXX": True,
                "NE": as_array(300),
                "NZZZ": 4,
            },
            "SITES": {"NL": 2},
            "XSITES": {"NR": 3, "FLAG": True, "FLOAT": 3.5},
        },
    )

    parameters = definition.read_from_file(
        io.StringIO(
            """ENERGY GRID={3}
                     NE={300,200}
                     NXXX


SITES NL=2

XSITES NR=3 FLAG
                     FLOAT=3.5
              """
        )
    )
    assert isinstance(parameters, input_parameters.InputParameters)
    assert isinstance(parameters["ENERGY"], Section)
    assert isinstance(parameters["ENERGY"]["NE"], Option)
    assert isinstance(parameters["ENERGY"]["NXXX"], CustomOption)
    assert isinstance(parameters.ENERGY.NXXX, CustomOption)

    parameters.ENERGY.NXXX.set(5)
    assert parameters.ENERGY.NXXX() == 5
    parameters.ENERGY.NXXX.remove()
    assert not hasattr(parameters.ENERGY, "NXXX")

    assert_equal(parameters["ENERGY"].NE(), np.array((300, 200)))
    copied = parameters.copy(copy_values=True)
    copied.ENERGY.NE[0] = 100
    assert_equal(parameters["ENERGY"].NE(), np.array((300, 200)))
    shared = copied.copy()
    shared.ENERGY.NE[0] = 150
    assert_equal(copied["ENERGY"].NE(), np.array((150, 200)))
    copied.ENERGY.NE[0] = 180
    assert_equal(shared["ENERGY"].NE(), np.array((180, 200)))

    assert parameters["SITES"].NL() == 2
    assert parameters.find("NL").get_path() == "SITES.NL"
    parameters.find("NL").set(3)
    assert parameters["SITES"].NL() == 3
    assert isinstance(parameters["XSITES"], CustomSection)

    output = io.StringIO()
    parameters.save_to_file(output)
    output.seek(0)
    reread = definition.read_from_file(output)
    assert str(parameters.as_dict()) == str(reread.as_dict())

    with pytest.raises(DataValidityError):
        parameters.ENERGY.NE = "sss"
    parameters.ENERGY.NE.set_dangerous("sss")
    output = io.StringIO()
    parameters.save_to_file(output)
    output.seek(0)
    with pytest.raises(pp.ParseBaseException):
        definition.read_from_file(output)
    output.seek(0)
    reread = definition.read_from_file(output, allow_dangerous=True)
    assert parameters.ENERGY.NE() == "sss"
    assert str(parameters.as_dict()) == str(reread.as_dict())

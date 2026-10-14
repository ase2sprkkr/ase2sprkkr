import io
import re
from functools import partial

import numpy as np
import pyparsing as pp
import pytest

from ...common import grammar_types as gt
from ...common.grammar import generate_grammar
from ...common.warnings import DataValidityError
from .. import input_parameters_definitions as cd
from ._pytest_helpers import as_array, assert_equal, assert_not_valid, assert_parse

V = cd.InputValueDefinition


def compact(value):
    return re.sub(r"\s+", " ", value).strip()


def test_set_values():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("GRID", gt.SetOf(int, length=1), fixed_value=3),
                V("NE", gt.SetOf(int, min_length=1)),
                V("Ime", float, 0.0),
                V("ENERGY", float, 0.0),
            ],
            "SITES": [V("NL", 1)],
        }
    )
    parameters = definition.create_object()

    parameters.set({"ENERGY.Ime": 5.0, "NE": 10, "ENERGY": 7.0})
    assert parameters.find("energy") is parameters.ENERGY.ENERGY
    assert parameters.find("energy.ime") is parameters.ENERGY.Ime
    assert parameters.get_member("ENERGY") is parameters.ENERGY
    assert parameters.get_member("energy.ime") is parameters.ENERGY.Ime
    assert parameters.get_member(
        "ENERGY",
        accept=lambda member: member._definition.accept_value(6.0),
    ) is parameters.ENERGY.ENERGY
    with pytest.raises(KeyError, match="MISSING"):
        parameters.find("MISSING")
    assert (
        next(
            parameters.get_members(
                "ENERGY",
                accept=lambda member: member._definition.accept_value(6.0),
            )
        )
        is parameters.ENERGY.ENERGY
    )

    assert parameters.ENERGY.Ime() == 5.0
    assert_equal(parameters.ENERGY.NE(), 10)
    assert parameters.ENERGY.ENERGY() == 7.0

    parameters.set({"ENERGY.ENERGY": 4.0})
    assert parameters.ENERGY.ENERGY() == 4.0
    parameters.set("ENERGY.ENERGY", 5.0)
    assert parameters.ENERGY.ENERGY() == 5.0
    parameters.set("ENERGY", 6.0)
    assert parameters.ENERGY.ENERGY() == 6.0

    parameters.set("ENERGY.ENERGYY", 7.0, unknown="ignore")
    assert "ENERGYY" not in parameters.ENERGY
    with pytest.raises(KeyError):
        parameters.set("ENERGY.ENERGYY", 7.0, unknown="fail")
    with pytest.raises(KeyError, match="dotted path"):
        parameters.set("ENERGY.ENERGYY", 7.0, unknown="add")

    parameters.ENERGY.set("ENERGYY", 7.0, unknown="add")
    assert parameters.ENERGY.ENERGYY() == 7.0


def test_repeated_value():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"ENERGY": [V("A", 1), V("B", 2, is_repeated=True), V("C", 3)]}
    )

    data = """ENERGY A=3
    B=2
    C=77"""
    assert_parse(data, {"ENERGY": {"A": 3, "B": [2], "C": 77}}, definition.grammar())
    parameters = definition.read_from_string(data)
    assert_equal({"ENERGY": {"A": 3, "B": [2], "C": 77}}, parameters.to_dict())
    assert compact(parameters.ENERGY.to_string()) == "ENERGY A=3 B=2 C=77"

    data = """ENERGY A=3
    B=2
    B=5
    B=8
    C=77"""
    grammar = definition.grammar()
    assert_parse(data, {"ENERGY": {"A": 3, "B": [2, 5, 8], "C": 77}}, grammar)
    parameters = definition.read_from_string(data)
    assert_equal(
        {"ENERGY": {"A": 3, "B": as_array([2, 5, 8]), "C": 77}},
        parameters.to_dict(),
    )
    assert compact(parameters.ENERGY.to_string()) == "ENERGY A=3 B=2 B=5 B=8 C=77"

    definition["ENERGY"].force_order = True
    assert_parse(data, {"ENERGY": {"A": 3, "B": [2, 5, 8], "C": 77}}, grammar)
    parameters = definition.read_from_string(data)
    assert_equal(
        {"ENERGY": {"A": 3, "B": as_array([2, 5, 8]), "C": 77}},
        parameters.to_dict(),
    )
    assert compact(parameters.ENERGY.to_string()) == "ENERGY A=3 B=2 B=5 B=8 C=77"

    definition["ENERGY"]["B"].is_repeated = V.Repeated.NUMBERED
    grammar = definition.grammar()
    data = """ENERGY A=3
    B1=2
    B2=5
    B3=8
    C=77"""
    pp.ParserElement.verbose_stacktrace = True
    assert_parse(data, {"ENERGY": {"A": 3, "B": [2, 5, 8], "C": 77}}, grammar)

    data = """ENERGY A=3
    B1=2
    B3=5
    B2=8
    C=77"""
    assert_parse(data, {"ENERGY": {"A": 3, "B": [2, 8, 5], "C": 77}}, grammar)

    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("A", 1),
                V("B", gt.Array(int, length=2), is_repeated=True),
                V("C", 3),
            ]
        }
    )
    grammar = definition.grammar()
    data = """ENERGY A=3
    B=2 3
    B=5 9
    B=8 16
    C=77"""
    assert_parse(
        data,
        {
            "ENERGY": {
                "A": 3,
                "B": [as_array([2, 3]), as_array([5, 9]), as_array([8, 16])],
                "C": 77,
            }
        },
        grammar,
    )
    parameters = definition.read_from_string(data)
    assert compact(parameters.ENERGY.to_string()) == compact(data)

    definition["ENERGY"]["B"].is_repeated = V.Repeated.NUMBERED
    grammar = definition.grammar()
    assert_not_valid(data, grammar)

    data = """ENERGY A=3
    B1=2 3
    B2=5 9
    B3=8 16
    C=77"""
    assert_parse(
        data,
        {
            "ENERGY": {
                "A": 3,
                "B": [as_array([2, 3]), as_array([5, 9]), as_array([8, 16])],
                "C": 77,
            }
        },
        grammar,
    )
    parameters = definition.read_from_string(data)
    assert compact(parameters.ENERGY.to_string()) == compact(data)


def test_numbered_if():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("MODE", 0),
                V(
                    "KA",
                    int,
                    is_repeated=V.Repeated.NUMBERED_IF(
                        lambda option: option._container["MODE"]() == 1
                    ),
                ),
            ]
        }
    )

    parameters = definition.read_from_string("ENERGY MODE=0 KA1=4")
    assert np.array_equal(parameters.ENERGY.KA(), [4])
    assert compact(parameters.to_string()) == "ENERGY MODE=0 KA=4"

    parameters = definition.read_from_string("ENERGY KA=4 MODE=1")
    assert np.array_equal(parameters.ENERGY.KA(), [4])
    assert compact(parameters.to_string()) == "ENERGY MODE=1 KA1=4"

    parameters = definition.read_from_string("ENERGY MODE=1 KA1=4 KA2=5")
    assert np.array_equal(parameters.ENERGY.KA(), [4, 5])
    assert compact(parameters.to_string()) == "ENERGY MODE=1 KA1=4 KA2=5"

    with pytest.warns(DataValidityError, match="can contain only one value"):
        unnumbered = definition.read_from_string("ENERGY KA1=4 KA2=5 MODE=0")
    assert np.array_equal(unnumbered.ENERGY.KA(), [4, 5])
    with pytest.raises(DataValidityError, match="can contain only one value"):
        unnumbered.to_string(validate=False)

    with pytest.raises(DataValidityError, match="can contain only one value"):
        parameters.ENERGY.MODE = 0
    assert parameters.ENERGY.MODE() == 1
    assert compact(parameters.to_string()) == "ENERGY MODE=1 KA1=4 KA2=5"

    parameters = definition.create_object()
    with pytest.raises(DataValidityError, match="can contain only one value"):
        parameters.ENERGY.KA = [7, 8]
    assert parameters.ENERGY.KA() is None

    with pytest.raises(pp.ParseBaseException):
        definition.read_from_string("ENERGY MODE=0 KA=4 KA1=5")

    copied = definition.copy()
    parameters = copied.read_from_string("ENERGY MODE=1 KA=3")
    assert compact(parameters.to_string()) == "ENERGY MODE=1 KA1=3"

    vector_definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("MODE", 0),
                V(
                    "KA",
                    gt.SetOf(int, length=3),
                    is_repeated=V.Repeated.NUMBERED_IF(
                        lambda option: option._container["MODE"]() == 1
                    ),
                ),
            ]
        }
    )
    vector = vector_definition.read_from_string("ENERGY MODE=0 KA={1,2,3}")
    assert vector.ENERGY.KA().shape == (1, 3)
    assert compact(vector.to_string()) == "ENERGY MODE=0 KA={1,2,3}"
    vector.ENERGY.MODE = 1
    assert compact(vector.to_string()) == "ENERGY MODE=1 KA1={1,2,3}"


def test_sparse_numbered():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("GRID", gt.SetOf(int, length=1), fixed_value=3),
                V("NE", gt.SetOf(int, min_length=1)),
                V("Ime", float, 0.0),
            ],
            "SITES": [V("NL", 1)],
        }
    )
    with generate_grammar():
        grammar = definition._grammar_of_values()

    check_parse = partial(assert_parse, grammar=lambda: grammar)
    check_invalid = partial(assert_not_valid, grammar=lambda: grammar)

    check_invalid("ENERGY NE=1 Ime=0.5 Ime=2.0")
    check_parse("ENERGY NE=1 Ime=0.5", {"ENERGY": {"NE": as_array(1), "Ime": 0.5}})

    definition.sections["ENERGY"]["Ime"].is_repeated = V.Repeated.DICT
    with generate_grammar():
        grammar = definition._grammar_of_values()
    check_invalid("ENERGY NE=1 Ime=0.5")
    check_invalid("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8")
    check_parse(
        "ENERGY NE=1 Ime1=0.4 Ime5=0.8",
        {"ENERGY": {"NE": as_array(1), "Ime": {1: 0.4, 5: 0.8}}},
    )

    definition.sections["ENERGY"]["Ime"].is_repeated = V.Repeated.DEFAULTDICT
    with generate_grammar():
        grammar = definition._grammar_of_values()
    check_parse(
        "ENERGY NE=1 Ime1=0.4 Ime5=0.8",
        {"ENERGY": {"NE": as_array(1), "Ime": {1: 0.4, 5: 0.8}}},
    )
    check_invalid("ENERGY NE=1 Ime=0.5 Ime=2.0")
    check_parse(
        "ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8",
        {
            "ENERGY": {
                "NE": as_array(1),
                "Ime": {"def": 0.5, 1: 0.4, 5: 0.8},
            }
        },
    )

    parameters = definition.read_from_file(
        io.StringIO("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8")
    )
    assert parameters.ENERGY.Ime() == 0.5
    assert parameters.ENERGY.Ime[5] == 0.8
    assert_equal(parameters.ENERGY.Ime[:], [0.4, 0.5, 0.5, 0.5, 0.8])
    assert parameters.ENERGY.Ime.as_dict() == {"def": 0.5, 1: 0.4, 5: 0.8}
    assert parameters.ENERGY.Ime(all_values=True) == {"def": 0.5, 1: 0.4, 5: 0.8}

    parameters.ENERGY.Ime[1] = 0.7
    parameters.ENERGY.Ime[9] = 0.9
    assert compact(parameters.ENERGY.to_string()) == (
        "ENERGY GRID={3} NE={1} Ime=0.5 Ime1=0.7 Ime5=0.8 Ime9=0.9"
    )

    parameters.ENERGY.Ime[[1, 4, 9]] = 0.2
    assert_equal(parameters.ENERGY.Ime[1, 5, 9], [0.2, 0.8, 0.2])
    parameters.ENERGY.Ime[[1, 4, 9]] = 0.3
    assert_equal(parameters.ENERGY.Ime[[1, 5, 9]], [0.3, 0.8, 0.3])
    parameters.ENERGY.Ime[5:9] = 0.2
    assert_equal(parameters.ENERGY.Ime[5:10], [0.2, 0.2, 0.2, 0.2, 0.3])

    with pytest.raises(KeyError):
        _ = parameters.ENERGY.Ime[5.0]
    with pytest.raises(KeyError):
        parameters.ENERGY.Ime[5.0] = 1
    with pytest.raises(KeyError):
        _ = parameters.ENERGY.Ime["5"]
    with pytest.raises(KeyError):
        parameters.ENERGY.Ime["5"] = 1

    assert_equal(parameters.ENERGY.Ime[[1, 5, 9]], [0.3, 0.2, 0.3])
    parameters.ENERGY.Ime = {3: 5.0, "def": 3.0}
    assert_equal(parameters.ENERGY.Ime[:], [3.0, 3.0, 5.0])
    parameters.ENERGY.Ime[2:5] = [2.0, 3.0, 8.0]
    assert_equal(parameters.ENERGY.Ime[:], [3.0, 2.0, 3.0, 8.0])
    parameters.ENERGY.Ime[1:3] = 7.0
    assert_equal(parameters.ENERGY.Ime[:], [7.0, 7.0, 3.0, 8.0])

    parameters = definition.read_from_file(
        io.StringIO("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8")
    )
    with pytest.raises(DataValidityError):
        parameters.ENERGY.Ime = "ss"
    parameters.ENERGY.Ime.set_dangerous("uu")
    assert parameters.ENERGY.Ime() == "uu"
    out = parameters.ENERGY.to_string()
    assert compact(out) == "ENERGY GRID={3} NE={1} Ime=uu Ime1=0.4 Ime5=0.8"
    with pytest.raises(pp.ParseBaseException):
        definition.read_from_file(io.StringIO(out))
    parameters = definition.read_from_file(io.StringIO(out), allow_dangerous=True)
    assert parameters.ENERGY.Ime(all_values=True) == {"def": "uu", 1: 0.4, 5: 0.8}

    parameters.ENERGY.Ime = 1.0
    with pytest.raises(DataValidityError):
        parameters.ENERGY.Ime[5] = "ss"
    parameters.ENERGY.Ime.set_dangerous("yy", index=5)
    out = parameters.ENERGY.to_string()
    assert compact(out) == "ENERGY GRID={3} NE={1} Ime=1.0 Ime1=0.4 Ime5=yy"
    with pytest.raises(pp.ParseBaseException):
        definition.read_from_file(io.StringIO(out))
    parameters = definition.read_from_file(io.StringIO(out), allow_dangerous=True)
    assert parameters.ENERGY.Ime(all_values=True) == {"def": 1.0, 1: 0.4, 5: "yy"}

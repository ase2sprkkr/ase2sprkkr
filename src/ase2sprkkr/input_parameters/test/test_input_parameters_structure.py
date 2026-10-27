import re
import warnings

import numpy as np
import pyparsing as pp
import pytest

from ...common import grammar_types as gt
from ...common.configuration_definitions import VirtualDefinition, gather, switch
from ...common.configuration_transaction import ConfigurationTransaction
from ...common.generated_configuration_definitions import GeneratedValueDefinition, Length
from ...common.options import DangerousValue
from ...common.warnings import DataValidityError
from .. import input_parameters_definitions as cd
from ._pytest_helpers import as_array, assert_equal, assert_not_valid, assert_parse

V = cd.InputValueDefinition


def compact(value):
    return re.sub(r"\s+", " ", value).strip()


def length_is_optional(option):
    """Require an input count whenever any of its source fields has a value."""
    return all(
        option._container[name]() is None
        for name in option._definition._length_of
    )


def test_gather():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("GRID", gt.SetOf(int, length=1), fixed_value=3),
                *gather(V("A", 1), V("B", 2)),
                V("C", 3),
            ]
        }
    )
    assert_parse(
        "ENERGY GRID={3} A B=1 2 C=3",
        {"ENERGY": {"GRID": as_array(3), "A": 1, "B": 2, "C": 3}},
        definition.grammar(),
    )
    parameters = definition.read_from_string("ENERGY GRID={3} A B=1 2 C=3")
    assert re.sub(r"[\s\t\n]+", " ", parameters.to_string()) == (
        "ENERGY GRID={3} A B=1 2 C=3 "
    )


def test_numpy_array():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"ENERGY": [V("C", gt.NumpyArray(lines=3), name_in_grammar=False)]}
    )
    assert_parse(
        """ENERGY
    1 1
    2 2
    3 3
    """,
        {"ENERGY": {"C": as_array([1, 1, 2, 2, 3, 3]).reshape(3, 2)}},
        definition.grammar(),
    )

    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("A", 1),
                V("B", 2),
                V("C", gt.NumpyArray(lines=3), name_in_grammar=False),
                V("D", gt.NumpyArray(lines=2), name_in_grammar=False),
                V("E", 3),
            ]
        }
    )
    data = """ENERGY A=3 B=2
    1 1
    2 2
    3 3
    4 5 6
    5 9 8
    E=77"""
    expected = {
        "ENERGY": {
            "A": 3,
            "B": 2,
            "C": as_array([1.0, 1, 2, 2, 3, 3]).reshape(3, 2),
            "D": as_array([4.0, 5, 6, 5, 9, 8]).reshape((2, 3)),
            "E": 77,
        }
    }
    assert_parse(data, expected, definition.grammar())
    definition["ENERGY"].force_order = True
    assert_parse(data, expected, definition.grammar())

    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("A", 1),
                V("B", 2),
                V("C", gt.NumpyArray(lines="A"), name_in_grammar=False),
                V(
                    "D",
                    gt.NumpyArray(lines=lambda B: int(B)),
                    name_in_grammar=False,
                ),
                V("E", 3),
            ]
        }
    )
    assert_parse(data, expected, definition.grammar())

    def check(array_type, value, expected_text):
        assert array_type.string(value) == expected_text
        with generate_grammar_context():
            assert_equal(np.asarray(value), array_type.parse(expected_text))

    check(
        gt.NumpyArray(item_format="%2.0f"),
        [[1, 2, 3], [4, 5, 6]],
        " 1  2  3\n 4  5  6",
    )
    check(
        gt.NumpyArray(
            shape=(5,),
            items_per_line=3,
            item_format="%2d",
            written_delimiter="",
            dtype=int,
        ),
        [1, 2, 3, 4, 5],
        " 1 2 3\n 4 5",
    )
    check(
        gt.NumpyArray(indented=2, item_format="%2.0f"),
        [[1, 2, 3], [4, 5, 6]],
        "   1  2  3\n   4  5  6",
    )
    check(
        gt.NumpyArray(indented=(8, 2), item_format="%2.0f"),
        [[1, 2, 3, 4, 5, 6], [4, 5, 6, 7, 8, 9]],
        " 1  2  3\n    4  5\n    6\n 4  5  6\n    7  8\n    9",
    )
    check(
        gt.NumpyArray(
            indented=(8, 2),
            no_newline_at_end=False,
            item_format="%2.0f",
        ),
        [[1, 2, 3, 4, 5, 6], [4, 5, 6, 7, 8, 9]],
        " 1  2  3\n    4  5\n    6\n 4  5  6\n    7  8\n    9\n",
    )
    check(
        gt.NumpyArray(line_length=8, item_format="%2.0f", shape=(2, -1)),
        [[1, 2, 3, 4, 5, 6], [4, 5, 6, 7, 8, 9]],
        " 1  2  3\n 4  5  6\n 4  5  6\n 7  8  9",
    )


def generate_grammar_context():
    from ...common.grammar import generate_grammar

    return generate_grammar()


def test_copy():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"ENERGY": [V("A", 1), V("B", 2, is_repeated=True), V("C", 3)]}
    )
    copied = definition.copy()
    copied["ENERGY"]["A"].default_value = 5
    assert definition["ENERGY"]["A"].default_value == 1


@pytest.mark.parametrize("remove_replaced", [False, True])
def test_copy_skips_discarded_members(monkeypatch, remove_replaced):
    section = cd.InputSectionDefinition(
        "CONTROL", [V("KEPT", 1), V("REPLACED", 2), V("LAST", 3), V("REMOVED", 4)]
    )

    def discarded_copy(**kwargs):
        raise AssertionError("Discarded definitions must not be copied")

    for name in ("REPLACED", "REMOVED"):
        monkeypatch.setattr(section[name], "copy", discarded_copy)
    replacement = V("REPLACED", 5)
    removed = ["REMOVED", "REPLACED"] if remove_replaced else ["REMOVED"]
    copied = section.copy(add=[replacement, V("NEW", 6)], remove=removed)

    expected = ["KEPT", "LAST", "REPLACED", "NEW"] if remove_replaced else ["KEPT", "REPLACED", "LAST", "NEW"]
    assert list(copied.names()) == expected
    assert copied["KEPT"] is not section["KEPT"]
    assert copied["REPLACED"] is replacement
    copied["KEPT"].default_value = 7
    assert section["KEPT"].default_value == 1
    assert list(section.names()) == ["KEPT", "REPLACED", "LAST", "REMOVED"]
    with pytest.raises(KeyError):
        section.copy(remove=["MISSING"])


@pytest.mark.parametrize("task", ["aipes", "band", "clpes", "spleed"])
def test_customized_spectroscopy_definitions_are_independent(task):
    from importlib import import_module

    factory = import_module(f"ase2sprkkr.input_parameters.definitions.{task}").input_parameters_definition
    first = factory()
    second = factory()
    for section in first.members():
        assert section is not second[section.name]
        assert section.container is first
        for value in section.members():
            assert value is not second[section.name][value.name]
            assert value.container is section
    first["ENERGY"]["NE"].default_value = [10]
    assert not np.array_equal(first["ENERGY"]["NE"].default_value, second["ENERGY"]["NE"].default_value)


def test_length():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("A", 1),
                V(
                    "NK", Length("K1", "K2"),
                    is_optional=length_is_optional, is_required=False,
                ),
                V("K1", gt.Array(int), is_optional=True),
                V("K2", gt.Array(int), is_optional=True),
            ]
        }
    )
    parameters = definition.create_object()
    parameters.ENERGY.set(K1=[1, 2, 3], K2=[2, 2, 3])
    assert compact(parameters.to_string()) == "ENERGY A=1 NK=3 K1=1 2 3 K2=2 2 3"
    parsed = definition.read_from_string(parameters.to_string())
    assert parsed.ENERGY.NK() == 3

    for sources in ("K1=1 2 3", "K2=2 2 3", "K1=1 2 3 K2=2 2 3"):
        with pytest.raises(pp.ParseBaseException, match="Missing non-optional input item ENERGY.NK"):
            definition.read_from_string(f"ENERGY A=1 {sources}")

    with warnings.catch_warnings():
        warnings.simplefilter("error", DataValidityError)
        absent = definition.read_from_string("ENERGY A=1")
        assert "NK=" not in absent.to_string()
    assert absent.ENERGY.NK() is None
    assert not hasattr(absent.ENERGY, "_parsed_values")

    with pytest.raises(DataValidityError):
        parsed.ENERGY.NK = 7
    with pytest.raises(DataValidityError):
        parameters.ENERGY.K2 = [2, 2]
    assert np.array_equal(parameters.ENERGY.K2(), [2, 2, 3])

    with pytest.warns(DataValidityError):
        invalid = definition.read_from_string("ENERGY A=1 NK=3 K1=1 2 3 K2=2 2")
    with pytest.raises(DataValidityError):
        invalid.to_string(validate="save")

    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("A", 1),
                V("NK", Length("K1", "K2", default_values=(1, 2))),
                V("K1", gt.Array(int), is_optional=True),
                V("K2", gt.Array(int), is_optional=True),
            ]
        }
    )
    parameters = definition.create_object()
    parameters.ENERGY.NK = 3
    assert np.array_equal(parameters.ENERGY.K1(), [1, 1, 1])
    assert np.array_equal(parameters.ENERGY.K2(), [2, 2, 2])
    parameters.ENERGY.NK = 5
    assert np.array_equal(parameters.ENERGY.K1(), [1, 1, 1, 1, 1])
    parameters.ENERGY.NK = 2
    assert np.array_equal(parameters.ENERGY.K1(), [1, 1])

    section = parameters.ENERGY
    old_k1 = section.K1().copy()
    old_k2 = section.K2().copy()
    with pytest.raises(DataValidityError):
        section.set(K1=[1, 1, 1], K2=[2, 2])
    assert np.array_equal(section.K1(), old_k1)
    assert np.array_equal(section.K2(), old_k2)

    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("A", 1),
                V("NK", Length("K1", default_values=[1, 2, 3])),
                V(
                    "K1",
                    gt.Array(int, length=3),
                    is_repeated="REPEATED",
                    is_optional=True,
                ),
            ]
        }
    )
    parameters = definition.create_object()
    parameters.ENERGY.NK = 2
    assert np.array_equal(parameters.ENERGY.K1(), [[1, 2, 3], [1, 2, 3]])
    parameters.ENERGY.K1[1] = [3, 3, 3]
    parameters.ENERGY.NK = 3
    assert np.array_equal(
        parameters.ENERGY.K1(),
        [[1, 2, 3], [3, 3, 3], [1, 2, 3]],
    )


def test_standard_container_access_during_transactions():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [V("VALUE", 1)]}
    )
    section = definition.create_object().CONTROL
    section_definition = definition["CONTROL"]
    assert not section_definition.is_option
    assert section_definition["VALUE"].is_option

    assert "VALUE" in section
    assert section["VALUE"]() == 1
    assert "EXTRA" not in section
    section["VALUE"].set_dangerous("not-an-integer")
    assert section["VALUE"].is_dangerous()

    dangerous = DangerousValue("unchecked")
    with ConfigurationTransaction(section) as transaction:
        with transaction.savepoint() as savepoint:
            section.stage_value(transaction, "EXTRA", dangerous, unknown="add")
            assert "VALUE" in section
            assert "EXTRA" in section
            assert "MISSING" not in section
            assert section["EXTRA"].is_dangerous()
            savepoint.rollback()
    assert "EXTRA" not in section

    generated_definition = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [V("VALUES", gt.Array(int)), V("NVALUES", Length("VALUES"))]}
    )
    generated_section = generated_definition.create_object().CONTROL
    with ConfigurationTransaction(generated_section) as transaction:
        with transaction.savepoint() as savepoint:
            generated_section["VALUES"].stage(transaction, [1, 2, 3])
            assert generated_section.NVALUES() == 3
            savepoint.rollback()

    copied = generated_definition.copy().create_object()
    copied.CONTROL.VALUES = [1, 2]
    assert copied.CONTROL.NVALUES() == 2


def test_definition_storage_metadata():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "CONTROL": [
                V("VALUES", gt.Array(int), is_optional=True),
                V(
                    "NVALUES", Length("VALUES"),
                    is_optional=length_is_optional, is_required=False,
                ),
                V("NOT_STORED", int, is_stored=False, is_optional=False, is_required=False),
                GeneratedValueDefinition("GENERATED", lambda section: 1),
            ]
        }
    )

    assert not VirtualDefinition("VIRTUAL").is_stored
    assert definition.is_stored
    assert definition["CONTROL"].is_stored
    assert definition["CONTROL"]["VALUES"].is_stored
    assert definition["CONTROL"]["NVALUES"].is_stored
    assert not definition["CONTROL"]["NOT_STORED"].is_stored
    assert not definition["CONTROL"]["GENERATED"].is_stored

    # Length is generated but explicitly stored; other generated/non-stored
    # options must not be required by the input grammar.
    with warnings.catch_warnings():
        warnings.simplefilter("error", DataValidityError)
        parsed = definition.read_from_string("CONTROL VALUES=1 2 3 NVALUES=3")
    assert parsed.CONTROL.NVALUES() == 3

    with pytest.raises(pp.ParseBaseException, match="Missing non-optional input item"):
        definition.read_from_string("CONTROL VALUES=1 2 3")


def test_generated_option_uses_its_runtime_section():
    def generated_getter(section, _key=None):
        return section["SOURCE"]()

    def generated_setter(section, value, _key=None):
        transaction = ConfigurationTransaction.current(section)
        section["SOURCE"].stage(transaction, value)

    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "CONTROL": [
                V("SOURCE", 1),
                GeneratedValueDefinition("GENERATED", generated_getter, generated_setter),
            ]
        }
    )
    parameters = definition.create_object()
    assert parameters.CONTROL.GENERATED() == 1
    parameters.CONTROL.GENERATED.set(2)
    assert parameters.CONTROL.SOURCE() == 2


def test_parse_validation_sees_keyword_assignments():
    def parsed_keywords(section, _values, why):
        if why == "parse":
            assert "MODE" in section._parsed_values

    definition = cd.InputParametersDefinition(
        [
            cd.InputSectionDefinition(
                "CONTROL",
                [V("MODE", 1)],
                validators=parsed_keywords,
            )
        ]
    )
    parameters = definition.create_object()
    parameters.CONTROL.set(validation_reason="parse", MODE=2)
    assert parameters.CONTROL.MODE() == 2
    assert not hasattr(parameters.CONTROL, "_parsed_values")


def test_adding_empty_custom_member_is_validated():
    def reject_extra(root, _values, _why):
        if "EXTRA" in root:
            return DataValidityError("EXTRA is not allowed")

    definition = cd.InputParametersDefinition([], validators=reject_extra)
    parameters = definition.create_object()

    with pytest.raises(DataValidityError, match="EXTRA is not allowed"):
        parameters.add("EXTRA")
    assert "EXTRA" not in parameters


def test_conditional_length():
    def mode_is_one(_definition, section):
        return section["MODE"]() == 1

    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "CONTROL": [
                V("VALUES", gt.Array(int), is_optional=True),
                V(
                    "NVALUES", Length("VALUES"), condition=mode_is_one,
                    is_optional=length_is_optional, is_required=False,
                ),
                V("MODE", 1),
            ]
        }
    )

    active = definition.read_from_string("CONTROL VALUES=1 2 3 NVALUES=3 MODE=1")
    assert active.CONTROL.NVALUES() == 3
    assert not hasattr(active.CONTROL, "_parsed_values")
    assert "NVALUES=3" in active.to_string()

    with pytest.warns(DataValidityError):
        definition.read_from_string("CONTROL VALUES=1 2 3 NVALUES=2 MODE=1")
    with pytest.raises(pp.ParseBaseException):
        definition.read_from_string("CONTROL VALUES=1 2 3 NVALUES=3 MODE=2")
    with pytest.raises(pp.ParseBaseException):
        definition.read_from_string("CONTROL VALUES=1 2 3 NVALUES=2 MODE=2")
    with pytest.raises(pp.ParseBaseException, match="Missing non-optional input item"):
        definition.read_from_string("CONTROL VALUES=1 2 3 MODE=1")

    for mode in (1, 2):
        with warnings.catch_warnings():
            warnings.simplefilter("error", DataValidityError)
            absent = definition.read_from_string(f"CONTROL MODE={mode}")
            assert "NVALUES=" not in absent.to_string()
        assert absent.CONTROL["NVALUES"]() is None
        assert not hasattr(absent.CONTROL, "_parsed_values")

    with warnings.catch_warnings():
        warnings.simplefilter("error", DataValidityError)
        inactive = definition.read_from_string("CONTROL VALUES=1 2 3 MODE=2")
    assert inactive.CONTROL["NVALUES"]() == 3
    assert "NVALUES" not in inactive.to_string()


def test_cross_section_conditional_length():
    def mode_is_one(_definition, section):
        return section._get_root_container()["SECTION_B"]["MODE"]() == 1

    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "SECTION_A": [
                V("OTHER", 1),
                V("VALUES", gt.Array(int), is_optional=True),
                V(
                    "NVALUES", Length("VALUES"), condition=mode_is_one,
                    is_optional=length_is_optional, is_required=False,
                ),
            ],
            "SECTION_B": [V("MODE", 1)],
        }
    )

    active = definition.read_from_string(
        "SECTION_A VALUES=1 2 3 NVALUES=3\nSECTION_B MODE=1"
    )
    assert active.SECTION_A.NVALUES() == 3
    assert "NVALUES=3" in active.to_string()

    with pytest.warns(DataValidityError):
        definition.read_from_string(
            "SECTION_A VALUES=1 2 3 NVALUES=2\nSECTION_B MODE=1"
        )
    with pytest.raises(pp.ParseBaseException):
        definition.read_from_string(
            "SECTION_A VALUES=1 2 3 NVALUES=3\nSECTION_B MODE=2"
        )
    with pytest.raises(pp.ParseBaseException):
        definition.read_from_string(
            "SECTION_A VALUES=1 2 3 NVALUES=2\nSECTION_B MODE=2"
        )
    with pytest.raises(pp.ParseBaseException, match="Missing non-optional input item"):
        definition.read_from_string("SECTION_A VALUES=1 2 3\nSECTION_B MODE=1")

    for mode in (1, 2):
        with warnings.catch_warnings():
            warnings.simplefilter("error", DataValidityError)
            absent = definition.read_from_string(f"SECTION_A OTHER=1\nSECTION_B MODE={mode}")
            assert "NVALUES=" not in absent.to_string()
        assert absent.SECTION_A["NVALUES"]() is None
        assert not hasattr(absent.SECTION_A, "_parsed_values")

    with warnings.catch_warnings():
        warnings.simplefilter("error", DataValidityError)
        inactive = definition.read_from_string(
            "SECTION_A VALUES=1 2 3\nSECTION_B MODE=2"
        )
    assert inactive.SECTION_A["NVALUES"]() == 3
    assert "NVALUES" not in inactive.to_string()


def test_switch():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("A", 1),
                *switch(
                    "A",
                    {
                        1: [V("B", 2, is_optional=True)],
                        2: V("C", 3, is_optional=True),
                    },
                ),
            ]
        }
    )
    definition.custom_class = None
    definition["ENERGY"].custom_class = None
    definition["ENERGY"].force_order = True
    grammar = definition.grammar()

    assert_parse("ENERGY A=2 C=4", {"ENERGY": {"A": 2, "C": 4}}, grammar)
    out = definition.read_from_string("ENERGY A=2 C=4")
    assert out["ENERGY"].to_string() == "ENERGY\n\tA=2\n\tC=4\n"
    assert_parse("ENERGY A=1 B=2", {"ENERGY": {"A": 1, "B": 2}}, grammar)
    out = definition.read_from_string("ENERGY A=1 B=4")
    assert out["ENERGY"].to_string() == "ENERGY\n\tA=1\n\tB=4\n"

    for text in (
        "ENERGY A=2 B=2",
        "ENERGY A=1 C=1",
        "ENERGY A=1 B=2 C=1",
        "ENERGY A=2 B=2 C=1",
        "ENERGY A=3 C=3",
    ):
        assert_not_valid(text, grammar)
    assert_parse("ENERGY A=2", {"ENERGY": {"A": 2}}, grammar)

    definition["ENERGY"]["C"].is_optional = False
    assert_not_valid("ENERGY A=2", definition.grammar())

    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("A", 1),
                *switch(
                    "A",
                    {
                        1: [
                            V("B", 2, is_optional=True),
                            V("D", 3, is_optional=True),
                            V("C", 1, is_optional=True),
                        ],
                        2: [
                            V("C", 3, is_optional=True),
                            V("E", 1, is_optional=True),
                            V("B", 2, is_optional=True),
                        ],
                    },
                ),
            ]
        }
    )
    definition["ENERGY"].force_order = True
    grammar = definition.grammar()

    cases = (
        ("ENERGY A=2 C=4 E=2 B=6", {"A": 2, "C": 4, "E": 2, "B": 6}),
        ("ENERGY A=1 B=5 D=3 C=5", {"A": 1, "B": 5, "D": 3, "C": 5}),
        ("ENERGY A=2 B=7", {"A": 2, "B": 7}),
        ("ENERGY A=1 B=7", {"A": 1, "B": 7}),
        ("ENERGY A=1 B=7 D=6", {"A": 1, "B": 7, "D": 6}),
    )
    for text, expected in cases:
        assert_parse(text, {"ENERGY": expected}, grammar)

    for text in (
        "ENERGY A=1 D=2 B=3",
        "ENERGY A=1 B=5 D=3 C=5 B=6",
        "ENERGY A=1 B=5 D=3 C=5 E=6",
        "ENERGY A=2 B=5 C=4 E=3",
    ):
        assert_not_valid(text, grammar)

    for name in ("B", "C", "D", "E"):
        definition["ENERGY"][name].is_optional = False
    grammar = definition.grammar()
    for text in (
        "ENERGY A=2 B=2",
        "ENERGY A=1 B=2",
        "ENERGY A=1 B=2 D=3",
        "ENERGY A=1 D=2 B=3",
        "ENERGY A=1 B=5 D=3 C=5 B=6",
        "ENERGY A=1 B=5 D=3 C=5 E=6",
        "ENERGY A=2 B=5 C=4 E=3",
    ):
        assert_not_valid(text, grammar)

    out = definition.read_from_string("ENERGY A=2 C=4 E=2 B=6")
    assert out["ENERGY"].to_dict() == {"A": 2, "C": 4, "E": 2, "B": 6}
    assert out["ENERGY"].to_string() == "ENERGY\n\tA=2\n\tC=4\n\tE=2\n\tB=6\n"

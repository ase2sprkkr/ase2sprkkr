import inspect
import io
import pickle
import warnings
from functools import partial

import numpy as np
import pyparsing as pp
import pytest

from ...common.grammar import generate_grammar
from ...common.warnings import DataValidityError, DataValidityWarning
from .. import input_parameters as input_parameters
from .. import input_parameters_definitions as cd
from ._pytest_helpers import assert_equal, assert_not_valid, assert_parse, parse
from ._test_callbacks import (
    root_validator_matching_values as _root_validator_matching_values,
    root_validator_reject_all as _root_validator_reject_all,
    root_validator_warning as _root_validator_warning,
    section_validator_warning as _section_validator_warning,
    validator_error_if_three as _validator_error_if_three,
    validator_warning_if_two as _validator_warning_if_two,
    warning_condition_if_two as _warning_condition_if_two,
)

V = cd.InputValueDefinition


@pytest.mark.parametrize("copy_values", [False, True])
def test_copy_preserves_requested_and_effective_task_views(copy_values):
    parameters = input_parameters.InputParameters.create("bsfek")

    copied = parameters.copy(copy_values=copy_values)

    assert copied.requested_task_name == "BSFEK"
    assert copied.task_name == "BSFEK"
    assert copied.task_definition.name == "BSFEK"
    assert copied.definition_name == "BSF"
    np.testing.assert_array_equal(copied.ENERGY.NE(), [200])


def test_section_delimiter_value():
    with generate_grammar():
        grammar = cd.InputParametersDefinition.delimiter()
        grammar = "a" + grammar + "b"

    invalid = ["a b", "a\n b", "a\n\n b", "a\n \n b", "a \n\n b", "a\n\n\n b"]
    for value in invalid:
        with pytest.raises(pp.ParseException):
            grammar.parse_string(value, True)

    valid = ["a\nb", "a \nb", "a\n  \nb", "a  \n \nb", "a \n\t\nb", "a\n\n\nb"]
    for value in valid:
        assert_equal(grammar.parse_string(value, True).as_list(), ["a", "b"])


def test_custom_value():
    with generate_grammar():
        grammar = cd.InputSectionDefinition.custom_member_grammar(
            lambda value: value[0] not in ["aaa"]
        )
    assert "\n" not in grammar.whiteChars

    check_parse = partial(assert_parse, grammar=lambda: grammar)
    check_invalid = partial(assert_not_valid, grammar=lambda: grammar)

    check_parse("bbb=1", ("bbb", 1))
    check_parse("bbb=1.3", ("bbb", 1.3))
    check_parse("bbb", ("bbb", True))
    check_invalid("aaa")
    check_invalid("aaa=1")

    with generate_grammar():
        grammar = grammar + "a"
    assert_equal(list(parse("bbb=1 a", grammar)), [("bbb", 1), "a"])
    assert_equal(list(parse("bbb=1a a", grammar)), [("bbb", "1a"), "a"])
    check_invalid("bbb=1\na")


def test_dangerous_value():
    with generate_grammar():
        grammar = cd.InputValueDefinition("aaa", int, is_required=True).grammar()

    check_parse = partial(assert_parse, grammar=lambda: grammar)
    check_invalid = partial(assert_not_valid, grammar=lambda: grammar)

    check_parse("aaa=1", ("aaa", 1))
    check_invalid("aaa")
    check_invalid("aaa=sss")
    check_invalid("aaa=1.0")

    with generate_grammar():
        grammar = cd.InputValueDefinition("aaa", int, is_required=True).grammar(
            allow_dangerous=True
        )

    check_parse("aaa=1", ("aaa", 1))

    def check_dangerous(text, expected):
        actual = grammar.parse_string(text, True)[0]
        assert actual[0] == expected[0]
        assert_equal(actual[1].value, expected[1])

    check_dangerous("aaa=1.0", ("aaa", 1.0))
    check_dangerous("aaa=sss", ("aaa", "sss"))
    check_dangerous("aaa={2,3}", ("aaa", np.array([2, 3])))

    definition = cd.InputValueDefinition("aaa", 2, is_required=True)
    option = definition.create_object()
    option.set_dangerous("AAA")
    assert option.to_string() == "\taaa=AAA"
    option.set_dangerous(None)
    assert option.to_string() == ""


def test_is_required():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "ENERGY": [
                V("E", 1),
                V(
                    "F",
                    int,
                    None,
                    is_required=lambda option: option._container["E"]() != 2,
                ),
            ]
        }
    )

    parameters = definition.create_object()
    with pytest.raises(DataValidityError):
        parameters.validate()
    energy = parameters.ENERGY
    energy.F = 2
    parameters.validate()
    with pytest.raises(DataValidityError):
        energy.F = None
    energy.E = 2
    parameters.validate()
    energy.F = None
    parameters.validate()

    assert definition.read_from_string("ENERGY E=2").ENERGY.to_dict() == {"E": 2}
    assert definition.read_from_string("ENERGY E=1 F=2").ENERGY.to_dict() == {
        "E": 1,
        "F": 2,
    }
    with pytest.warns(DataValidityError):
        definition.read_from_string("ENERGY E=1")


def test_write_condition():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"ENERGY": [V("E", 1)]}
    )
    parameters = definition.read_from_file(io.StringIO("ENERGY E=1"))
    assert parameters.to_string() == "ENERGY\n\tE=1\n"

    definition["ENERGY"]["E"].write_condition = lambda option: True
    assert parameters.to_string() == "ENERGY\n\tE=1\n"
    definition["ENERGY"]["E"].write_condition = lambda option: False
    assert parameters.to_string() == "ENERGY\n"
    definition["ENERGY"].write_condition = lambda section: True
    assert parameters.to_string() == "ENERGY\n"
    definition["ENERGY"].write_condition = lambda section: False
    assert parameters.to_string() == ""


def test_validators():
    definition = V(
        "VALUE",
        1,
        validators=(_validator_warning_if_two, _validator_error_if_three),
    )
    copied = definition.copy()
    assert copied.validators == (_validator_warning_if_two, _validator_error_if_three)
    pickle.dumps(copied)

    root = cd.InputParametersDefinition.definition_from_dict({"CONTROL": [copied]})

    parameters = root.create_object()
    with pytest.warns(DataValidityWarning, match="set: VALUE should not be two") as caught:
        parameters.CONTROL.VALUE = 2
    assert len(caught) == 1

    parameters = root.create_object()
    with pytest.warns(DataValidityWarning, match="set: VALUE should not be two") as caught:
        parameters.set(VALUE=2)
    assert len(caught) == 1

    with pytest.warns(DataValidityWarning, match="parse: VALUE should not be two") as caught:
        parameters = root.read_from_string("CONTROL\n\tVALUE=2\n")
    assert len(caught) == 1

    with pytest.warns(DataValidityWarning, match="save: VALUE should not be two") as caught:
        parameters.to_string(validate="save")
    assert len(caught) == 1

    parameters = root.create_object()
    with pytest.raises(DataValidityError, match="set: VALUE must not be three"):
        parameters.CONTROL.VALUE = 3
    assert parameters.CONTROL.VALUE() == 1

    with pytest.raises(DataValidityError, match="set: VALUE must not be three"):
        parameters.set(VALUE=3)
    assert parameters.CONTROL.VALUE() == 1

    with pytest.warns(DataValidityError, match="set: VALUE must not be three"):
        parameters.set(
            VALUE=3,
            retain_invalid="typed",
            report_invalid="warn",
        )
    assert parameters.CONTROL.VALUE() == 3

    parameters = root.create_object()
    with pytest.warns(DataValidityError, match="set: VALUE must not be three"):
        parameters.CONTROL.VALUE.set(
            3,
            retain_invalid="typed",
            report_invalid="warn",
        )
    assert parameters.CONTROL.VALUE() == 3

    with pytest.warns(DataValidityError, match="parse: VALUE must not be three"):
        invalid = root.read_from_string("CONTROL\n\tVALUE=3\n")
    with pytest.warns(DataValidityError, match="save: VALUE must not be three"):
        invalid.validate("save", report_invalid="warn")
    with pytest.warns(DataValidityError, match="save: VALUE must not be three"):
        found = invalid.check_for_errors()
    assert len(found) == 1
    assert isinstance(found[0].message, DataValidityError)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        found = invalid.check_for_errors(print=False)
    assert not caught
    assert len(found) == 1
    assert isinstance(found[0].message, DataValidityError)

    with pytest.raises(DataValidityError, match="parse: VALUE must not be three"):
        invalid.validate("parse", report_invalid="raise")
    with pytest.raises(ValueError, match="Unknown validation reason"):
        invalid.validate("warning")

    section = cd.InputSectionDefinition(
        "CONTROL",
        [V("VALUE", 1)],
        validators=_section_validator_warning,
    )
    root = cd.InputParametersDefinition(
        [section],
        validators=_root_validator_warning,
        is_optional=True,
        info="root info",
        description="root description",
    )
    copied_root = root.copy()
    signature = inspect.signature(cd.InputParametersDefinition)
    for argument in ("info", "description", "validators", "is_optional"):
        assert argument in signature.parameters
    assert copied_root.validators == (_root_validator_warning,)
    assert copied_root.is_optional is True
    assert copied_root._info == "root info"
    assert copied_root._description == "root description"
    pickle.dumps(copied_root)

    parameters = root.create_object()
    with pytest.warns(
        DataValidityWarning, match="set: section received its runtime object"
    ) as caught:
        parameters.CONTROL.VALUE = 4
    assert len(caught) == 1
    with pytest.warns(
        DataValidityWarning, match="set: root received its runtime object"
    ) as caught:
        parameters.CONTROL.VALUE = 5
    assert len(caught) == 1

    warning_definition = V(
        "VALUE",
        1,
        warning_condition=_warning_condition_if_two,
    )
    warning_copy = warning_definition.copy()
    assert warning_copy.warning_condition is _warning_condition_if_two
    pickle.dumps(warning_copy)
    warning_parameters = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [warning_copy]}
    ).create_object()
    with pytest.warns(DataValidityWarning, match="VALUE should not be two"):
        warning_parameters.CONTROL.VALUE = 2


def test_cross_section_validator_for_proposed_values():
    root = cd.InputParametersDefinition(
        [
            cd.InputSectionDefinition("SECTION_A", [V("VALUE", 0)]),
            cd.InputSectionDefinition("SECTION_B", [V("MODE", 0)]),
        ],
        validators=_root_validator_matching_values,
    )
    parameters = root.create_object()

    parameters.set({"SECTION_A.VALUE": 1, "SECTION_B.MODE": 1})
    assert parameters.SECTION_A.VALUE() == 1
    assert parameters.SECTION_B.MODE() == 1

    with pytest.raises(
        DataValidityError, match="values in SECTION_A and SECTION_B must match"
    ):
        parameters.set({"SECTION_A.VALUE": 2, "SECTION_B.MODE": 3})
    assert parameters.SECTION_A.VALUE() == 1
    assert parameters.SECTION_B.MODE() == 1

    parameters.set(SECTION_A={"VALUE": 4}, SECTION_B={"MODE": 4})
    assert parameters.SECTION_A.VALUE() == 4
    assert parameters.SECTION_B.MODE() == 4


def test_set_runs_unconditional_root_validator():
    root = cd.InputParametersDefinition(
        [cd.InputSectionDefinition("CONTROL", [V("VALUE", 1)])],
        validators=_root_validator_reject_all,
    )
    parameters = root.create_object()

    with pytest.raises(DataValidityError, match="root validator rejects"):
        parameters.CONTROL.VALUE = 2

    assert parameters.CONTROL.VALUE() == 1

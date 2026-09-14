import pyparsing as pp
import inspect
import io
import pickle
import re
import warnings
import numpy as np
import pytest
from functools import partial

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:  # Just a linter worshiping
    from ...common.grammar import generate_grammar
    from ...common import grammar_types as gt
    from .. import input_parameters_definitions as cd
    from .. import input_parameters as input_parameters
    from ...common.configuration_containers import Section, CustomSection
    from ...common.configuration_transaction import (
        ConfigurationTransaction,
        PostChangeHookError,
    )
    from ...common.backward_compatibility import ExceptionGroup
    from ...common.options import Option, CustomOption, DangerousValue
    from ...common.configuration_definitions import gather, switch
    from ...common.generated_configuration_definitions import GeneratedValueDefinition, Length, NumpyViewDefinition
    from ...common.warnings import (
        DataValidityError,
        DataValidityErrors,
        DataValidityWarning,
    )

V = cd.InputValueDefinition


def ar(x):
    return np.atleast_1d(x)


def _validator_warning_if_two(option, values, why):
    assert values is option._container
    if option() == 2:
        return DataValidityWarning(f"{why}: VALUE should not be two")


def _validator_error_if_three(option, values, why):
    if option() == 3:
        return DataValidityError(f"{why}: VALUE must not be three")


def _warning_condition_if_two(value):
    if value == 2:
        return "VALUE should not be two"


def _section_validator_warning(section, values, why):
    assert values is section
    if values["VALUE"]() == 4:
        return DataValidityWarning(f"{why}: section received its runtime object")


def _root_validator_warning(root, values, why):
    assert values is root
    if values["CONTROL"]["VALUE"]() == 5:
        return DataValidityWarning(f"{why}: root received its runtime object")


def _root_validator_matching_values(root, values, why):
    if values["SECTION_A"]["VALUE"]() != values["SECTION_B"]["MODE"]():
        return DataValidityError(f"{why}: values in SECTION_A and SECTION_B must match")


def _root_validator_reject_all(root, values, why):
    return DataValidityError(f"{why}: root validator rejects the configuration")


def _repeated_validator_reject_two(section, values, why):
    if values["VALUE"]() == 2:
        return DataValidityError(f"{why}: repeated VALUE must not be two")


def _generated_source_length_validator(root, values, why):
    source = values["CONTROL"]["VALUES"]()
    if source is not None and len(source) > 1:
        return DataValidityError(f"{why}: generated setter produced too many values")


def _get_cross_section_generated(section, key=None):
    target = section._get_root_container().find("SECTION_B.VALUE")
    return target()


def _set_cross_section_generated(section, value, key=None):
    transaction = ConfigurationTransaction.current(section)
    target = section._get_root_container().find("SECTION_B.VALUE")
    target.stage(transaction, value)


def _reject_cross_section_two(root, values, why):
    if values["SECTION_B"]["VALUE"]() == 2:
        return DataValidityError(f"{why}: cross-section generated value must not be two")


def _reject_negative_raw_data(root, values, why):
    raw = values["DATA"]["RAW"]()
    if raw is not None and np.any(raw < 0):
        return DataValidityError(f"{why}: raw data must not be negative")


def _required_in_mode_one(option):
    return option._container["MODE"]() == 1


class TestInputParameters(TestCase):
    def assertParse(self, text, value, grammar):
        out = self.parse(text, grammar).asList()
        self.assertEqual(len(out), 1)
        out = out[0]
        if hasattr(out, "to_dict"):
            out = out.to_dict()
        self.assertEqual(out, value)

    def parse(self, text, grammar):
        if not isinstance(grammar, pp.ParserElement):
            grammar = grammar()
        return grammar.parse_string(text, True)

    def assertNotValid(self, text, grammar=None):
        self.assertRaises(pp.ParseBaseException, lambda: self.parse(text, grammar))

    def test_section_delimiter_value(self):
        with generate_grammar():
            grammar = cd.InputParametersDefinition.delimiter()
            grammar = "a" + grammar + "b"
        for w in ["a b", "a\n b", "a\n\n b", "a\n \n b", "a \n\n b", "a\n\n\n b"]:
            self.assertRaises(pp.ParseException, lambda: grammar.parse_string(w, True))
        for w in ["a\nb", "a \nb", "a\n  \nb", "a  \n \nb", "a \n\t\nb", "a\n\n\nb"]:
            self.assertEqual(["a", "b"], grammar.parse_string(w, True).asList())

    #

    def test_custom_value(self):
        with generate_grammar():
            cv = cd.InputSectionDefinition.custom_member_grammar(lambda x: x[0] not in ["aaa"])
        self.assertTrue("\n" not in cv.whiteChars)

        assertParse = partial(self.assertParse, grammar=lambda: cv)
        assertNotValid = partial(self.assertNotValid, grammar=lambda: cv)

        assertParse("bbb=1", ("bbb", 1))
        assertParse("bbb=1.3", ("bbb", 1.3))
        assertParse("bbb", ("bbb", True))
        assertNotValid("aaa")
        assertNotValid("aaa=1")
        with generate_grammar():
            cv = cv + "a"
        out = self.parse("bbb=1 a", cv)
        self.assertEqual(list(out), [("bbb", 1), "a"])
        out = self.parse("bbb=1a a", cv)
        self.assertEqual(list(out), [("bbb", "1a"), "a"])
        assertNotValid("bbb=1\na")

    def test_dangerous_value(self):
        with generate_grammar():
            cv = cd.InputValueDefinition("aaa", int, is_required=True).grammar()

        assertParse = partial(self.assertParse, grammar=lambda: cv)
        assertNotValid = partial(self.assertNotValid, grammar=lambda: cv)

        def assertParseDangerous(text, result):
            out = cv.parse_string(text, True)[0]
            assert out[0] == result[0]
            self.assertEqual(out[1].value, result[1])

        assertParse("aaa=1", ("aaa", 1))
        assertNotValid("aaa")
        assertNotValid("aaa=sss")
        assertNotValid("aaa=1.0")

        with generate_grammar():
            cv = cd.InputValueDefinition("aaa", int, is_required=True).grammar(allow_dangerous=True)

        assertParse("aaa=1", ("aaa", 1))
        assertParseDangerous("aaa=1.0", ("aaa", 1.0))
        assertParseDangerous("aaa=sss", ("aaa", "sss"))
        assertParseDangerous("aaa={2,3}", ("aaa", np.array([2, 3])))

        vd = cd.InputValueDefinition("aaa", 2, is_required=True)
        o = vd.create_object()
        o.set_dangerous("AAA")
        self.assertEqual(o.to_string(), "\taaa=AAA")
        o.set_dangerous(None)
        self.assertEqual(o.to_string(), "")
        #

    def test_is_required(self):
        ipd = cd.InputParametersDefinition.definition_from_dict(
            {"ENERGY": [V("E", 1), V("F", int, None, is_required=lambda option: option._container["E"]() != 2)]}
        )

        ip = ipd.create_object()
        with pytest.raises(DataValidityError):
            ip.validate()
        e = ip.ENERGY
        e.F = 2
        ip.validate()
        with pytest.raises(DataValidityError):
            e.F = None
        e.E = 2
        ip.validate()
        e.F = None
        ip.validate()

        assert ipd.read_from_string("ENERGY E=2").ENERGY.to_dict() == {"E": 2}
        assert ipd.read_from_string("ENERGY E=1 F=2").ENERGY.to_dict() == {"E": 1, "F": 2}
        with pytest.warns(DataValidityError):
            ipd.read_from_string("ENERGY E=1")

    def test_write_condition(self):
        input_parameters_def = cd.InputParametersDefinition.definition_from_dict({"ENERGY": [V("E", 1)]})
        # del input_parameters_def._members['TASK']
        id = input_parameters_def.read_from_file(io.StringIO("ENERGY E=1"))
        self.assertEqual(id.to_string(), "ENERGY\n\tE=1\n")
        input_parameters_def["ENERGY"]["E"].write_condition = lambda o: True
        self.assertEqual(id.to_string(), "ENERGY\n\tE=1\n")
        input_parameters_def["ENERGY"]["E"].write_condition = lambda o: False
        self.assertEqual(id.to_string(), "ENERGY\n")
        input_parameters_def["ENERGY"].write_condition = lambda o: True
        self.assertEqual(id.to_string(), "ENERGY\n")
        input_parameters_def["ENERGY"].write_condition = lambda o: False
        self.assertEqual(id.to_string(), "")
        #

    def test_validators(self):
        definition = V(
            "VALUE",
            1,
            validators=(_validator_warning_if_two, _validator_error_if_three),
        )
        copied = definition.copy()
        assert copied.validators == (_validator_warning_if_two, _validator_error_if_three)
        pickle.dumps(copied)

        ipd = cd.InputParametersDefinition.definition_from_dict({"CONTROL": [copied]})

        parameters = ipd.create_object()
        with pytest.warns(DataValidityWarning, match="set: VALUE should not be two") as caught:
            parameters.CONTROL.VALUE = 2
        assert len(caught) == 1

        parameters = ipd.create_object()
        with pytest.warns(DataValidityWarning, match="set: VALUE should not be two") as caught:
            parameters.set(VALUE=2)
        assert len(caught) == 1

        with pytest.warns(DataValidityWarning, match="parse: VALUE should not be two") as caught:
            parameters = ipd.read_from_string("CONTROL\n\tVALUE=2\n")
        assert len(caught) == 1

        with pytest.warns(DataValidityWarning, match="save: VALUE should not be two") as caught:
            parameters.to_string(validate="save")
        assert len(caught) == 1

        parameters = ipd.create_object()
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

        parameters = ipd.create_object()
        with pytest.warns(DataValidityError, match="set: VALUE must not be three"):
            parameters.CONTROL.VALUE.set(
                3,
                retain_invalid="typed",
                report_invalid="warn",
            )
        assert parameters.CONTROL.VALUE() == 3

        with pytest.warns(DataValidityError, match="parse: VALUE must not be three"):
            invalid = ipd.read_from_string("CONTROL\n\tVALUE=3\n")
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
        with pytest.warns(DataValidityWarning, match="set: section received its runtime object") as caught:
            parameters.CONTROL.VALUE = 4
        assert len(caught) == 1
        with pytest.warns(DataValidityWarning, match="set: root received its runtime object") as caught:
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

    def test_cross_section_validator_for_proposed_values(self):
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

        with pytest.raises(DataValidityError, match="values in SECTION_A and SECTION_B must match"):
            parameters.set({"SECTION_A.VALUE": 2, "SECTION_B.MODE": 3})
        assert parameters.SECTION_A.VALUE() == 1
        assert parameters.SECTION_B.MODE() == 1

        parameters.set(SECTION_A={"VALUE": 4}, SECTION_B={"MODE": 4})
        assert parameters.SECTION_A.VALUE() == 4
        assert parameters.SECTION_B.MODE() == 4

    def test_set_runs_unconditional_root_validator(self):
        root = cd.InputParametersDefinition(
            [cd.InputSectionDefinition("CONTROL", [V("VALUE", 1)])],
            validators=_root_validator_reject_all,
        )
        parameters = root.create_object()

        with pytest.raises(DataValidityError, match="root validator rejects"):
            parameters.CONTROL.VALUE = 2

        assert parameters.CONTROL.VALUE() == 1

    def test_transaction_converts_each_value_once(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", 1)]}
        )
        parameters = definition.create_object()
        option = parameters.CONTROL.VALUE
        value_definition = option._definition
        original_convert = value_definition.convert_value
        packed = []

        def count_conversion(option, value, item=False):
            packed.append(value)
            return original_convert(option, value, item)

        value_definition.convert_value = count_conversion
        option.set(2)

        assert packed == [2]
        assert option() == 2

    def test_transaction_commit_uses_the_validated_values(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("FIRST", 1), V("SECOND", 1)]}
        )
        parameters = definition.create_object()
        option = parameters.CONTROL.SECOND
        value_definition = option._definition
        original_convert = value_definition.convert_value
        calls = 0

        def fail_on_second_conversion(option, value, item=False):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise RuntimeError("value was converted again during commit")
            return original_convert(option, value, item)

        value_definition.convert_value = fail_on_second_conversion
        parameters.CONTROL.set(FIRST=2, SECOND=2)

        assert calls == 1
        assert parameters.CONTROL.FIRST() == 2
        assert parameters.CONTROL.SECOND() == 2

    def test_numbered_item_assignment_converts_new_value_once(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("VALUES", int, is_repeated=V.Repeated.DICT),
                ]
            }
        )
        option = definition.create_object().CONTROL.VALUES
        value_definition = option._definition
        original_convert = value_definition.convert_value
        converted = []

        def count_conversion(option, value, item=False):
            converted.append(value)
            return original_convert(option, value, item)

        value_definition.convert_value = count_conversion
        option[1] = 7

        assert converted == [7]
        assert option[1] == 7

    def test_post_change_hook_errors_are_aggregated_after_commit(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("FIRST", 1), V("SECOND", 1)]}
        )
        parameters = definition.create_object()
        observed = []

        def reject_first(option):
            observed.append((option.name, ConfigurationTransaction.current(option)))
            raise RuntimeError("first hook failed")

        def reject_second(option):
            observed.append((option.name, ConfigurationTransaction.current(option)))
            raise ValueError("second hook failed")

        parameters.CONTROL.FIRST.add_hook(reject_first)
        parameters.CONTROL.SECOND.add_hook(reject_second)

        with pytest.raises(PostChangeHookError, match="already committed") as caught:
            parameters.CONTROL.set(FIRST=2, SECOND=2)

        assert observed == [("FIRST", None), ("SECOND", None)]
        assert [type(error) for error in caught.value.exceptions] == [
            RuntimeError,
            ValueError,
        ]
        assert isinstance(caught.value, ExceptionGroup)
        assert parameters.CONTROL.FIRST() == 2
        assert parameters.CONTROL.SECOND() == 2

    def test_validation_and_post_change_hook_errors_are_reported_together(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("VALUE", 1, validators=_validator_error_if_three),
                ]
            }
        )
        parameters = definition.create_object()

        def reject(_option):
            raise RuntimeError("hook failed")

        parameters.CONTROL.VALUE.add_hook(reject)

        with pytest.raises(DataValidityErrors) as caught:
            parameters.CONTROL.VALUE.set(
                3,
                retain_invalid="typed",
                report_invalid="raise",
            )

        validation_error, hook_error = caught.value.exceptions
        assert isinstance(validation_error, DataValidityError)
        assert isinstance(hook_error, PostChangeHookError)
        assert isinstance(hook_error.exceptions[0], RuntimeError)
        assert isinstance(
            caught.value.subgroup(DataValidityError),
            DataValidityErrors,
        )
        assert not isinstance(
            caught.value.subgroup(PostChangeHookError),
            DataValidityErrors,
        )
        assert parameters.CONTROL.VALUE() == 3

    def test_failed_indexed_array_hook_does_not_roll_back_value(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUES", gt.Array(int), [1, 2])]}
        )
        parameters = definition.create_object()

        def reject(_option):
            raise RuntimeError("hook rejected indexed assignment")

        parameters.CONTROL.VALUES.add_hook(reject)

        with pytest.raises(PostChangeHookError, match="already committed"):
            parameters.CONTROL.VALUES[0] = 9

        assert np.array_equal(parameters.CONTROL.VALUES(), [9, 2])
        assert parameters.CONTROL.VALUES.is_set()

    def test_indexed_array_materializes_default_without_modifying_definition(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUES", gt.Array(int), [1, 2])]}
        )
        first = definition.create_object()
        second = definition.create_object()

        first.CONTROL.VALUES[0] = 9

        assert np.array_equal(first.CONTROL.VALUES(), [9, 2])
        assert np.array_equal(second.CONTROL.VALUES(), [1, 2])
        assert np.array_equal(
            definition["CONTROL"]["VALUES"].default_value,
            [1, 2],
        )

    def test_indexed_array_materializes_result(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUES", gt.Array(int), [1, 2])]}
        )
        option = definition.create_object().CONTROL.VALUES
        option.result = np.array([3, 4])

        option[0] = 8

        assert np.array_equal(option(), [8, 4])
        assert option._value is not None
        assert not hasattr(option, "_result")

    def test_transaction_context_manager_commits_or_rolls_back(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", 1)]}
        )
        parameters = definition.create_object()
        option = parameters.CONTROL.VALUE

        with ConfigurationTransaction(parameters) as transaction:
            option.stage(transaction, 2)
            assert option() == 2
        assert option() == 2

        with pytest.raises(RuntimeError):
            with ConfigurationTransaction(parameters) as transaction:
                option.stage(transaction, 3)
                assert option() == 3
                raise RuntimeError("abort")
        assert option() == 2

    def test_direct_stage_uses_transaction_validation_policy(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", int, 1)]}
        )
        parameters = definition.create_object()
        option = parameters.CONTROL.VALUE

        with pytest.raises(DataValidityError):
            with ConfigurationTransaction(parameters) as transaction:
                option.stage(transaction, "not an integer")

        assert option() == 1

    def test_transaction_mutation_requires_context(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", 1)]}
        )
        parameters = definition.create_object()
        transaction = ConfigurationTransaction(parameters)

        with pytest.raises(RuntimeError, match="context manager"):
            parameters.CONTROL.VALUE.stage(transaction, 2)
        assert parameters.CONTROL.VALUE() == 1

        with transaction:
            parameters.CONTROL.VALUE.stage(transaction, 2)
        assert parameters.CONTROL.VALUE() == 2

        with pytest.raises(RuntimeError, match="context manager"):
            parameters.CONTROL.VALUE.stage(transaction, 3)
        assert parameters.CONTROL.VALUE() == 2

        with transaction:
            parameters.CONTROL.VALUE.stage(transaction, 3)
        assert parameters.CONTROL.VALUE() == 3
        assert transaction._changes == []

        with transaction:
            parameters.CONTROL.VALUE.stage(transaction, 4)
            transaction.abort()
        assert parameters.CONTROL.VALUE() == 3

        assert not hasattr(transaction, "commit")
        assert not hasattr(transaction, "rollback")

    def test_transaction_use_joins_active_transaction(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", 1)]}
        )
        parameters = definition.create_object()

        with ConfigurationTransaction.use(parameters) as transaction:
            assert parameters._active_transaction is transaction
            with ConfigurationTransaction.use(parameters.CONTROL) as nested:
                assert nested is transaction
                parameters.CONTROL.VALUE.set(2)
            boundaries = [
                change
                for change in transaction._changes
                if isinstance(change, ConfigurationTransaction._Boundary)
            ]
            assert len(boundaries) > 1
            assert boundaries[0]._active
            assert not boundaries[-1]._active

        assert parameters.CONTROL.VALUE() == 2
        assert not hasattr(parameters, "_active_transaction")

    def test_caught_nested_set_rolls_back_only_nested_changes(self):
        def reject_two(option, values, _why):
            if option() == 2:
                return DataValidityError("two rejected")

        definition = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("FIRST", 1),
                    V("SECOND", 1, validators=reject_two),
                ]
            }
        )
        parameters = definition.create_object()

        with ConfigurationTransaction(parameters):
            parameters.CONTROL.FIRST.set(2)
            with pytest.raises(DataValidityError, match="two rejected"):
                parameters.CONTROL.SECOND.set(2)

        assert parameters.CONTROL.FIRST() == 2
        assert parameters.CONTROL.SECOND() == 1

    def test_uncaught_nested_set_rolls_back_the_whole_transaction(self):
        def reject_two(option, values, _why):
            if option() == 2:
                return DataValidityError("two rejected")

        definition = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("FIRST", 1),
                    V("SECOND", 1, validators=reject_two),
                ]
            }
        )
        parameters = definition.create_object()

        with pytest.raises(DataValidityError, match="two rejected"):
            with ConfigurationTransaction(parameters):
                parameters.CONTROL.FIRST.set(2)
                parameters.CONTROL.SECOND.set(2)

        assert parameters.CONTROL.FIRST() == 1
        assert parameters.CONTROL.SECOND() == 1

    def test_transaction_savepoint_rolls_back_only_its_changes(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("FIRST", 1), V("SECOND", 1)]}
        )
        parameters = definition.create_object()

        with ConfigurationTransaction(parameters) as transaction:
            parameters.CONTROL.FIRST.stage(transaction, 2)
            with pytest.raises(RuntimeError, match="abort savepoint"):
                with transaction.savepoint():
                    parameters.CONTROL.SECOND.stage(transaction, 3)
                    raise RuntimeError("abort savepoint")

            assert parameters.CONTROL.FIRST() == 2
            assert parameters.CONTROL.SECOND() == 1

        assert parameters.CONTROL.FIRST() == 2
        assert parameters.CONTROL.SECOND() == 1

    def test_transaction_savepoint_rolls_back_repeated_option_change(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", 1)]}
        )
        parameters = definition.create_object()

        with ConfigurationTransaction(parameters) as transaction:
            parameters.CONTROL.VALUE.stage(transaction, 2)
            with pytest.raises(RuntimeError, match="abort savepoint"):
                with transaction.savepoint():
                    parameters.CONTROL.VALUE.stage(transaction, 3)
                    raise RuntimeError("abort savepoint")

            assert parameters.CONTROL.VALUE() == 2

        assert parameters.CONTROL.VALUE() == 2

    def test_option_hook_is_coalesced_and_runs_outside_transaction(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", 1)]}
        )
        parameters = definition.create_object()
        calls = []
        parameters.CONTROL.VALUE.add_hook(
            lambda option: calls.append(
                (option(), ConfigurationTransaction.current(option))
            )
        )

        with ConfigurationTransaction(parameters) as transaction:
            parameters.CONTROL.VALUE.stage(transaction, 2)
            parameters.CONTROL.VALUE.stage(transaction, 3)

        assert calls == [(3, None)]

    def test_option_stage_registers_its_change(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", 1)]}
        )
        option = definition.create_object().CONTROL.VALUE

        with ConfigurationTransaction(option) as transaction:
            with transaction.savepoint() as savepoint:
                assert option.stage(transaction, 2)
                assert isinstance(transaction._changes[-1], Option.Change)
                option.stage(transaction, 3)
                savepoint.rollback()
        assert option() == 1

    def test_unknown_member_is_converted_before_transaction_commit(self):
        section = cd.InputSectionDefinition("CONTROL", [V("KNOWN", 1)])
        section.custom_class = CustomOption.factory(V, int)
        parameters = cd.InputParametersDefinition([section]).create_object()

        with pytest.raises(DataValidityError):
            parameters.CONTROL.set(
                {"KNOWN": 2, "EXTRA": "not-an-integer"},
                unknown="add",
            )

        assert parameters.CONTROL.KNOWN() == 1
        assert "EXTRA" not in parameters.CONTROL

    def test_add_with_invalid_value_rolls_back_custom_member(self):
        section = cd.InputSectionDefinition("CONTROL", [V("KNOWN", 1)])
        section.custom_class = CustomOption.factory(V, int)
        definition = cd.InputParametersDefinition([section])

        def reject_extra(_root, values, _why):
            if "EXTRA" in values["CONTROL"]:
                return DataValidityError("custom member rejected")

        definition.validators = (reject_extra,)
        parameters = definition.create_object()

        with pytest.raises(DataValidityError, match="custom member rejected"):
            parameters.CONTROL.add("EXTRA", 2)

        assert "EXTRA" not in parameters.CONTROL

    def test_discard_and_ignore_unconvertible_staged_values(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", 1)]}
        )
        parameters = definition.create_object()

        parameters.set(
            VALUE="not-an-integer",
            retain_invalid="none",
            report_invalid="ignore",
        )

        assert parameters.CONTROL.VALUE() == 1

    def test_discarded_invalid_value_does_not_add_custom_member(self):
        section = cd.InputSectionDefinition("CONTROL", [V("KNOWN", 1)])
        section.custom_class = CustomOption.factory(V, int)
        parameters = cd.InputParametersDefinition([section]).create_object()

        parameters.CONTROL.set(
            {"EXTRA": "not-an-integer"},
            unknown="add",
            retain_invalid="none",
            report_invalid="ignore",
        )

        assert "EXTRA" not in parameters.CONTROL

    def test_typed_and_unconvertible_invalid_values_are_distinguished(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", gt.Unsigned(), 1)]}
        )

        parameters = definition.create_object()
        with pytest.warns(DataValidityError, match="positive"):
            parameters.set(
                VALUE=-1,
                retain_invalid="none",
                report_invalid="warn",
            )
        assert parameters.CONTROL.VALUE() == 1

        with pytest.warns(DataValidityError, match="positive"):
            parameters.set(
                VALUE=-1,
                retain_invalid="typed",
                report_invalid="warn",
            )
        assert parameters.CONTROL.VALUE() == -1
        assert not parameters.CONTROL.VALUE.is_dangerous()

        parameters = definition.create_object()
        with pytest.warns(DataValidityError, match="unsigned integer"):
            parameters.set(
                VALUE="not-an-integer",
                retain_invalid="typed",
                report_invalid="warn",
            )
        assert parameters.CONTROL.VALUE() == 1

        with pytest.warns(DataValidityError, match="unsigned integer"):
            parameters.set(
                VALUE="not-an-integer",
                retain_invalid="all",
                report_invalid="warn",
            )
        assert parameters.CONTROL.VALUE() == "not-an-integer"
        assert parameters.CONTROL.VALUE.is_dangerous()
        assert "VALUE=not-an-integer" in parameters.to_string(validate=False)

        array_definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUES", gt.Array(gt.Unsigned()), [1])]}
        )
        array = array_definition.create_object()
        with pytest.warns(DataValidityError, match="positive"):
            array.set(
                VALUES=[-1],
                retain_invalid="typed",
                report_invalid="warn",
            )
        assert np.array_equal(array.CONTROL.VALUES(), [-1])

        array = array_definition.create_object()
        with pytest.warns(DataValidityError):
            array.set(
                VALUES=["not-an-integer"],
                retain_invalid="typed",
                report_invalid="warn",
            )
        assert np.array_equal(array.CONTROL.VALUES(), [1])

        batch = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("VALID", 1),
                    V("INVALID", gt.Unsigned(), 1),
                ]
            }
        ).create_object()
        with pytest.warns(DataValidityError):
            batch.set(
                VALID=2,
                INVALID="not-an-integer",
                retain_invalid="none",
                report_invalid="warn",
            )
        assert batch.CONTROL.VALID() == 2
        assert batch.CONTROL.INVALID() == 1

    def test_retention_and_reporting_are_independent(self):
        typed = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", gt.Unsigned(), 1)]}
        ).create_object()

        with pytest.raises(DataValidityError, match="positive"):
            typed.set(
                VALUE=-1,
                retain_invalid="typed",
                report_invalid="raise",
            )
        assert typed.CONTROL.VALUE() == -1

        semantic = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("VALUE", 1, validators=_validator_error_if_three)
                ]
            }
        ).create_object()
        with pytest.warns(DataValidityError, match="must not be three"):
            semantic.set(
                VALUE=3,
                retain_invalid="none",
                report_invalid="warn",
            )
        assert semantic.CONTROL.VALUE() == 1

        with pytest.raises(DataValidityError, match="must not be three"):
            semantic.set(
                VALUE=3,
                retain_invalid="typed",
                report_invalid="raise",
            )
        assert semantic.CONTROL.VALUE() == 3

        required = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUE", int, None)]}
        ).create_object()
        required.CONTROL.VALUE = 1
        with pytest.warns(DataValidityError, match="must have a value"):
            required.CONTROL.VALUE.set(
                None,
                retain_invalid="none",
                report_invalid="warn",
            )
        assert required.CONTROL.VALUE() == 1

        with pytest.raises(DataValidityError, match="must have a value"):
            required.CONTROL.VALUE.set(
                None,
                retain_invalid="typed",
                report_invalid="raise",
            )
        assert required.CONTROL.VALUE() is None

    def test_strict_batch_collects_all_value_errors_before_raising(self):
        parameters = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("FIRST", gt.Unsigned(), 1),
                    V("SECOND", gt.Unsigned(), 2),
                    V("THIRD", int, 3),
                    V("FOURTH", int, 4),
                ]
            }
        ).create_object()

        with pytest.raises(DataValidityErrors) as caught:
            parameters.CONTROL.set(
                FIRST=-1,
                SECOND=-2,
                THIRD="not-an-integer",
                FOURTH="also-not-an-integer",
            )

        assert len(caught.value.exceptions) == 4
        assert parameters.CONTROL.FIRST() == 1
        assert parameters.CONTROL.SECOND() == 2
        assert parameters.CONTROL.THIRD() == 3
        assert parameters.CONTROL.FOURTH() == 4

    def test_repeated_dictionary_collects_all_errors_atomically(self):
        parameters = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V(
                        "VALUES",
                        gt.Unsigned(),
                        is_optional=True,
                        is_repeated=V.Repeated.DICT,
                    )
                ]
            }
        ).create_object()

        with pytest.raises(DataValidityErrors) as caught:
            parameters.CONTROL.VALUES = {1: -1, 2: -2}

        assert len(caught.value.exceptions) == 2
        assert parameters.CONTROL.VALUES(all_values=True) is None

    def test_repeated_values_use_invalid_value_policy(self):
        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "ROWS",
                    [V("VALUE", gt.Unsigned(), 1)],
                    is_repeated=True,
                )
            ]
        )
        parameters = definition.create_object()

        with pytest.warns(DataValidityError, match="positive"):
            parameters.ROWS.set(
                [{"VALUE": -1}],
                retain_invalid="typed",
                report_invalid="warn",
            )
        assert parameters.ROWS[0].VALUE() == -1

        with pytest.warns(DataValidityError, match="unsigned integer"):
            parameters.ROWS.set(
                [{"VALUE": "not-an-integer"}],
                retain_invalid="all",
                report_invalid="warn",
            )
        assert parameters.ROWS[0].VALUE() == "not-an-integer"
        assert parameters.ROWS[0].VALUE.is_dangerous()

    def test_dotted_unknown_member_is_not_added(self):
        section = cd.InputSectionDefinition("CONTROL", [V("KNOWN", 1)])
        section.custom_class = CustomOption.factory(V, int)
        parameters = cd.InputParametersDefinition([section]).create_object()

        with pytest.raises(KeyError, match="dotted path"):
            parameters.CONTROL.set(
                {"EXTRA.CHILD": 1},
                unknown="add",
            )

        assert "EXTRA" not in parameters.CONTROL

    def test_unknown_add_accepts_only_non_dotted_names(self):
        parameters = cd.InputParametersDefinition(
            [cd.InputSectionDefinition("CONTROL", [V("KNOWN", 1)])]
        ).create_object()

        with pytest.raises(KeyError, match="EXTRA.VALUE"):
            parameters.get_member("EXTRA.VALUE")
        assert list(parameters.get_members("EXTRA.VALUE")) == []
        with pytest.raises(KeyError, match="dotted path"):
            parameters.set("EXTRA.VALUE", 2, unknown="add")
        assert "EXTRA" not in parameters

        parameters.set({"EXTRA": {"VALUE": 2}}, unknown="add")
        assert parameters.EXTRA.VALUE() == 2

    def test_failed_read_rolls_back_clear_and_parsed_values_together(self):
        def reject_parsed_two(_root, values, why):
            if why == "parse" and values["CONTROL"]["VALUE"]() == 2:
                raise RuntimeError("parsed value rejected")

        definition = cd.InputParametersDefinition(
            [cd.InputSectionDefinition("CONTROL", [V("VALUE", 0)])],
            validators=reject_parsed_two,
        )
        parameters = definition.create_object()
        parameters.CONTROL.VALUE = 1

        with pytest.raises(RuntimeError, match="parsed value rejected"):
            parameters.read_from_file(io.StringIO("CONTROL\n\tVALUE=2\n"))

        assert parameters.CONTROL.VALUE() == 1
        assert not hasattr(parameters, "_parsed_values")
        assert not hasattr(parameters.CONTROL, "_parsed_values")

    def test_exact_resolution_stops_at_option(self):
        parameters = cd.InputParametersDefinition(
            [cd.InputSectionDefinition("CONTROL", [V("KNOWN", 1)])]
        ).create_object()

        with pytest.raises(KeyError, match="CONTROL.KNOWN.CHILD"):
            parameters.get_member("CONTROL.KNOWN.CHILD")
        assert list(parameters.get_members("CONTROL.KNOWN.CHILD")) == []
        with pytest.raises(KeyError, match="CONTROL.KNOWN.CHILD"):
            parameters.set("CONTROL.KNOWN.CHILD", 2, unknown="add")

    def test_unknown_add_does_not_replace_existing_member(self):
        definition = cd.InputParametersDefinition(
            [cd.InputSectionDefinition("CONTROL", [V("VALUE", 1)])]
        )
        parameters = definition.create_object()
        original = parameters.CONTROL

        with pytest.raises(TypeError, match="already in"):
            parameters.set({"CONTROL": 2}, unknown="add")

        assert parameters.CONTROL is original

    def test_transaction_uses_standard_container_access(self):
        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "OUTER",
                    [cd.InputSectionDefinition("INNER", [V("VALUE", 1)])],
                )
            ]
        )
        parameters = definition.create_object()

        with ConfigurationTransaction(parameters):
            assert parameters["OUTER"]["INNER"]["VALUE"]() == 1

    def test_post_change_hook_cannot_mutate_its_configuration(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("FIRST", 1), V("SECOND", 1)]}
        )
        parameters = definition.create_object()

        def mutate(_option):
            parameters.CONTROL.SECOND.set(2)

        parameters.CONTROL.FIRST.add_hook(mutate)

        with pytest.raises(PostChangeHookError, match="already committed") as caught:
            parameters.CONTROL.FIRST.set(2)

        assert len(caught.value.exceptions) == 1
        assert "post-change hooks" in str(caught.value.exceptions[0])
        assert parameters.CONTROL.FIRST() == 2
        assert parameters.CONTROL.SECOND() == 1

    def test_rolled_back_savepoint_does_not_register_hook(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("FIRST", 1),
                    V("SECOND", 2),
                ]
            }
        )
        parameters = definition.create_object()

        calls = []
        parameters.CONTROL.FIRST.add_hook(lambda option: calls.append(option()))
        parameters.CONTROL.SECOND.add_hook(lambda option: calls.append(option()))

        with ConfigurationTransaction(parameters) as transaction:
            parameters.CONTROL.FIRST.stage(transaction, 3)
            with transaction.savepoint() as savepoint:
                parameters.CONTROL.SECOND.stage(transaction, 4)
                savepoint.rollback()

        assert calls == [3]
        assert parameters.CONTROL.FIRST() == 3
        assert parameters.CONTROL.SECOND() == 2

    def test_clear_is_validated_and_rolled_back(self):
        def reject_missing(option, values, _why):
            if option() is None:
                return DataValidityError("missing value rejected")

        definition = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V(
                        "FIRST",
                        int,
                        None,
                        is_optional=True,
                        is_required=False,
                        validators=reject_missing,
                    ),
                    V("SECOND", 1),
                ]
            }
        )
        parameters = definition.create_object()
        parameters.CONTROL.FIRST = 1

        with pytest.raises(DataValidityError, match="missing value rejected"):
            parameters.CONTROL.FIRST.clear()
        assert parameters.CONTROL.FIRST() == 1

        with pytest.raises(DataValidityError, match="missing value rejected"):
            parameters.CONTROL.clear()
        assert parameters.CONTROL.FIRST() == 1
        assert parameters.CONTROL.SECOND() == 1

        required_definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("REQUIRED", int, None, is_required=True)]}
        )
        required = required_definition.create_object()
        required.CONTROL.REQUIRED = 1

        with pytest.raises(DataValidityError, match="must have a value"):
            required.CONTROL.clear()

        assert required.CONTROL.REQUIRED() == 1

    def test_setting_none_uses_required_clear_semantics(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("REQUIRED", int, None, is_required=True)]}
        )
        parameters = definition.create_object()
        parameters.CONTROL.REQUIRED.set(2)

        with pytest.raises(DataValidityError, match="must have a value"):
            parameters.CONTROL.REQUIRED.set(None)

        assert parameters.CONTROL.REQUIRED() == 2

    def test_container_clear_skips_generated_setters(self):
        generated_clears = []

        def get_generated(section, key=None):
            return section["SOURCE"]()

        def set_generated(section, value, key=None):
            generated_clears.append(value)
            transaction = ConfigurationTransaction.current(section)
            section["SOURCE"].stage(transaction, value)

        definition = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("SOURCE", int, None, is_optional=True),
                    GeneratedValueDefinition(
                        "GENERATED", get_generated, set_generated
                    ),
                ]
            }
        )
        parameters = definition.create_object()
        parameters.CONTROL.SOURCE = 1

        parameters.CONTROL.clear()

        assert parameters.CONTROL.SOURCE() is None
        assert generated_clears == []

        parameters.CONTROL.SOURCE = 2
        parameters.CONTROL.GENERATED.clear()

        assert parameters.CONTROL.SOURCE() is None
        assert generated_clears == [None]

    def test_repeated_clear_is_validated_and_rolled_back(self):
        def reject_empty(_root, values, _why):
            if len(values["ROWS"]) == 0:
                return DataValidityError("empty rows rejected")

        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "ROWS",
                    [V("VALUE", 1)],
                    is_repeated=True,
                )
            ],
            validators=reject_empty,
        )
        parameters = definition.create_object()
        parameters.ROWS.set([{"VALUE": 2}])

        with pytest.raises(DataValidityError, match="empty rows rejected"):
            parameters.ROWS.clear()

        assert len(parameters.ROWS) == 1
        assert parameters.ROWS[0].VALUE() == 2

    def test_empty_repeated_section_as_dict(self):
        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "ROWS",
                    [V("VALUE", 1)],
                    is_repeated=True,
                )
            ]
        )

        assert definition.create_object().ROWS.as_dict(False) is None

    def test_option_validation_uses_its_container_values(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("MODE", 1),
                    V(
                        "DEPENDENT",
                        int,
                        2,
                        is_required=_required_in_mode_one,
                    ),
                ]
            }
        )
        parameters = definition.create_object()

        parameters.CONTROL.DEPENDENT.validate("save")

    def test_setting_selector_validates_dependent_required_values(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("MODE", 2),
                    V(
                        "DEPENDENT",
                        int,
                        None,
                        is_optional=True,
                        is_required=_required_in_mode_one,
                    ),
                ]
            }
        )
        parameters = definition.create_object()

        with pytest.raises(DataValidityError, match="DEPENDENT"):
            parameters.CONTROL.MODE = 1

        assert parameters.CONTROL.MODE() == 2

    def test_unconditional_required_values_are_save_time_completeness_checks(self):
        definition = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("REQUIRED", int, None), V("OTHER", 1)]}
        )
        parameters = definition.create_object()

        parameters.CONTROL.OTHER = 2

        assert parameters.CONTROL.OTHER() == 2
        with pytest.raises(DataValidityError, match="REQUIRED"):
            parameters.validate("save")

    def test_disabled_section_is_not_validated(self):
        disabled = cd.InputSectionDefinition(
            "DISABLED", [V("REQUIRED", int, None)]
        )
        disabled.condition = lambda _definition, _container: False
        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition("ACTIVE", [V("VALUE", 1)]),
                disabled,
            ]
        )

        definition.create_object().validate("save")

    def test_copy_does_not_revalidate_existing_values(self):
        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "CONTROL",
                    [V("VALUE", 1, validators=_validator_error_if_three)],
                )
            ],
        )
        with pytest.warns(DataValidityError, match="VALUE must not be three"):
            parameters = definition.read_from_string("CONTROL\n\tVALUE=3\n")

        copied = parameters.copy()

        assert copied.CONTROL.VALUE() == 3

        with pytest.warns(DataValidityError):
            parameters.CONTROL.VALUE.set(
                "not-an-integer",
                retain_invalid="all",
                report_invalid="warn",
            )

        copied = parameters.copy(copy_values=True)

        assert copied.CONTROL.VALUE() == "not-an-integer"
        assert copied.CONTROL.VALUE.is_dangerous()

        custom = cd.InputParametersDefinition([]).create_object()
        custom.set({"EXTRA": {"VALUE": 2}}, unknown="add")

        copied_custom = custom.copy()

        assert type(copied_custom.EXTRA) is type(custom.EXTRA)
        assert type(copied_custom.EXTRA.VALUE) is type(custom.EXTRA.VALUE)

    def test_repeated_and_generated_sets_are_transactional(self):
        repeated = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "ROWS",
                    [V("VALUE", 0)],
                    is_repeated=True,
                    validators=_repeated_validator_reject_two,
                )
            ]
        ).create_object()
        repeated.ROWS.set([{"VALUE": 1}])
        with pytest.raises(DataValidityError, match="repeated VALUE must not be two"):
            repeated.ROWS.set([{"VALUE": 2}])
        assert repeated.ROWS[0].VALUE() == 1

        with pytest.raises(DataValidityError, match="repeated VALUE must not be two"):
            repeated.set(ROWS=[{"VALUE": 2}])
        assert repeated.ROWS[0].VALUE() == 1

        generated = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("NVALUES", Length("VALUES", default_values=1)),
                    V("VALUES", gt.Array(int), is_optional=True),
                ]
            }
        )
        generated.validators = (_generated_source_length_validator,)
        parameters = generated.create_object()
        parameters.CONTROL.NVALUES = 1
        with pytest.raises(DataValidityError, match="generated setter produced too many values"):
            parameters.CONTROL.NVALUES = 2
        assert np.array_equal(parameters.CONTROL.VALUES(), [1])

    def test_dictionary_repeated_set_replaces_or_merges_values(self):
        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "ROWS",
                    [V("VALUE", 0), V("OTHER", 0)],
                    is_repeated=V.Repeated.DICT_SECTION,
                )
            ]
        )
        parameters = definition.create_object()
        parameters.ROWS.set(
            {1: {"VALUE": 1, "OTHER": 10}, 2: {"VALUE": 2}}
        )

        parameters.ROWS.set({1: {"VALUE": 3}})

        assert list(parameters.ROWS) == [1]
        assert parameters.ROWS[1].VALUE() == 3
        assert parameters.ROWS[1].OTHER() == 0

        parameters.ROWS.set(
            {1: {"VALUE": 1, "OTHER": 10}, 2: {"VALUE": 2}}
        )
        parameters.ROWS.set({1: {"VALUE": 3}}, merge=True)

        assert list(parameters.ROWS) == [1, 2]
        assert parameters.ROWS[1].VALUE() == 3
        assert parameters.ROWS[1].OTHER() == 10
        assert parameters.ROWS[2].VALUE() == 2

        parameters.ROWS._definition.validators = (_repeated_validator_reject_two,)
        with pytest.raises(DataValidityError, match="repeated VALUE must not be two"):
            parameters.ROWS.set({1: {"VALUE": 2}}, merge=True)

        assert list(parameters.ROWS) == [1, 2]
        assert parameters.ROWS[1].VALUE() == 3
        assert parameters.ROWS[2].VALUE() == 2

    def test_list_repeated_set_replaces_or_merges_values(self):
        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "ROWS",
                    [V("VALUE", 0)],
                    is_repeated=True,
                )
            ]
        )
        parameters = definition.create_object()
        parameters.ROWS.set([{"VALUE": 1}, {"VALUE": 2}])

        parameters.ROWS.set([{"VALUE": 3}])
        assert [parameters.ROWS[i].VALUE() for i in parameters.ROWS] == [3]

        parameters.ROWS.set([{"VALUE": 1}, {"VALUE": 2}])
        parameters.ROWS.set([{"VALUE": 3}], merge=True)
        assert [parameters.ROWS[i].VALUE() for i in parameters.ROWS] == [
            1,
            2,
            3,
        ]

    def test_list_repeated_parse_merge_attaches_metadata_to_new_rows(self):
        def validate_parsed(option):
            parsed = getattr(option._container, "_parsed_values", None)
            if option() == 2 and parsed is None:
                raise RuntimeError("new row has no parse metadata")
            if parsed is not None and parsed["VALUE"] != option():
                raise RuntimeError("parse metadata belongs to another row")

        value = V("VALUE", 0)
        value.validate_parsed = validate_parsed
        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "ROWS", [value], is_repeated=True
                )
            ]
        )
        parameters = definition.create_object()
        parameters.ROWS.set([{"VALUE": 1}])

        parameters.ROWS.set(
            [{"VALUE": 2}], merge=True, validation_reason="parse"
        )

        assert [parameters.ROWS[i].VALUE() for i in parameters.ROWS] == [1, 2]

    def test_generated_setter_uses_mutable_cross_section_proposal(self):
        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "SECTION_A",
                    [
                        GeneratedValueDefinition(
                            "GENERATED",
                            _get_cross_section_generated,
                            _set_cross_section_generated,
                        )
                    ],
                ),
                cd.InputSectionDefinition("SECTION_B", [V("VALUE", 0)]),
            ],
            validators=_reject_cross_section_two,
        )
        parameters = definition.create_object()

        parameters.SECTION_A.GENERATED = 1
        assert parameters.SECTION_A.GENERATED() == 1
        assert parameters.SECTION_B.VALUE() == 1

        with pytest.raises(DataValidityError, match="cross-section generated value must not be two"):
            parameters.SECTION_A.GENERATED = 2
        assert parameters.SECTION_B.VALUE() == 1

        repeated_definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "ROWS",
                    [
                        GeneratedValueDefinition(
                            "GENERATED",
                            _get_cross_section_generated,
                            _set_cross_section_generated,
                        )
                    ],
                    is_repeated=True,
                ),
                cd.InputSectionDefinition("SECTION_B", [V("VALUE", 0)]),
            ],
            validators=_reject_cross_section_two,
        )
        repeated = repeated_definition.create_object()
        repeated.ROWS.set([{"GENERATED": 1}])
        assert repeated.SECTION_B.VALUE() == 1
        assert len(repeated.ROWS) == 1

        with pytest.raises(DataValidityError, match="cross-section generated value must not be two"):
            repeated.ROWS.set([{"GENERATED": 2}, {"GENERATED": 2}])
        assert repeated.SECTION_B.VALUE() == 1
        assert len(repeated.ROWS) == 1

    def test_generated_numpy_view_mutates_without_rollback(self):
        definition = cd.InputParametersDefinition(
            [
                cd.InputSectionDefinition(
                    "DATA",
                    [
                        V("RAW", gt.Array(float), is_optional=True),
                        NumpyViewDefinition("VIEW", "RAW"),
                    ],
                )
            ],
            validators=_reject_negative_raw_data,
        )
        parameters = definition.create_object()
        parameters.DATA.RAW = [1.0, 2.0, 3.0]

        with pytest.raises(ValueError):
            parameters.DATA.VIEW[:] = [4.0, 5.0]
        assert np.array_equal(parameters.DATA.RAW(), [1.0, 2.0, 3.0])

        with pytest.raises(DataValidityError, match="raw data must not be negative"):
            parameters.DATA.VIEW[1] = -2.0
        assert np.array_equal(parameters.DATA.RAW(), [1.0, -2.0, 3.0])

    def test_missing_sections_are_validation_results(self):
        definition = cd.InputParametersDefinition(
            [cd.InputSectionDefinition("REQUIRED", [V("VALUE", int, None)])]
        )
        found = definition.create_object().check_for_errors(print=False)
        assert any(
            isinstance(issue.message, DataValidityError)
            and "Non-optional section REQUIRED" in str(issue.message)
            for issue in found
        )

    @pytest.mark.slow
    def test_input_parameters_definition(self):
        input_parameters_def = cd.InputParametersDefinition.definition_from_dict(
            {
                "ENERGY": [
                    V("GRID", gt.SetOf(int, length=1), fixed_value=3),
                    V("NE", gt.SetOf(int, min_length=1)),
                    V("Ime", float, 0.0),
                ],
                "SITES": [V("NL", int)],
            }
        )
        grammar = input_parameters_def.sections["SITES"].values["NL"].grammar()
        assertParse = partial(self.assertParse, grammar=lambda: grammar)
        assertNotValid = partial(self.assertNotValid, grammar=lambda: grammar)

        assertParse("NL=3", ("NL", 3))
        grammar = input_parameters_def.sections["ENERGY"].values["NE"].grammar()
        assertParse("NE={3}", ("NE", ar(3)))
        grammar = input_parameters_def.sections["ENERGY"].values["Ime"].grammar()
        assertParse("Ime= 0.5", ("Ime", 0.5))
        grammar = input_parameters_def.sections["ENERGY"].values["GRID"].grammar()
        assertParse("GRID={3}", ("GRID", ar(3)))
        grammar = input_parameters_def.sections["ENERGY"].grammar()
        assertParse("ENERGY Ime= 0.5", ("ENERGY", {"Ime": 0.5}))
        assertParse("ENERGY Ime= 0.5 NE={5}", ("ENERGY", {"Ime": 0.5, "NE": ar(5)}))
        assertParse(
            """ENERGY Ime= 0.5
                                   NE={5}""",
            ("ENERGY", {"Ime": 0.5, "NE": ar(5)}),
        )
        assertParse(
            """ENERGY Ime= 0.5

                                   NE={5}""",
            ("ENERGY", {"Ime": 0.5, "NE": ar(5)}),
        )
        assertNotValid("""ENERGY Ime= 0.5

NE={5}""")

        # fixed value = 3
        assertNotValid("""ENERGY GRID={1}""")
        # no space before a section name
        assertNotValid(""" ENERGY GRID={3}""")
        assertParse("""ENERGY GRID={3}""", ("ENERGY", {"GRID": ar(3)}))

        grammar = input_parameters_def.grammar()
        assertParse("""ENERGY GRID={3}""", {"ENERGY": {"GRID": ar(3)}})
        assertParse(
            """ENERGY GRID={3}
                     NE={300}


          """,
            {"ENERGY": {"GRID": 3, "NE": 300}},
        )
        assertParse(
            """ENERGY
                     NE={300}
                     GRID={3}


SITES NL=2""",
            {"ENERGY": {"NE": 300, "GRID": 3}, "SITES": {"NL": 2}},
        )

        # custom values
        with generate_grammar():
            grammar = input_parameters_def["ENERGY"]._grammar_of_values()
        assertParse(
            """GRID={3}
                     NE={300}
                     """,
            {"GRID": ar(3), "NE": ar(300)},
        )

        assertParse(
            """GRID={3}
                   NE={300}
                   NXXX=5""",
            {"GRID": ar(3), "NE": ar(300), "NXXX": 5},
        )

        grammar = input_parameters_def.grammar()
        assertParse(
            """ENERGY GRID={3}
                     NE={300}
                     NXXX=5


SITES NL=2""",
            {"ENERGY": {"GRID": ar(3), "NE": ar(300), "NXXX": 5}, "SITES": {"NL": 2}},
        )

        assertParse(
            """ENERGY GRID={3}
                     NE={300}
                     NXXX

SITES NL=2

              """,
            {"ENERGY": {"GRID": ar(3), "NE": ar(300), "NXXX": True}, "SITES": {"NL": 2}},
        )

        # SITES do not start on the begin of the line, so it is not the start of the section
        assertParse(
            """ENERGY GRID={3}
                     NE={300}
                     NXXX

     SITES NL=2

              """,
            {"ENERGY": {"GRID": ar(3), "NE": ar(300), "NXXX": True, "SITES": True, "NL": 2}},
        )

        assertParse(
            """ENERGY GRID={3}
                     NE={300}
                     NXXX

  SITES NL=2

              """,
            {"ENERGY": {"GRID": ar(3), "NE": ar(300), "NXXX": True, "SITES": True, "NL": 2}},
        )

        assertNotValid("""ENERGY GRID={3}
                     NE={300}
                     NXXX


SITES NL=2

SITES NL=3
              """)

        # custom section
        assertParse(
            """ENERGY GRID={3}
                     NE={300}
                     NXXX

    #numbered_arrays

SITES NL=2

XSITES NR=3
              """,
            {"ENERGY": {"GRID": ar(3), "NE": ar(300), "NXXX": True}, "SITES": {"NL": 2}, "XSITES": {"NR": 3}},
        )

        # multiline custom section
        assertParse(
            """ENERGY GRID={3}
                     NE={300}
                     NXXX


SITES NL=2

XSITES NR=3 NF=1
                     NZ=5.5
              """,
            {
                "ENERGY": {"GRID": ar(3), "NE": ar(300), "NXXX": True},
                "SITES": {"NL": 2},
                "XSITES": {"NR": 3, "NF": 1, "NZ": 5.5},
            },
        )

        # a custom section with a flag
        assertParse(
            """ENERGY GRID={3}
                     NXXX
                     NE={300}
                     NZZZ=4


SITES NL=2

XSITES NR=3 FLAG
                     FLOAT=3.5
              """,
            {
                "ENERGY": {"GRID": ar(3), "NXXX": True, "NE": ar(300), "NZZZ": 4},
                "SITES": {"NL": 2},
                "XSITES": {"NR": 3, "FLAG": True, "FLOAT": 3.5},
            },
        )

        ips = input_parameters_def.read_from_file(
            io.StringIO("""ENERGY GRID={3}
                     NE={300,200}
                     NXXX


SITES NL=2

XSITES NR=3 FLAG
                     FLOAT=3.5
              """)
        )
        self.assertTrue(isinstance(ips, input_parameters.InputParameters))
        self.assertTrue(isinstance(ips["ENERGY"], Section))
        self.assertTrue(isinstance(ips["ENERGY"]["NE"], Option))
        self.assertTrue(isinstance(ips["ENERGY"]["NXXX"], CustomOption))
        self.assertTrue(isinstance(ips.ENERGY.NXXX, CustomOption))
        ips.ENERGY.NXXX.set(5)
        self.assertEqual(ips.ENERGY.NXXX(), 5)
        ips.ENERGY.NXXX.remove()
        self.assertFalse(hasattr(ips.ENERGY, "NXXX"))

        self.assertEqual(ips["ENERGY"].NE(), np.array((300, 200)))
        cps = ips.copy(copy_values=True)
        cps.ENERGY.NE[0] = 100
        self.assertEqual(ips["ENERGY"].NE(), np.array((300, 200)))
        dps = cps.copy()
        dps.ENERGY.NE[0] = 150
        self.assertEqual(cps["ENERGY"].NE(), np.array((150, 200)))
        cps.ENERGY.NE[0] = 180
        self.assertEqual(dps["ENERGY"].NE(), np.array((180, 200)))

        self.assertEqual(ips["SITES"].NL(), 2)
        self.assertEqual(ips.find("NL").get_path(), "SITES.NL")
        ips.find("NL").set(3)
        self.assertEqual(ips["SITES"].NL(), 3)
        self.assertTrue(isinstance(ips["XSITES"], CustomSection))

        output = io.StringIO()
        ips.save_to_file(output)
        output.seek(0)
        ips2 = input_parameters_def.read_from_file(output)
        self.assertEqual(str(ips.as_dict()), str(ips2.as_dict()))

        with pytest.raises(DataValidityError):
            ips.ENERGY.NE = "sss"
        ips.ENERGY.NE.set_dangerous("sss")
        output = io.StringIO()
        ips.save_to_file(output)
        output.seek(0)
        with pytest.raises(pp.ParseBaseException):
            ips2 = input_parameters_def.read_from_file(output)
        output.seek(0)
        ips2 = input_parameters_def.read_from_file(output, allow_dangerous=True)
        self.assertEqual("sss", ips.ENERGY.NE())
        self.assertEqual(str(ips.as_dict()), str(ips2.as_dict()))

    #

    def test_set_values(self):
        input_parameters_def = cd.InputParametersDefinition.definition_from_dict(
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
        ips = input_parameters_def.create_object()

        ips.set({"ENERGY.Ime": 5.0, "NE": 10, "ENERGY": 7.0})
        assert ips.find("energy") is ips.ENERGY.ENERGY
        assert ips.find("energy.ime") is ips.ENERGY.Ime
        assert ips.get_member("ENERGY") is ips.ENERGY
        assert ips.get_member("energy.ime") is ips.ENERGY.Ime
        assert ips.get_member(
            "ENERGY",
            accept=lambda member: member._definition.accept_value(6.0),
        ) is ips.ENERGY.ENERGY
        with pytest.raises(KeyError, match="MISSING"):
            ips.find("MISSING")
        assert (
            next(
                ips.get_members(
                    "ENERGY",
                    accept=lambda member: member._definition.accept_value(6.0),
                )
            )
            is ips.ENERGY.ENERGY
        )
        self.assertEqual(ips.ENERGY.Ime(), 5.0)
        self.assertEqual(ips.ENERGY.NE(), 10)
        self.assertEqual(ips.ENERGY.ENERGY(), 7.0)
        ips.set({"ENERGY.ENERGY": 4.0})
        self.assertEqual(ips.ENERGY.ENERGY(), 4.0)
        ips.set("ENERGY.ENERGY", 5.0)
        self.assertEqual(ips.ENERGY.ENERGY(), 5.0)
        ips.set("ENERGY", 6.0)
        self.assertEqual(ips.ENERGY.ENERGY(), 6.0)
        ips.set("ENERGY.ENERGYY", 7.0, unknown="ignore")
        self.assertFalse("ENERGYY" in ips.ENERGY)
        self.assertRaises(KeyError, lambda: ips.set("ENERGY.ENERGYY", 7.0, unknown="fail"))
        with pytest.raises(KeyError, match="dotted path"):
            ips.set("ENERGY.ENERGYY", 7.0, unknown="add")
        ips.ENERGY.set("ENERGYY", 7.0, unknown="add")
        self.assertEqual(ips.ENERGY.ENERGYY(), 7.0)

    #

    def test_repeated_value(self):

        t = lambda x: re.sub(r"\s+", " ", x).strip()  # noqa E731
        assertParse = self.assertParse
        assertNotValid = partial(self.assertNotValid, grammar=lambda: grammar)

        ipd = cd.InputParametersDefinition.definition_from_dict(
            {"ENERGY": [V("A", 1), V("B", 2, is_repeated=True), V("C", 3)]}
        )

        data = """ENERGY A=3
    B=2
    C=77"""
        assertParse(data, {"ENERGY": {"A": 3, "B": [2], "C": 77}}, ipd.grammar())
        ip = ipd.read_from_string(data)
        self.assertEqual({"ENERGY": {"A": 3, "B": [2], "C": 77}}, ip.to_dict())
        self.assertEqual("ENERGY A=3 B=2 C=77", t(ip.ENERGY.to_string()))

        data = """ENERGY A=3
    B=2
    B=5
    B=8
    C=77"""
        grammar = ipd.grammar()
        assertParse(data, {"ENERGY": {"A": 3, "B": [2, 5, 8], "C": 77}}, grammar)
        ip = ipd.read_from_string(data)
        self.assertEqual({"ENERGY": {"A": 3, "B": ar([2, 5, 8]), "C": 77}}, ip.to_dict())
        self.assertEqual("ENERGY A=3 B=2 B=5 B=8 C=77", t(ip.ENERGY.to_string()))
        ipd["ENERGY"].force_order = True
        assertParse(data, {"ENERGY": {"A": 3, "B": [2, 5, 8], "C": 77}}, grammar)
        ip = ipd.read_from_string(data)
        self.assertEqual({"ENERGY": {"A": 3, "B": ar([2, 5, 8]), "C": 77}}, ip.to_dict())
        self.assertEqual("ENERGY A=3 B=2 B=5 B=8 C=77", t(ip.ENERGY.to_string()))

        ipd["ENERGY"]["B"].is_repeated = V.Repeated.NUMBERED
        grammar = ipd.grammar()
        data = """ENERGY A=3
    B1=2
    B2=5
    B3=8
    C=77"""
        pp.ParserElement.verbose_stacktrace = True
        assertParse(data, {"ENERGY": {"A": 3, "B": [2, 5, 8], "C": 77}}, grammar)
        self.assertEqual(t(data), t(ip.ENERGY.to_string()))
        data = """ENERGY A=3
    B1=2
    B3=5
    B2=8
    C=77"""
        assertParse(data, {"ENERGY": {"A": 3, "B": [2, 8, 5], "C": 77}}, grammar)

        ipd = cd.InputParametersDefinition.definition_from_dict(
            {"ENERGY": [V("A", 1), V("B", gt.Array(int, length=2), is_repeated=True), V("C", 3)]}
        )
        grammar = ipd.grammar()

        data = """ENERGY A=3
    B=2 3
    B=5 9
    B=8 16
    C=77"""
        assertParse(data, {"ENERGY": {"A": 3, "B": [ar([2, 3]), ar([5, 9]), ar([8, 16])], "C": 77}}, grammar)
        ip = ipd.read_from_string(data)
        self.assertEqual(t(data), t(ip.ENERGY.to_string()))

        ipd["ENERGY"]["B"].is_repeated = V.Repeated.NUMBERED
        grammar = ipd.grammar()
        assertNotValid(data)
        data = """ENERGY A=3
    B1=2 3
    B2=5 9
    B3=8 16
    C=77"""
        assertParse(data, {"ENERGY": {"A": 3, "B": [ar([2, 3]), ar([5, 9]), ar([8, 16])], "C": 77}}, grammar)
        ip = ipd.read_from_string(data)
        self.assertEqual(t(data), t(ip.ENERGY.to_string()))

    def test_numbered_if(self):
        compact = lambda value: re.sub(r"\s+", " ", value).strip()  # noqa E731
        ipd = cd.InputParametersDefinition.definition_from_dict(
            {
                "ENERGY": [
                    V("MODE", 0),
                    V(
                        "KA",
                        int,
                        is_repeated=V.Repeated.NUMBERED_IF(lambda option: option._container["MODE"]() == 1),
                    ),
                ]
            }
        )

        # Parsing accepts both spellings independent of the current mode and
        # normalizes the spelling when the complete configuration is written.
        ip = ipd.read_from_string("ENERGY MODE=0 KA1=4")
        assert np.array_equal(ip.ENERGY.KA(), [4])
        assert compact(ip.to_string()) == "ENERGY MODE=0 KA=4"

        ip = ipd.read_from_string("ENERGY KA=4 MODE=1")
        assert np.array_equal(ip.ENERGY.KA(), [4])
        assert compact(ip.to_string()) == "ENERGY MODE=1 KA1=4"

        ip = ipd.read_from_string("ENERGY MODE=1 KA1=4 KA2=5")
        assert np.array_equal(ip.ENERGY.KA(), [4, 5])
        assert compact(ip.to_string()) == "ENERGY MODE=1 KA1=4 KA2=5"

        with pytest.warns(DataValidityError, match="can contain only one value"):
            unnumbered = ipd.read_from_string("ENERGY KA1=4 KA2=5 MODE=0")
        assert np.array_equal(unnumbered.ENERGY.KA(), [4, 5])
        with pytest.raises(DataValidityError, match="can contain only one value"):
            unnumbered.to_string(validate=False)

        # A change that would make the stored value invalid is rejected before
        # modifying the configuration.
        with pytest.raises(DataValidityError, match="can contain only one value"):
            ip.ENERGY.MODE = 0
        assert ip.ENERGY.MODE() == 1
        assert compact(ip.to_string()) == "ENERGY MODE=1 KA1=4 KA2=5"

        ip = ipd.create_object()
        # Setting an overlong value while it is unnumbered is rejected.
        with pytest.raises(DataValidityError, match="can contain only one value"):
            ip.ENERGY.KA = [7, 8]
        assert ip.ENERGY.KA() is None

        # The two spellings denote the same first item, not two different ones.
        with pytest.raises(pp.ParseBaseException):
            ipd.read_from_string("ENERGY MODE=0 KA=4 KA1=5")

        copied = ipd.copy()
        ip = copied.read_from_string("ENERGY MODE=1 KA=3")
        assert compact(ip.to_string()) == "ENERGY MODE=1 KA1=3"

        # KA-like vector items retain the outer repetition dimension as well.
        vector_ipd = cd.InputParametersDefinition.definition_from_dict(
            {
                "ENERGY": [
                    V("MODE", 0),
                    V(
                        "KA",
                        gt.SetOf(int, length=3),
                        is_repeated=V.Repeated.NUMBERED_IF(lambda option: option._container["MODE"]() == 1),
                    ),
                ]
            }
        )
        vector = vector_ipd.read_from_string("ENERGY MODE=0 KA={1,2,3}")
        assert vector.ENERGY.KA().shape == (1, 3)
        assert compact(vector.to_string()) == "ENERGY MODE=0 KA={1,2,3}"
        vector.ENERGY.MODE = 1
        assert compact(vector.to_string()) == "ENERGY MODE=1 KA1={1,2,3}"

    #
    def test_sparse_numbered(self):
        input_parameters_def = cd.InputParametersDefinition.definition_from_dict(
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
            grammar = input_parameters_def._grammar_of_values()

        assertParse = partial(self.assertParse, grammar=lambda: grammar)
        assertNotValid = partial(self.assertNotValid, grammar=lambda: grammar)

        assertNotValid("ENERGY NE=1 Ime=0.5 Ime=2.0")
        assertParse("ENERGY NE=1 Ime=0.5", {"ENERGY": {"NE": ar(1), "Ime": 0.5}})
        input_parameters_def.sections["ENERGY"]["Ime"].is_repeated = V.Repeated["DICT"]

        with generate_grammar():
            grammar = input_parameters_def._grammar_of_values()
        assertNotValid("ENERGY NE=1 Ime=0.5")
        assertNotValid("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8")
        assertParse("ENERGY NE=1 Ime1=0.4 Ime5=0.8", {"ENERGY": {"NE": ar(1), "Ime": {1: 0.4, 5: 0.8}}})

        input_parameters_def.sections["ENERGY"]["Ime"].is_repeated = V.Repeated["DEFAULTDICT"]
        with generate_grammar():
            grammar = input_parameters_def._grammar_of_values()

        assertParse("ENERGY NE=1 Ime1=0.4 Ime5=0.8", {"ENERGY": {"NE": ar(1), "Ime": {1: 0.4, 5: 0.8}}})
        assertNotValid("ENERGY NE=1 Ime=0.5 Ime=2.0")
        assertParse(
            "ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8", {"ENERGY": {"NE": ar(1), "Ime": {"def": 0.5, 1: 0.4, 5: 0.8}}}
        )

        ip = input_parameters_def.read_from_file(io.StringIO("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8"))
        self.assertEqual(0.5, ip.ENERGY.Ime())
        self.assertEqual(0.8, ip.ENERGY.Ime[5])
        self.assertEqual([0.4, 0.5, 0.5, 0.5, 0.8], ip.ENERGY.Ime[:])
        self.assertEqual({"def": 0.5, 1: 0.4, 5: 0.8}, ip.ENERGY.Ime.as_dict())
        self.assertEqual({"def": 0.5, 1: 0.4, 5: 0.8}, ip.ENERGY.Ime(all_values=True))
        ip.ENERGY.Ime[1] = 0.7
        ip.ENERGY.Ime[9] = 0.9
        self.assertEqual(
            "ENERGY GRID={3} NE={1} Ime=0.5 Ime1=0.7 Ime5=0.8 Ime9=0.9",
            re.sub(r"\s+", " ", ip.ENERGY.to_string()).strip(),
        )
        ip.ENERGY.Ime[[1, 4, 9]] = 0.2
        self.assertEqual([0.2, 0.8, 0.2], ip.ENERGY.Ime[1, 5, 9])
        ip.ENERGY.Ime[[1, 4, 9]] = 0.3
        self.assertEqual([0.3, 0.8, 0.3], ip.ENERGY.Ime[[1, 5, 9]])
        ip.ENERGY.Ime[5:9] = 0.2
        self.assertEqual([0.2, 0.2, 0.2, 0.2, 0.3], ip.ENERGY.Ime[5:10])
        self.assertRaises(KeyError, lambda: ip.ENERGY.Ime[5.0])

        def e():
            ip.ENERGY.Ime[5.0] = 1

        self.assertRaises(KeyError, e)
        self.assertRaises(KeyError, lambda: ip.ENERGY.Ime["5"])

        def e():
            ip.ENERGY.Ime["5"] = 1

        self.assertRaises(KeyError, e)
        self.assertEqual([0.3, 0.2, 0.3], ip.ENERGY.Ime[[1, 5, 9]])
        ip.ENERGY.Ime = {3: 5.0, "def": 3.0}
        self.assertEqual([3.0, 3.0, 5.0], ip.ENERGY.Ime[:])
        ip.ENERGY.Ime[2:5] = [2.0, 3.0, 8.0]
        self.assertEqual([3.0, 2.0, 3.0, 8.0], ip.ENERGY.Ime[:])
        ip.ENERGY.Ime[1:3] = 7.0
        self.assertEqual([7.0, 7.0, 3.0, 8.0], ip.ENERGY.Ime[:])

        ip = input_parameters_def.read_from_file(io.StringIO("ENERGY NE=1 Ime=0.5 Ime1=0.4 Ime5=0.8"))
        with pytest.raises(DataValidityError):
            ip.ENERGY.Ime = "ss"
        ip.ENERGY.Ime.set_dangerous("uu")
        self.assertEqual(ip.ENERGY.Ime(), "uu")
        out = ip.ENERGY.to_string()
        self.assertEqual("ENERGY GRID={3} NE={1} Ime=uu Ime1=0.4 Ime5=0.8", re.sub(r"\s+", " ", out).strip())
        with pytest.raises(pp.ParseBaseException):
            input_parameters_def.read_from_file(io.StringIO(out))
        ip = input_parameters_def.read_from_file(io.StringIO(out), allow_dangerous=True)
        self.assertEqual(ip.ENERGY.Ime(all_values=True), {"def": "uu", 1: 0.4, 5: 0.8})
        ip.ENERGY.Ime = 1.0

        with pytest.raises(DataValidityError):
            ip.ENERGY.Ime[5] = "ss"
        ip.ENERGY.Ime.set_dangerous("yy", index=5)
        out = ip.ENERGY.to_string()
        self.assertEqual(
            "ENERGY GRID={3} NE={1} Ime=1.0 Ime1=0.4 Ime5=yy", re.sub(r"\s+", " ", ip.ENERGY.to_string()).strip()
        )
        with pytest.raises(pp.ParseBaseException):
            input_parameters_def.read_from_file(io.StringIO(out))
        ip = input_parameters_def.read_from_file(io.StringIO(out), allow_dangerous=True)
        self.assertEqual(ip.ENERGY.Ime(all_values=True), {"def": 1.0, 1: 0.4, 5: "yy"})

    def test_gather(self):
        assertParse = self.assertParse
        ipd = cd.InputParametersDefinition.definition_from_dict(
            {"ENERGY": [V("GRID", gt.SetOf(int, length=1), fixed_value=3), *gather(V("A", 1), V("B", 2)), V("C", 3)]}
        )
        assertParse("ENERGY GRID={3} A B=1 2 C=3", {"ENERGY": {"GRID": ar(3), "A": 1, "B": 2, "C": 3}}, ipd.grammar())
        out = ipd.read_from_string("ENERGY GRID={3} A B=1 2 C=3")
        # self.assertEqual("ENERGY GRID={3} A B=1 2 C=3 TASK INPUTPARAMETERSDEFINITION ", re.sub(r'[\s\t\n]+',' ', out.to_string()))
        self.assertEqual("ENERGY GRID={3} A B=1 2 C=3 ", re.sub(r"[\s\t\n]+", " ", out.to_string()))

    def test_numpy_array(self):
        assertParse = self.assertParse

        ipd = cd.InputParametersDefinition.definition_from_dict(
            {"ENERGY": [V("C", gt.NumpyArray(lines=3), name_in_grammar=False)]}
        )
        assertParse(
            """ENERGY
    1 1
    2 2
    3 3
    """,
            {"ENERGY": {"C": ar([1, 1, 2, 2, 3, 3]).reshape(3, 2)}},
            ipd.grammar(),
        )

        ipd = cd.InputParametersDefinition.definition_from_dict(
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

        g = ipd.grammar()
        assertParse(
            data,
            {
                "ENERGY": {
                    "A": 3,
                    "B": 2,
                    "C": ar([1.0, 1, 2, 2, 3, 3]).reshape(3, 2),
                    "D": ar([4.0, 5, 6, 5, 9, 8]).reshape((2, 3)),
                    "E": 77,
                }
            },
            g,
        )
        ipd["ENERGY"].force_order = True
        g = ipd.grammar()
        assertParse(
            data,
            {
                "ENERGY": {
                    "A": 3,
                    "B": 2,
                    "C": ar([1.0, 1, 2, 2, 3, 3]).reshape(3, 2),
                    "D": ar([4.0, 5, 6, 5, 9, 8]).reshape((2, 3)),
                    "E": 77,
                }
            },
            g,
        )

        ipd = cd.InputParametersDefinition.definition_from_dict(
            {
                "ENERGY": [
                    V("A", 1),
                    V("B", 2),
                    V("C", gt.NumpyArray(lines="A"), name_in_grammar=False),
                    V("D", gt.NumpyArray(lines="B"), name_in_grammar=False),
                    V("E", 3),
                ]
            }
        )
        assertParse(
            """ENERGY A=3 B=2
    1 1
    2 2
    3 3
    4 5 6
    5 9 8
    E=77""",
            {
                "ENERGY": {
                    "A": 3,
                    "B": 2,
                    "C": ar([1.0, 1, 2, 2, 3, 3]).reshape(3, 2),
                    "D": ar([4.0, 5, 6, 5, 9, 8]).reshape((2, 3)),
                    "E": 77,
                }
            },
            ipd.grammar(),
        )

        def test(na, val, str):
            self.assertEqual(na.string(val), str)
            with generate_grammar():
                self.assertEqual(np.asarray(val), na.parse(str))

        test(gt.NumpyArray(item_format="%2.0f"), [[1, 2, 3], [4, 5, 6]], " 1  2  3\n 4  5  6")
        test(gt.NumpyArray(indented=2, item_format="%2.0f"), [[1, 2, 3], [4, 5, 6]], "   1  2  3\n   4  5  6")
        test(
            gt.NumpyArray(indented=(8, 2), item_format="%2.0f"),
            [[1, 2, 3, 4, 5, 6], [4, 5, 6, 7, 8, 9]],
            " 1  2  3\n    4  5\n    6\n 4  5  6\n    7  8\n    9",
        )
        test(
            gt.NumpyArray(indented=(8, 2), no_newline_at_end=False, item_format="%2.0f"),
            [[1, 2, 3, 4, 5, 6], [4, 5, 6, 7, 8, 9]],
            " 1  2  3\n    4  5\n    6\n 4  5  6\n    7  8\n    9\n",
        )
        test(
            gt.NumpyArray(line_length=8, item_format="%2.0f", shape=(2, -1)),
            [[1, 2, 3, 4, 5, 6], [4, 5, 6, 7, 8, 9]],
            " 1  2  3\n 4  5  6\n 4  5  6\n 7  8  9",
        )

    def test_copy(self):
        ipd = cd.InputParametersDefinition.definition_from_dict(
            {"ENERGY": [V("A", 1), V("B", 2, is_repeated=True), V("C", 3)]}
        )
        ipd2 = ipd.copy()
        ipd2["ENERGY"]["A"].default_value = 5
        assert ipd["ENERGY"]["A"].default_value == 1

    def test_length(self):
        t = lambda x: re.sub(r"\s+", " ", x).strip()  # noqa E731

        ipd = cd.InputParametersDefinition.definition_from_dict(
            {
                "ENERGY": [
                    V("A", 1),
                    V("NK", Length("K1", "K2")),
                    V("K1", gt.Array(int), is_optional=True),
                    V("K2", gt.Array(int), is_optional=True),
                ]
            }
        )
        ip = ipd.create_object()
        ip.ENERGY.set(K1=[1, 2, 3], K2=[2, 2, 3])
        assert t(ip.to_string()) == "ENERGY A=1 NK=3 K1=1 2 3 K2=2 2 3"
        ip2 = ipd.read_from_string(ip.to_string())
        assert ip2.ENERGY.NK() == 3

        with pytest.raises(DataValidityError):
            ip2.ENERGY.NK = 7

        with pytest.raises(DataValidityError):
            ip.ENERGY.K2 = [2, 2]
        assert np.array_equal(ip.ENERGY.K2(), [2, 2, 3])

        with pytest.warns(DataValidityError):
            invalid = ipd.read_from_string("ENERGY A=1 NK=3 K1=1 2 3 K2=2 2")
        with pytest.raises(DataValidityError):
            invalid.to_string(validate="save")

        ipd = cd.InputParametersDefinition.definition_from_dict(
            {
                "ENERGY": [
                    V("A", 1),
                    V("NK", Length("K1", "K2", default_values=(1, 2))),
                    V("K1", gt.Array(int), is_optional=True),
                    V("K2", gt.Array(int), is_optional=True),
                ]
            }
        )
        ip = ipd.create_object()
        ip.ENERGY.NK = 3
        assert (ip.ENERGY.K1() == [1, 1, 1]).all()
        assert (ip.ENERGY.K2() == [2, 2, 2]).all()
        ip.ENERGY.NK = 5
        assert (ip.ENERGY.K1() == [1, 1, 1, 1, 1]).all()
        ip.ENERGY.NK = 2
        assert (ip.ENERGY.K1() == [1, 1]).all()
        ip.ENERGY.NK = 2

        section = ip.ENERGY
        old_k1 = section.K1().copy()
        old_k2 = section.K2().copy()
        with pytest.raises(DataValidityError):
            section.set(K1=[1, 1, 1], K2=[2, 2])
        assert (section.K1() == old_k1).all()
        assert (section.K2() == old_k2).all()

        ipd = cd.InputParametersDefinition.definition_from_dict(
            {
                "ENERGY": [
                    V("A", 1),
                    V("NK", Length("K1", default_values=[1, 2, 3])),
                    V("K1", gt.Array(int, length=3), is_repeated="REPEATED", is_optional=True),
                ]
            }
        )
        ip = ipd.create_object()
        ip.ENERGY.NK = 2
        assert (ip.ENERGY.K1() == [[1, 2, 3], [1, 2, 3]]).all()
        ip.ENERGY.K1[1] = [3, 3, 3]
        ip.ENERGY.NK = 3
        assert (ip.ENERGY.K1() == [[1, 2, 3], [3, 3, 3], [1, 2, 3]]).all()

    def test_standard_container_access_during_transactions(self):
        ipd = cd.InputParametersDefinition.definition_from_dict({"CONTROL": [V("VALUE", 1)]})
        section = ipd.create_object().CONTROL
        definition = ipd["CONTROL"]
        assert not definition.is_option
        assert definition["VALUE"].is_option

        assert "VALUE" in section
        assert section["VALUE"]() == 1
        assert "EXTRA" not in section
        section["VALUE"].set_dangerous("not-an-integer")
        assert section["VALUE"].is_dangerous()

        dangerous = DangerousValue("unchecked")
        with ConfigurationTransaction(section) as transaction:
            with transaction.savepoint() as savepoint:
                section.stage_value(
                    transaction, "EXTRA", dangerous, unknown="add"
                )
                assert "VALUE" in section
                assert "EXTRA" in section
                assert "MISSING" not in section
                assert section["EXTRA"].is_dangerous()
                savepoint.rollback()
        assert "EXTRA" not in section

        generated_ipd = cd.InputParametersDefinition.definition_from_dict(
            {"CONTROL": [V("VALUES", gt.Array(int)), V("NVALUES", Length("VALUES"))]}
        )
        generated_section = generated_ipd.create_object().CONTROL
        with ConfigurationTransaction(generated_section) as transaction:
            with transaction.savepoint() as savepoint:
                generated_section["VALUES"].stage(transaction, [1, 2, 3])
                assert generated_section.NVALUES() == 3
                savepoint.rollback()

        copied_generated = generated_ipd.copy().create_object()
        copied_generated.CONTROL.VALUES = [1, 2]
        assert copied_generated.CONTROL.NVALUES() == 2

    def test_generated_option_uses_its_runtime_section(self):
        def generated_getter(section, _key=None):
            return section["SOURCE"]()

        def generated_setter(section, value, _key=None):
            transaction = ConfigurationTransaction.current(section)
            section["SOURCE"].stage(transaction, value)

        definition = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("SOURCE", 1),
                    GeneratedValueDefinition(
                        "GENERATED", generated_getter, generated_setter
                    ),
                ]
            }
        )
        parameters = definition.create_object()
        assert parameters.CONTROL.GENERATED() == 1
        parameters.CONTROL.GENERATED.set(2)
        assert parameters.CONTROL.SOURCE() == 2

    def test_parse_validation_sees_keyword_assignments(self):
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

    def test_adding_empty_custom_member_is_validated(self):
        def reject_extra(root, _values, _why):
            if "EXTRA" in root:
                return DataValidityError("EXTRA is not allowed")

        definition = cd.InputParametersDefinition([], validators=reject_extra)
        parameters = definition.create_object()

        with pytest.raises(DataValidityError, match="EXTRA is not allowed"):
            parameters.add("EXTRA")
        assert "EXTRA" not in parameters

    def test_conditional_length(self):
        mode_is_one = lambda definition, section: section["MODE"]() == 1  # noqa E731
        ipd = cd.InputParametersDefinition.definition_from_dict(
            {
                "CONTROL": [
                    V("VALUES", gt.Array(int), is_optional=True),
                    V("NVALUES", Length("VALUES"), condition=mode_is_one),
                    V("MODE", 1),
                ]
            }
        )

        active = ipd.read_from_string("CONTROL VALUES=1 2 3 NVALUES=3 MODE=1")
        assert active.CONTROL.NVALUES() == 3
        assert not hasattr(active.CONTROL, "_parsed_values")
        assert "NVALUES=3" in active.to_string()

        with pytest.warns(DataValidityError):
            ipd.read_from_string("CONTROL VALUES=1 2 3 NVALUES=2 MODE=1")

        with pytest.raises(pp.ParseBaseException):
            ipd.read_from_string("CONTROL VALUES=1 2 3 NVALUES=3 MODE=2")

        with pytest.raises(pp.ParseBaseException):
            ipd.read_from_string("CONTROL VALUES=1 2 3 NVALUES=2 MODE=2")

        with pytest.raises(pp.ParseBaseException):
            ipd.read_from_string("CONTROL VALUES=1 2 3 MODE=1")

        with warnings.catch_warnings():
            warnings.simplefilter("error", DataValidityError)
            inactive = ipd.read_from_string("CONTROL VALUES=1 2 3 MODE=2")
        assert inactive.CONTROL["NVALUES"]() == 3
        assert "NVALUES" not in inactive.to_string()

    def test_cross_section_conditional_length(self):
        def mode_is_one(definition, section):
            return section._get_root_container()["SECTION_B"]["MODE"]() == 1

        ipd = cd.InputParametersDefinition.definition_from_dict(
            {
                "SECTION_A": [
                    V("VALUES", gt.Array(int), is_optional=True),
                    V("NVALUES", Length("VALUES"), condition=mode_is_one),
                ],
                "SECTION_B": [V("MODE", 1)],
            }
        )

        active = ipd.read_from_string("SECTION_A VALUES=1 2 3 NVALUES=3\nSECTION_B MODE=1")
        assert active.SECTION_A.NVALUES() == 3
        assert "NVALUES=3" in active.to_string()

        with pytest.warns(DataValidityError):
            ipd.read_from_string("SECTION_A VALUES=1 2 3 NVALUES=2\nSECTION_B MODE=1")

        with pytest.raises(pp.ParseBaseException):
            ipd.read_from_string("SECTION_A VALUES=1 2 3 NVALUES=3\nSECTION_B MODE=2")

        with pytest.raises(pp.ParseBaseException):
            ipd.read_from_string("SECTION_A VALUES=1 2 3 NVALUES=2\nSECTION_B MODE=2")

        with pytest.raises(pp.ParseBaseException):
            ipd.read_from_string("SECTION_A VALUES=1 2 3\nSECTION_B MODE=1")

        with warnings.catch_warnings():
            warnings.simplefilter("error", DataValidityError)
            inactive = ipd.read_from_string("SECTION_A VALUES=1 2 3\nSECTION_B MODE=2")
        assert inactive.SECTION_A["NVALUES"]() == 3
        assert "NVALUES" not in inactive.to_string()

    #
    def test_switch(self):
        assertParse = self.assertParse
        assertNotValid = self.assertNotValid
        ipd = cd.InputParametersDefinition.definition_from_dict(
            {"ENERGY": [V("A", 1), *switch("A", {1: [V("B", 2, is_optional=True)], 2: V("C", 3, is_optional=True)})]}
        )

        ipd.custom_class = None
        ipd["ENERGY"].custom_class = None
        ipd["ENERGY"].force_order = True
        grammar = ipd.grammar()
        assertParse("ENERGY A=2 C=4", {"ENERGY": {"A": 2, "C": 4}}, grammar)
        assertParse("ENERGY A=2 C=4", {"ENERGY": {"A": 2, "C": 4}}, grammar)
        out = ipd.read_from_string("ENERGY A=2 C=4")
        self.assertEqual(out["ENERGY"].to_string(), "ENERGY\n\tA=2\n\tC=4\n")
        self.assertEqual(out["ENERGY"].to_string(), "ENERGY\n\tA=2\n\tC=4\n")
        assertParse("ENERGY A=1 B=2", {"ENERGY": {"A": 1, "B": 2}}, grammar)
        out = ipd.read_from_string("ENERGY A=1 B=4")
        self.assertEqual(out["ENERGY"].to_string(), "ENERGY\n\tA=1\n\tB=4\n")

        assertNotValid("ENERGY A=2 B=2", grammar)
        assertNotValid("ENERGY A=1 C=1", grammar)
        assertNotValid("ENERGY A=1 B=2 C=1", grammar)
        assertNotValid("ENERGY A=2 B=2 C=1", grammar)
        assertNotValid("ENERGY A=3 C=3", grammar)
        assertParse("ENERGY A=2", {"ENERGY": {"A": 2}}, grammar)
        ipd["ENERGY"]["C"].is_optional = False
        grammar = ipd.grammar()
        assertNotValid("ENERGY A=2", grammar)
        ipd = cd.InputParametersDefinition.definition_from_dict(
            {
                "ENERGY": [
                    V("A", 1),
                    *switch(
                        "A",
                        {
                            1: [V("B", 2, is_optional=True), V("D", 3, is_optional=True), V("C", 1, is_optional=True)],
                            2: [V("C", 3, is_optional=True), V("E", 1, is_optional=True), V("B", 2, is_optional=True)],
                        },
                    ),
                ]
            }
        )
        ipd["ENERGY"].force_order = True
        grammar = ipd.grammar()
        assertParse("ENERGY A=2 C=4 E=2 B=6", {"ENERGY": {"A": 2, "C": 4, "E": 2, "B": 6}}, grammar)
        assertParse("ENERGY A=1 B=5 D=3 C=5", {"ENERGY": {"A": 1, "B": 5, "D": 3, "C": 5}}, grammar)
        assertParse("ENERGY A=2 B=7", {"ENERGY": {"A": 2, "B": 7}}, grammar)
        assertParse("ENERGY A=1 B=7", {"ENERGY": {"A": 1, "B": 7}}, grammar)
        assertParse("ENERGY A=1 B=7 D=6", {"ENERGY": {"A": 1, "B": 7, "D": 6}}, grammar)
        assertNotValid("ENERGY A=1 D=2 B=3", grammar)
        assertNotValid("ENERGY A=1 B=5 D=3 C=5 B=6", grammar)
        assertNotValid("ENERGY A=1 B=5 D=3 C=5 E=6", grammar)
        assertNotValid("ENERGY A=2 B=5 C=4 E=3", grammar)

        ipd["ENERGY"]["B"].is_optional = False
        ipd["ENERGY"]["C"].is_optional = False
        ipd["ENERGY"]["D"].is_optional = False
        ipd["ENERGY"]["E"].is_optional = False
        grammar = ipd.grammar()
        assertNotValid("ENERGY A=2 B=2", grammar)
        assertNotValid("ENERGY A=1 B=2", grammar)
        assertNotValid("ENERGY A=1 B=2 D=3", grammar)
        assertNotValid("ENERGY A=1 D=2 B=3", grammar)
        assertNotValid("ENERGY A=1 B=5 D=3 C=5 B=6", grammar)
        assertNotValid("ENERGY A=1 B=5 D=3 C=5 E=6", grammar)
        assertNotValid("ENERGY A=2 B=5 C=4 E=3", grammar)
        out = ipd.read_from_string("ENERGY A=2 C=4 E=2 B=6")
        self.assertEqual(out["ENERGY"].to_dict(), {"A": 2, "C": 4, "E": 2, "B": 6})
        self.assertEqual(out["ENERGY"].to_string(), "ENERGY\n\tA=2\n\tC=4\n\tE=2\n\tB=6\n")

import numpy as np
import pytest

from ...common import grammar_types as gt
from ...common.backward_compatibility import ExceptionGroup
from ...common.configuration_transaction import (
    ConfigurationTransaction,
    PostChangeHookError,
)
from ...common.options import CustomOption, Option
from ...common.warnings import DataValidityError, DataValidityErrors
from .. import input_parameters_definitions as cd
from ._test_callbacks import validator_error_if_three as _validator_error_if_three

V = cd.InputValueDefinition


def test_transaction_converts_each_value_once():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [V("VALUE", 1)]}
    )
    parameters = definition.create_object()
    option = parameters.CONTROL.VALUE
    value_definition = option._definition
    original_convert = value_definition.convert_value
    converted = []

    def count_conversion(option, value, item=False):
        converted.append(value)
        return original_convert(option, value, item)

    value_definition.convert_value = count_conversion
    option.set(2)

    assert converted == [2]
    assert option() == 2


def test_transaction_commit_uses_the_validated_values():
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


def test_numbered_item_assignment_converts_new_value_once():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [V("VALUES", int, is_repeated=V.Repeated.DICT)]}
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


def test_post_change_hook_errors_are_aggregated_after_commit():
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


def test_validation_and_post_change_hook_errors_are_reported_together():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [V("VALUE", 1, validators=_validator_error_if_three)]}
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
    assert isinstance(caught.value.subgroup(DataValidityError), DataValidityErrors)
    assert not isinstance(caught.value.subgroup(PostChangeHookError), DataValidityErrors)
    assert parameters.CONTROL.VALUE() == 3


def test_failed_indexed_array_hook_does_not_roll_back_value():
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


def test_indexed_array_materializes_default_without_modifying_definition():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [V("VALUES", gt.Array(int), [1, 2])]}
    )
    first = definition.create_object()
    second = definition.create_object()

    first.CONTROL.VALUES[0] = 9

    assert np.array_equal(first.CONTROL.VALUES(), [9, 2])
    assert np.array_equal(second.CONTROL.VALUES(), [1, 2])
    assert np.array_equal(definition["CONTROL"]["VALUES"].default_value, [1, 2])


def test_indexed_array_materializes_result():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [V("VALUES", gt.Array(int), [1, 2])]}
    )
    option = definition.create_object().CONTROL.VALUES
    option.result = np.array([3, 4])

    option[0] = 8

    assert np.array_equal(option(), [8, 4])
    assert option._value is not None
    assert not hasattr(option, "_result")


def test_transaction_context_manager_commits_or_rolls_back():
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


def test_direct_stage_uses_transaction_validation_policy():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [V("VALUE", int, 1)]}
    )
    parameters = definition.create_object()
    option = parameters.CONTROL.VALUE

    with pytest.raises(DataValidityError):
        with ConfigurationTransaction(parameters) as transaction:
            option.stage(transaction, "not an integer")

    assert option() == 1


def test_transaction_mutation_requires_context():
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


def test_transaction_use_joins_active_transaction():
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


def test_caught_nested_set_rolls_back_only_nested_changes():
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


def test_uncaught_nested_set_rolls_back_the_whole_transaction():
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


def test_transaction_savepoint_rolls_back_only_its_changes():
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


def test_transaction_savepoint_rolls_back_repeated_option_change():
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


def test_option_hook_is_coalesced_and_runs_outside_transaction():
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


def test_option_stage_registers_its_change():
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


def test_unknown_member_is_converted_before_transaction_commit():
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


def test_add_with_invalid_value_rolls_back_custom_member():
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


def test_discard_and_ignore_unconvertible_staged_values():
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


def test_discarded_invalid_value_does_not_add_custom_member():
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

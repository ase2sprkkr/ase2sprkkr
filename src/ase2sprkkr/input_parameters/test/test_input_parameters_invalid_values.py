import io

import numpy as np
import pytest

from ...common import grammar_types as gt
from ...common.configuration_transaction import (
    ConfigurationTransaction,
    PostChangeHookError,
)
from ...common.generated_configuration_definitions import (
    GeneratedValueDefinition,
    Length,
    NumpyViewDefinition,
)
from ...common.options import CustomOption
from ...common.warnings import DataValidityError, DataValidityErrors
from .. import input_parameters_definitions as cd
from ._test_callbacks import (
    generated_source_length_validator,
    get_cross_section_generated,
    reject_cross_section_two,
    reject_negative_raw_data,
    repeated_validator_reject_two,
    required_in_mode_one,
    set_cross_section_generated,
    validator_error_if_three,
)

V = cd.InputValueDefinition


def test_typed_and_unconvertible_invalid_values_are_distinguished():
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


def test_retention_and_reporting_are_independent():
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
        {"CONTROL": [V("VALUE", 1, validators=validator_error_if_three)]}
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


def test_strict_batch_collects_all_value_errors_before_raising():
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


def test_repeated_dictionary_collects_all_errors_atomically():
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


def test_repeated_values_use_invalid_value_policy():
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


def test_invalid_dotted_paths_are_atomic_and_simple_custom_paths_remain_supported():
    section = cd.InputSectionDefinition("CONTROL", [V("KNOWN", 1)])
    section.custom_class = CustomOption.factory(V, int)
    parameters = cd.InputParametersDefinition([section]).create_object()

    # An unknown dotted member must not be installed inside a section.
    with pytest.raises(KeyError, match="dotted path"):
        parameters.CONTROL.set({"EXTRA.CHILD": 1}, unknown="add")
    assert "EXTRA" not in parameters.CONTROL

    # The same rule holds for unknown root paths.
    with pytest.raises(KeyError, match="EXTRA.VALUE"):
        parameters.get_member("EXTRA.VALUE")
    assert list(parameters.get_members("EXTRA.VALUE")) == []
    with pytest.raises(KeyError, match="dotted path"):
        parameters.set("EXTRA.VALUE", 2, unknown="add")
    assert "EXTRA" not in parameters

    # A path beyond an existing option must not resolve or create a child.
    with pytest.raises(KeyError, match="CONTROL.KNOWN.CHILD"):
        parameters.get_member("CONTROL.KNOWN.CHILD")
    assert list(parameters.get_members("CONTROL.KNOWN.CHILD")) == []
    with pytest.raises(KeyError, match="CONTROL.KNOWN.CHILD"):
        parameters.set("CONTROL.KNOWN.CHILD", 2, unknown="add")
    assert parameters.CONTROL.KNOWN() == 1

    # A non-dotted custom member remains supported after the rejected attempts.
    parameters.set({"EXTRA": {"VALUE": 2}}, unknown="add")
    assert parameters.EXTRA.VALUE() == 2


def test_failed_read_rolls_back_clear_and_parsed_values_together():
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


def test_unknown_add_does_not_replace_existing_member():
    definition = cd.InputParametersDefinition(
        [cd.InputSectionDefinition("CONTROL", [V("VALUE", 1)])]
    )
    parameters = definition.create_object()
    original = parameters.CONTROL

    with pytest.raises(TypeError, match="already in"):
        parameters.set({"CONTROL": 2}, unknown="add")

    assert parameters.CONTROL is original


def test_transaction_uses_standard_container_access():
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


def test_post_change_hook_cannot_mutate_its_configuration():
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


def test_rolled_back_savepoint_does_not_register_hook():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [V("FIRST", 1), V("SECOND", 2)]}
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


def test_clear_is_validated_and_rolled_back():
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

    # Verify set(None) on a different nonempty value, as in the old test.
    required.CONTROL.REQUIRED.set(2)
    with pytest.raises(DataValidityError, match="must have a value"):
        required.CONTROL.REQUIRED.set(None)
    assert required.CONTROL.REQUIRED() == 2


def test_container_clear_skips_generated_setters():
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
                GeneratedValueDefinition("GENERATED", get_generated, set_generated),
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


def test_repeated_clear_is_validated_and_rolled_back():
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
    assert parameters.ROWS.as_dict(False) is None
    parameters.ROWS.set([{"VALUE": 2}])

    with pytest.raises(DataValidityError, match="empty rows rejected"):
        parameters.ROWS.clear()

    assert len(parameters.ROWS) == 1
    assert parameters.ROWS[0].VALUE() == 2


def test_option_validation_uses_its_container_values():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "CONTROL": [
                V("MODE", 1),
                V("DEPENDENT", int, 2, is_required=required_in_mode_one),
            ]
        }
    )
    definition.create_object().CONTROL.DEPENDENT.validate("save")


def test_setting_selector_validates_dependent_required_values():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {
            "CONTROL": [
                V("MODE", 2),
                V(
                    "DEPENDENT",
                    int,
                    None,
                    is_optional=True,
                    is_required=required_in_mode_one,
                ),
            ]
        }
    )
    parameters = definition.create_object()

    with pytest.raises(DataValidityError, match="DEPENDENT"):
        parameters.CONTROL.MODE = 1

    assert parameters.CONTROL.MODE() == 2


def test_unconditional_required_values_are_save_time_completeness_checks():
    definition = cd.InputParametersDefinition.definition_from_dict(
        {"CONTROL": [V("REQUIRED", int, None), V("OTHER", 1)]}
    )
    parameters = definition.create_object()

    parameters.CONTROL.OTHER = 2
    assert parameters.CONTROL.OTHER() == 2
    with pytest.raises(DataValidityError, match="REQUIRED"):
        parameters.validate("save")


def test_disabled_section_is_not_validated():
    disabled = cd.InputSectionDefinition("DISABLED", [V("REQUIRED", int, None)])
    disabled.condition = lambda _definition, _container: False
    definition = cd.InputParametersDefinition(
        [
            cd.InputSectionDefinition("ACTIVE", [V("VALUE", 1)]),
            disabled,
        ]
    )

    definition.create_object().validate("save")


def test_copy_does_not_revalidate_existing_values():
    definition = cd.InputParametersDefinition(
        [
            cd.InputSectionDefinition(
                "CONTROL",
                [V("VALUE", 1, validators=validator_error_if_three)],
            )
        ]
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


def test_repeated_and_generated_sets_are_transactional():
    repeated = cd.InputParametersDefinition(
        [
            cd.InputSectionDefinition(
                "ROWS",
                [V("VALUE", 0)],
                is_repeated=True,
                validators=repeated_validator_reject_two,
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
    generated.validators = (generated_source_length_validator,)
    parameters = generated.create_object()
    parameters.CONTROL.NVALUES = 1
    with pytest.raises(DataValidityError, match="generated setter produced too many values"):
        parameters.CONTROL.NVALUES = 2
    assert np.array_equal(parameters.CONTROL.VALUES(), [1])


def test_dictionary_repeated_set_replaces_or_merges_values():
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
    parameters.ROWS.set({1: {"VALUE": 1, "OTHER": 10}, 2: {"VALUE": 2}})

    parameters.ROWS.set({1: {"VALUE": 3}})
    assert list(parameters.ROWS) == [1]
    assert parameters.ROWS[1].VALUE() == 3
    assert parameters.ROWS[1].OTHER() == 0

    parameters.ROWS.set({1: {"VALUE": 1, "OTHER": 10}, 2: {"VALUE": 2}})
    parameters.ROWS.set({1: {"VALUE": 3}}, merge=True)
    assert list(parameters.ROWS) == [1, 2]
    assert parameters.ROWS[1].VALUE() == 3
    assert parameters.ROWS[1].OTHER() == 10
    assert parameters.ROWS[2].VALUE() == 2

    parameters.ROWS._definition.validators = (repeated_validator_reject_two,)
    with pytest.raises(DataValidityError, match="repeated VALUE must not be two"):
        parameters.ROWS.set({1: {"VALUE": 2}}, merge=True)

    assert list(parameters.ROWS) == [1, 2]
    assert parameters.ROWS[1].VALUE() == 3
    assert parameters.ROWS[2].VALUE() == 2


def test_list_repeated_set_replaces_or_merges_values():
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
    assert [parameters.ROWS[i].VALUE() for i in parameters.ROWS] == [1, 2, 3]


def test_list_repeated_parse_merge_attaches_metadata_to_new_rows():
    def validate_parsed(option):
        parsed = getattr(option._container, "_parsed_values", None)
        if option() == 2 and parsed is None:
            raise RuntimeError("new row has no parse metadata")
        if parsed is not None and parsed["VALUE"] != option():
            raise RuntimeError("parse metadata belongs to another row")

    value = V("VALUE", 0)
    value.validate_parsed = validate_parsed
    definition = cd.InputParametersDefinition(
        [cd.InputSectionDefinition("ROWS", [value], is_repeated=True)]
    )
    parameters = definition.create_object()
    parameters.ROWS.set([{"VALUE": 1}])

    parameters.ROWS.set([{"VALUE": 2}], merge=True, validation_reason="parse")

    assert [parameters.ROWS[i].VALUE() for i in parameters.ROWS] == [1, 2]


def test_generated_setter_uses_mutable_cross_section_proposal():
    definition = cd.InputParametersDefinition(
        [
            cd.InputSectionDefinition(
                "SECTION_A",
                [
                    GeneratedValueDefinition(
                        "GENERATED",
                        get_cross_section_generated,
                        set_cross_section_generated,
                    )
                ],
            ),
            cd.InputSectionDefinition("SECTION_B", [V("VALUE", 0)]),
        ],
        validators=reject_cross_section_two,
    )
    parameters = definition.create_object()

    parameters.SECTION_A.GENERATED = 1
    assert parameters.SECTION_A.GENERATED() == 1
    assert parameters.SECTION_B.VALUE() == 1

    with pytest.raises(
        DataValidityError, match="cross-section generated value must not be two"
    ):
        parameters.SECTION_A.GENERATED = 2
    assert parameters.SECTION_B.VALUE() == 1

    repeated_definition = cd.InputParametersDefinition(
        [
            cd.InputSectionDefinition(
                "ROWS",
                [
                    GeneratedValueDefinition(
                        "GENERATED",
                        get_cross_section_generated,
                        set_cross_section_generated,
                    )
                ],
                is_repeated=True,
            ),
            cd.InputSectionDefinition("SECTION_B", [V("VALUE", 0)]),
        ],
        validators=reject_cross_section_two,
    )
    repeated = repeated_definition.create_object()
    repeated.ROWS.set([{"GENERATED": 1}])
    assert repeated.SECTION_B.VALUE() == 1
    assert len(repeated.ROWS) == 1

    with pytest.raises(
        DataValidityError, match="cross-section generated value must not be two"
    ):
        repeated.ROWS.set([{"GENERATED": 2}, {"GENERATED": 2}])
    assert repeated.SECTION_B.VALUE() == 1
    assert len(repeated.ROWS) == 1


def test_generated_numpy_view_mutates_without_rollback():
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
        validators=reject_negative_raw_data,
    )
    parameters = definition.create_object()
    parameters.DATA.RAW = [1.0, 2.0, 3.0]

    with pytest.raises(ValueError):
        parameters.DATA.VIEW[:] = [4.0, 5.0]
    assert np.array_equal(parameters.DATA.RAW(), [1.0, 2.0, 3.0])

    with pytest.raises(DataValidityError, match="raw data must not be negative"):
        parameters.DATA.VIEW[1] = -2.0
    assert np.array_equal(parameters.DATA.RAW(), [1.0, -2.0, 3.0])


def test_missing_sections_are_validation_results():
    definition = cd.InputParametersDefinition(
        [cd.InputSectionDefinition("REQUIRED", [V("VALUE", int, None)], is_optional=False)]
    )
    found = definition.create_object().check_for_errors(print=False)
    assert any(
        isinstance(issue.message, DataValidityError)
        and "Non-optional section REQUIRED" in str(issue.message)
        for issue in found
    )

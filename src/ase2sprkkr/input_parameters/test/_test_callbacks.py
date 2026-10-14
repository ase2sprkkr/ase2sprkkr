import numpy as np

from ...common.configuration_transaction import ConfigurationTransaction
from ...common.warnings import DataValidityError, DataValidityWarning


def validator_warning_if_two(option, values, why):
    assert values is option._container
    if option() == 2:
        return DataValidityWarning(f"{why}: VALUE should not be two")


def validator_error_if_three(option, values, why):
    if option() == 3:
        return DataValidityError(f"{why}: VALUE must not be three")


def warning_condition_if_two(value):
    if value == 2:
        return "VALUE should not be two"


def section_validator_warning(section, values, why):
    assert values is section
    if values["VALUE"]() == 4:
        return DataValidityWarning(f"{why}: section received its runtime object")


def root_validator_warning(root, values, why):
    assert values is root
    if values["CONTROL"]["VALUE"]() == 5:
        return DataValidityWarning(f"{why}: root received its runtime object")


def root_validator_matching_values(root, values, why):
    if values["SECTION_A"]["VALUE"]() != values["SECTION_B"]["MODE"]():
        return DataValidityError(
            f"{why}: values in SECTION_A and SECTION_B must match"
        )


def root_validator_reject_all(root, values, why):
    return DataValidityError(f"{why}: root validator rejects the configuration")


def repeated_validator_reject_two(section, values, why):
    if values["VALUE"]() == 2:
        return DataValidityError(f"{why}: repeated VALUE must not be two")


def generated_source_length_validator(root, values, why):
    source = values["CONTROL"]["VALUES"]()
    if source is not None and len(source) > 1:
        return DataValidityError(f"{why}: generated setter produced too many values")


def get_cross_section_generated(section, key=None):
    target = section._get_root_container().find("SECTION_B.VALUE")
    return target()


def set_cross_section_generated(section, value, key=None):
    transaction = ConfigurationTransaction.current(section)
    target = section._get_root_container().find("SECTION_B.VALUE")
    target.stage(transaction, value)


def reject_cross_section_two(root, values, why):
    if values["SECTION_B"]["VALUE"]() == 2:
        return DataValidityError(f"{why}: cross-section generated value must not be two")


def reject_negative_raw_data(root, values, why):
    raw = values["DATA"]["RAW"]()
    if raw is not None and np.any(raw < 0):
        return DataValidityError(f"{why}: raw data must not be negative")


def required_in_mode_one(option):
    return option._container["MODE"]() == 1

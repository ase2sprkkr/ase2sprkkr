"""Experimental task and schema metadata regression tests."""

import argparse

import pytest

from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.input_parameters.input_parameters_definitions import (
    InputSectionDefinition,
    InputValueDefinition,
)
from ase2sprkkr.input_parameters.task_definition import TaskDefinition
from ase2sprkkr.tools.commands import info as info_command


LEGACY_TASKS = {
    "ARPES", "BSFEK", "BSFKK", "DOS", "GILBERT",
    "JXC", "PHAGEN", "SCF", "TORQUE", "XAS",
}


def test_experimental_values_sections_copy_and_help():
    value = InputValueDefinition("TEST", 1, experimental=True)
    section = InputSectionDefinition("SECTION", [value], experimental=True)
    assert value.experimental is True
    assert section.experimental is True
    assert value.copy().experimental is True
    assert section.copy().experimental is True
    assert section.copy()["TEST"].experimental is True
    assert "experimental" in value.data_description()
    assert "experimental" in section.data_description()
    assert InputValueDefinition("STABLE", 1).experimental is False


def test_task_experimental_metadata():
    assert TaskDefinition("SCF").experimental is False
    assert TaskDefinition("SCF", experimental=True).experimental is True
    with pytest.raises(TypeError, match="experimental must be a bool"):
        TaskDefinition("SCF", experimental=1)


def test_legacy_tasks_remain_stable():
    tasks = InputParameters.task_registry.values()
    stable = {task.name for task in tasks if not task.experimental}
    assert stable == LEGACY_TASKS
    assert InputParameters.task_registry["DMI"].experimental
    assert InputParameters.task_registry["AES"].experimental
    assert InputParameters.task_registry["PLOTSFN"].experimental


def test_experimental_schema_metadata_propagates():
    registry = InputParameters.input_parameters_registry
    tasks = InputParameters.task_registry
    assert registry["AES"].experimental
    assert registry["SCF"].experimental is False
    assert registry["PLOTSFN"].experimental
    assert tasks["AES"] is registry["AES"]._task_definition


def test_info_lists_stable_tasks_before_experimental(capsys):
    parser = argparse.ArgumentParser()
    info_command.parser(parser)
    info_command.run(parser.parse_args(["--task"]), None)
    output = capsys.readouterr().out
    heading = "  Experimental tasks\n  ~~~~~~~~~~~~~~~~~~"
    assert output.count(heading) == 1
    regular, experimental = output.split(heading)
    assert "SCF" in regular
    assert "BSFEK" in regular
    assert "BSFKK" in regular
    assert "AES" in experimental
    assert "DMI" in experimental
    assert "PLOTSFN" in experimental
    assert " (experimental)" not in output

    info_command.run(parser.parse_args(["--task=dmi"]), None)
    filtered = capsys.readouterr().out
    assert heading in filtered
    assert "DMI" in filtered
    assert " (experimental)" not in filtered

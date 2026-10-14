from ..input_parameters import InputParameters


def test_task_registry_discovers_only_task_modules_without_loading_them():
    registry = InputParameters.task_registry

    names = registry.names(include_aliases=False)

    assert "SCF" in names
    assert "SECTIONS" not in names
    assert "BSFEK" not in names
    assert "BSFKK" not in names


def test_task_registry_aliases_are_task_definitions():
    canonical = InputParameters.task_registry["BSF"]
    bsfek = InputParameters.task_registry["BSFEK"]
    bsfkk = InputParameters.task_registry["BSFKK"]

    assert bsfek is not canonical
    assert bsfkk is not canonical
    assert bsfek.input_parameters_definition is canonical.input_parameters_definition
    assert bsfkk.input_parameters_definition is canonical.input_parameters_definition
    assert bsfek.name == "BSFEK"
    assert bsfkk.name == "BSFKK"

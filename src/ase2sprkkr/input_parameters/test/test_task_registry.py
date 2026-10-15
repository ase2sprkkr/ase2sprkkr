from ..input_parameters import InputParameters


def test_standard_task_comes_from_input_parameters_definition():
    schema = InputParameters.input_parameters_registry["SCF"]

    assert schema._task_definition is InputParameters.task_registry["SCF"]
    assert "SCF" not in InputParameters.task_registry.extra_tasks.names()


def test_spec_tasks_share_one_schema_and_resolve_from_values():
    schema = InputParameters.input_parameters_registry["ARPES"]
    spleed = InputParameters.create("SPLEED")

    assert schema._task_definition is None
    assert spleed._definition is schema
    assert spleed.task_definition.name == "SPLEED"

    spleed.TASK.TASK.set("ARPES")
    assert spleed.task_definition.name == "ARPES"


def test_bsf_views_share_schema_but_defaults_stay_in_schema():
    schema = InputParameters.input_parameters_registry["BSF"]
    ek = InputParameters.create("BSFEK")
    kk = InputParameters.create("BSFKK")

    assert schema._task_definition is None
    assert ek._definition is schema
    assert kk._definition is schema
    assert ek.ENERGY.NE()[0] > 1
    assert kk.ENERGY.NE()[0] == 1
    assert ek.task_definition.name == "BSFEK"
    assert kk.task_definition.name == "BSFKK"


def test_input_and_task_registry_namespaces_are_distinct():
    input_names = set(InputParameters.input_parameters_registry.names())
    task_names = set(InputParameters.task_registry.names())

    assert {"ARPES", "BSF"} <= input_names
    assert "SPLEED" not in input_names
    assert "BSFEK" not in input_names
    assert {"ARPES", "SPLEED", "BSF", "BSFEK", "BSFKK"} <= task_names

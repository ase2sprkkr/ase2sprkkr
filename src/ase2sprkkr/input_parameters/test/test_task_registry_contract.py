from ..input_parameters import InputParameters


def test_task_registry_unifies_colocated_and_shared_task_declarations():
    task_registry = InputParameters.task_registry
    input_registry = InputParameters.input_parameters_registry

    scf = task_registry["SCF"]
    scf_schema = input_registry.definition_for_task(scf)
    assert scf_schema.name.upper() == "SCF"
    assert scf_schema._task_definition is scf

    arpes = task_registry["ARPES"]
    spleed = task_registry["SPLEED"]
    arpes_schema = input_registry.definition_for_task(arpes)
    assert arpes_schema is not input_registry.definition_for_task(spleed)
    assert arpes_schema._task_definition is None
    assert arpes.input_parameters_definition_name == arpes_schema.name.upper()
    assert spleed.input_parameters_definition_name == "SPLEED"
    assert input_registry.definition_for_task(spleed)._task_definition is spleed

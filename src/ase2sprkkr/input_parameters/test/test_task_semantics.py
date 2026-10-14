from ase2sprkkr.input_parameters.input_parameters import InputParameters
from ase2sprkkr.outputs.task_result import TaskResult


def test_task_definition_declares_potential_production():
    assert InputParameters.definition("SCF").produces_potential
    assert not InputParameters.definition("DOS").produces_potential


def test_task_result_exposes_definition_potential_semantics(tmp_path):
    scf_result = TaskResult(InputParameters.create("SCF"), None, tmp_path)
    dos_result = TaskResult(InputParameters.create("DOS"), None, tmp_path)

    assert scf_result.produces_potential
    assert not dos_result.produces_potential

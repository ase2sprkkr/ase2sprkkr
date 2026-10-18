from ..task_definition import TaskDefinition
from ._spec_files import SPEC_INPUT_FILES, SPEC_OUTPUT_FILES


def task_definition():
    return TaskDefinition("ARPES", name="SPLEED", executable="kkrspec", mpi=False,
        input_file_templates=SPEC_INPUT_FILES,
        output_file_templates=SPEC_OUTPUT_FILES)

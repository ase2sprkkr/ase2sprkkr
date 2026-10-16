from ..definitions._file_manifests import kkrspec_task_file_templates
from ..task_definition import TaskDefinition


def task_definition():
    input_files, output_files = kkrspec_task_file_templates()
    return TaskDefinition(
        "ARPES",
        name="BAND",
        executable="kkrspec",
        input_file_templates=input_files,
        output_file_templates=output_files,
    )

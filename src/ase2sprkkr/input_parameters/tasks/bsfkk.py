from ..task_definition import TaskDefinition


def task_definition():
    return TaskDefinition(
        "BSF",
        name="BSFKK",
        executable="kkrgen",
        mpi=True,
        output_reader="BSF",
        output_file_templates=("{dataset}_BLOCHSF.bsf", "{dataset}_BLOCHSF_spol.bsf"),
    )

from ..task_definition import TaskDefinition, DeclaredFile


def task_definition():
    return TaskDefinition(
        "BSF",
        name="BSFKK",
        executable="kkrgen",
        mpi=True,
        output_reader="BSF",
        output_file_templates=("{dataset}_BLOCHSF.bsf", DeclaredFile("{dataset}_BLOCHSF_spol.bsf", optional=True)),
    )

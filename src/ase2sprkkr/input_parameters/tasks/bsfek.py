from ..definitions._file_manifests import tau_file
from ..task_definition import FileSpec, TaskDefinition, TaskFiles


def _bsf_output(parameters, _directory):
    dataset = parameters.CONTROL.DATASET()
    adsi = parameters.CONTROL.ADSI()
    if not dataset or not adsi:
        return None
    return f"{dataset}_{adsi}.bsf"


def task_definition():
    return TaskDefinition(
        "BSF",
        name="BSFEK",
        executable="kkrgen",
        mpi=True,
        output_reader="BSF",
        files=TaskFiles(
            output=(FileSpec(_bsf_output),),
            cache=(FileSpec(tau_file),),
        ),
    )

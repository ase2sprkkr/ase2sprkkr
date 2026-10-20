from ase2sprkkr.execution import ExecutionSpec, LaunchSpec
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor


class _Transport:
    host = "cluster"


def test_remote_project_runner_uses_manager_run_id(tmp_path):
    executor = RemoteEndpointExecutor(_Transport())
    spec = ExecutionSpec(
        argv=("true",),
        input_file=None,
        output_file=None,
        directory=str(tmp_path),
        launch=LaunchSpec(local_root=tmp_path, remote_root="/scratch/work"),
    )

    runner = executor.create_project_run(
        spec, run_id="ase2sprkkr-run-id"
    )

    assert runner._run_id == "ase2sprkkr-run-id"

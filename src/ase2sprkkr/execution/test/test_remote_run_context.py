from ase2sprkkr.execution import Executor, LaunchSpec
from ase2sprkkr.execution.remote_endpoint import RemoteEndpointExecutor


class _Transport:
    host = "cluster"



def test_run_context_is_staging_policy_not_executor(tmp_path):
    endpoint = RemoteEndpointExecutor(_Transport())
    context = endpoint._context(
        endpoint,
        LaunchSpec(
            local_root=tmp_path,
            remote_root="/scratch/a",
            input_sync="always",
            output_sync="missing",
        ),
    )

    assert not isinstance(context, Executor)
    assert context.endpoint is endpoint
    assert context.transport is endpoint.transport
    assert context.local_root == tmp_path.resolve()
    assert context.remote_root == "/scratch/a"



def test_different_run_staging_policies_share_one_endpoint(tmp_path):
    endpoint = RemoteEndpointExecutor(_Transport())
    first = endpoint._context(
        endpoint,
        LaunchSpec(local_root=tmp_path, remote_root="/scratch/one"),
    )
    second = endpoint._context(
        endpoint,
        LaunchSpec(local_root=tmp_path, remote_root="/scratch/two"),
    )

    assert first is not second
    assert first.endpoint is endpoint
    assert second.endpoint is endpoint
    assert first.remote_root == "/scratch/one"
    assert second.remote_root == "/scratch/two"

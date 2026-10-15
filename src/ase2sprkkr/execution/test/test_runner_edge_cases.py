from types import SimpleNamespace

import pytest

from ase2sprkkr.execution import ExecutionSpec
from ase2sprkkr.execution.remote import RemoteRunner


class _Reader:
    print_output = False

    def create_result(self, output_file, input_file):
        return SimpleNamespace(return_code=None)


def test_remote_cancel_before_submit_is_rejected(tmp_path):
    spec = ExecutionSpec(
        argv=(),
        input_file=None,
        output_file=None,
        directory=str(tmp_path),
    )
    runner = RemoteRunner(SimpleNamespace(), _Reader(), spec)

    with pytest.raises(RuntimeError, match="not been submitted"):
        runner.cancel()

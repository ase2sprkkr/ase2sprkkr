import multiprocessing as mp
import subprocess as stdlib_subprocess
import sys

import pytest

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from .. import subprocess as subprocess_module


class TestCommonSubprocess(TestCase):

    def test_import_does_not_change_start_method(self):
        methods = [
            method
            for method in ("fork", "forkserver")
            if method in mp.get_all_start_methods()
        ]
        if not methods:
            pytest.skip("No non-spawn multiprocessing start method is available")

        code = f"""
import multiprocessing as mp
mp.set_start_method({methods[0]!r}, force=True)
before = mp.get_start_method()
import ase2sprkkr.common.subprocess
assert mp.get_start_method() == before
"""
        result = stdlib_subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)

    def test_executor_uses_spawn(self):
        try:
            method = subprocess_module.in_subprocess(
                "get_start_method",
                "multiprocessing",
                (),
                {},
            )
            self.assertEqual(method, "spawn")
        finally:
            if subprocess_module._pool:
                subprocess_module._pool.shutdown()
                subprocess_module._pool = None

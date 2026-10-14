import multiprocessing as mp
import signal
import subprocess as stdlib_subprocess
import sys

import pytest

from .. import subprocess as subprocess_module


@pytest.fixture(autouse=True)
def shutdown_subprocess_pool():
    yield
    pool = subprocess_module._pool
    if pool is not None:
        pool.shutdown()
        subprocess_module._pool = None


def _get_sigchld_handler():
    return signal.getsignal(signal.SIGCHLD)


def test_import_does_not_change_start_method():
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
    assert result.returncode == 0, result.stderr


def test_executor_uses_spawn():
    method = subprocess_module.in_subprocess(
        "get_start_method",
        "multiprocessing",
        (),
        {},
    )
    assert method == "spawn"


def test_sigchld_is_ignored_only_in_worker(monkeypatch):
    def fail_if_called_in_parent(_num):
        raise AssertionError("SIGCHLD must not be changed in the parent process")

    monkeypatch.setattr(subprocess_module, "ignore_signal", fail_if_called_in_parent)

    handler = subprocess_module.in_subprocess(
        "getsignal",
        "signal",
        (signal.SIGCHLD,),
        {},
    )
    assert handler is signal.SIG_IGN

    handler = subprocess_module.in_subprocess(_get_sigchld_handler)
    assert handler is signal.SIG_IGN

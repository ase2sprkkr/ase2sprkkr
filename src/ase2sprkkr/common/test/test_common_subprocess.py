import multiprocessing as mp
import signal
import subprocess as stdlib_subprocess
import sys

import pytest

from .. import subprocess as subprocess_module


@pytest.fixture(autouse=True)
def shutdown_subprocess_pool():
    yield
    subprocess_module.close_subprocess_pool()


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


def test_spawned_worker_uses_spawn_and_ignores_sigchld(monkeypatch):
    def fail_if_called_in_parent(_num):
        raise AssertionError("SIGCHLD must not be changed in the parent process")

    monkeypatch.setattr(subprocess_module, "ignore_signal", fail_if_called_in_parent)

    # These independent checks run in the same worker. Spawning a new
    # process for the signal checks would only repeat expensive setup.
    method = subprocess_module.in_subprocess(
        "get_start_method",
        "multiprocessing",
        (),
        {},
    )
    assert method == "spawn"

    handler = subprocess_module.in_subprocess(
        "getsignal",
        "signal",
        (signal.SIGCHLD,),
        {},
    )
    assert handler is signal.SIG_IGN

    handler = subprocess_module.in_subprocess(_get_sigchld_handler)
    assert handler is signal.SIG_IGN


def test_failed_call_shuts_down_pool():
    with pytest.raises(RuntimeError, match="subprocess"):
        subprocess_module.in_subprocess(lambda: None)

    assert subprocess_module._pool is None


def test_failed_old_call_does_not_discard_a_new_pool(monkeypatch):
    closed = []
    replacement = object()

    class Future:
        def result(self):
            subprocess_module._pool = replacement
            raise ValueError("old worker failed")

    class Pool:
        def submit(self, *args, **kwargs):
            return Future()

        def shutdown(self, **kwargs):
            closed.append(self)

    previous = Pool()
    monkeypatch.setattr(subprocess_module, "_pool", previous)
    try:
        with pytest.raises(RuntimeError, match="old worker failed"):
            subprocess_module.in_subprocess(int)
        assert subprocess_module._pool is replacement
        assert closed == [previous]
    finally:
        subprocess_module._pool = None

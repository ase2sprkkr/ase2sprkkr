from concurrent.futures import ProcessPoolExecutor
import importlib
import signal
import contextlib

_pool = None

import multiprocessing as mp


@contextlib.contextmanager
def ignore_signal(num):
    old = signal.getsignal(num)
    signal.signal(num, signal.SIG_IGN)
    try:
        yield
    finally:
        signal.signal(num, old)


def _call(lmbda):
    with ignore_signal(signal.SIGCHLD):
        return lmbda()


def _call_function(name, module, args, kwargs):
    mod = importlib.import_module(module)
    with ignore_signal(signal.SIGCHLD):
        return getattr(mod, name)(*args, **kwargs)


def in_subprocess(lmbda, module=None, *args, **kwargs):
    global _pool

    if not _pool:
        _pool = ProcessPoolExecutor(
            max_workers=1,
            mp_context=mp.get_context("spawn"),
        )
    try:
        if module:
            future = _pool.submit(_call_function, lmbda, module, *args, **kwargs)
        else:
            future = _pool.submit(_call, lmbda)
        return future.result()
    except Exception as e:
        pool = _pool
        _pool = None
        pool.shutdown(wait=True)
        raise RuntimeError(f"Error when running the subprocess: {e}") from e

import contextlib
import importlib
import multiprocessing as mp
import signal
from concurrent.futures import ProcessPoolExecutor
from threading import Lock

_pool = None
_pool_lock = Lock()


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


def close_subprocess_pool():
    """Release the shared worker, waiting for its already submitted calls."""
    global _pool
    with _pool_lock:
        pool, _pool = _pool, None
    if pool is not None:
        pool.shutdown()


def in_subprocess(lmbda, module=None, *args, **kwargs):
    global _pool

    with _pool_lock:
        if _pool is None:
            _pool = ProcessPoolExecutor(
                max_workers=1,
                mp_context=mp.get_context("spawn"),
            )
        pool = _pool
    try:
        if module:
            future = pool.submit(_call_function, lmbda, module, *args, **kwargs)
        else:
            future = pool.submit(_call, lmbda)
        return future.result()
    except Exception as e:
        with _pool_lock:
            # Another caller may already have replaced the failed pool.
            if _pool is pool:
                _pool = None
        pool.shutdown(wait=True)
        raise RuntimeError(f"Error when running the subprocess: {e}") from e

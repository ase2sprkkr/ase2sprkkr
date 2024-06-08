from concurrent.futures import ProcessPoolExecutor
import importlib
import signal
import contextlib

_pool=None


@contextlib.contextmanager
def ignore_signal(num):
    old = signal.getsignal(num)
    signal.signal(num, signal.SIG_IGN)
    try:
      yield
    finally:
      signal.signal(num, old)


def _call_function(name, module, args, kwargs):
    mod = importlib.import_module(module)
    out = getattr(mod, name)(*args, **kwargs)
    return out


def in_subprocess(lmbda, module=None, *args, **kwargs):
    global _pool

    if not _pool:
       _pool = ProcessPoolExecutor(max_workers=1)
    try:
        with ignore_signal(signal.SIGCHLD):
            if module:
                future = _pool.submit(_call_function, lmbda, module, *args, **kwargs)
            else:
                future = _pool.submit(lmbda)
    except Exception as e:
        _pool.shutdown(False)
        _pool = None
        raise RuntimeError("Error when running the supbprocess") from e

    return future.result()

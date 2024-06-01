from concurrent.futures import ProcessPoolExecutor
import importlib
_pool=None


def _call_function(name, module, args, kwargs):
    mod = importlib.import_module(module)
    return getattr(mod, name)(*args, **kwargs)


def in_subprocess(lmbda, module=None, *args, **kwargs):
    global _pool
    if not _pool:
       _pool = ProcessPoolExecutor(max_workers=1)
    if module:
        future = _pool.submit(_call_function, lmbda, module, *args, **kwargs)
    else:
        future = _pool.submit(lmbda)

    return future.result()

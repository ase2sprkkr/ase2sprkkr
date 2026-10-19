"""Run async operations from synchronous APIs, including inside event loops."""

import asyncio
from concurrent.futures import Future
from threading import Thread


def run_coro_sync(coro):
    """Run a coroutine synchronously without nesting the caller's event loop.

    A running loop requires a separate thread. Propagate its return value or
    *any* exception to the original caller instead of losing worker failures.
    """
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coro)

    result = Future()

    def worker():
        try:
            result.set_result(asyncio.run(coro))
        except BaseException as exc:
            result.set_exception(exc)

    thread = Thread(target=worker)
    thread.start()
    thread.join()
    return result.result()

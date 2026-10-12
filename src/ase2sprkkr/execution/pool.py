"""Reusable execution backend pooling independent of run persistence."""

from __future__ import annotations

from threading import RLock


class ExecutorPool:
    """Process-local cache of reusable execution backends."""

    def __init__(self):
        self._executors = {}
        self._lock = RLock()

    @staticmethod
    def key_for_profile(profile):
        """Return the backend identity used for executor reuse.

        Local and PSI/J executors depend only on the endpoint. Remote execution
        still keeps staging policy on ``RemoteExecutor`` today, so SSH profiles
        temporarily include launch parameters in the key until that policy is
        moved into the per-run execution specification.
        """
        if profile.executor == "ssh":
            return profile.endpoint, profile.launch_parameters
        return profile.endpoint

    def get(self, key, factory):
        with self._lock:
            executor = self._executors.get(key)
            if executor is None:
                executor = factory()
                self._executors[key] = executor
            return executor

    def for_profile(self, profile):
        """Return the shared backend suitable for *profile*."""
        from ..sprkkr.profile import executor_from_profile

        key = self.key_for_profile(profile)
        return self.get(key, lambda: executor_from_profile(profile))

    def discard(self, key):
        """Forget one pooled backend and close it when supported."""
        with self._lock:
            executor = self._executors.pop(key, None)
        if executor is not None:
            close = getattr(executor, "close", None)
            if close is not None:
                close()
        return executor

    def discard_profile(self, profile):
        return self.discard(self.key_for_profile(profile))

    def close(self):
        """Close and remove all currently pooled executors."""
        with self._lock:
            items = tuple(self._executors.items())
            self._executors.clear()
        for _, executor in items:
            close = getattr(executor, "close", None)
            if close is not None:
                close()


_DEFAULT_EXECUTOR_POOL = ExecutorPool()


def default_executor_pool():
    """Return the process-wide pool used by ordinary calculations."""
    return _DEFAULT_EXECUTOR_POOL

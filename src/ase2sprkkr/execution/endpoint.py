"""Serializable execution endpoint identity and reconstruction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class ExecutorEndpoint:
    """Identity of a reusable execution backend.

    The endpoint contains only connection/backend properties. Launch and
    staging policy deliberately live elsewhere so one executor can be shared
    by profiles that differ only in how individual calculations are launched.
    The value is serializable and hashable, so it can also be persisted with a
    run and used as an :class:`ExecutorPool` key after process restart.
    """

    executor: str = "local"
    host: Any = None
    ssh_options: Any = ()
    remote_command: Any = ()
    psij_backend: Any = None
    psij_url: Any = None

    def __post_init__(self):
        object.__setattr__(self, "executor", str(self.executor or "local"))
        for name in ("ssh_options", "remote_command"):
            value = getattr(self, name)
            if value is None:
                value = ()
            elif isinstance(value, str):
                value = (value,)
            else:
                value = tuple(value)
            object.__setattr__(self, name, value)


def endpoint_for_executor(executor) -> ExecutorEndpoint:
    """Describe a live executor without persisting launch/staging policy."""
    from .local import LocalExecutor
    from .psij import PsijExecutor
    from .remote_endpoint import RemoteEndpointExecutor
    from .transport import SshTransport

    if isinstance(executor, LocalExecutor):
        return ExecutorEndpoint()
    if isinstance(executor, PsijExecutor):
        return ExecutorEndpoint(
            executor="psij",
            psij_backend=executor.backend,
            psij_url=executor.url,
        )
    if isinstance(executor, RemoteEndpointExecutor):
        transport = executor.transport
        psij_backend = None
        psij_url = None
        if executor.backend == "psij":
            options = executor.backend_options or {}
            psij_backend = options.get("backend") or "slurm"
            psij_url = options.get("url")

        if not isinstance(transport, SshTransport):
            return ExecutorEndpoint(
                executor=(
                    "custom:{}.{}:{}".format(
                        type(executor).__module__,
                        type(executor).__qualname__,
                        executor.name or "remote",
                    )
                ),
                host=getattr(transport, "host", None),
                remote_command=executor.command,
                psij_backend=psij_backend,
                psij_url=psij_url,
            )
        return ExecutorEndpoint(
            executor="ssh",
            host=transport.host,
            ssh_options=transport.options,
            remote_command=executor.command,
            psij_backend=psij_backend,
            psij_url=psij_url,
        )

    return ExecutorEndpoint(
        executor="custom:{}.{}:{}".format(
            type(executor).__module__,
            type(executor).__qualname__,
            executor.name or type(executor).__qualname__,
        )
    )


def create_executor_from_endpoint(endpoint: ExecutorEndpoint):
    """Reconstruct a standard executor from persisted endpoint identity."""
    from .local import LocalExecutor
    from .psij import PsijExecutor
    from .remote_endpoint import RemoteEndpointExecutor
    from .transport import SshTransport

    if isinstance(endpoint, dict):
        endpoint = ExecutorEndpoint(**endpoint)

    if endpoint.executor == "local":
        return LocalExecutor()
    if endpoint.executor == "psij":
        return PsijExecutor(
            endpoint.psij_backend or "slurm",
            url=endpoint.psij_url,
        )
    if endpoint.executor == "ssh":
        if not endpoint.host:
            raise ValueError("SSH execution endpoint requires 'host'")
        kwargs = {}
        if endpoint.remote_command:
            kwargs["command"] = endpoint.remote_command
        if endpoint.psij_backend is not None or endpoint.psij_url is not None:
            kwargs["backend"] = "psij"
            kwargs["backend_options"] = {
                "backend": endpoint.psij_backend or "slurm",
                "url": endpoint.psij_url,
            }
        return RemoteEndpointExecutor(
            SshTransport(endpoint.host, options=endpoint.ssh_options),
            **kwargs,
        )
    raise ValueError(
        f"Cannot reconstruct unsupported execution endpoint {endpoint.executor!r}"
    )

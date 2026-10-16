"""Remote endpoint public compatibility import.

The implementation lives in :mod:`remote_endpoint_client`; remote lifecycle
state is owned exclusively by the client.
"""

from .remote_endpoint_client import RemoteEndpointExecutor

__all__ = ["RemoteEndpointExecutor"]

"""Temporary compatibility shim for the historical ``OutputFile`` name."""

from .result_files import ResultFile

# TODO: Remove this compatibility alias in a future release after downstream
# users have migrated to ``ResultFile``.
OutputFile = ResultFile

__all__ = ["OutputFile"]

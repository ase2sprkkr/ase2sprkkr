"""Deprecated compatibility package for the former structured result-file API."""

import warnings

from ..result_files import ResultFile as OutputFile

warnings.warn(
    "ase2sprkkr.output_files is deprecated; use ase2sprkkr.result_files instead.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ["OutputFile"]

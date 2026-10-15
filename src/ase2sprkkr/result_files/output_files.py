"""Internal compatibility imports for result-file format definitions."""

from .result_files import Arithmetic, CommonOutputFile, ResultFile

# Concrete format definitions historically use this internal base name.
# Public compatibility lives exclusively in ase2sprkkr.output_files.
OutputFile = ResultFile

__all__ = ["Arithmetic", "CommonOutputFile"]

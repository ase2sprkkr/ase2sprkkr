"""
Module version
==============

Just the version of the ASE2SPRKKR package.
"""
try:
  from importlib.metadata import version, PackageNotFoundError
except ImportError:
  from importlib_metadata import version, PackageNotFoundError

#: Version number of the ASE2SPRKKR package
try:
    __version__ = version('ase2sprkkr')
except PackageNotFoundError:
    __version__ = None

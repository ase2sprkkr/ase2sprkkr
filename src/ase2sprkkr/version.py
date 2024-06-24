"""
Module version
==============

Just the version of the ASE2SPRKKR package.
"""
#: Version number of the ASE2SPRKKR package
import pkg_resources

try:
    __version__ = pkg_resources.get_distribution('ase2sprkkr').version
except pkg_resources.DistributionNotFound:
    __version__ = None

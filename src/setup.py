#vim: set et ts=4 sw=4:
# coding: utf-8

from setuptools import setup, find_packages
import sprkkr

setup(name='sprkkr',
      version=sprkkr.__version__,
      include_package_data=True,
      packages=find_packages(include='sprkkr.*'),
      install_requires=[
	'numpy',
        'ase',
	'spglib',
      ])

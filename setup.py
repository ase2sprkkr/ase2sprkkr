import setuptools
import os

version_file = os.path.join(os.path.dirname(__file__), "src/ase2sprkkr/version.py")
version = {}
with open(version_file) as ver_file:
    exec(ver_file.read(), version)

setuptools.setup(
    version = version['__version__']
)

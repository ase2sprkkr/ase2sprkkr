""" This module contains configuration, that could be changed, preferrably
by .config/ase2sprkkr/__init__.py file"""

import os

""" These values overrides calculators defaults """
calculator_parameters = {
  'empty_spheres' : 'auto',
  'mpi' : True,
  'print_output': 'info'
}

""" This suffix is appended (if not stated otherwise) the SPRKKR executable names. """
sprkkr_executable_suffix = os.getenv('SPRKKR_EXECUTABLE_SUFFIX', '')

""" The executable to run MPI, given as list (with possible arguments).
In the most cases, the ``None`` value - autodetection - is sufficient """
mpi_runner = None

del os

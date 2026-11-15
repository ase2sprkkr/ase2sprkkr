Ase2sprkkr is a modern, comprehensive replacement for
`XBAND <https://sprkkr.org/xband/>`_, providing a unified framework for
`SPR-KKR <https://sprkkr.org/>`_ calculations, ranging from easy
interactive use to efficient, highly automated workflows on HPC systems.
It integrates with the
`Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/>`_
and supports preparing, running, and analyzing calculations through Python.

Why use it?
-----------

The underlying SPR-KKR package uses the Korringa-Kohn-Rostoker (KKR)
multiple-scattering, Green's-function approach. Its applications include
density of states, Bloch spectral functions, magnetism, chemically
disordered alloys (CPA), relativistic calculations and spectroscopy.
The physics, approximations and numerical choices belong to
`the SPR-KKR manual <https://sprkkr.org/manual/>`_;
this documentation focuses on working with the software.

**ASE2SPRKKR** creates native SPR-KKR inputs and potentials from ASE
structures, selects and configures tasks, starts calculations, and parses
and visualizes supported results. The independent
`KKRunner <https://ase2sprkkr.github.io/kkrunner/>`_ package provides
local, MPI, SSH and scheduler execution and optional persistent history.
A separate `ASE2SPRKKR GUI <https://github.com/ase2sprkkr/ase2sprkkr-gui>`_
provides a graphical interface for preparing and managing calculations.

SPR-KKR executables are **installed separately**; ASE2SPRKKR does not
include those binaries. Calculations are task-oriented (SCF, DOS,
BSF, ARPES and others), and not every standard ASE calculator property
is available.

Begin with :doc:`quickstart` to run an SCF-to-DOS example, then
:doc:`installation` for configuration and upgrades, or
:doc:`usage` for the Python API.

Useful links
------------

* `SPR-KKR download/access <https://sprkkr.org/download/>`_
* `SPR-KKR manual <https://sprkkr.org/manual/>`_
* `Source code and issue tracker <https://github.com/ase2sprkkr/ase2sprkkr>`_
* `KKRunner documentation <https://ase2sprkkr.github.io/kkrunner/>`_
* `KKRunner GitHub repository <https://github.com/ase2sprkkr/kkrunner>`_
* `ASE2SPRKKR GUI <https://github.com/ase2sprkkr/ase2sprkkr-gui>`_
* `ASE documentation <https://wiki.fysik.dtu.dk/ase/>`_
* `ASE2SPRKKR discussion group <https://groups.google.com/g/ase2sprkkr/>`_

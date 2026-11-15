ASE2SPRKKR — a Python interface for SPR-KKR
===========================================

Ase2sprkkr is a modern, comprehensive replacement for
`XBAND <https://sprkkr.org/xband/>`_, providing a unified framework for
`SPR-KKR <https://sprkkr.org/>`_ calculations, ranging from easy
interactive use to efficient, highly automated workflows on HPC systems.
It integrates with the
`Atomic Simulation Environment (ASE) <https://wiki.fysik.dtu.dk/ase/>`_
and supports preparing, running, and analyzing calculations through Python.

SPR-KKR uses the multiple-scattering KKR Green's-function method and
supports electronic structure, magnetism, chemically disordered alloys
(CPA), Bloch spectral functions and spectroscopy. See the
`SPR-KKR manual <https://sprkkr.org/manual/>`_ for the physics and
numerical choices; this documentation focuses on the Python interface.

ASE2SPRKKR creates input files and starting potentials from ASE
structures, runs native SPR-KKR tasks and parses supported results.
`KKRunner <https://ase2sprkkr.github.io/kkrunner/>`_ provides local,
MPI, SSH and scheduler execution. **SPR-KKR executables are installed
separately**; they are not included in the Python distribution.

Start here
----------

* `First calculation: SCF to DOS <doc_src/quickstart.rst>`_
* `Installation, upgrading, virtual environments and optional dependencies
  <doc_src/installation.rst>`_
* `Python usage and task parameters <doc_src/usage.rst>`_
* `Examples generated from Python sources <doc_src/examples.rst>`_
* `Command-line interface <doc_src/cli.rst>`_

Install or upgrade the latest stable release::

   python -m pip install --upgrade ase2sprkkr
   ase2sprkkr --version

To install an available **beta** release (including the 4.0 series),
use ``python -m pip install --upgrade --pre ase2sprkkr``. Some
features described in the development documentation may not yet be
available in a published version.

For a separate test environment or a distribution-managed Python installation,
see `virtual environments <doc_src/installation.rst#virtual-environments>`_.

The native executables (e.g. ``kkrscf`` and ``kkrgen``) can be found
on ``PATH`` or configured through the active execution profile::

   ase2sprkkr config --set profiles.default.executable_dir /opt/sprkkr/bin

The ``profiles.default`` alias uses an explicitly selected named profile
or otherwise ASE's ``[sprkkr]`` configuration. To edit the ASE profile
directly, use ``profiles.ase``. A suffix can also be configured::

   ase2sprkkr config --set profiles.ase.executable_suffix _v8

A short SCF calculation starts with ASE atoms:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import SPRKKR

   result = SPRKKR(atoms=bulk("Al", "fcc", a=4.05)).calculate()
   print("Energy (eV):", result.energy)
   print("SCF converged:", result.converged)
   print("Potential:", result.potential_filename)

Always check convergence before reusing an SCF potential for other
tasks. The `first-calculation tutorial <doc_src/quickstart.rst>`_
continues with a DOS calculation.

Further documentation
---------------------

* `Result files, plotting and NeXus <doc_src/result-files.rst>`_
* `Remote running, MPI and execution profiles <doc_src/remote.rst>`_
* `Persistent project calculations <doc_src/execution.rst>`_
* `Crystal symmetry <doc_src/symmetry.rst>`_
* `Architecture for developers and TaskDefinition <doc_src/philosophy.rst>`_
* `Changelog <doc_src/changelog.rst>`_

For upgrades from ASE2SPRKKR 3.x, configuration migration and
optional features such as plotting or NeXus export, see the
`installation guide <doc_src/installation.rst>`_.

Contributing
------------

Report problems and propose changes on
`GitHub <https://github.com/ase2sprkkr/ase2sprkkr>`_.

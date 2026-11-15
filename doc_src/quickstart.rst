First calculation: aluminium SCF and DOS
========================================

This tutorial shows how to run a first SPR-KKR calculation using ASE2SPRKKR
and reuse the converged potential to calculate the density of states.
For physics, numerical settings and convergence studies, consult the
`SPR-KKR manual <https://sprkkr.org/manual/>`_.

Install the Python interface and SPR-KKR
----------------------------------------

Install the native SPR-KKR programs separately, following the
`official download instructions <https://sprkkr.org/download/>`_.
Then install the Python interface::

   python -m pip install --upgrade ase2sprkkr
   ase2sprkkr --version

This tutorial covers the 4.x interface. If only a beta of version 4 is
available, install it with ``python -m pip install --upgrade --pre ase2sprkkr``.
For more on versions see :doc:`installation`.

For a separate Python environment (for example, to avoid dependency
conflicts), see :ref:`installation-venv`.

Check whether SPR-KKR is available on ``PATH``, for example with
``command -v kkrscf`` on Unix-like systems. If not, configure
where its executables are installed::

   ase2sprkkr config --set profiles.default.executable_dir /opt/sprkkr/bin

This modifies the active execution profile, which uses ASE's
``[sprkkr]`` configuration by default.
If your binaries use a suffix such as ``kkrscf_v8``, configure it::

   ase2sprkkr config --set profiles.default.executable_suffix _v8

MPI selection is automatic by default; see :doc:`remote` for
execution profiles and MPI configuration.

Run an SCF calculation
----------------------

Save the following in ``first_calculation.py`` in a writable directory:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import SPRKKR

   atoms = bulk("Al", "fcc", a=4.05)
   calculator = SPRKKR(
       atoms=atoms,
       input_parameters="SCF",
       directory="al-scf",
   )
   scf = calculator.calculate()

   print("Total energy (eV):", scf.energy)
   print("SCF converged:", scf.converged)
   print("Iterations:", len(scf.iterations))
   print("Resulting potential:", scf.potential_filename)

   if not scf.converged:
       raise RuntimeError(
           "SCF did not converge; inspect the native output "
           "and convergence parameters before reusing its potential."
       )

Run with ``python first_calculation.py``. ``calculate()`` prepares
input files, runs the native executable and returns a parsed result.
Successful process termination does not by itself establish physical
convergence; inspect ``scf.converged`` and the native output.

SPR-KKR potential files (``.pot``) contain atomic-site and electronic
potential data, including information needed for further calculations.
They are not the equivalent of VASP's POTCAR or a PAW dataset.

Run a DOS calculation
---------------------

Add the following after the SCF convergence check:

.. code-block:: python

   dos = SPRKKR(
       input_parameters="DOS",
       directory="al-dos",
   ).calculate(potential=scf.potential_filename)
   print("DOS file:", dos.dos_filename)

This runs a different native task with the potential produced by SCF.
For plotting, install the :ref:`plot extra <installation-optional>`::

   python -m pip install 'ase2sprkkr[plot]'

Plot the DOS, replacing the path with the ``dos.dos_filename``
printed above::

   ase2sprkkr plot path/to/result.dos
   ase2sprkkr plot path/to/result.dos --show-values
   ase2sprkkr plot path/to/result.dos -o dos.png

For structured result readers and formats see :doc:`result-files`.
This introductory script is not a recommended set of production numerical
parameters; see the SPR-KKR manual for physical assumptions and convergence.

Troubleshooting
---------------

* **Executable not found:** check ``PATH``, the installation path and
  ``ase2sprkkr config --show``.
* **MPI version was not compiled:** if automatic selection cannot use
  the available serial executables, explicitly turn MPI off::

     ase2sprkkr config --set profiles.default.mpi false

* **SCF did not converge:** inspect the output and adjust relevant
  numerical settings based on the SPR-KKR manual.
* **Migrating from older ASE2SPRKKR:** see
  :ref:`upgrading from 3.x <installation-upgrade>` and run
  ``ase2sprkkr migrate`` if you have legacy user configuration.

Next, see :doc:`usage` for task parameters, or :doc:`examples`
for more applications.

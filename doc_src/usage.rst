Using ASE2SPRKKR
================

The ``SPRKKR`` calculator prepares and runs native SPR-KKR tasks from
ASE structures or existing potentials. If you have not yet run a
calculation, follow the :doc:`quickstart` first. For physics and
convergence settings, consult the
`SPR-KKR manual <https://sprkkr.org/manual/>`_.

Tasks and input parameters
--------------------------

The default task is ``SCF``. Other tasks include ``DOS``, ``ARPES``,
``SPLEED``, ``BSFEK``, ``BSFKK`` and ``JXC``; some less established
tasks are marked experimental. Inspect the installed version to find
task names and available options::

   ase2sprkkr info --task
   ase2sprkkr info --task=SCF
   ase2sprkkr info --task=SCF.ENERGY

A calculator can receive a task name, an input-file path or an
``InputParameters`` object. For example:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import InputParameters, SPRKKR

   parameters = InputParameters.for_task("SCF")
   parameters.ENERGY.NE = 5

   atoms = bulk("Al", "fcc", a=4.05)
   calculator = SPRKKR(atoms=atoms, input_parameters=parameters)
   result = calculator.calculate()

   print("SCF energy (eV):", result.energy)
   print("Converged:", result.converged)

You can change an option on the calculator with
``calculator.set({"ENERGY.NE": 6})``, or override it for one run with
``calculator.calculate(options={"ENERGY.NE": 6})``.
A section-qualified name avoids ambiguity when different sections have
options with the same name.

To inspect an existing SPR-KKR input file without running a calculation:

.. code-block:: python

   from ase2sprkkr import InputParameters

   parameters = InputParameters.from_file("existing.inp")
   print(parameters.task_name)

Incomplete files can be inspected with
``InputParameters.from_file("existing.inp", allow_incomplete=True)``.
Input files do not contain all the physical data needed to run
a calculation; a potential or atoms must also be provided.
For task-definition internals and shared parameter schemas, see
:doc:`philosophy`.

Site occupations
----------------

ASE does not represent all SPR-KKR site properties. Convert an ASE
structure to ``SPRKKRAtoms`` to edit mixed site occupations:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import promote_ase_atoms

   atoms = promote_ase_atoms(bulk("LiCl", "rocksalt", a=5.64))
   atoms.sites[1].occupation.set({"Cl": 0.5, "I": 0.5})

Symmetry-equivalent atoms may share a site type. See :doc:`symmetry`
before modifying site equivalence or transforming a disordered cell.

Potentials and consecutive tasks
--------------------------------

An SPR-KKR ``.pot`` file contains atomic-site information and the
electronic potential used by native calculations. It is not a
pseudopotential library such as VASP's POTCAR.

Create a starting potential from an ASE structure, or load an existing
one:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import Potential

   potential = Potential.from_atoms(bulk("Al", "fcc", a=4.05))
   previous = Potential.from_file("existing.pot")

A converged SCF potential can be used for a subsequent DOS calculation:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import SPRKKR

   scf = SPRKKR(atoms=bulk("Al", "fcc", a=4.05)).calculate()
   if not scf.converged:
       raise RuntimeError("SCF did not converge")

   dos = SPRKKR(input_parameters="DOS").calculate(
       potential=scf.potential_filename
   )
   print(dos.dos_filename)

Passing the potential explicitly identifies which SCF run supplies it.
For the full first-calculation walkthrough, see :doc:`quickstart`.

To generate native input files **without executing SPR-KKR**:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import SPRKKR

   SPRKKR(atoms=bulk("Al", "fcc", a=4.05)).save_input(
       input_file="example.inp", potential_file="example.pot"
   )

Results and further workflows
-----------------------------

``calculate()`` returns a task-specific ``TaskResult``, including
parsed quantities and paths to output files. You can also read
recognized result files independently:

.. code-block:: python

   from ase2sprkkr import ResultFile

   dos_data = ResultFile.from_file("calculation.dos")

See :doc:`result-files` for plotting, Grace files and NeXus export,
:doc:`remote` for MPI and execution profiles, and :doc:`execution`
for persistent calculations and reconnection. Command-line commands
are collected in :doc:`cli`, with runnable examples in :doc:`examples`.

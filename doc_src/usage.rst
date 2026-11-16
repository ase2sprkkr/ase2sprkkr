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

Interactive parameter editing
-----------------------------

Input-parameter defaults belong to the calculator and can be inspected
and changed before running a task. For example:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import SPRKKR

   calculator = SPRKKR(atoms=bulk("Al", "fcc", a=4.05),
                       input_parameters="SCF")
   calculator.set({"ENERGY.NE": 12})
   print(calculator.get("ENERGY.NE"))

   # The same option can be accessed through its section:
   calculator.input_parameters.ENERGY.NE = 16
   print(calculator.input_parameters.ENERGY.NE())

   # This changes the option for just one calculation:
   result = calculator.calculate(options={"ENERGY.NE": 24})

The calculator retains its default value of 16 after this one-time
override. To select a different task or replace all parameter defaults,
assign to ``calculator.input_parameters`` (for example, ``"DOS"``,
an ``InputParameters`` instance, or an input-file path). The
``input_parameters`` argument to ``calculate()`` can also select
a task for a single calculation without replacing the defaults.
Use ``ase2sprkkr info --task=TASK`` to discover valid settings.

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

Site occupations and symmetry
-----------------------------

ASE structures can describe geometry, but SPR-KKR needs additional
site properties, including chemical occupations and atomic-type data.
``promote_ase_atoms()`` adds the ``sites`` interface to an ASE
``Atoms`` object. Promotion also happens when atoms are supplied to
the ``SPRKKR`` calculator.

Symmetry-equivalent sites can share a ``SiteType``. To change
**only one** site, first call ``break_symmetry()``; otherwise editing
a shared occupation raises an error:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import promote_ase_atoms

   atoms = promote_ase_atoms(
       bulk("LiCl", "rocksalt", a=5.64) * (2, 1, 1)
   )
   site = atoms.sites[3]  # one of the chlorine sites
   site.break_symmetry()
   site.occupation.set({"Cl": 0.5, "I": 0.5})

   # Adding Br at 20% reduces the existing fractions proportionally:
   site.occupation["Br"] = 0.2
   print(site.occupation["Br"])
   print(site.occupation.as_dict)

To modify **all** equivalent sites together, edit their shared
``site.site_type.occupation`` instead of breaking symmetry.
Occupation values are normalized; atomic types can be specified by
their chemical symbols or atomic numbers. Changes are synchronized to
ASE's ``atoms.info["occupancy"]`` representation using the
``spacegroup_kinds`` mapping, and existing ASE occupancy metadata is
recognized when atoms are promoted.

Different occupation fractions represent chemical disorder (such as
a substitutional alloy), not a larger explicitly ordered supercell.
The consequences of changing equivalence or converting a disordered
cell are discussed in :doc:`symmetry`.

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

When starting from atoms, ASE2SPRKKR prepares a potential containing
the site information needed for the native calculation. When starting
from a previous ``.pot`` file, the stored electronic state can be
reused. Input parameters specify the task and its controls; they do
**not** replace the structure and electronic data in the potential.
For that reason, loading an ``.inp`` file alone is insufficient
unless its referenced potential can also be found.

A calculator owns default atoms or a potential, along with default
input parameters. Explicit ``atoms``, ``potential`` and
``input_parameters`` arguments to ``calculate()`` select the
physical state and task for that run. Do not supply both an ``Atoms``
object and a ``Potential`` object to the calculator constructor;
the latter already contains the atomic structure. During preparation,
input parameters are copied before per-run options and generated paths
are inserted, but atoms and potential objects are mutable working
state: promoting atoms, adding empty spheres, or synchronizing the
potential may modify them. Once prepared, a run uses its own
materialized files.

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
Alternatively, the SCF result provides a calculator associated with
the new potential, suitable for the next task:

.. code-block:: python

   dos = scf.calculator.calculate(input_parameters="DOS")

A successful synchronous SCF calculation can also update the
original calculator's potential when it uses that calculator's
default physical state. If a run overrides its atoms or potential,
do not assume those one-time changes have become calculator defaults.
Use the result's potential or ``scf.calculator`` to choose the
state explicitly.

For the full first-calculation walkthrough, see :doc:`quickstart`.

To generate native input files **without executing SPR-KKR**:

.. code-block:: python

   from ase.build import bulk
   from ase2sprkkr import SPRKKR

   SPRKKR(atoms=bulk("Al", "fcc", a=4.05)).save_input(
       input_file="example.inp", potential_file="example.pot"
   )

``save_input()`` returns the prepared input parameters and the paths
of the input, potential and output files. The destination directory
must be persistent; ``directory=False`` is only for actual runs
using ``create_run()`` or ``calculate()``.

Filename templates
------------------

The calculator's ``input_file``, ``potential_file`` and
``output_file`` names can contain placeholders that are resolved
while preparing a calculation:

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Placeholder
     - Replacement
   * - ``%a``
     - Chemical symbols / structure label derived from the atoms
   * - ``%t``, ``%T``
     - Task name, as supplied and in uppercase respectively
   * - ``%l``
     - Calculator label (empty if none was supplied)
   * - ``%d``
     - Local date and time (format ``YYYY-MM-DD_HH:MM``)
   * - ``%c``
     - Incrementing calculation counter
   * - ``%_``
     - Underscore separator: omitted at the start of a filename;
       consecutive ``%_`` separators are collapsed

For example:

.. code-block:: python

   calculator = SPRKKR(
       atoms=bulk("Al", "fcc", a=4.05),
       label="test",
       input_file="%l%_%a_%T.inp",
       potential_file="%a.pot",
   )
   _, input_name, potential_name, _ = calculator.save_input()
   print(input_name, potential_name)

The values are resolved consistently for each prepared calculation,
including when multiple filenames use the same placeholder.
The ``%d`` format contains a colon, which is not portable in
Windows filenames; choose another template there.

Inspecting SCF results
----------------------

The SCF result contains parsed iterations, convergence information
and access to the newly written potential. After a calculation:

.. code-block:: python

   scf = SPRKKR(atoms=bulk("Al", "fcc", a=4.05)).calculate()
   print("Iterations:", len(scf.iterations))
   print("Converged:", scf.converged)
   print("Final convergence error:", scf.last_iteration.error())
   print("Final spin moment:", scf.last_iteration.moment.spin())
   print("Total energies:", scf.iterations.values_of(("energy", "ETOT")))

   print("Output potential:", scf.potential_filename)
   print("Charge data:", scf.potential.CHARGE())

The iteration data can be used to monitor convergence or to
postprocess many calculations. The final error and magnetic
moments come from the native output; whether the result is physically
converged must be checked separately from successful process
completion. Plotting the error history (with the optional ``plot``
dependencies installed) is also possible through
``scf.plot("error")``.

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

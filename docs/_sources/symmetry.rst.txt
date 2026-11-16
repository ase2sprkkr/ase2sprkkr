Crystal symmetry and cell transformations
==========================================

The ``symmetry`` command uses spglib to provide three distinct operations:

* Analysis identifies the space group and leaves the structure unchanged.
* Cell conversion finds a primitive or standardized conventional cell.
* Idealization projects lattice metrics and atom positions onto the detected
  symmetry, preserving the original basis and atom ordering.

These operations can be requested independently or combined.

Examples
--------

.. code-block:: console

   ase2sprkkr symmetry input.cif
   ase2sprkkr symmetry input.cif --primitive -o primitive.cif
   ase2sprkkr symmetry input.cif --standardize -o conventional.cif
   ase2sprkkr symmetry input.cif --idealize -o idealized.cif
   ase2sprkkr symmetry input.cif --primitive --idealize -o primitive_ideal.cif

Use ``--equivalence occupation`` when chemically equivalent sites
should be compared by occupation rather than their current site types.

CIF-declared symmetry
---------------------

For CIF input, the command reports both the declared space group
(if present) and the symmetry calculated from the atomic positions.
Different space-group numbers are highlighted. Use ``--operations``
to list symmetry operations, or save a JSON report::

   ase2sprkkr symmetry input.cif --report report.json

Python API:

.. code-block:: python

   from ase.io import read
   from ase2sprkkr.sprkkr.symmetry import analyze_symmetry

   atoms = read("input.cif")
   analysis = analyze_symmetry(atoms)
   primitive = analysis.find_primitive(atoms)
   conventional = analysis.standardize(atoms)
   idealized_in_original_cell = analysis.idealize(atoms)
   idealized_primitive = analysis.find_primitive(atoms, idealize=True)

The default ``site-types`` equivalence distinguishes SPR-KKR site types
and magnetic information. Region identity is preserved when finding
equivalent sites. The ``occupation`` mode can merge sites with identical
chemical occupations even when their electronic site-type data differ;
recheck the resulting electronic data before further calculations.

Extra ASE arrays use one of three transfer policies when converting cells:

* require-equal (default): Values of atoms in the same primitive translation class must agree.
* discard: Do not transfer additional arrays.
* choose-first: Select the first source representative of each primitive site, and warn
  when values differ.

Extra arrays are checked when a cell is converted, and tensor/vector
frame transformations are not performed automatically. Idealization by
itself retains atom ordering and electronic data. If the geometry changes,
a new SCF calculation is recommended before relying on an old potential.

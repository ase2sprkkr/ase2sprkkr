Crystal symmetry and cell transformations
==========================================

New crystal symmetry workflows use spglib to distinguish three operations:

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

Python API (analysis retains no reference to the input Atoms):

.. code-block:: python

   from ase2sprkkr.sprkkr.symmetry import analyze_symmetry

   analysis = analyze_symmetry(atoms)
   primitive = analysis.find_primitive(atoms)
   conventional = analysis.standardize(atoms)
   idealized_in_original_cell = analysis.idealize(atoms)
   idealized_primitive = analysis.find_primitive(atoms, idealize=True)

By default, site-types are distinguished by SPR-KKR SiteType identity,
existing site magnetisation, and ASE initial magnetic moments.
Region membership is part of the site identity for SpacegroupInfo:
atoms in different regions are not merged even when otherwise equivalent.
Occupation equivalence explicitly permits merging different SiteTypes that
have the same chemical composition; unselected electronic data may be lost.

Extra ASE arrays use one of three transfer policies when converting cells:

* require-equal (default): Values of atoms in the same primitive translation class must agree.
* discard: Do not transfer additional arrays.
* choose-first: Select the first source representative of each primitive site, and warn
  when values differ.

Different primitive sites remain distinct even when they have the same
species tag or SiteType. Tensor and vector frame transformations are not
handled by this first version.

Idealization alone preserves additional arrays and site order. Radial
potentials, charges, and moments are also preserved. If idealization changes
geometry with electronic data present, a warning recommends a new SCF run.
Cell conversion retains unrelated ASE info metadata while rebuilding
spacegroup information and occupancy records for the new geometry.

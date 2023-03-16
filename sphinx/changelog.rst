Changelog
=========

Version 1.0.6
-------------

* Sites data moved to ASE arrays to allow merging two Atoms structures
* Testing switched to pytest
* ARPES task fixies

Version 1.0.7
-------------
* a2s_visualise_in_struct script fixed
* a2s_visualise_in_struct accepts scale-radii argument to control the size of visualised atoms
* make now by default install the ase2sprkkr even if the version number have not been changed

Version 2.0.0-beta1
-------------------

User interface
~~~~~~~~~~~~~~
* Runtime documentation available.
* Runtime documentation is added to the docstring and to the generated documentation.
  (so far for input parameters).
* Class names refactored - abuse of 'BaseSomething' names solved.
* Dangerous values (that do not pass the validity checks) are allowed.

Architecture changes
~~~~~~~~~~~~~~~~~~~~
* Allow the Keywords arguments to accept descirptions of the keywords.
* Complex GrammarType were added.
* Option has the result attribute, that can hold the processed value of an user input
* Possibility to add 'expert' values to a configuration definition. The expert
  values are outputed only if they differ from the defaults.
* Expert sections have been introduced. They are printed out only if there is any changes (from defaults).
* The ARPES task have been documented and more options have been added.
* The SCF task have been documented and more options have been added.
* Default GrammarType for bool default values in InputParameters is now Flag.
* Numbered arrays have been introduced to allow options like MDIR, MDIR1, MDIR2 etc...
* Python 3.11 support added.


Version 2.0.0-beta2
-------------------

User interface
~~~~~~~~~~~~~~
* es_finder integration for empty spheres finding
* Support for 2D problems
* Routines for building 2D problem
* change_task method for InputParameters
* calculate(..., directory=False) runs the calculation in a temporary directory

Internals
~~~~~~~~~
* Sections validation
* LatticeData class refactored
* Brackets in value names are allowed

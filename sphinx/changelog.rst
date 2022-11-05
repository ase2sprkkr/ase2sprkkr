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

Unpublished
-----------

User interface
~~~~~~~~~~~~~~
* Runtime documentation available.
* Runtime documentation is added to the docstring and to the generated documentation.
  (so far for input parameters).
* Class names refactored - abuse of 'BaseSomething' names solved.

Architecture changes
~~~~~~~~~~~~~~~~~~~~
* Allow Keywords arguments to accept descirptions of keywords.
* Complex GrammarType were added.
* Option has the result attribute, that can hold the processed value of user input
* Possibility to add 'expert' values to configuration definition. Expert
  values are outputed only if they differ from the defaults.
* Expert sections have been introduced. They are printed out only if there is any changed (from default)
* ARPES task have been documented and more options have been added
* SCF task have been documented and more options have been added
* Default GrammarType for bool default values in InputParameters is now Flag
* Numbered arrays have been introduced to allow options like MDIR, MDIR1, MDIR2 etc...

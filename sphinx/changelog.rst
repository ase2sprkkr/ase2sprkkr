Changelog
=========

User interface
~~~~~~~~~~~~~~
* The JXC task have been added
* ARPES configuration improved
* The possibility to define own tasks or modify the defaults of the tasks for the user have been added
* Parsing of the task outputs improved
* Energies in SCF output refactored
* Parsing of POTENTIAL and CHARGE sections of potential file
* ProcessOutputReader.read_from_file can be called as classmethod
* Occupancies can now have string indices (since some other codes produce such occupancies)
* Partial occupancy is not upnormalized, vacum pseudoatom is added instead

Internals
~~~~~~~~~
* NumpyArray grammar type improved
* Some fixies and code lininting
* Repeated configuration allowed - for parsing output files (e.g. iterations of SCF cycle)
* The possibility to emit a warning to the options was added

Refactoring
~~~~~~~~~~~
* Output files moved to the separate directory
* Configuration definitions splitted to several files

Version 2.2.1
-------------
* Fix of ConfigurationContainer.set_values

Version 2.2.0-beta
------------------
* a2s_plot_output can handle DOS and BSF output files
* ARPES task fix
* Arithmetic can be done with DOS result
* executable_postfix argument of calculator renamed to executable_suffix to make it consistent with the name of the environment variable
* User-defined input parameters for repeatedly used task

Internals
~~~~~~~~~
* Generated type improvement
* Switch grammar element: format of a parsed file can depend on the previously parsed values
* Gather grammar element for ``NAME1 NAME2 = VALUE1 VALUE2`` syntax
* Routines for plotting the results are (i hope) stabilized
* Various small improvements and fixies

Version 2.1.1
-------------

User interface
~~~~~~~~~~~~~~
* ARPES and SCF task definition improved
* ARPES SPC results can be parsed and plotted
* DOS results parsing
* FULLPOT mode for SCF calculation
* a2s_plot_output script to plot SPC results
* Better naming of input and potential files
* Gilbert TASK added (experimental, not tested)
* input_parameters.change_task method fixed
* [] access to array options/values (no need for VARIABLE()[] notation)
* Numbered arrays (e.g. CONTROL.MDIR) can be set using arrays
* Better formating of input parameters
* Some minor tweaks and corrections of input parameters

Internals
~~~~~~~~~
* GrammmarTypes refactorized (splitted to more files)
* NumpyArray and RestOfTheFile grammar types for output files
* Generated grammar types and values for easy access to output files
* Calculator.save_input refactored

Version 2.0.4
-------------
* Some fixies
* ASR repcipies available as ase2sprkkr.asr subpackage

Version 2.0.1
-------------
* Fix of sys-file generation
* Hastily written ASE2SPRKKR slides included

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


Version 1.0.7
-------------
* a2s_visualise_in_struct script fixed
* a2s_visualise_in_struct accepts scale-radii argument to control the size of visualised atoms
* make now by default install the ase2sprkkr even if the version number have not been changed


Version 1.0.6
-------------

* Sites data moved to ASE arrays to allow merging two Atoms structures
* Testing switched to pytest
* ARPES task fixies





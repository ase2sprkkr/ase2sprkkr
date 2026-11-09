Changelog
=========


Version 4.0.0-beta1
-------------------


Execution, remote calculations and projects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Managed local and remote calculations with asynchronous execution,
  status monitoring, output following, cancellation and persistent run history.
* The independent ``kkrunner`` package provides the execution layer.
  ASE2SPRKKR prepares SPR-KKR inputs and interprets scientific outputs;
  KKRunner manages processes, MPI, SSH, schedulers, file transfers and run state.
* ``SPRKKR.create_run()`` separates preparation, submission and waiting
  (``submit()`` and ``wait()``) while ``SPRKKR.calculate()`` remains available
  for synchronous calculations.
* Host-aware MPI launchers, including ``srun`` in Slurm allocations
  and ``mpirun``/``mpiexec`` elsewhere.
* Remote execution over SSH, with remote working directories, input
  staging, output retrieval, configurable synchronization policies,
  and remote executable and MPI configuration.
* Optional PSI/J scheduler integration, including Slurm, via the
  ``psij`` extra, with scheduler-specific resources and backend behavior.
* Persistent project-based calculations with durable run identities,
  statuses, execution records and calculation history; tracking is opt-in.
* Reconnection to active jobs after client restart, resumption of
  interrupted downloads and explicit postprocessing retries without
  resubmitting SPR-KKR. Recovery depends on backend capabilities.
* Historical ``TaskResult`` reconstruction using ``TaskResult.from_run()``
  or ``load_project_result()``; ``attach_project_run()`` reconnects active
  persistent project runs.
* Declarative task input/output/cache file manifests, including
  optional files and wildcard patterns, synchronization policies and
  integrity checks to detect missing or overwritten historical results.
* Managed stdout/stderr and output listeners. ``ase2sprkkr remote`` provides
  user-facing status, follow and cancel commands; the execution protocol
  uses ``kkrunner remote``.

Configuration and profiles
~~~~~~~~~~~~~~~~~~~~~~~~~~
* YAML user configuration and named execution profiles for executable
  locations, MPI, local/SSH/scheduler backends, remote directories,
  filesystem identities and file synchronization settings.
* Migration from legacy Python configuration with
  ``ase2sprkkr config --migrate``; extended CLI to list, inspect and
  edit profiles and their settings.
* Profiles and per-calculation overrides can select execution
  backends without changing task input parameters.

New and extended SPR-KKR tasks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* More than 40 new native task input definitions, many experimental.
  ``ase2sprkkr info`` lists stable tasks before experimental ones;
  definitions do not imply numerical validation of every advanced mode.
* Photoemission, spectroscopy and related tasks: ``AES``, ``APS``,
  ``NRAES``, ``AIPES``, ``BAND``, ``SPLEED``, ``CLPES``, ``CLXPS``,
  ``VBPES``, ``XES``, ``XRS`` and ``CHIXAS``.
* Magnetic response, optical and transport tasks: ``CHI``,
  ``CHILANDAU``, ``SIGMA``, ``OPTABS``, ``T1``, ``NUCSSCPL``,
  ``STONER``, ``SOCPAR``, ``FMAG`` and ``COMPTON``.
* Electronic-structure, scattering and current-density tasks:
  ``BONDING``, ``BSDOSIJ``, ``FSCAT``, ``EKREL``, ``JALF``,
  ``JORB`` and ``TDOBS``.
* Positron workflows: ``POSANIPREP``, ``POSANITAU`` and ``POSANI``.
* Plotting and diagnostic tasks: ``PLOTPRS``, ``PLOTPOT``,
  ``PLOTRHO``, ``PLOTSFN``, ``PLOT2DPR``, ``PLOT2DPOT``,
  ``PLOT2DRHO``, ``WFPLOT``, ``MECHECK``, ``MEPLOT``,
  ``PSHIFT`` and ``TUTORIAL``.
* BSF modes have separate ``BSFEK`` and ``BSFKK`` task identities
  sharing the BSF input schema; BSF validation and k-paths improved.
* Existing ``SCF``, ``DOS``, ``ARPES``, ``PHAGEN``, ``JXC``, ``XAS``,
  ``GILBERT`` and other definitions received improved options,
  default values, validation or declared files.
* Task identity and executable metadata separated from reusable
  input parameter schemas; user-defined presets preserve canonical
  task identities instead of introducing unwanted new tasks.
* Task-specific input, output and cache declarations, including
  conditional filenames, ``DATASET``/``ADSI`` and
  ``TAUDIR``/``TAUPATH`` handling.
* Support for reading incomplete input files with warnings,
  conditional requiredness, and cross-section validation of
  dependent options and array lengths.

Symmetry and structure tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* New ``ase2sprkkr symmetry`` command and ``analyze_symmetry()`` API
  for spglib-based space-group and site-symmetry analysis without
  modifying or retaining the input Atoms object.
* Primitive and standardized conventional cell transformations,
  optional idealization, and policies for preserving SPR-KKR site
  metadata and extra ASE arrays.
* Equivalent sites, Wyckoff positions and symmetry operations
  can be displayed, with optional JSON report output.
* CIF-declared symmetry is reported separately from calculated
  symmetry, including detection of differing space-group numbers.

Results and plotting
~~~~~~~~~~~~~~~~~~~~
* Structured result files exposed as ``ResultFile`` from ``ase2sprkkr``
  and ``ase2sprkkr.result_files``. The top-level ``OutputFile``
  remains available as a deprecated compatibility alias.
* SFN file reading, writing and plotting, explicit ``SFNFIL``
  selection, and information about SFN file generation.
* SCF parsing supports runs without ``EMIN`` (thanks to Chriss Woodgate).
* ARPES results provide all polarization spectra through
  ``ArpesResult.spectra``; NeXus/HDF5 export with metadata
  and interactive ARPES plotting controls.
* ``print_output`` supports informational output, regular expressions
  and custom output listeners.
* Corrected DOS, BSF and ARPES arithmetic, DOS energy handling,
  and Gaussian and Lorentzian spectral broadening.

Calculator and command-line interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* ``SPRKKR.prepare()`` generates scientific calculation inputs without
  submitting a run.
* Explicit ``restart_scf`` option to control SCF restarts.
* Unsuccessful synchronous calculations raise ASE's
  ``CalculationFailed`` instead of updating calculator state as if
  the calculation had succeeded.
* Lightweight ``ase2sprkkr_cli`` package now implements the CLI entry point,
  with lazy imports of command-specific functionality.

Compatibility
~~~~~~~~~~~~~
* Python 3.8 support dropped; Python 3.9 or newer required.
* ``kkrunner`` is a required dependency; scheduler integration uses
  the optional ``psij`` extra.
* ``ArpesResult.spc`` deprecated in favor of ``ArpesResult.spectra``.
* Legacy execution-related configuration options are deprecated
  in favor of named profiles.

Internals
~~~~~~~~~
* Transactional configuration updates and cross-section validation.
* Refactored dependent, repeated and generated configuration values,
  including NumPy arrays with specified or dependent line counts.
* Optional fast parsing of simple numerical arrays via ``numpy.fromstring``,
  with fallback for incompatible input formats.
* Reorganized task and result-reader registries, with SPR-KKR
  file metadata and result reconstruction in ASE2SPRKKR and
  generic execution state and persistence delegated to KKRunner.
* Improved symmetry/site-equivalence handling and regression tests.
* Various bug fixes and compatibility improvements.

Version 3.6.0
-------------
* Added support for ANDERSON mixing (experimental), BSF task improved
* A better debugging in `ase2sprkkr shell`
* `promote_ase_atoms` can be directly imported from the `ase2sprkkr` root package
* Changing properties via `site` is disallowed in cases where `site_type` is shared
* JXC output file parsing and plotting
* XAS output file parsing and plotting
* Beyond-DFT parammeters (DFTM, DFT+U,...)
* Hyperfine-field switches
* Semi-infinite structure building improved
* Uppasd improvements: direct calling of `sd`, input file stub writing
* Improved layout for multiple plots
* Additional values parsed for SFC task
* Tool: Python-like output for `ase2sprkkr` k-path tool
* Tool: `ase2sprkkr` example search (--regex) and example display
* Tool: `ase2sprkkr shell --pdb` for running examples in the Python debugger
* Tool: `ase2sprkkr plot -V <filename>` for showing plottable values in a file
* Improved Hall number determination from Cell
* Docstrings added for examples
* Various bug fixes
* SLURM work manager auto-detection
* Examples from ASE2SPRKKR article were added
* Runtime dependencies were revisited

Internals
~~~~~~~~~
* Adoption of Ruff: code formatting and cleanup
* Picklability of Atom objects (for better ase.view compatibility)
* Refactoring of ProcessOutput and ProcessReaders
* Actions for Outputs and OutputFiles (for GUI support)
* BaseCalculator is used as a base class for SPRKKRCalculator
* Improvements for GUI:
  * ResultValues are implemented as Options, to reuse existing functionality in GUI
  * Better result files handling
  * Files in results provide both relative and absolute path

Version 3.4
-----------

* Magnetization moments parse
* Units uses unyts for better units support (retaining the unit user used)
* Empty spheres uses spglib for symmetry (which solves some cases that have not been computed correctly )
* ASE read/write compatibility
* Allow to run computation inside Jupyter notebook
* promote_ase_atoms exported
* BSF and SCP outputs improved
* Calculator do not longer overwrite the potential file, if it is given in constructor
* Symmetry handling improved
* Some fixies and enhancements for plotting
* Allow async running of the process
* Python 3.14 compatibility
* MDIR option has no default value
* SPEC task parser
* Dij output file parser

Version 3.3.4
------------------
* Sites are update according to the symbols during save
* XAS task
* Torque Task
* Empty-spheres finding improved
* MALF, MBET and MGAM options
* BSF-SPN file added
* K-path generating tool for BSF task
* Layer and fermi argument to BSF plotting
* --separate_plots argument for plotting routines
* Uppasd task for generating Uppasd input files
* Various fixies


Version 3.3.3
------------------
* Added mesh parameter to empty spheres finding
* Spacegroup from ASE is no longer used, since it have been deprecated

Version 3.3.2
------------------
* Buildsystem changed to scikit-build-core (mesonpy do not support editable build at least in some setups)
* Atoms rotating, shifting and flipping routines
* Config file printing, value-setting etc. routines from cmdline

Version 3.2.5-beta
------------------

User interface
~~~~~~~~~~~~~~
* Better data validation, warnings on data validity
* Bash ase2sprkkr commandline completing
* Faster package loading, lazy loading of some modules/data definitions
* Task BSF improved
* Python 3.7, pyparsing 2.x, and numpy < 1.23 support dropped
* Nomad uploading proof of concept
* Configuration file improved, help available
* Occupancy and symmetry handling fixed
* show_structure tools a bit improved (it uses temp files by default)
* Radial charge and potentials for multiple occupancy and FULLPOT
* EMIN/EMAX can be always specified as EMINEV/EMAXEV

For developers
~~~~~~~~~~~~~~
* Switch element improved to allow implementation of if_defined element
* Repeated arrays
* Length element
* Profiling

Version 3.1.0-beta
------------------
* BSF task and result
* Examples in package and documentation

Version 3.0.2-beta
------------------
* Host Madelung potential and Charge moments parsing
* Sites ASE array
* Various fixies
* Test tool to run tests
* has_symmetry and break_symmetry methods of Site and SiteType

Version 3.0.0-beta
------------------

User interface
~~~~~~~~~~~~~~
* Empty spheres finding has been added
* SiteType class added: sites are now always distinct, site_types can be shared between sites
  due to symmetry
* The JXC task have been added
* Total DOS for DOS plotting added
* DOS improved (arithmetics, access by atom type)
* ARPES configuration improved
* The possibility to define own tasks or modify the defaults of the tasks for the user have been added
* Parsing of the task outputs improved
* Energies in SCF output refactored
* Parsing of POTENTIAL and CHARGE sections of potential file
* ProcessOutputReader.read_from_file can be called as classmethod
* Occupancies can now have string indices (since some other codes produce such occupancies)
* Partial occupancy is not upnormalized, vacum pseudoatom is added instead
* Energies in SCF result are returned in more uniform way
* Writing of modified output files is possible
* KRMT and KRWS configuration values added

Tools
~~~~~
* Tools (``a2s_...``) have been integrated into one ase2sprkkr tool
* Test tool introduced to test the ase2sprkkr installation
* Run tool for running the calculations using prepared .pot files
* Various tools and plotting fixies

Internals
~~~~~~~~~
* NumpyArray grammar type improved
* Some fixies and code lininting
* Repeated configuration allowed - for parsing output files (e.g. iterations of SCF cycle)
* The possibility to emit a warning to the options was added
* Build backend changed to Meson to allow build cython empty spheres
* Testing fully switched to pytest, the unittest dependency has been dropped

Refactoring
~~~~~~~~~~~
* Output files moved to a separate directory
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
* Plotting improved
* Better help values for output files

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

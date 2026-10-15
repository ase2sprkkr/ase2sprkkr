ASE2SPR-KKR's philosophy for the developers
============================================

The base philosophy of the package
----------------------------------
The whole package is built around configuration objects.
A configuration object is an object, that can both read and write configuration files.
Each configuration object has the capabilities to both parse the configuration file 
(e.g. input file for SPR-KKR task or a potential file) and to write the file.

The configuration object has a tree structure: it is divided by section, while each section
can contain values. (However, some (mostly custom: i.e. e.g. sections of a potential file, that
are not currently parsed and their value is given as is) sections are just containers for a single value.

The leaves of the configuration object are options
(:class:`ase2sprkkr.common.options.Option`)
- single value containers. They are described by
configuration definitions
(see :class:`ase2sprkkr.common.configuration_definitions.ValueDefinition`),
which hold the required type of the value (see :mod:`ase2sprkkr.common.grammar_types`) and its few other properties, such as the default value or parameters that determine the format of the value in the parsed/writed configuration file.

The non-leaves nodes of a configuration trees are containers (:class:`ase2sprkkr.common.configuration_containers.ConfigurationContainer`). Their content and the way how they are parsed from file and written to file is again determined by the corresponding definitions (:class:`ase2sprkkr.common.configuration_definitions.ContainerDefinition`), thus, the definition form a tree with the same structure as the configuration object (so, there is many to one relationship between configuration (tree) object ad its definition, as for each "definition tree" can be instantiated many containers for various contents of a configuration file).

Structure of the package
--------------------------------------
There are the following subpackages in the ase2sprkkr package:

:mod:`common <ase2sprkkr.common>`
"""""""""""""""""""""""""""""""""
There are the common shared modules: the options, containers, their definitions and grammar types and miscellaneous stuff.

:mod:`input_parameters <ase2sprkkr.input_parameters>`
""""""""""""""""""""""""""""""""""""""""""""""""""""""
There are the specializations of the configuration objects and definitions for the input files.
In the (:mod:`ase2sprkkr.input_parameters.definitions`) subpackage are the definitions of the
possible tasks, that can be run by SPR-KKR.

:mod:`potentials<ase2sprkkr.potentials>`
""""""""""""""""""""""""""""""""""""""""""""""""""""
There are specializations for the potential file definition.
There is currently just one definition of a potential
in (:mod:`ase2sprkkr.potentials.definitions.potential`)
(there is just
one format of the potential file, there is no support for creating
different versions of potential in the 1.0 version of the package).
The definitions of the potential sections are in
(:mod:`ase2sprkkr.potentials.definitions.sections`).

:mod:`outputs<ase2sprkkr.outputs>`
""""""""""""""""""""""""""""""""""""
Each input file also prescribes the type of the reader: object, that
reads and possibly parses the output of the SPR-KKR. These readers
(located at :mod:`ase2sprkkr.outputs.readers`) are
just simply objects, that read the output line by line, searching
for a regex pattern and parsing selected lines of the file (skipping
the rest) to obtain some values (e.g. convergence results from iterations).
The parsed values are then returned as the result of the calculations.

In order to be capable of parsing both standard output and error
output stream of SPR-KKR, asyncio is employed for the readers.

In version 1.0 of the package, only the output reader of SCF task
is implemented, results of the other tasks are read by the dummy :class:`default reader<ase2sprkkr.outputs.readers.default.DefaultOutputReader>`, which does not parse anything
from the output.

Note, that some parts of the output can and are parsed using the generated
grammar, as it is described above, see e.g. :func:`ase2sprkkr.outputs.readers.scf.ScfOutputReader.read_output`.

:mod:`sprkkr<ase2sprkkr.sprkkr>`
"""""""""""""""""""""""""""""""""""
ASE calculator for SPR-KKR, and data structures for holding the properties,
that are special for SPR-KKR.

To allow ASE Atoms object (of the ASE Atoms class) to hold the SPR-KKR specific
properties, it has to change its (OOP) ancestor to the class
(:class:`ase2sprkkr.sprkkr.sprkkr_atoms.SPRKKRAtoms`).
This is done by :func:`ase2sprkkr.sprkkr_atoms.SPRKKRAtoms.promote_ase_atoms`
method, which is either called automatically when any Atoms object is passed to
ase2sprkkr routines, or can be called manually if it is needed.

:mod:`tools<ase2sprkkr.tools>`
"""""""""""""""""""""""""""""""
Scripts for postprocessing of the SPRKRR results.


Reading the configuration files
-------------------------------------------
Configuration definition utilizes `pyparsing<https://pyparsing-docs.readthedocs.io/en/latest/>`
for creating the grammar, that parses the configuration file.

Thus, the reading of a file has two phases. First, the definition parses the content of the file and created the tree of the python values (dictionaries,
numpy arrays, etc.). In the second step, the tree is transformed, according to the definitions, to the configuration object.
If a potential file is parsed, then the third phase occurs: the ASE Atoms object is set up accordingly to the values in the configuration object.

The grammar is defined accordingly to the constants and/or class methods of the configuration definitions descendants, that customizes the shared base objects, see e.g.
:mod:`ase2sprkkr.input_parameters.input_prameters_definitions` or
:mod:`ase2sprkkr.potentials.potential_definitions`

Writing the configuration files
--------------------------------------------
The writing to the file is performed by save_to_file :func:`ase2sprkkr.common.configuration_containers.ConfigurationContainer.save_to_file`, which is actually implemented in its descendants. The Options are written accordingly to their definitions and types, see the methods :func:`ase2sprkkr.common.configuration_definitions.ValueDefinition.write` and :func:`ase2sprkkr.common.grammar_types.GrammarType.write`.

Again, the exact way how to write the file - i.e. how to separate name-value pairs, how to separate values each from others etc. is given in the configuration_definition descendants.


Running the program
---------------------------------------------
When SPR-KKR calculation is to be invoked, the calculator:
 * creates or updates the potential from the Atoms object (however, it is possible to run calculation just
   according to the given potential, without an Atoms object)
 * (if the task is given by its name) creates the InputParameters object
 * saves the potential file and the input file
 * runs the executable (given by the definition object of the InputParameters)
 * creates the OutputReader (again, the type is given by the definition object of the InputParameters)
 * let the OutputReader read and parse the results from the output of the run executable

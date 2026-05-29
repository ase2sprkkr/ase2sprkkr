How to use ASE2SPRKKR
**********************


Installation
------------

To install the package using either pip:

::

   pip install ase2sprkkr

or conda package manager:

::

   conda install -c ase2sprkkr ase2sprkkr

| 
| If you want to use source codes from GitHub repository

::

   git clone https://github.com/ase2sprkkr/ase2sprkkr.git

you can install the package using

::

   install.sh

You are encouraged to run the test after installation, using

::

   make test


The package requires, that the SPR-KKR is installed (and thus
their executables are located somewhere in the path). If you
use a build of SPR-KKR with an executables suffix, please set
the environment variable ``SPRKKR_EXECUTABLE_SUFFIX`` accordingly,
e.g. if you use bash as your shell and the name of the "main"
SPR-KKR executable is ``kkrscf8.6``,
place into your ``.bash.rc`` the following line:

::

    export SPRKKR_EXECUTABLE_SUFFIX=8.6

For more verbose documentation about the installation the package, please see
the topic :ref:`More about the package installation<Installation>`


Usage
-----

The package allows us to define the problem to be computed using Python
and ASE. The current version of the package is aimed primarily at
defining the problem. The capabilities of obtaining the results and
their postprocessing are still limited: it is assumed, that the user of
the package is experienced with SPR-KKR package and is capable to
analyze the resulting output of SPR-KKR. Bellow, you can see several
basic examples of how to use the package.

Computing a bulk material
~~~~~~~~~~~~~~~~~~~~~~~~~

The basic usage of the package is to use ASE to define the computed
structure and, using SPR-KKR calculator, run the calculation:

::

   from ase.build import bulk
   from ase2sprkkr.sprkkr.calculator import SPRKKR

   atoms = bulk('Li')
   calculator = SPRKKR(atoms=atoms)
   calculator.calculate()

Calculation requires to have SPR-KKR in the PATH. In some cases, SPR-KKR
is compiled with a version (and/or hostname) suffix. In this case,
please define environmental variable ``SPRKKR_COMMAND_SUFFIX`` if the
SPR-KKR executable for SCF task has the name ``kkrscf8.6``, please run
(the example is for the bash shell):

::

   export SPRKKR_EXECUTABLE_SUFFIX=8.6

Defining the material properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the properties of the material can be defined as is described in
the `ASE documentation. <https://wiki.fysik.dtu.dk/ase/>`__ However, ASE
do not support occupation of the sites. To define it, you most promote
the atoms object to the ASE2SPRKKR one. It is done automatically, if you
pass the atoms object into the SPRKKR calculator constructor, however,
you can do it manually, if needed.

::

   from ase.build import bulk
   from ase2sprkkr.sprkkr.sprkkr_atoms import SPRKKRAtoms
   atoms=bulk('LiCl', 'rocksalt', a=5.64) * (2, 1, 1)
   atoms = SPRKKRAtoms.promote_ase_atoms(atoms)

Promoted atoms have the sites attribute, which is an array of the sites
in the primitive cell. If two sites are indistinguishable, they share
the sites object --- this is computed automatically according to the
symmetry of the structure (see the SPRKKRAtoms.compute_sites_symmetry
method). Sites object in the sites array holds all the SPRKKR specific
properties of the sites. You can change the occupation of a given site
by calling the :meth:`Site.occupation.set<ase2sprkkr.sprkkr.occupations.Occupation.set>` method:

::

   atoms.sites[3].occupation.set({'Cl':0.5, 'I' : 0.5 })

Atomic type can be given both by a chemical symbol or atomic number. You
can also just add an atomic type to the site:

::

   atoms.sites[3].occupation['Br'] = 0.2

-- the occupations of the currently presented atomic types are lowered
accordingly --- or just query the occupation of a given atomic type:

::

   atoms.sites[3].occupation['Br']

Some codes stores the occupation information in
``atoms.info['occupancy']``. To allow interoperability with these codes
the occupancy is readed from this array when
meth:`SPRKKRAtoms.promote_ase_atoms<ase2sprkkr.sprkkr.sprkkr_atoms.SPRKKRAtoms.promote_ase_atoms>` is called and the
``atoms.info['occupancy']`` is updated each time the
``atoms.sites[x].occupation`` is changed.

Reading the results
~~~~~~~~~~~~~~~~~~~

Some of the results of the computation are returned by the calculate
method (see e.g.
:class:`ase2sprkkr.outputs.readers.scf.ScfResult`.

::

   out = calculator.calculate()
   print(f"Energy is {out.energy}")
   print(f"There were {len(out.iterations}) iterations")
   print(f"The convergence error is {out.iterations[-1]['error']}")
   print(f"The moments are {out.last_iteration['moment']}")

For the SCF task, you can request the newly-created potential file,
either by filename, or by a Potential object, from which you can read
the other results in a text format (you are encouraged to contribute to
ase2sprkkr to define the format of not-yet-implemented sections of
potential).

::

   print(f"Potential has been saved to: {out.potential_filename}")
   potential = out.potential
   print(potential.CHARGE())

Setting the task type and input parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SPRKKR allows users to compute different tasks and each task can receive
many input parameters. The default task (if none is specified) is SCF:
to do a self-consistent cycle to compute the wavefunctions. You can
either choose the task and its parameters either during creating the
calculator via ``input_parameters`` argument:

::

   calculator = SPRKKR(atoms=atoms,
                       input_parameters='PHAGEN',
                       options={'NE': 5})

or one-time in the ``calculate`` methods:

::

   calculator.calculate(input_parameters='PHAGEN', options={'NE': 5})

The ``input_parameters`` argument accepts either a name of one of the
predefined task (SCF, PHAGEN, ARPES, DOS), which uses a predefined set
of parameters for the task, the filename (containing either dot or
slash) from where the parameters will be loaded or the
:class:`InputParameters<ase2sprkkr.input_parameters.input_parameters.InputParameters>` object (see the later example).

::

   calculator = SPRKKR(atoms=atoms, input_parameters='PHAGEN',
                       options={'NE': 5})
   calculator = SPRKKR(atoms=atoms, input_parameters='./input.inp',
                       options={'NE': 5})

You can see and modify the input parameters of a calculator using its
:meth:`set<ase2sprkkr.sprkkr.sprkkr.calculator.SPRKKRCalculator.set>`
and
:meth:`set<ase2sprkkr.sprkkr.sprkkr.calculator.SPRKKRCalculator.get>`
methods:

::

   calculator.set(NE = 5)
   print(calculator.get('NE'))

We can see, that both the
:meth:`set<ase2sprkkr.sprkkr.sprkkr.calculator.SPRKKRCalculator.set>`
and
:meth:`set<ase2sprkkr.sprkkr.sprkkr.calculator.SPRKKRCalculator.get>`
methods and the ``options``
parameter do not require to specify the section names (see the SPR-KKR
manual). They get/set the first (known, according to the definition of
the input parameters) parameter with the given name. However, if it is
necessary to avoid a name conflict, you can either use ``SECTION.VALUE``
notation

::

   calculator.set({'ENERGY.NE':6})
   print(calculator.get('ENERGY.NE'))

...or use the :meth:`input_parameters<ase2sprkkr.sprkkr.sprkkr_calculator.input_parameters>` property of the calculator. Through
this property (which contains :class:`InputParameters<ase2sprkkr.input_parameters.input_parameters.InputParameters>` object) the sections
and their parameters are directly acessible. During interactive work,
you can use tab-completion to see the sections and their parameters.

::

   calculator = SPRKKR(input_paramters = 'SCF')
   calculator.input_parameters.ENERGY.NE = 5
   print(calculator.input_parameters.ENERGY.NE())
   calculator.input_parameters.ENERGY.<tab>

The sections and their parameters have their names in uppercase (at
least partialy, e.g. ``ImE``).

Working with InputParameters object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can also directly create the :class:`InputParameters<ase2sprkkr.input_parameters.input_parameters.InputParameters>` object and pass it
into the
:meth:`InputParameters<ase2sprkkr.sprkkr.calculator.SPRKKR.input_parameters>`
property,
:meth:`calculator constructor<ase2sprkkr.sprkkr.calculator.SPRKKR.__init__>`
or to the
:meth:`InputParameters<ase2sprkkr.sprkkr.calculator.SPRKKR.calculate>`
method

::

   from ase2sprkkr.input_parameters.input_parameters import InputParameters
   input_parameters = InputParameters.create('SCF')
   input_parameters.ENERGY.NE = 5
   input_parameters.TAU.NKTAB = 13
   calculator = SPRKKR(atoms=atoms, input_parameters=input_parameters)

You can also add your own custom (not-predefined) parameters and
sections, if there is a need.

::

   input_parameters.add('MY_CUSTOM_SECTION')
   input_parameters.MY_CUSTOM_SECTION.add('MY_CUSTOM_VALUE', 17)

If you work with input parameters readed from an already created input
file, you can use the
:meth:`calculate<ase2sprkkr.input_parameters.input_parameters.InputParameters.calculate>`
method to avoid the necessity to
create the calculator (manually):

::

   input_parameters = InputParameters.from_file('input.inp')
   input_parameters.calculate()

The input parameters can be set or reset (to the task-predefined values,
to the values provided in a :class:`InputParameters<ase2sprkkr.sprkkr.input_parameters.input_parameters.InputParameters>` object or to the values contained
in an input file) via the calculator :meth:`input_parameters<ase2sprkkr.sprkkr.calculator.SPRKKR.input_parameters>` property.

::

   calculator.input_parameters = 'SCF'
   calculator.input_parameters = input_parameters
   calculator.input_parameters = './input_file.inp'

Running more subsequent tasks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A common usecase is to run more subseuent task with the same atomic
structure. To do so, you can either modify the
:meth:`input_parameters<ase2sprkkr.sprkkr.calculator.SPRKKR.input_parameters>`
and
:class:`potential_file<ase2sprkkr.sprkkr.calculator.SPRKKR.potential_file>`
properties of a calculator (the
:meth:`out.potential_filename<ase2sprkkr.outputs.readers.scf.ScfResult.potential_filename>`
property contains the name of the \`converged potential')

::

   calculator = SPRKKR(atoms=atoms, input_parameters='SCF')
   out = calculator.calculate()
   calculator.calculate(input_parameters='PHAGEN', potential=out.potential_filename)

...or you can use calculator associated with the converged potential,
which is available in the :meth:`output.calculator<ase2sprkkr.outputs.readers.scf.ScfResult.calculator>` property of a SCF task
result:

::

   out = calculator.calculate(input_parameters='SCF')
   out.calculator.calculate(input_parameters='PHAGEN')

Reading the input file
^^^^^^^^^^^^^^^^^^^^^^

If you want to repeat a calculation, you can read (and modify) the input
parameters from an existing input file. In this case, it may be useful
to use the method
:meth:`calculate<ase2sprkkr.input_parameters.input_parameters.InputParameters.calculate>`
of InputParameters object to avoid the
necessity to create the calculator.

::

   input_parameters = InputParameters.from_file('an_input_file.inp')
   input_parameters.calculate(potential='a_potential_file.pot')

or you can of course just set the calculator
:meth:`input_parameters<ase2sprkkr.sprkkr.calculator.SPRKKR.input_parameters>`
property to the filename, or pass the filename as ``input_parameters``
argument to its
:meth:`constructor<ase2sprkkr.sprkkr.calculator.SPRKKR.__init__>`
or
:meth:`calculate<ase2sprkkr.sprkkr.calculator.SPRKKR.calculate>`
method.

Working with potential files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The most of the properties that determine the computed problem in
SPR-KKR are given in potential file. The potential file is created from
the atoms automatically by the calculator. However, it is possible to
pass your own manually (or earlier) created potential file to the
calculation:

::

   calculator.calculate(potential = 'my_potential_file')

or to create the :class:`Potential<ase2sprkkr.potentials.potentials.Potential>` object manually from ASE atoms object and to
alter/check its properties before the computation

::

   from ase2sprkkr.potential.potentials import Potential
   potential = Potential.from_atoms(atoms)
   print(potential.GLOBAL_SYSTEM_PARAMETER.IREL())
   calculator.calculate(potential = potential)

Generating the input files
~~~~~~~~~~~~~~~~~~~~~~~~~~

If you just want to generate the input files (the input and potential ones)
and not to run the calculation, you can use method
:meth:`save_input<ase2sprkkr.sprkkr.calculator.SPRKKR.save_input>`.
It acceps, same as the method
:meth:`calculate<ase2sprkkr.sprkkr.calculator.SPRKKR.calculate>`
the arguments to specify the filenames of the task and potential (and,
for the
:meth:`calculate<ase2sprkkr.sprkkr.calculator.SPRKKR.calculate>`
method, output) file.

::

   calculator.calculate(input_file = ..., potential_file = ..., output_file = ....)

If you pass ``None`` to any of the arguments above, the temp file will
be used. String argument will be interpreted as filename where to store
the input and output files. However, these strings are interpereted as
templates, where the following placeholders can be used:

+----+----------------------------------------------------------------+
| %d | Current date and time                                          |
+----+----------------------------------------------------------------+
| %t | InputParameters name (SCF, PHAGEN, ...)                        |
+----+----------------------------------------------------------------+
| %a | Chemical structure (e.g. CF4)                                  |
+----+----------------------------------------------------------------+
| %c | Counter, starting from one, advanced each function call (that  |
|    | uses the counter)                                              |
+----+----------------------------------------------------------------+

A note about potentials and atoms and input_parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the calculator, one can (in various methods) specify either atoms
object or/and potential, which can lead to the duplicity of the
information provided to the calculator. Therefore, there is the
following logic behind the scenes:

If the potential is not provided, it is created according to the atoms
object (which in this case has to be supplied). If the potential is
given by a filename, it is used "as is", while if it is given by an
object (of a
:class:`Potential<ase2sprkkr.potentials.potentials.Potential>`
class), this object is updated according to
the atoms object (if the latter is supplied). The special value False
for the potential means to use the potential file specified by the task:
in this case, the InputParameters object or task file (with specified
potential) has to be supplied (as the argument).

The input parameters can be given either by a filename (containing the
"input file"), by object (of class
:class:`InputParameters<ase2sprkkr.input_parameters.input_parameters.InputParameters>`
), or by string
(which creates the
:class:`InputParameters<ase2sprkkr.input_parameters.input_parameters.InputParameters>`
object containing the predefined
values). If the parameters are given by object, the input file is
created. If it is given by filename, the file is readed and modified
only if it is required: that is either if the potential argument of
:meth:`save_input<ase2sprkkr.sprkkr.calculator.SPRKKR.save_input>`,
:meth:`run<ase2sprkkr.sprkkr.calculator.SPRKKR.run>` or,
:meth:`calculate<ase2sprkkr.sprkkr.calculator.SPRKKR.calculate>`,
methods is not False - which indicates that the content of the file should
be replaced by the new potential - or if the non-empty options argument is specified -
which leads to the modification of the potential in the file according to the options.

Note, that:

  -  If the potential is set to an SPRKKR calculator (either by the
       :meth:`potential<ase2sprkkr.sprkkr.calculator.SPRKKR.potential>`,
       property setter or in the
       :meth:`constructor<ase2sprkkr.sprkkr.calculator.SPRKKR.__init__>`,
       ), the atoms
       object is created from the potential (and stored in the
       :meth:`atoms<ase2sprkkr.sprkkr.calculator.SPRKKR.atoms>`,
       property). However, this object does not reflect
       the changes made to the potential thereafter.
  -  After the computation of the SCF task, the result provide a new
       potential (in the
       :meth:`result.potential<ase2sprkkr.outputs.readers.scf.ScfResult.potential>`
       property), the old one is not updated, nor the calculator potential property.

MPI calculations
~~~~~~~~~~~~~~~~

You can use MPI in your calculations using the
``mpi`` parameter (of
both the
:meth:`calculator constructor<ase2sprkkr.sprkkr.calculator.SPRKKR.__init__>`,
and 
:meth:`calculate<ase2sprkkr.sprkkr.calculator.SPRKKR.calculate>`,
method). On clusters,
where the number of processes is determined by the batch system, you can
just pass True to the argument, otherwise supply an integer that denotes
the number of wanted processes. The mpi runner is detected
automatically. If the detection failed, the message is printed and
non-MPI calculation is runned. In this case, or if you want to use a
different MPI implementation, you have to pass the correct mpi runner
and its arguments to the calculator:

::

   calculator = SPRKKR(atoms=atoms, mpi = 'path/to/my/mpirunner')
   calculator = SPRKKR(atoms=atoms, mpi = ['path/to/my/mpirunner','-np','4'])

In the case you use common mpi runners (``mpirun``, ``mpirun.openmpi``
or ``mpirun.mpich``), you can either pass just ``True`` to the ``mpi``
parameter - if the runner can detect the number of processes, e.g. in a
cluster environment - or integer that denotes the number of requested
processes.

::

   calculator = SPRKKR(atoms=atoms, mpi=True)
   calculator = SPRKKR(atoms=atoms, mpi=4)

Bundled tools
~~~~~~~~~~~~~

In
:mod:`tools<ase2sprkkr.tools>`
subdirectory of the package, there is a tool
visualise_in_struct for visualisation of the surfaces computed by the
SPR-KKR package. To view a surface, you can run

::

   a2s_visualise_in_struct.py -i in_structure.inp -p potential.pot

where ``potential.pot`` is a potential file and ``in_structure.inp`` is
a file containing the structure of the surface. Run

::

   a2s_visualise_in_struct.py --help

| for the further options of this tool.
| Note: if you have not installed the package using pip or conda, you
  must run the script using its actual location, e.g.

::

   ./tools/a2s_visualise_in_struct.py --help

Contributing and error reporting
--------------------------------

Please feel free to make a pullrequest or post an issue at our `GitHub
repository <https://github.com/ase2sprkkr/ase2sprkkr>`__. Any form of
contributing, error reporting and/or feature request posting will be
highly appreciated.


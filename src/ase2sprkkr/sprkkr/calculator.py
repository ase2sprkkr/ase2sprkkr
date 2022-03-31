"""
Module calculator
=================

This module contains different classes used to define a new calculator for
specific spectroscopies understood by MsSpec.

For more informations about SPR-KKR, follow this
`link <http://olymp.cup.uni-muenchen.de/index.php?option=com_content&
view=article&id=8%3Asprkkr&catid=4%3Asoftware&Itemid=7&lang=en>`

"""
import os, sys
from ase.calculators.calculator import Calculator, all_changes

from .sprkkr_atoms import SPRKKRAtoms
from ..potentials.potentials import Potential
from ..common.misc import add_to_signature
import shutil
import copy
import subprocess
from typing import Union, Any, Dict

class SPRKKR(Calculator):
    """
    ASE calculator for SPR-KKR.

    :cvar ase2sprkkr.input_parameters.input_parameters.InputParameters InputParameters: (just for an easier access to the class)
    :cvar ase2sprkkr.potentials.potentials.Potential Potential: (dto.)

    The parameters of performed calculations are determined mostly by ``input_parameters``, given either
    in :meth:`constructor<ase2sprkkr.sprkkr.calculator.SPRKKR.__init__>` or in the :meth:`calculate<ase2sprkkr.sprkkr.calculator.SPRKKR.calculate>` method.
    The following attributes of the calculator are mostly the default values of the arguments of the :meth:`calculate<ase2sprkkr.sprkkr.calculator.SPRKKR.calculate>` or :meth:`calculate<ase2sprkkr.sprkkr.calculator.SPRKKR.calculate>` methods.

        """
    implemented_properties = ['energy']

    def __init__(self, restart=None,
                 label=None, atoms=None, directory='.',
                 input_file=True, output_file=True, potential_file=True,
                 print_output='info',
                 mpi=False,
                 input_parameters=None, options={}, potential=True,
                 executable_postfix=True,
                 **kwargs):
        """
        Parameters
        ----------
        label: str or None
          Path and begin of the file names used for default templates
          for input, output and potential files.

        directory: str or None
          Directory, where the files will be created.

        input_file: str or file or bool
          Template according to which the input file name will be created.
          True means the default template.
          False means to use a temporary file.
          The placeholders for the task name, date, etc... can be used.

        output_file: str or file or bool (default).
          The template for the output file name.
          True means use the input file name, with replaced .in[p] by .out
          or appended .out suffix
          False means not to write the output file at all.
          See the input_file parameter for the meaning of other values.

        potential_file: str or file or bool
          The template for the potential file name (see the input_file parameter).

        print_output: Union[bool,str]
          Write the output of runned executables to stdout (in addition to the output file)?
          (Default value for the calculator)
          If print_output = 'info', only a few info lines per iteration will be printed.

        mpi: Union[list,string,int,bool]
          Runner for mpi to run a mpi calculation. True and int means autodetect: use True
          for a cluster where mpi is able to autodetect the number of the processes, otherwise
          use integer to specify the number of processes.
          E.g. mpi = [ 'mpirun', '-np', '4' ], mpi = 4

        input_parameters: sprkkr.input_parameters.input_parameters.InputParameters or str or None
          The default input parameters, according to which the input file for the calculation
          will be created. None means that the parameters will be specified in the calculate
          (or save_input) method.
          If str is given, it is interpreted as a name of the task (if no dot and no slash
          are contained in it), for which (see :py:mod:`ase2sprkkr.input_parameters.definitions`) the default
          input parameters will be used, or as a filename contained the input file (which will
          be readed)

        options: dict
          Further parameters for input parameters

        potential: sprkkr.potential.Potential or str or bool or None
          The (default) potential to be used in calculations.
          If a string is given, the potential will be read from the given filename.
          True means to generate the potential from atoms (the default value).
          False means that the potential will be given by the input_parameters.
          None means that the potential is required to be supplied, when calculate or save_input
          methods are called.

        executable_postfix: str or boolean
          String to be added to the runned executable. In some environments, the version
          and the hostname is added to the name of sprkkr executables.
          True: use SPRKKR_EXECUTABLE_SUFFIX environment variable
          False: do not append anything (the same as '')
       """
        if potential and not isinstance(potential, bool):
           if isinstance(potential, str):
               potential = Potential.from_file(potential, atoms = atoms)
           atoms = potential.atoms
        elif atoms:
           atoms = SPRKKRAtoms.promote_ase_atoms(atoms)

        self._potential = None
        self.atoms = atoms
        self.potential = potential
        self.executable_postfix = executable_postfix
        """ The default postfix for the names of executables. False means no postfix, True means
        to use ``SKRKKR_EXECUTABLE_POSTFIX`` environmental variable. """

        super().__init__(restart,
                         label=label, atoms=atoms, directory=directory)
        self.atoms = atoms
        self.potential = potential
        self.mpi = mpi
        """ Default value for mpi parameter of calculate - it determine the way whether and how
        the mpi is employed. """
        self.input_file = input_file
        """ A pathname or an open named file for the input file (if not specified otherwise in the
        called method) """
        self.output_file = output_file
        """ A pathname or an open file for the output of the runned SPR-KKR task (if not specified
        otherwise in the called method) """
        self.potential_file = potential_file
        """ A pathname or an open named file for the potential file to be used (if not specified otherwise in the called method) """
        self.input_parameters = input_parameters
        if options:
          self.input_parameters.set(options, unknown = 'find')
        self.print_output = print_output
        """ The default parameter for the print_output arguments of
        :meth:ase2sprkkr.sprkkr.calculator.SPRKKR.calculate` method - whether ouptut of the
        runned SPR-KKR should be written to the standard ouput (in addition to the output file). """

        self._counter = 0
        #For %c template in file names

    @property
    def input_parameters(self):
        """
        The parameters for the runned task, given either by the
        :class:`ase2sprkkr.input_parameters.input_parameters.InputParameters`
        object, or by string resulting in the default parameters for a taks with the given name.

        If a string is set to the property, the
        :class:`InputPamaters<ase2sprkkr.input_parameters.input_parameters.InputParameters>`
        object is created with the default values for the task with a given name.

        The parameters can be either overriden by ``input_parameters`` argument of the
        :meth:`calculate<ase2sprkkr.sprkkr.calculator.SPRKKR.calculate>` method,
        or (one-time) modified by its ``options`` argument.
        """
        if self._input_parameters is None:
           self.input_parameters = 'SCF'
        return self._input_parameters

    @input_parameters.setter
    def input_parameters(self, value):
        if value is not None:
            value = InputParameters.create_input_parameters(value)
        self._input_parameters=value

    @property
    def potential(self) -> Potential:
       """ The potential associated with the calculator. It will be used in
       calculate and save_input methods, if it is not explicitly overriden by
       the potential argument.
       By default, the potential is created from atoms object (if it is set).

       Return
       ------
       potential: ase2sprkkr.potential.potentials.Potential

       """
       if self._potential is None:
          if not self._atoms:
             return None
          self._potential = Potential.from_atoms(self._atoms)
       return self._potential

    def set(self, options:Union[Dict[str,Any],str,None]={}, value=None, *, unknown='find', **kwargs):
       """
       ASE method to set the parameters.

       Beware, currently, the method sets only input file parameters, not the potential parameters (see the SPRKKR documentation), and do not call the reset() function at all.

      Parameters
      ----------

      options:
        Dictionary of values to be set, or the name of the value, if the value is given.

      value:
        Value to be set. Setting this argument require to pass string name to the options argument.

      unknown: 'add', 'find' or None
        How to handle unknown (not known by the definition) parameters.
        If 'find', try to find the values in descendant containers, throw and exception if none is found.
        If 'add', add unknown values as custom values (use SECTION.OPTION_NAME notation)
        If None, throw an exception.
        Keyword only argument.

      **kwargs: dict
        The values to be set (an alternative syntax as syntactical sugar)


       """
       if options.__class__ is str:
          options = { options : value }
       elif value is not None:
          raise ValueError("If value argument is given to SPRKKR.set, the options"
              " argument have to be string name of the value")
       if kwargs or options:
          self.input_parameters.set(options, **kwargs, unknown=unknown)

    def get(self, name):
       """ Get the value of an input parameter.
           The parameter of the given name will be sought in sections
           and the value of the first such one will be returned. If there
           is an ambiguity, the name of the parameter can be given in
           SECTION.NAME notation.
       """
       return self.input_parameters.get(name)

    @potential.setter
    def potential(self, pot):
       self._potential = pot
       if pot and not isinstance(pot, (bool, str)):
         atoms = pot.atoms
         if atoms:
            self._atoms = atoms
         else:
           pot._atoms = self._atoms

    @property
    def atoms(self) -> SPRKKRAtoms:
       """ Atoms object, associate with the calculator.

       Return
       ------
       atoms: ase2sprkkr.sprkkr.sprkkr_atoms.SPRKKRAtoms
       """
       if self._atoms is None:
          if not self._potential or isinstance(self._potential, (bool, str)):
             return None
          self._atoms = self._potential.atoms
       return self._atoms

    @atoms.setter
    def atoms(self, atoms):
       atoms = SPRKKRAtoms.promote_ase_atoms(atoms)
       self._atoms = atoms
       if self._potential and not self._potential is True:
          self._potential.atoms = atoms

    def _advance_counter(self):
        """ Advance counter for generating filenames with %c counter placeholder """
        self._counter+=1
        return self.counter

    def _open_file(self, filename, templator=None, named=False, mode='w+',
                         allow_temporary=True, create_subdirs=False):
        """
        Open a file given a 'template' filename file name

        Parameters
        ----------
        filename: str or file or None
          if False, temporary object is used
          if str is given, use it as a template
          if file is given, return it unchanged

        templator: callable
          If it is given, the filemae is processed using this function,
          see FilenameTemplator

        named: bool
          Whether to use named temporary file.

        mode: str
          Mode to open the file

        allow_temporary: bool
          If False, throw an exception if the resulting filename should be temporary,
          i.e. if None as filename is given

        create_subdirs: bool
          If true, create all non-existent directories in the given path
        """

        if filename is False:
           if not allow_temporary:
              raise ValueError("Creation of temporary files is not allowed in this context")
           f = tempfile.NamedTemporaryFile if named else tempfile.TemporaryFile
           return f(mode = mode)
        if not isinstance(filename, str):
           return filename

        if templator:
           filename = templator(filename)

        if '/' not in filename and self._directory:
           filename = os.path.join(self._directory, filename)

        if create_subdirs:
           dirs = os.path.dirname(filename)
           try:
             os.makedirs(dirs)
           except FileExistsError:
             pass

        return open(filename, mode)

    def save_input(self, atoms=None, input_parameters=None,
                  potential=None, input_file=None, potential_file=None, output_file=False,
                  create_subdirs=False,
                  options={},return_files=False,mpi=None):
        """
        Save input and potential files for a calculation.

        Parameters
        ----------
        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.

        input_parameters: sprkkr.input_parameters.input_parameters.InputParameters or str or None
            If None, the task specified in __init__ is used, or the default parameters for 'SCF' task
               will be created.
            If a string is given then if it is a task name, an InputParameters object will be created
               (with the default options for a given task) and saved to the input_file.
               Otherwise the task is interpreted as a task file name from which the input parameters will be readed.
            The input parameters is written into the input file,
               whose name is given by input_file argument.

        potential: sprkkr.potential.potentials.Potential or str or None or False
            If it is None, the self.potential value is used instead (see the later options)
            If False, the potential filename is specified in the input_parameters.
            If True,  create the potential according to the atoms and write to the resulting <potential_file> file.
            If str, then it is the filename of the potential, that is left unchanged on input
            If None, use the potential value of the calcualtor.

        input_file: str or bool or None
            Filename or template (see FilenameTemplator class) where to save the input file,
            None means to use the default value from the Calculator.
            See the __init__ method for the meaning of other options.

        potential_file: str or bool or None
            Filename or template (see FilenameTemplator class) where to save the potential_file
            None means to use the default value from the Calculator.
            See the __init__ method for the meaning of other options.

        output_file: str or bool or None
            A filename or a filename template (see FilenameTemplator class) where to save the output
            None means to use the default value from the Calculator.
            See the __init__ method for the meaning of other options.

        create_subdirs: boolean
            If true, create directories if they don't exists.

        options: dict
            Options to set to the input_parameters. If input_parameters are given by a filename,
            they are readed from the file, altered by the options and then the modified input file
            (in a temporary file) will be used.

        return_files: boolean
            Return open files object instead of just string filenames.

        mpi: bool or None
            Save input for a mpi calculation. None means to use the mpi value specified in the
            constructor. Actually, the only difference is that the temporary input file has
            to have filename.

        Returns
        -------
        input_parameters: InputParameters
            Created task object to be run.

        input_filename: str or file
            Input/task file. See return_files

        potential_filename: str or file
            Potential file. See return_files

        output_filename: str or file
            Output file. See return_files

        """

        def makepath(path, must_exist):
           if '/' not in path and self._directory:
              path = os.path.join(self._directory, path)
           if must_exist and not os.path.isfile(path):
              raise ValueError(must_exist.format(path=path))
           return path

        def from_input_name(template, replacement, default):
            if template is not True:
               return template
            if input_file.name:
               if input_file.name.endswith('.inp'):
                  return input_file.name[:-4] + replacement
               if input_file.name.endswith('.in'):
                  return input_file.name[:-3] + replacement
            return default

        templator = FilenameTemplator(self)

        """ Resolve potential argument - whether is in fact given or not"""
        if potential is None:
           potential = self.potential

        """ Get the input file """
        save_input = True
        if input_parameters:
           if isinstance(input_parameters, str):
              if InputParameters.is_it_a_input_parameters_name(input_parameters):
                 input_parameters = InputParameters.create_input_parameters(input_parameters)
              else:
                 if not options and not potential and not input_file:
                   save_input = False
                   input_file = makepath(input_parameters, "'{path}' is not a task file nor a known name of input_parameters.")
                 input_parameters = InputParameters.from_file(input_parameters)
        else:
           input_parameters = self.input_parameters
        templator.input_parameters = input_parameters

        input_file = input_file or self.input_file
        if input_file is True:
            input_file = (self.label or '') + '%a_%t.inp'

        input_file = self._open_file(input_file, templator, input_parameters.is_mpi(self.mpi if mpi is None else mpi),
                                    allow_temporary=return_files,
                                    create_subdirs=create_subdirs,
                                    mode = 'w+' if save_input else 'r'
                                    )


        if potential and potential is not True and atoms:
            raise ValueError("You can not provide both a potential and atoms object to the SPRKKR calculate method")

        """ Get the potential file """
        # potential_file - the file containing the potential (possibly a template)
        # potential - potential object (to be updated by atoms)

        potential = potential if potential is not None else self.potential
        atoms = atoms or self._atoms
        if isinstance(potential, str):
              if potential_file:
                 potential = Potential.from_file(potential)
              else:
                 potential_file = makepath(potential, "'{path}' is not a potential file")
                 potential = None
        elif potential is True:
             if not atoms:
                 raise ValueError("Potential set to <True> which means to generate the potential from the ASE-atoms object."
                                  "However, this object has not been supplied")
             potential=Potential.from_atoms(atoms)
        elif potential is False:
              if potential_file is not None:
                  raise ValueError("When potential is True, the value of POTENTIAL option from the input (parameters) file will be used "
                                   " as the potential file. Thus, specifing the <potential_file> argument is not allowed.")
              potential = None
        elif potential is None:
           raise ValueError("The potential can not be <None>. However, consider supplying either <True> as the potential to generate"
                             " it from the atoms object, or False to use the POTENTIAL option from the input (parameters) file")

        if potential:
           if potential_file is None:
               potential_file = self.potential_file
           pf = from_input_name(potential_file, '.pot', '%a.pot')
           pf = self._open_file(pf, templator, True,
                                            allow_temporary=return_files,
                                            create_subdirs=create_subdirs)
           potential.save_to_file(pf, atoms)
           pf.close()
           potential_file = pf.name

        """ update input_parameters by the potential """
        if save_input:
          if options:
            input_parameters.set(options, unknown = 'find')
          if potential_file:
            #use the relative potential file name to avoid too-long-
            #-potential-file-name problem
            dirr = os.path.dirname( input_file.name ) if input_file.name else self.directory
            input_parameters.CONTROL.POTFIL = os.path.relpath( potential_file, dirr )
          else:
            potential_file = input_parameters.CONTROL.POTFIL()
          input_parameters.save_to_file(input_file)
          input_file.seek(0)
        #This branch can occur if the potential have been explicitly set to False,
        #which means not to create the potential (and take it from the input_parameters)
        else:
            potential_file = input_parameters.CONTROL.POTFIL

        if output_file is not False:
          """ Get the output file """
          output_file = output_file or self.output_file
          output_file = from_input_name(output_file, '.out', '%a_%t.out')
          output_file = self._open_file(output_file, templator, False, mode='wb',
                                            allow_temporary=return_files,
                                            create_subdirs=create_subdirs)

        if not return_files:
          for f in input_file, output_file:
              if hasattr(f,'close'):
                 f.close()
          input_file = input_file.name
          if output_file:
             output_file = output_file.name

        if output_file:
           return input_parameters, input_file, potential_file, output_file
        else:
           return input_parameters, input_file, potential_file



    def run(self, atoms=None, input_parameters=None,
                  potential=None, input_file=None, potential_file=None, output_file=None,
                  create_subdirs=False,
                  options={},
                  print_output=None, executable_postfix=None, mpi=None):
        """
        Do the calculation, return various results.

        Parameters
        ----------

        print_output: bool or str or None
            Print output to stdout, too.
            If print_output=='info' only a few lines per iteration will be printed.
            None means to use a default value (specified in constructor)

        executable_postfix: str or bool or None
            If not None, it overrides the executable_postifx, that have been specified when the
            calculator have been created.

        mpi: bool or None
            Runner for mpi to run a mpi calculation. True and int means autodetect: use True
            for a cluster where mpi is able to autodetect the number of the processes, otherwise
            use integer to specify the number of processes.
            None means use the default value (specified in the calculator's contructor).
            E.g. mpi = [ 'mpirun', '-np', '4' ], mpi = 4

        ---------
        For the other parameters, see the save_input() method
        """

        if print_output is None:
           print_output = self.print_output

        if output_file==False:
           output_file = None

        input_parameters, input_file, _, output_file = self.save_input(
                            atoms=atoms, input_parameters=input_parameters,
                            potential=potential, output_file=output_file, options=options,
                            return_files=True, mpi=mpi)

        return input_parameters.run_process(self, input_file, output_file, print_output if print_output is not None else self.print_output,
                                     executable_postfix = executable_postfix,
                                     mpi=self.mpi if mpi is None else mpi
                                    )

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes,
                  input_parameters=None, potential=None,
                  input_file=None, potential_file=None, output_file=None,
                  create_subdirs=False,
                  options={},
                  print_output=None, executable_postfix=None, mpi=None):
        """
        ASE-interface method for the calculation.  This method runs the appropriate task(s)
        for the requested properties (currently always the SCF one) and updates the
        calculator object with the results.

        See
        https://wiki.fysik.dtu.dk/ase/development/calculators.html#ase.calculators.calculator.Calculator.calculate
        for the documentation of ASE interface.

        Parameters
        ----------
        properties: list or str
            List of what needs to be calculated. Can be any combination of ‘energy’, ‘forces’, ‘stress’, ‘dipole’, ‘charges’, ‘magmom’ and ‘magmoms’.

        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these six: 'positions', 'numbers', 'cell',
            'pbc', 'initial_charges' and 'initial_magmoms'.

            This calculator ignore properties and system_changes.

        ---------
        For the other parameters, see the run() and save_input() methods

        ---------

        From ASE documentation: Calculated properties should be
        inserted into results dictionary like shown in this dummy example::

            self.results = {'energy': 0.0,
                            'forces': np.zeros((len(atoms), 3)),
                            'stress': np.zeros(6),
                            'dipole': np.zeros(3),
                            'charges': np.zeros(len(atoms)),
                            'magmom': 0.0,
                            'magmoms': np.zeros(len(atoms))}

        """
        super().calculate(atoms, properties, system_changes)
        out = self.run(
                  atoms, input_parameters,
                  potential, input_file, potential_file, output_file,
                  create_subdirs,
                  options,
                  print_output, executable_postfix, mpi=mpi
        )
        if hasattr(out, 'energy'):
            self.results.update({
              'energy' : out.energy,
            })
        return out


    def scf(self, *args, **kwargs):
         """A shortcut for calculating a SCF task"""
         self.calculate(input_parameters='SCF', *args, **kwargs)

    def phagen(self, *args, **kwargs):
         """A shortcut for calculating a PHAGEN task"""
         self.calculate(input_parameters='PHAGEN', *args, **kwargs)

    def kkrgen(self, *args, **kwargs):
         """A shortcut for calculating a KKRGEN task"""
         self.calculate(input_parameters='KKRGEN', *args, **kwargs)

    def kkrspec(self, *args, **kwargs):
         """A shortcut for calculating a KKRSPEC task"""
         self.calculate(input_parameters='KKRSPEC', *args, **kwargs)

    def kkrch(self, *args, **kwargs):
         """A shortcut for calculating a KKRCH task"""
         self.calculate(input_parameters='KKRCH', *args, **kwargs)

    def run_xband(self, executable='xband'):
        subprocess.Popen(executable, stdout=subprocess.DEVNULL)

    def copy_with_potential(self, potential):
        """ Return copy of self, with the potential variable set.
            Use the method to create a new calculator with the given "result potential"
        """
        out = copy.copy(self)
        out.potential = potential
        out.results = out.results.copy()
        return out

    @property
    def potential_object(self):
        """ Convert self.potential to a Potential object """
        if isinstance(self.potential, str):
          self.potential = Potential.from_file(self.potential)
        return self.potential


class FilenameTemplator:
    """ Class that replaces the placeholders in filenames with propper values,
        see the replacements property of the class.
        The used values are remembered so all the processed filenames use
        the same value (e.g. counter, datetime etc...)
    """

    def __init__(self, calculator):
        self.data = {}
        self.calculator = calculator
        self.input_parameters = None

    def _get(self, name, default, calculator):
        if not name in self.data:
            self.data[name] = default(self, calculator)
        return self.data[name]

    replacements = {
          "%d" : lambda self, calc: datetime.today().strftime('%Y-%m-%d_%H:%M'),
          "%t" : lambda self, calc: self.input_parameters.task_name if self.input_parameters else 'SPRKKR',
          "%a" : lambda self, calc: str(calc.atoms.symbols) if calc.atoms else 'custom',
          "%c" : lambda self, calc: calc._advance_counter()
        }

    def __call__(self, template):
        for i,v in self.replacements.items():
          if i in template:
             template = template.replace(i, self._get(i, v, self.calculator))
        return template

#at least, to avoid a circular import
from ..input_parameters.input_parameters import InputParameters

#is there a better way to not document an inner classes?
if os.path.basename(sys.argv[0]) != 'sphinx-build':
    SPRKKR.InputParameters = InputParameters
    SPRKKR.Potential = Potential

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

from .sprkkr_atoms import SprKkrAtoms
from ..task.tasks import Task
from ..potential.potentials import Potential
from ..common.misc import add_to_signature
import shutil

class SprKkr(Calculator):
    """
    ASE calculator for SPR-KKR.

    :cvar ase2sprkkr.task.tasks.Task Task: (just for an easier access to the class)
    :cvar ase2sprkkr.potential.potentials.Potential Potential: (dto.)

    """
    implemented_properties = ['energy']

    def __init__(self, restart=None,
                 label=None, atoms=None, directory='.',
                 input_file=True, output_file=True, potential_file=True,
                 print_output=False,
                 mpi=True, task=None, potential=None,
                 command_postfix=True,
                 **kwargs):
        """
        Parameters
        ----------
        label: str or None
          Path and begin of the file names used for default templates
          for input, output and potential files.

        directory: str or None
          Directory, where the files will be created.

        input_file: str or file or None or True
          Template according to the task (input) file name will be created.
          None means to use a temporary file.
          The placeholders for task name, date, etc... can be used.
          True means the default template.

        output_file: str or file or None or True (default).
          The template for the output file name (see the input_file parameter).
          True means use the input file name, with replaced .in[p] by .out
          or appended .out suffix

        potential_file: str or file or None
          The template for the potential file name (see the input_file parameter).

        print_output: bool
          Write the output of runned executables to stdout (in addition to output file).
          (Default value for are calculations)

        mpi: list or string or True
          Runner for mpi to run mpi calculation. True means autodetect.
          E.g. 'mpirun'

        potential: sprkkr.potential.Potential or str or None
          Potential to be used.
          If string is given, the potential will be read from the file.
          The atoms will be created (or updated, if it is given) from this potential.

        command_postfix: str or boolean
          String to be added to the runned command. In some environments, the version
          and the hostname is added to the name of sprkkr executables.
          True: use SPRKKR_COMMAND_SUFFIX environment variable
          False: do not append anything (the same as '')

        task: sprkkr.task.tasks.Task or str or None
          Task to compute. Can be None, if no task is specified before calculate,
          SCF task is created.
          If str is given, it is interpreted as name of task (if no dot and no slash is contained
          in it), or filename contained the task
        """
        if potential:
           if isinstance(potential, str):
               potential = Potential.read_from_file(potential, atoms = atoms)
           atoms = potential.atoms
        elif atoms:
           atoms = SprKkrAtoms.promote_ase_atoms(atoms)

        self._potential = None
        self.atoms = atoms
        self.potential = potential
        self.command_postfix = command_postfix

        super().__init__(restart,
                         label=label, atoms=atoms, directory=directory)
        self.atoms = atoms
        self.potential = potential

        if mpi is True:
           self.mpi = [ 'mpirun' ] if shutil.which('mpirun') else []
        else:
           self.mpi = [ mpi ] if isinstance(mpi, str) else mpi

        self.input_file = (self.label or '') + '%a_%t.inp' if input_file is True else input_file
        self.output_file = output_file
        self.potential_file = potential_file
        self.task = Task.create_task(task)

        #For %c template in file names
        self._counter = 0

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

    @potential.setter
    def potential(self, pot):
       self._potential = pot
       if pot:
         atoms = pot.atoms
         if atoms:
            self._atoms = atoms
         else:
           pot._atoms = self._atoms

    @property
    def atoms(self) -> SprKkrAtoms:
       """ Atoms object, associate with the calculator.

       Return
       ------
       atoms: ase2sprkkr.sprkkr.sprkkr_atoms.SprKkrAtoms
       """
       if self._atoms is None:
          if not self._potential:
             return None
          self._atoms = self._potential.atoms
       return self._atoms

    @atoms.setter
    def atoms(self, at):
       self._atoms = at
       if self._potential:
          self._potential.atoms = at

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
          if None, temporary object is used
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

        if filename is None:
           if not allow_temporary:
              raise ValueException("Creation of temporary files is not allowed in this context")
           f = tempfile.NamedTemporaryFile if named else tempfile.TemporaryFile
           return f(mode = mode)
        if not isinstance(filename, str):
           return tem

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

    def save_input(self, atoms=None, task=None, system_changes=all_changes,
                  potential=None, task_file=None, potential_file=None, output_file=False,
                  create_subdirs=False,
                  options = {}, return_files=False):
        """
        Save input (task) and potential files for a calculation.

        Paramters
        ---------
        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.

        task: sprkkr.task.tasks.Task or str or None
            If None, task specified in __init__ is used, or the default 'SCF' task
               will be created.
            If string is given then if it is a task name, the Task object will be created
               and the task file created using it.
               Otherwise the task is interpreted as task file, which will be left unchanged.
            The task is written into the location given by self.input_file.

        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these six: 'positions', 'numbers', 'cell',
            'pbc', 'initial_charges' and 'initial_magmoms'.

            This calculator ignore properties and system_changes.

        potential: sprkkr.potential.potentials.Potential or str or None or False
            If False, the potential file is specified in the task.
            If None, create the potential according to the atoms and self.potential_file property.
            If str, then it is the filename of the potential, that is left unchanged on input

        task_file: str or None
            Filename or template (see FilenameTemplator class) where to save the task file.
            None means to use the default value from the Calculator.

        potential_file: str or None
            Filename or template (see FilenameTemplator class) where to save the potential_file
            None means to use the default value from the Calculator.
            If potential is given as filename, the potential is NOT readed and stored again,
            so in this case setting potential_file has no effect.

        output_file: str or None or False
             Filename or template (see FilenameTemplator class) where to save the output
             None means to use the default value from the Calculator.
             False means not to create the file.

        create_subdirs: boolean
             If true, create directories if they don't exists.

        options: dict
            Options to set to the task. If task is given by filename, it is readed, altered by the
            options and then the modified task (in a temporary file) will be used.

        return_files: boolean
            Return open files object instead of just string filenames.

        Returns
        -------
        task: Task
            Created task object to be run.

        task_filename: str or file
            Input/task file. See return_files

        potential_filename: str or file
            Potential file. See return_files

        output_filename: str or file
            Output file. See return_files

        """

        def makepath(path, must_exist):
           if '/' not in path and self._directory:
              path = os.path.join(self._directory, path)
           if must_exist and not os.path.isfile(task):
              raise ValueError(must_exist.format(path))
           return path

        def from_task_name(template, replacement, default):
            if template is not True:
               return template
            if task_file.name:
               if task_file.name.endswith('.inp'):
                  return task_file.name[:-4] + replacement
               if task_file.name.endswith('.in'):
                  return task_file.name[:-3] + replacement
            return default

        templator = FilenameTemplator(self)

        """ Get the task file """
        save_task = True
        if task:
           if isinstance(task, str):
              if Task.is_it_a_task_name(task):
                 task = task.create_task(task)
              else:
                 task_file = makepath(task, "'{path}' is not a task file nor a known name of task.")
                 task = task.from_file(task_file)
                 if not options and not potname and not potential:
                   save_task = False
        else:
           task = self.task or Task.default_task()
        templator.task = task

        task_file = self._open_file(task_file or self.input_file, templator, False,
                                    allow_temporary=return_files,
                                    create_subdirs=create_subdirs
                                    )

        """ Get the potential file """
        # potential - potntial object
        if potential:
           if isinstance(potential, str):
              potential_file = makepath(potential, "'{path}' is not a potential file")
              potential = None
              potname = potential_file
        elif potential is not False:
           potential = self.potential
        else:
           potname = potential_file = None

        if potential:
           pf = potential_file or self.potential_file
           pf = from_task_name(pf, '.pot', '%a.pot')
           potential_file = self._open_file(pf, templator, True,
                                            allow_temporary=return_files,
                                            create_subdirs=create_subdirs)
           potential.save_to_file(potential_file, self.atoms)
           potential_file.close()
           potname = potential_file.name

        """ update task by the potential """
        if save_task:
          if options:
            task.set(options, unknown = 'find')
          if potname:
            #use the relative potential file name to avoid too-long-
            #-potential-file-name problem
            dirr = os.path.dirname( task_file.name ) if task_file.name else self.directory
            task.CONTROL.POTFIL = os.path.relpath( potname, dirr )
            task.save_to_file(task_file)
            task_file.seek(0)
        #This branch can occures only if the potential have been explicitly set to False,
        #which means not to create the potential (and take it from the task file)
        else:
            potname = task.CONTROL.POTFIL

        if output_file is not False:
          """ Get the output file """
          output_file = output_file or self.output_file
          output_file = from_task_name(output_file, '.out', '%a_%t.out')
          output_file = self._open_file(output_file, templator, False, mode='wb',
                                            allow_temporary=return_files,
                                            create_subdirs=create_subdirs)

        if not return_files:
          task_file.close()
          for f in task_file, output_file:
              if hasattr(f,'close'):
                 f.close()
          task_file = task_file.name
          if output_file:
             output_file = output_file.name

        if output_file:
           return task, task_file, potname, output_file
        else:
           return task, task_file, potname



    def calculate(self, atoms=None, task=None, system_changes=all_changes,
                  potential=None, task_file=None, potential_file=None, output_file=None,
                  create_subdirs=False,
                  options={},
                  print_output=False, command_postfix=None):
        """
        Do the calculation.

        From ASE documentation: Calculated properties should be
        inserted into results dictionary like shown in this dummy example::

            self.results = {'energy': 0.0,
                            'forces': np.zeros((len(atoms), 3)),
                            'stress': np.zeros(6),
                            'dipole': np.zeros(3),
                            'charges': np.zeros(len(atoms)),
                            'magmom': 0.0,
                            'magmoms': np.zeros(len(atoms))}
        Parameters
        ----------

        print_output: bool
            Print output to stdout, too.

        command_postfix: str or bool or None
            If not None, it overrides the command_postifx, that have been specified when the
            calculator have been created.



        ---------
        For other parameters, see save_input() method
        """

        super().calculate(atoms, system_changes)

        if output_file==False:
           output_file = None

        task, task_file, _, output_file = self.save_input(atoms=atoms, task=task,
                            potential=potential, output_file=output_file, options=options,
                            return_files=True)

        return task.run_task_process(self, task_file, output_file, print_output,
                                    command_postfix = command_postfix)


    def scf(self, *args, **kwargs):
         """A shortcut for calculating a SCF task"""
         self.calculate(task='SCF', *args, **kwargs)

    def phagen(self, *args, **kwargs):
         """A shortcut for calculating a PHAGEN task"""
         self.calculate(task='PHAGEN', *args, **kwargs)

    def kkrgen(self, *args, **kwargs):
         """A shortcut for calculating a KKRGEN task"""
         self.calculate(task='KKRGEN', *args, **kwargs)

    def kkrspec(self, *args, **kwargs):
         """A shortcut for calculating a KKRSPEC task"""
         self.calculate(task='KKRSPEC', *args, **kwargs)

    def kkrch(self, *args, **kwargs):
         """A shortcut for calculating a KKRCH task"""
         self.calculate(task='KKRCH', *args, **kwargs)

#is there a better way to not document an inner classes?
if os.path.basename(sys.argv[0]) != 'sphinx-build':
    SprKkr.Task = Task
    SprKkr.Potential = Potential


class FilenameTemplator:
    """ Class that replaces the placeholders in filenames with propper values,
        see the replacements property of the class.
        The used values are remembered so all the processed filenames use
        the same value (e.g. counter, datetime etc...)
    """

    def __init__(self, calculator):
        self.data = {}
        self.calculator = calculator
        self.task = None

    def _get(self, name, default, calculator):
        if not name in self.data:
            self.data[name] = default(self, calculator)
        return self.data[name]

    replacements = {
          "%d" : lambda self, calc: datetime.today().strftime('%Y-%m-%d_%H:%M'),
          "%t" : lambda self, calc: self.task.task_name if self.task else 'SPRKKR',
          "%a" : lambda self, calc: str(calc.atoms.symbols) if calc.atoms else 'custom',
          "%c" : lambda self, calc: calc._advance_counter()
        }

    def __call__(self, template):
        for i,v in self.replacements.items():
          if i in template:
             template = template.replace(i, self._get(i, v, self.calculator))
        return template

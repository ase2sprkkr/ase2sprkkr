"""
Module calculator
=================

This module contains different classes used to define a new calculator for
specific spectroscopies understood by MsSpec.

For more informations about SPR-KKR, follow this
`link <http://olymp.cup.uni-muenchen.de/index.php?option=com_content&
view=article&id=8%3Asprkkr&catid=4%3Asoftware&Itemid=7&lang=en>`__

"""
import os
from ase.calculators.calculator import Calculator, all_changes

from .sprkkr_atoms import SprKkrAtoms
from ..task.tasks import Task
from ..potential.potentials import Potential
from ..common.misc import add_to_signature
import shutil

class SprKkr(Calculator):
    implemented_properties = ['energy']

    def __init__(self, restart=None,
                 label=None, atoms=None, directory='.',
                 input_file=True, output_file=True, potential_file=True,
                 print_output=False,
                 mpi=True, task=None, potential=None,
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
          Input file. None means use temporary file. Can use placeholders
          for task name, date, etc....
          True means default template form the label.

        output_file: str or file or None or True (default)
          Same as the input file. Moreover, True means use the input file name,
          with replaced .in[p] by .out or appended .out suffix.

        print_output: bool
          Write the output of runned executables to stdout (in addition to output file).
          (Default value for are calculations)

        potential_file: str or file or None
          Same as the input_file

        mpi: list or string or True
          Runner for mpi to run mpi calculation. True means autodetect.
          E.g. 'mpirun'


        task: sprkkr.task.tasks.Task or str or None
          Task to compute. Can be None, if no task is specified before calculate,
          SCF task is created.
          If str is given, it is interpreted as name of task (if no dot and no slash is contained
          in it), or filename contained the task
        """

        if atoms:
           atoms = SprKkrAtoms.promote_ase_atoms(atoms)
        elif potential:
           atoms = potential.atoms
        self._potential = None
        self.atoms = atoms
        self.potential = potential

        super().__init__(restart,
                         label=label, atoms=atoms, directory=directory)
        self.atoms = atoms
        self.potential = potential

        if mpi is True:
           self.mpi = [ 'mpirun' ] if shutil.which('mpirun') else []
        else:
           self.mpi = [ mpi ] if isinstance(mpi, str) else mpi

        self.input_file = (self.label or '') + 'a_%t.inp' if input_file is True else input_file
        self.output_file = output_file
        self.potential_file = (self.label or '') + '%a.pot' if potential_file is True else potential_file
        self.task = Task.create_task(task)

        """ For %c template in file names """
        self._counter = 0

    @property
    def potential(self):
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
    def atoms(self):
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

    def advance_counter(self):
        self._counter+=1
        return self.counter

    def _open_file(self, filename, templator=None, named=False, mode='w+'):
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
        """

        if filename is None:
           f = tempfile.NamedTemporaryFile if named else tempfile.TemporaryFile
           return f(mode = mode)
        if not isinstance(filename, str):
           return tem

        if templator:
           filename = templator(filename)

        if '/' not in filename and self._directory:
           filename = os.path.join(self._directory, filename)

        return open(filename, mode)


    def calculate(self, atoms=None, task=None, system_changes=all_changes,
                  potential = None, output_file=False, print_output=True,
                  options = {}):
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


        Paramters
        ---------
        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these six: 'positions', 'numbers', 'cell',
            'pbc', 'initial_charges' and 'initial_magmoms'.

            This calculator ignore properties and system_changes.a

        task: sprkkr.task.tasks.Task or str or None
            If None, task specified in __init__ is used, or the default 'SCF' task
               will be created.
            If string is given then if it is a task name, the Task object will be created
               and the task file created using it.
               Otherwise the task is interpreted as task file, which will be left unchanged.
            The task is written into the location given by self.input_file.

        potential: sprkkr.potential.potentials.Potential or str or None or False
            If False, the potential file is specified in the task.
            If None, create the potential according to the atoms and self.potential_file property.
            If str, then it is the filename of the potential, that is left unchanged on input
                (however, it will be modified by SPRKKR itself)

        output_file: str or False or True
            Filename to store output of the SPRKKR. False mean do not store (just read), True
            mean to use the default template. An open binary file can be also given.

        print_output: bool
            Print output to stdout, too.

        options: dict
            Options to set to the task. If task is given by filename, it is readed, altered by the
            options and then the modified task (in a temporary file) will be used.

        """

        def makepath(path, must_exist):
           if '/' not in path and self._directory:
              path = os.path.join(self._directory, path)
           if must_exist and not os.path.isfile(task):
              raise ValueError(must_exist.format(path))
           return path

        super().calculate(atoms, task, system_changes)

        templator = FilenameTemplator(self)

        """ Get the potential file """
        if potential:
           if isinstance(name, str):
              potential_file = makepath(task, "'{path}' is not a task file nor a known name of task.")
              potential = None
        elif potential is not False:
           potential = Potential.from_atoms(self.atoms)
        if potential:
           potential_file = self._open_file(self.potential_file, templator, True)
           potential.save_to_file(potential_file, self.atoms)

        """ Get the task file """
        if task:
           if isinstance(name, str):
              if Task.is_it_a_task_name(task):
                 task = task.create_task(task)
              else:
                 task_file = makepath(task, "'{path}' is not a task file nor a known name of task.")
                 task = None
        else:
           task = self.task or Task.default_task()

        if options and not task:
           task = task.from_file(task)
           task_file = None #use a temporary file

        if task:
           if options:
               task.set(options, unknown = 'find')
           task_file = self._open_file(self.input_file, templator, False)
           if potential_file.name:
              task.CONTROL.POTFIL = potential_file.name
           task.save_to_file(task_file)
           task_file.seek(0)

        """ Get the output file """
        output_file = output_file or self.output_file
        if output_file is True:
           if getattr(task_file, 'name', None):
              output_file = re.sub('(.in(p?))?', '.out', task_file, name)
              output_file = open(output_file, 'wb')
           else:
              output_file = self._open_file('%a_%t.out', templator, False, mode='wb')
        elif output_file:
              output_file = self._open_file(output_file, templator, False, mode='wb')

        task.run_task_process(self, task_file, output_file, print_output=print_output if print_output is not None else self.print_output)


    def phagen(self):
         self.calculate(self.atoms, None, None)

    def kkrgen(self):
         self.calculate(self.atoms, None, None)

    def kkrspec(self):
         self.calculate(self.atoms, None, None)

    def kkrchi(self):
         self.calculate(self.atoms, None, None)

    def read(self):
        raise NotImplementedError



class FilenameTemplator:
    """ Class that replaces the placeholders in filenames with propper values,
        see the replacements property of the class.
        The used values are remembered so all the processed filenames use
        the same value (e.g. counter, datetime etc...)
    """

    def __init__(self, calculator):
        self.data = {}
        self.calculator = calculator

    def _get(self, name, default, calculator):
        if not name in self.data:
            self.data[name] = default(calculator)
        return self.data[name]

    replacements = {
          "%d" : lambda calc: datetime.today().strftime('%Y-%m-%d_%H:%M'),
          "%t" : lambda calc: calc.task.task_name if calc.task else 'sprkkr',
          "%a" : lambda calc: str(self.atoms.symbols) if calc.atoms else 'custom',
          "%c" : lambda calc: calc.advance_counter()
        }

    def __call__(self, template):
        for i,v in self.replacements.items():
          if i in template:
             template.replace(i, self._get(i, v, self.calculator))
        return template

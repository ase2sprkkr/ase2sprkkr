from ..output_definitions import OutputSectionDefinition as Section, \
                                 OutputValueDefinition as V, \
                                 OutputValueEqualDefinition as VE, \
                                 OutputNonameValueDefinition as VN
from ...common.grammar_types import Table, integer, string, Real, RealWithUnits,String, Sequence, Array
import pyparsing as pp
import numpy as np
from ase.units import Rydberg
from ..task_result import TaskResult, TaskResultReader
from ...common.misc import cached_property

class ScfResult(TaskResult):

  def __init__(self, task, iterations, error, return_code):
      self.iterations = iterations
      super().__init__(task, return_code)

  @cached_property
  def potential_filename(self):
      """ New (output) potential file name """
      potfil = self.task.CONTROL.POTFIL()
      if not potfil:
         raise ValueError("Please set CONTROL.POTFIL of the task to read the potential")
      fname = self.task.CONTROL.POTFIL() + '_new'
      if self.directory:
         fname = os.path.join(self.directory, fname)
      return fname

  @cached_property
  def potential(self):
      return Potential.from_file(self.potential_filename)

  @property
  def energy(self):
      return self.iterations[-1]['ETOT']

  @property
  def converged(self):
      return self.iterations[-1]['converged']

  def iteration_values(self, name):
      if not name in self.iterations[0]:
         raise KeyError(f"No such iteration value: {name}")

      return np.fromiter(
              (i[name] for i in self.iterations),
              count = len(self.iterations),
              dtype = self.iterations[0][name].__class__
            )

  def plot(self, what=['error', 'ETOT'], filename=None, logscale = set(['err']), **kwargs):
      import matplotlib.pyplot as pyplot

      fig, axs = pyplot.subplots(len(what), sharex=True)
      iterations = self.iteration_values('iteration')
      for name, ax in zip(what, axs):
        data = self.iteration_values(name)
        if name in logscale:
           ax.set_yscale("log")
        kwargs.update({
            'ls' : None,
            'mfc': None,
            'mew': 2,
            'ms' : 6
        })
        ax.plot(iterations, data, '+', **kwargs)
        ax.set_ylabel(name)

      ax.set_xlabel('Iteration')
      fig.tight_layout()
      if filename:
        fig.savefig(filename)
      else:
        pyplot.show()

class KkrScfProcessOutputReader(TaskResultReader):

  atoms_conf_type = Section('atoms', [
    VN('IECURR', integer),
    VE('E', float),
    VN('L', float),
    VE('IT', integer),
    VN('atom', string),
    VN('orbitals', Table({
        'l' : string,
        'DOS': float,
        'NOS': float,
        'P_spin' : float,
        'm_spin' : float,
        'P_orb' : float,
        'm_orb' : float,
        'B_val' : float,
        'B_val_desc': String(default_value = ''),
        'B_core' : Real(default_value = float('NaN'))
        }, free_header=True, default_values=True)),
    V('E_band', RealWithUnits(units = {'[Ry]' : Rydberg }), is_optional=True),
    V('dipole moment', Sequence(int, Array(float, length=3)), is_optional=True)
  ])

  async def read_output(self, stdout):

        iterations = []

        async def readline():
          line = await stdout.readline()
          if not line:
             raise EOFError()
          return line.decode('utf8')

        async def readlinecond(cond, canend=True):
          while True:
            line = await stdout.readline()
            if not line:
               if canend:
                  return ''
               raise EOFError()
            if cond(line):
              return line.decode('utf8')
        try:
          while True:
            out = {}
            line = await readlinecond(lambda line: b'SPRKKR-run for: ' in line)
            if not line:
                return iterations
            run = line.replace('SPRKKR-run for:', '').strip()
            out['run'] = run

            line = await readlinecond(lambda line: b' E=' in line)
            atoms = []
            while True:
              atoms.append(await self.atoms_conf_type.parse_from_stream(stdout,
                up_to=b'\n -------------------------------------------------------------------------------',
                start=line
              ))
              line = await readlinecond(lambda line: line!=b'\n')
              if not 'E=' in line:
                break
            out['atoms'] = atoms

            line = await readlinecond(lambda line: b' ERR' in line and b'EF' in line)
            items = line.split()
            out['iteration'] = int(items[0])
            out['error']=float(items[2])
            out['EF']=float(items[5])
            out['M']= float(items[10]), float(items[11])

            line = (await readline()).split()
            out['ETOT'] = float(line[1]) * Rydberg
            out['converged'] = line[5] == 'converged'
            iterations.append(out)
            out = {}

        except EOFError:
          raise Exception('The output ends unexpectedly')

  result_class = ScfResult


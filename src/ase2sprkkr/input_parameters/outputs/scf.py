from ..output_definitions import OutputSectionDefinition as Section, \
                                 OutputValueDefinition as V, \
                                 OutputValueEqualDefinition as VE, \
                                 OutputNonameValueDefinition as VN
from ...common.grammar_types import Table, integer, string, Real, RealWithUnits,String, Sequence, Array
import pyparsing as pp
import numpy as np
from ase.units import Rydberg
from ..task_result import TaskResult, OutputReader
from ...common.misc import cached_property
from ...potential.potentials import Potential
import os
import copy
from ...sprkkr.calculator import SPRKKR
from ...common.formats import fortran_format

class ScfResult(TaskResult):

  @property
  def iterations(self):
      return self.result

  @cached_property
  def potential_filename(self):
      """ New (output) potential file name """
      potfil = self.input_parameters.CONTROL.POTFIL()
      if not potfil:
         raise ValueError("Please set CONTROL.POTFIL of the input_parameters to read the potential")
      fname = self.input_parameters.CONTROL.POTFIL() + '_new'
      if self.directory:
         fname = os.path.join(self.directory, fname)
      return fname

  @cached_property
  def potential(self):
      return Potential.from_file(self.potential_filename)

  @cached_property
  def calculator(self):
      if self._calculator:
         return self._calculator.copy_with_potential(self.potential_filename)
      else:
         SPRKKR(potential = self.potential_filename)

  @property
  def energy(self):
      return self.last_iteration['ETOT']

  @property
  def converged(self):
      return self.last_iteration['converged']

  def iteration_values(self, name):
      if not name in self.iterations[0]:
         raise KeyError(f"No such iteration value: {name}")

      return np.fromiter(
              (i[name] for i in self.iterations),
              count = len(self.iterations),
              dtype = self.iterations[0][name].__class__
            )
  @property
  def last_iteration(self):
      if not self.iterations:
          raise AttributeError('No iteration has been finished')
      return self.iterations[-1]

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

class ScfOutputReader(OutputReader):

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
          first = True
          while True:
            out = {}
            line = await readlinecond(lambda line: b'SPRKKR-run for: ' in line)
            if not line:
                break
            line=line.strip()
            if first and self.print_output == 'info':
               print(line)
               first = False
            run = line.replace('SPRKKR-run for:', '')
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
            out['moment'] = {'spin' : float(items[10]),
                             'orbital' : float(items[11]) }
            line = (await readline()).split()
            out['ETOT'] = float(line[1]) * Rydberg
            out['converged'] = line[5] == 'converged'
            iterations.append(out)
            if self.print_output == 'info':
               error = fortran_format(out['error'], ":>12e")
               print(f"Iteration {out['iteration']:>5} error {error} "
                     f"spin moment: {out['moment']['spin']:>13.6e} "
                     f"orbital moment: {out['moment']['orbital']:>13.6e} ")

          if self.print_output == 'info':
            if not iterations:
                print("ERROR: No iteration has been finished. There is probably an error in the input files, please, examine the SPR-KKR output file.")
            elif iterations[-1]['converged']:
                print("OK: The computation converged.")
            else:
                print("WARNING: The computation does not converged!!!")

          return iterations

        except EOFError:
          raise Exception('The output ends unexpectedly')

  result_class = ScfResult


""" SCF (selfconsisten cycle) reader and result. """

import pyparsing as pp
import numpy as np
from ase.units import Rydberg

from ..output_definitions import OutputSectionDefinition as Section, \
                                 OutputValueDefinition as V, \
                                 OutputValueEqualDefinition as VE, \
                                 OutputNonameValueDefinition as VN
from ...common.grammar_types import Table, integer, string, Real, RealWithUnits,String, Sequence, Array
from ..task_result import TaskResult
from ...common.process_output_reader import readline, readline_until
from ..sprkkr_output_reader import SprKkrOutputReader
from ...common.decorators import cached_property
from ...sprkkr.calculator import SPRKKR
from ...common.formats import fortran_format
from ...common.grammar import replace_whitechars
from ..task_result import KkrProcess
from ...potentials.potentials import Potential


class RealOrStars(Real):
  """ A real value, where ``****`` means ``NaN`` """
  _grammar = Real._grammar | replace_whitechars(pp.Word('*')).setParseAction(lambda x: float('NaN'))


class ScfResult(TaskResult):
  """ Objects of this class holds the results of a self-consistent cycle and its iterations.
  It also allows to plot the convergence of values during iterations.
  """

  @cached_property
  def calculator(self):
      """ The calculator that has the new (output) potential assigned - i.e. that
      can be used for further calculations with the (hopefully) converged wavefunctions etc. """
      if self._calculator:
         return self._calculator.copy_with_potential(self.potential_filename)
      else:
         SPRKKR(potential = self.potential_filename)

  @property
  def potential_filename(self):
      fname = super().potential_filename + '_new'
      return fname

  @property
  def start_potential_filename(self):
      return super().potential_filename

  @property
  def start_potential(self):
      return Potential.from_file(self.potential_filename)

  @property
  def energy(self):
      """ Total energy of the last iteration """
      return self.last_iteration.energy.ETOT()

  @property
  def converged(self):
      """ The calculation coverged or not? """
      return self.last_iteration['converged']()

  def iteration_values(self, name):
      """ Return the array of values of given name from the
      iterations (i.e. "the column of the table of the values")

      E.g. result.iteration_values(('energy', 'ETOT')) return list of
      floats - the total energies computed in each iteration.

      Parameters
      ----------
      name: str
        Name of the parameter

      Return
      ------
      values: list
        The values.
      """
      if not name in self.iterations[0]:
         raise KeyError(f"No such iteration value: {name}")

      return np.fromiter(
              (i[name]() for i in self.iterations),
              count = len(self.iterations),
              dtype = self.iterations[0][name].__class__
            )

  @property
  def last_iteration(self):
      """ Return the data of the last iteration """
      if not self.iterations:
          raise AttributeError('No iteration has been finished')
      return self.iterations[-1]

  @property
  def energies(self):
      return self.last_iteration.energy()

  def plot(self, what=['error', ('energy','ETOT'), ('energy', 'EMIN')], filename=None, logscale = set(['err']), **kwargs):
      """ Plot the development of the given value(s) during iterations.

      Parameters
      ----------
      what: list[str]
        Names to be plotted

      filename: str or None
        If filename is given, the plot is rendered to the file, it is show on the screen
        otherwise.

      logscale: collections.abc.Set[str]
        Which values should be rendered in logscale

      **kwargs: dict
        All other arguments are passed to the matplotlib Axes.plot function
      """
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


"""
This definition parses the part of the output file that contains the datas
about atoms and computed data associated with them.
"""
atomic_types_definition = Section('atoms', [
  VN('IECURR', integer),
  VE('E', RealOrStars()),
  VN('L', float),
  VE('IT', integer),
  VN('symbol', string),
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
      'B_core' : Real(default_value = float('NaN'), nan=r'\*+')
    }, free_header=True, default_values=True)),
  V('E_band', RealWithUnits(units = {'[Ry]' : Rydberg }), is_required=False),
  V('dipole moment', Sequence(int, Array(float, length=3)), is_required=False)
])


"""
This definition is not used for parsing, but just for the data returned
from the reader.
"""
scf_section = Section('iteration', [
  V('system_name', str),
  V('iteration', int, info='Number of the iteration.'),
  V('error', float),
  V('converged', bool, info='True, if the SCF cycle converged this iteration.'),
  Section('moment', [
    V('spin', float),
    V('orbital', float)
  ]),
  Section('energy' , [
    V('EF', float, info='Fermi energy', alternative_names='fermi'),
    V('ETOT', float, info='Total energy', alternative_names='total'),
    V('EMIN', float, info='Bottom of energy contour for band states', alternative_names='band_states_min'),
    V('ESCBOT', float, info='Lower limit for semi-core states', alternative_names='semi_core_min', is_required=False),
    V('ECTOP', float, info='Upper limit for core states', alternative_names='core_max', is_required=False)
  ]),
  Section('atomic_types', atomic_types_definition.members(), is_repeated=True)
])


class ScfOutputReader(SprKkrOutputReader):
  """
  This class reads and parses the output of the SCF task of the SPR-KKR.
  """

  async def read_output(self, stdout, result):
        await self.read_commons(stdout, result)
        result.files['converged'] = result.files['potential'] + '_new'
        iterations = []
        try:
          first = True
          while True:
            out = {}
            line = await readline_until(stdout,lambda line: b'EMIN   = ' in line)
            if not line:
                break

            out['energy'] = {'EMIN' : float(line.split('=')[1]) }
            line = await stdout.readline()
            if b'ESCBOT' in line:
                out['energy']['ESCBOT'] = float(line.split(b'=')[1])
                line = await stdout.readline()
            if b'ECTOP' in line:
                out['energy']['ECTOP'] = float(line.split(b'=')[1])

            line = await readline_until(stdout,lambda line: b'SPRKKR-run for: ' in line)
            line=line.strip()
            if first and self.print_output == 'info':
               print(line)
               first = False
            run = line.replace('SPRKKR-run for:', '')
            out['system_name'] = run

            line = await readline_until(stdout,lambda line: b' E=' in line)
            atoms = []
            while True:
              atoms.append(await atomic_types_definition.parse_from_stream(stdout,
                up_to=b'\n -------------------------------------------------------------------------------',
                start=line
              ))
              line = await readline_until(stdout,lambda line: line!=b'\n')
              if not 'E=' in line:
                break
            out['atomic_types'] = atoms

            line = await readline_until(stdout,lambda line: b' ERR' in line and b'EF' in line)
            items = line.split()
            out['iteration'] = int(items[0])
            out['error']=float(items[2])
            out['energy']['EF']=float(items[5])
            out['moment'] = {'spin' : float(items[10]),
                             'orbital' : float(items[11]) }
            line = (await readline(stdout)).split()
            out['energy']['ETOT'] = float(line[1]) * Rydberg
            out['converged'] = line[5] == 'converged'

            iterations.append(scf_section.read_from_dict(out))
            if self.print_output == 'info':
               error = fortran_format(out['error'], ":>12e")
               print(f"Iteration {out['iteration']:>5} error {error} "
                     f"spin moment: {out['moment']['spin']:>13.6e} "
                     f"orbital moment: {out['moment']['orbital']:>13.6e} ")

          if self.print_output == 'info':
            if not iterations:
                print("ERROR: No iteration has been finished. There is probably an error in the input files, please, examine the SPR-KKR output file.")
            elif iterations[-1]['converged']():
                print("OK: The computation converged.")
            else:
                print("WARNING: The computation does not converged!!!")

          result.iterations = iterations

        except EOFError:
          raise Exception('The output ends unexpectedly')


class ScfProcess(KkrProcess):

  result_class = ScfResult
  reader_class = ScfOutputReader

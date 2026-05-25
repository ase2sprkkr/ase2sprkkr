""" SCF (selfconsisten cycle) reader and result. """

import re
import pyparsing as pp
import numpy as np
import unyt

from ..output_definitions import OutputSectionDefinition as Section, \
                                 OutputValueDefinition as V, \
                                 OutputValueEqualDefinition as VE, \
                                 OutputNonameValueDefinition as VN
from ...common.grammar_types import Table, integer, string, Real, RealWithUnits,String, Sequence, Array
from ..task_result import TaskResult, ResultValue
from ...common.process_output_reader import readline, readline_until
from ..sprkkr_output_reader import SprKkrOutputParser
from ...common.decorators import cached_property
from ...sprkkr.calculator import SPRKKR
from ...common.formats import fortran_format
from ...common.grammar import replace_whitechars
from ..task_result import KkrOutputReader
from ...potentials.potentials import Potential


class IterationValue(ResultValue):

  def __init__(self, name, value, plot):
      super().__init__(name, value)
      self._plot = plot

  def plot(self):
      self._plot()

  def actions(self):
      return [ 'plot' ]


class DataValue(ResultValue):

  def value_label(self):
      return ("<data...>")

  def actions(self):
      return [ 'data' ]

  def data(self):
      return self._value

def plot_iterations(iterations, what=['error', ('energy','ETOT'), ('energy', 'EMIN')], filename=None, logscale = None, **kwargs):

    import matplotlib.pyplot as pyplot

    if logscale is None:
        logscale = { 'error' }
    if not isinstance(what, list):
        what = [what]

    fig, axs = pyplot.subplots(len(what), sharex=True, squeeze=False)
    iterations_ids = iterations.values_of('iteration')
    for name, ax in zip(what, axs.ravel()):
      data = iterations.values_of(name)
      if name in logscale:
         ax.set_yscale("log")
      kwargs.update({
          'ls' : None,
          'mfc': None,
          'mew': 2,
          'ms' : 6
      })
      ax.plot(iterations_ids, data, '+', **kwargs)
      ax.set_ylabel(iterations[0][name]._definition.name)
    axs[-1,0].set_xlabel('Iteration')

    fig.tight_layout()
    if filename:
      fig.savefig(filename)
    else:
      pyplot.show()




class RealOrStars(Real):
  """ A real value, where ``****`` means ``NaN`` """
  _grammar = Real._grammar | replace_whitechars(pp.Word('*')).set_parse_action(lambda x: float('NaN'))


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
      try:
          return self.last_iteration.energy.ETOT()
      except Exception as e:
          raise AttributeError("No iteration has been finished") from e

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

      return self.iterations.values_of(name)

  @property
  def last_iteration(self):
      """ Return the data of the last iteration """
      if not self.iterations:
          raise AttributeError('No iteration has been finished')
      return self.iterations[-1]

  @property
  def energies(self):
      return self.last_iteration.energy()


  @cached_property
  def output_values(self):
      return {
          'converged' : ResultValue('Converged', str(self.last_iteration.converged()) ),
          'iterations': ResultValue('Number of iterations', len(self.iterations) ),
          'fermi energy': IterationValue('Fermi energy', self.last_iteration.energy.EF(), lambda: self.plot(('energy','EF'))),
          'error': IterationValue('Error', self.last_iteration.error(), lambda: self.plot('error')),
          'data': DataValue('Data', self),
      }

  def plot(self, what=['error', ('energy','ETOT'), ('energy', 'EMIN')], filename=None, logscale = None, **kwargs):
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
      return plot_iterations(self.iterations, what, filename, logscale, **kwargs)

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
  V('E_band', RealWithUnits(units = {'[Ry]' : unyt.Ry }), is_required=False),
  V('dipole moment', Sequence(int, Array(float, length=3)), is_required=False)
])


class PlottableValue(V):

    def enrich(self, option):
        def plot(options, **kwargs):
            c = option._container._container
            plot_iterations(option._container._container, self.get_path[len(c.get_path):], **kwargs)
        option.plot = lambda **kwargs: plot(option, **kwargs)
        super().enrich(option)

PV = PlottableValue

"""
This definition is not used for parsing, but just for the data returned
from the reader.
"""
scf_section = Section('iteration', [
  V('system_name', str),
  V('iteration', int, info='Number of the iteration.'),
  PV('error', float),
  V('b_error', float, is_required=False, info='RMS error of the exchange-correlation B-field at this iteration.'),
  V('duration', float, is_required=False, info='Execution time of this iteration in seconds.'),
  V('converged', bool, info='True, if the SCF cycle converged this iteration.'),
  Section('moment', [
    PV('spin', float),
    PV('orbital', float)
  ]),
  Section('energy' , [
    PV(('EF', 'fermi'), float, info='Fermi energy'),
    PV(('ETOT','total'), float, info='Total energy'),
    PV(('EMIN', 'band_states_min'), float, info='Bottom of energy contour for band states'),
    PV(('ESCBOT','semi_core_min'), float, info='Lower limit for semi-core states', is_required=False),
    PV(('ECTOP', 'core_max'), float, info='Upper limit for core states', is_required=False)
  ]),
  Section('atomic_types', atomic_types_definition.members(), is_repeated=True)
], is_repeated=True)

del PV

class ScfOutputParser(SprKkrOutputParser):
  """
  This class reads and parses the output of the SCF task of the SPR-KKR.
  """
  def set_print_output(self, print_output):
      self.print_info = print_output == 'info'
      self.print_output = print_output and not self.print_info

  async def read_output(self, stdout, result):
        await self.read_commons(stdout, result)
        result.files['converged'] = result.files['potential'] + '_new'
        iterations = scf_section.create_object()
        try:
          first = True
          _err_or_time = re.compile(rb'execution time for last iteration| ERR.*EF')

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
            if first and self.print_info:
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
            duration = None
            while True:
                line = await readline_until(stdout, _err_or_time.search)
                if 'execution time for last iteration' in line:
                    try:
                        duration = float(line.split()[-2])
                    except (ValueError, IndexError):
                        pass
                else:
                    break
            items = line.split()
            out['iteration'] = int(items[0])
            out['error']=float(items[2])
            out['b_error'] = float(items[3])
            out['energy']['EF']=float(items[5])
            out['moment'] = {'spin' : float(items[10]),
                             'orbital' : float(items[11]) }
            if duration is not None:
                out['duration'] = duration
            line = (await readline(stdout)).split()
            out['energy']['ETOT'] = float((float(line[1]) * unyt.Ry).to(unyt.eV))
            out['converged'] = line[5] == 'converged'
            section = iterations.add()
            section.set(out)

            if self.print_info:
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


class ScfOutputReader(KkrOutputReader):

  result_class = ScfResult
  parser_class = ScfOutputParser

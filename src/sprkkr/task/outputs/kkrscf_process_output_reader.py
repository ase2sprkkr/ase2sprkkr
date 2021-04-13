from ..output_definitions import OutputSectionDefinition as Section, \
                                 OutputValueDefinition as V, \
                                 OutputValueEqualDefinition as VE, \
                                 OutputNonameValueDefinition as VN
from ...common.grammar_types import Table, integer, string, Real, RealWithUnits,String, Sequence, Array
import pyparsing as pp
from ase.units import Rydberg
from ..task_result import TaskResult, TaskResultReader

class ScfResult(TaskResult):

  def __init__(self, task, iterations, error, return_code):
      self.iterations = iterations
      super().__init__(task, return_code)

  @property
  def energy(self):
      return self.iterations[-1]['ETOT']

  @property
  def converged(self):
      return self.iterations[-1]['converged']

  """
  def read_results(self):
       outstrg=self.read_output(os.path.join(self.directory,self.outfile))
       lastiter=len(outstrg['it'])
       self.niter=lastiter
       self.converged = outstrg['converged'][lastiter-1]
       if not self.converged:
           raise RuntimeError('SPRKKR did not converge! Check ' + self.outfile)

       self.results['raw_outfile'] = outstrg
       self.results['energy']=outstrg['ETOT'][lastiter-1]*Rydberg
  """

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
    V('dipole moment', Sequence(int, Array(float, length=3)))
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


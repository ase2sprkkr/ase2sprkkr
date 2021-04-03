from ...common.process_output_reader import BaseProcessOutputReader

class KkrScfProcessOutputReader(BaseProcessOutputReader):

  def __init__(self, cmd, outfile, **kwargs):
      super().__init__(cmd, outfile, **kwargs)
      self.iterations = []

      #TODO - remove
      self.data = {
            'it' : [],
            'ERR' : [],
            'ETOT' : [],
            'EF' : [],
            'M' : [],
            'converged' : [],
            'atom_confs' : [],
      }

  async def read_output(self, stdin):

      try:
        while True:
          line = await stdin.readline()
          if not line:
             break

          if 'ERR' in line and 'EF' in line:
              items = line.split()
              self.iterations.append(
                  {'iteration' : int(items[0]),
                    'error'    : float(items[2]),
                    'EF'       : float(items[5]),
                    'M'        : (float(items[10]), float(items[11])),
                  }
              )
      except StopIteration:
        pass

  def not_considered(self):
                if conitn:
                    items = _skip_lines(fd, 1).split()
                    out['ETOT'].append(float(items[1]))
                    flag = items[5] == 'converged'
                    out['converged'].append(flag)

                    out['atom_confs'].append(atom_confs)
                    atom_confs = {}

                    _skip_lines(fd, 1)

                elif 'SPRKKR-run for:' in line:
                    run = line.replace('SPRKKR-run for:', '').strip()
                    out['run'] = run

                elif ' E= ' in line:
                    atom = line.split()[-1]
                    akeys = _skip_lines(fd, 1).split()
                    line = _skip_lines_to(fd, 'sum').split()
                    avals = list(map(float, line[1:8])) + [float(line[9])]
                    line = _skip_lines_to(fd, 'E_band').split()
                    akeys.append(line[0])
                    if len(line) >= 2:
                       avals.append(line[1])
                    atom_conf = dict(zip(akeys, avals))
                    atom_confs[atom] = atom_conf

  def read_results(self):
       outstrg=self.read_output(os.path.join(self.directory,self.outfile))
       lastiter=len(outstrg['it'])
       self.niter=lastiter
       self.converged = outstrg['converged'][lastiter-1]
       if not self.converged:
           raise RuntimeError('SPRKKR did not converge! Check ' + self.outfile)

       self.results['raw_outfile'] = outstrg
       self.results['energy']=outstrg['ETOT'][lastiter-1]*Rydberg



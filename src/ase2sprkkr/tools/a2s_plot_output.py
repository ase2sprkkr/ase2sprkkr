#!/usr/bin/env python
"""
This is a sctipt to visualise in_struct.inp files. Run it to see the doc.
"""
import argparse
from pathlib import Path

if not __package__:
  __package__ = 'ase2sprkkr.tools'
import sys
sys.path.append(str(Path(__file__).resolve().parents[2]))
from ..common.tools import parse_tuple_function
from ..outputs.output_files import OutputFile

def plot():
  #READ IN COMAND LINE PARAMETERS
  description='This is a script for plotting SPR-KKR output files'
  parser = argparse.ArgumentParser(description=description)
  parser.add_argument('output', help='SPR-KKR output file, e.g. ARPES output (*.spc).')
  parser.add_argument('-v','--value', dest='value', type=str, help='Only the value of the given name will be plotted (can be repeated)', action='append', default=[], required=False)
  parser.add_argument('-o','--output_filename', dest='filename', type=str, help='The plot will be saved to a file with given name, instead of showing it on the screen', default=None, required=False)
  parser.add_argument('-L','--do_not_use_latex', dest='latex', action='store_false', default=True, help='Do not use LaTex for captions generating', required=False)
  parser.add_argument('-s','--plot_size', dest='figsize', default=(6,4), type=parse_tuple_function(float,2),   help='The plot size', required=False)
  parser.add_argument('-c','--colormap', dest='colormap', type=str, help='Matplotlib colormap', required=False)
  parser.add_argument('-n','--norm', dest='norm', choices=['lin', 'log'], help='Matplotlib colormap will use linear or logarithmic scale (the default behavior depends on the plotted data)', required=False)

  args = parser.parse_args()
  kwargs = vars(args)

  of = OutputFile.from_file(args.output)
  del kwargs['output']

  value = args.value
  del kwargs['value']

  if value:
    fn = kwargs['filename']
    for name in value:
      if fn:
         kwargs['filename'] = append_id_to_filename(fn, name)
      try:
         val = of[name.upper()]
      except KeyError:
         raise ValueError(f"There is no value named '{name.upper()}' in the output file.")
      if not hasattr(val, 'plot'):
         raise ValueError(f"Value '{name.upper()}' do not know, how it should be plotted.")
      val.plot(**kwargs)
  else:
    of.plot(**kwargs)

if __name__ == "__main__":
   plot()

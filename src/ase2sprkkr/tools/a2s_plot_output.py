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
from ..common.tools import parse_tuple_function, parse_named_option
from ..outputs.output_files import OutputFile

def plot():
  #READ IN COMAND LINE PARAMETERS
  description='This is a script for plotting SPR-KKR output files. The type of the file is guessed from the content of the file and from the extension. The currently supported files are: \n' + \
   '\n'.join(map(lambda x: f"    {x[0].upper()}: {x[1].definition.info()}", OutputFile.definitions.items()))

  parser = argparse.ArgumentParser(
      description=description,
      formatter_class=argparse.RawDescriptionHelpFormatter
      )
  parser.add_argument('output', help='SPR-KKR output file name (see the supported files above).')
  parser.add_argument('-o','--output_filename', dest='filename', type=str, help='The plot will be saved to a file with given name, instead of showing it on the screen', default=None, required=False)
  parser.add_argument('-v','--value', dest='value', type=str, help='Only the value of the given name will be plotted (option can be repeated)', action='append', default=[], required=False)

  parser.add_argument('-s','--plot_size', dest='figsize', default=(6,4), type=parse_tuple_function(float,2),   help='The plot size', required=False)
  parser.add_argument('-c','--colormap', dest='colormap', type=str, help='Matplotlib colormap', required=False)
  parser.add_argument('-n','--norm', dest='norm', choices=['lin', 'log'], help='Matplotlib colormap will use linear or logarithmic scale (the default behavior depends on the plotted data)', required=False)

  parser.add_argument('-L','--do_not_use_latex', dest='latex', action='store_false', default=True, help='Do not use LaTex for captions generating', required=False)
  parser.add_argument('-S','--set', dest='args', type=parse_named_option, help='Given a value of the format name=value, pass the value to the plotting function. You can so override various options that matplotlib plotting functions accept (e.g. vmin or vmax for pcolormesh), or the values that can be set using set_<something> functions (e.g. title or (x|y)label). This option can be repeated.', action='append', default=[], required=False)
  args = parser.parse_args()
  kwargs = vars(args)

  of = OutputFile.from_file(args.output)
  del kwargs['output']

  value = args.value
  del kwargs['value']

  kwargs = { k:v for k,v in kwargs.items() if v is not None }
  kwargs.update(dict(kwargs['args']))
  del kwargs['args']

  if value:
    fn = kwargs.get('filename', None)
    for name in value:
      if fn is not None:
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

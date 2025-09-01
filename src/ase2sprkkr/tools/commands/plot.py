#!/usr/bin/env python
"""
Plotting the values in SPRKKR output files
"""
from pathlib import Path
import sys
import argparse

if not __package__:
  __package__ = 'ase2sprkkr.tools.commands'
sys.path.append(str(Path(__file__).resolve().parents[3]))

from ...common.tools import parse_tuple_function, parse_named_option, \
        append_id_to_filename, parse_inches, main # NOQA
from ...common.lazy_string import LazyString      # NOQA


@LazyString
def description():
    from ...output_files.output_files import OutputFile
    out = 'The type of the file is guessed from the content of the file and from the extension. The currently supported files are: \n'
    defs = OutputFile.definitions.items()
    out += '\n'.join(map(lambda x: f"    {x[0].upper()}: {x[1].definition.info()}", defs))
    return out


help='Plot SPR-KKR output files.'


def parser(parser):
    def parse_layer(value):
        out=parse_tuple_function(int, 1, 2)(value)
        if len(out) == 1:
            return out[0]-1
        else:
            return slice(out[0]-1, out[1])

    parser.add_argument('output', help='SPR-KKR output file name (see the supported files above).')
    parser.add_argument('-o','--output_filename', dest='filename', type=str, help='The plot will be saved to a file with given name, instead of showing it on the screen. If more values are given (see -v), their name will be added to the filename.', default=None, required=False)
    parser.add_argument('-v','--value', dest='value', type=str, help='Only the value of the given name will be plotted (option can be repeated)', action='append', default=[], required=False)

    parser.add_argument('-s','--plot_size', dest='figsize', default=argparse.SUPPRESS, type=parse_tuple_function(parse_inches,2), help='The plot size. Example: "5cm,5cm", "6,4". The default units are inches.', required=False)
    parser.add_argument('-c','--colormap', dest='colormap', type=str, help='Matplotlib colormap', required=False)
    parser.add_argument('-d','--dpi', dest='dpi', type=float, help='DPI of the resulting image (default 600)', required=False)
    parser.add_argument('-n','--norm', dest='norm', choices=['lin', 'log'], help='Matplotlib colormap will use linear or logarithmic scale (the default behavior depends on the plotted data)', required=False)

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-l','--use_latex', dest='latex', action='store_true', help='Force use LaTex for generating captions. Default is to use the default matplotlib settings (can be configured per user).', required=False)
    group.add_argument('-L','--do_not_use_latex', dest='latex', action='store_false', help='Do not use LaTex for generating captions', required=False)
    parser.add_argument('-S','--set', dest='args', type=lambda x: parse_named_option(x,True), help='Given a value of the format name=value, pass the value to the plotting function. You can so override various options that matplotlib plotting functions accept (e.g. vmin or vmax for pcolormesh), or the values that can be set using set_<something> functions (e.g. title or (x|y)label). This option can be repeated.', action='append', default=[], required=False)
    parser.set_defaults(latex=None)  # or True/False if you want a default
    group.add_argument('--separate_plots', help='Plot each value from file in a separate window or/and file.', action='store_true', required=False)

    group = parser.add_argument_group('BSF specific options')
    group.add_argument('--layer', help='Select a layer for plotting. Either number or two comma delimited numbers from,to. Numbering starts from 1. ', type=parse_layer,
                       default=argparse.SUPPRESS, required=False)
    group.add_argument('--fermi', help='Draw a line at Fermi energy. Optional float specifies line width.', nargs='?', const=True, type=float,
                       default=argparse.SUPPRESS, required=False)

def run(args):
  from ...output_files.output_files import OutputFile
  kwargs = vars(args)

  of = OutputFile.from_file(args.output)
  for x in [ 'layer' ]:
      if x in kwargs and not x in of.plot_parameters:
          raise ValueError(f"Argument '--{x}' is not valid for {of}")

  del kwargs['output']

  value = args.value
  del kwargs['value']

  kwargs = { k:v for k,v in kwargs.items() if v is not None }
  kwargs.update(dict(kwargs['args']))
  del kwargs['args']

  if value:
    fn = kwargs.get('filename', None)
    for name in value:
      if fn is not None and len(value) > 1:
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
    main( globals() )

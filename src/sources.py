import os
import os.path
import re


def sources(dir, ext):

  if isinstance(ext, str):
      fce = lambda x: x.endswith(ext)
  elif isinstance(ext, re.Pattern):
      fce = lambda x: ext.match(x)
  elif ext is True:
      fce = lambda x: True
  else:
      fce = ext

  for dirpath, dirnames, filenames in os.walk(dir):
        fns = [f for f in filenames if fce(f) ]

        if not fns:
            continue
        print(f"\n{dirpath}\n")
        for f in fns:
            print(f"{dirpath}/{f}")


sources('ase2sprkkr', '.py')
sources('ase2sprkkr/input_parameters/examples', '.in')
sources('ase2sprkkr/output_files/examples', True)
sources('ase2sprkkr/outputs/examples', '.out')
sources('ase2sprkkr/bindings/xband/tests/', '.pot')
sources('ase2sprkkr/output/examples', '.out')
sources('ase2sprkkr/sprkkr/test', '.cif')
sources('ase2sprkkr/potentials/examples', re.compile(r'.*[.](pot|new)'))

import os
import glob

root = os.path.dirname(os.path.dirname(__file__))
os.chdir(root)

with open(os.path.join('sphinx', 'examples_generated.rst'), "w") as f:
  for i in glob.glob(os.path.join('src', 'ase2sprkkr', 'examples', '*', '*.py')):
      name, filename = i.rsplit('/', 2)[-2:]
      if filename.startswith('_'):
          continue
      with open(i, "r") as script:
          doc=script.read().split('"""', 3)[1]
          doc=doc.strip()
          doc=doc.split('\n\n',1)[0]
      modname = filename[:-3]
      header = f":py:mod:`{name} <ase2sprkkr.examples.{name}.{modname}.main>` "
      f.write(f"{header}\n{'-' * len(header)}\n\n")
      f.write(doc)
      f.write(f"  `[show source] <./_modules/ase2sprkkr/examples/{name}/{modname}.html#main>`_")
      f.write("\n\n\n")

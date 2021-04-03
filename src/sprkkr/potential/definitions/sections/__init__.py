import pkgutil
import os
import importlib

def f():
  names = (i for i in pkgutil.iter_modules([os.path.dirname(__file__)]))
  im = importlib.import_module
  modules = (im('.' + i.name, __package__) for i in names )
  return { m.section.__name__ : m.section for m in modules }

sections = f()
locals().update(sections)
__all__ = list(sections.keys())

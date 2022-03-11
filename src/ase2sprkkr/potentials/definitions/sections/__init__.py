""" This package contains the predefined sections to be used in
definitions of input_parameters """

import pkgutil
import os
import importlib

def _sections():
  """ Import all my modules and get the sections from them """
  names = (i for i in pkgutil.iter_modules([os.path.dirname(__file__)]))
  im = importlib.import_module
  modules = (im('.' + i.name, __package__) for i in names )
  return { m.section.__name__ : m.section for m in modules }

sections = _sections()
locals().update(sections)
__all__ = list(sections.keys())

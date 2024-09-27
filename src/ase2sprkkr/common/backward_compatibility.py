""" In this module, backward compatibility issues are solved by mocking """
import functools
import enum
import pyparsing as pp
from . import decorators

if not hasattr(pp.ParserElement, "setName"):
    pp.ParserElement.setName = getattr(pp.ParserElement, "set_name", lambda x: None)

if not hasattr(functools,'cache'):
    functools.cache = functools.lru_cache(maxsize=None)
    """ Functools.cache. Python 3.8 and earlier does not have this method, so it is mocked for this version of python. """

if not hasattr(enum, 'nonmember'):
    """ I hope it will be working """
    def _nonmember(cls):
        def fn(self):
             return cls

        class Nonmember(decorators.cached_class_property):

            def __getattr__(self, name):
                return getattr(cls, name)

        return Nonmember(lambda: cls)
    enum.nonmember = _nonmember

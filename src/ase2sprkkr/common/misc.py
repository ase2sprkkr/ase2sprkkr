""" Various classes and routines used thorough the package """

import functools
import inspect
import heapq
import copy
import itertools
import asyncio
import numpy as np
_x = object()
from collections import OrderedDict as _OrderedDict
from typing import Set

try:
  from numba import njit
except ImportError:
  def njit(fce):
      """ Mock the numba JIT compiler, if it is not available """
      return fce

def lazy_value(fce):
    """ The decorator for only once computed value. Same a functools.cache,
        but there is no need to take care of arguments.
    """
    x = []

    """ Hack for staticmethod decorator, which is in fact binded by the descriptor protocol """
    if isinstance(fce, staticmethod):
        fce = fce.__func__

    @functools.wraps(fce)
    def cached_fce():
        if not x:
           x.append(fce())
        return x[0]
    def clear():
        x = []

    cached_fce.clear = clear
    return cached_fce

""" Python 3.8 does not have functools.cache """
if hasattr(functools,'cache'):
   cache = functools.cache
else:
   cache = functools.lru_cache(maxsize=None)

def add_to_signature(func, prepend=False, self_name='self'):
  """
  Add the arguments in the <func> function to the list of arguments
  of the resulting function (as keyword_only arguments)
  The modified function has to have its arguments defined as
  in the following example:

  .. code-block::

      def parent(.....):
          ....

      @add_to_signature(parent)
      def child(new_param, *args, **kwargs):
          print(f'New param is {new_param})
          parent(*args, **kwargs)


  Parameters
  ----------
  func
    A function, whose parameters will be added to the resulting function signature.

  prepend
    If true, the new positional arguments will be the first.
    If false, they will be after the parents positional arguments.
    The same-named arguments just alter the properties of the "old"
    arguments, retaining their position.

  self_name
    Name of the self-parameter (for object method), that should be the first, even
    if ``prepend`` is True.
    Default self, set to None if you do not want special handling of self parameter
    (i.e. if your function has the first argument named *self* and it is not a method).

  """

  signature = inspect.signature(func)
  pars = list( signature.parameters.values())
  P = inspect.Parameter

  arg_ar = arg_kw = None
  for i in pars:
      if i.kind == P.VAR_KEYWORD:
         arg_kw = i.name
      elif i.kind == P.VAR_POSITIONAL:
         arg_ar = i.name

  def modify(mod_func):
      new_sig = inspect.signature(mod_func)
      new_pars = new_sig.parameters.values()

      #remove args/kwargs from the childs argument
      new_pars = { i.name:i for i in new_pars if i.kind != P.VAR_KEYWORD and i.kind != P.VAR_POSITIONAL }

      used = set()
      def use(i):
          nonlocal used
          used.add(i.name)
          return i

      old_pars = [ use(new_pars[i.name]) if i.name in new_pars else i for i in pars ]
      add_pars = [ i for i in new_pars.values() if i.name not in used ]

      #the easy and fast way - only keyword arguments are present, the function can be called as is,
      #just make its signature pretty
      for i in add_pars:
        if i.kind != P.KEYWORD_ONLY:
           break
      else:
        result = list(heapq.merge(old_pars, add_pars, key = lambda x: x.kind))
        result = new_sig.replace(parameters = result)
        mod_func.__signature__ = result
        return mod_func

      #the hard way, the arguments have to be rearranged
      if prepend:
         if old_pars and old_pars[0].name == self_name:
            self = old_pars[0]
            del old_pars[0]
         else:
            self = None
         result = heapq.merge(add_pars, old_pars, key = lambda x: x.kind)
         if self:
            result = itertools.chain((self,), result)
      else:
         result = heapq.merge(old_pars, add_pars, key = lambda x: x.kind)
         self = None

      result = list(result)
      result=new_sig.replace(parameters = result)

      @functools.wraps(mod_func)
      def wrapper(*args, **kwargs):
          ba=result.bind(*args, **kwargs)
          ba.apply_defaults()
          nargs = []
          nkwargs = {}
          for k,v in ba.arguments.items():
              if k==arg_kw:
                 nkwargs.update(v)
              elif k==arg_ar:
                 nargs.extend(v)
              elif k in new_pars and new_pars[k].kind != P.KEYWORD_ONLY:
                 nargs.append(v)
              else:
                 nkwargs[k]=v
          return mod_func(*nargs, **nkwargs)

      wrapper.__signature__ = result
      return wrapper

  return modify

def copy_list(src):
    """ Copy list of objects. Each object is copied just once (so
    the number of unique objects in the list is retained """

    mp = {}
    def cpy(i):
        if i in mp:
           return mp[i]
        mp[i] = copy.copy(i)
        return mp[i]

    return [ cpy(i) for i in src ]


class class_property:
    """
    Decorator that converts a method with a single cls argument into a property
    that can be accessed directly from the class.

    Example
    -------

    .. code-block::

       class Cls:
          @class_property
          def cls_property():
              return some_possibly_cached_value

          x = Cls.cls_property
    """
    def __init__(self, method=None):
        self.fget = method

    def __get__(self, instance, cls=None):
        return self.fget(cls)

class OrderedDict(_OrderedDict):

    def index(self, key):
        for i,k in enumerate(self.keys()):
            if k == key: return i
        raise KeyError(f"No suych key {key}");

    def first_item(self):
        return self[next(iter(self))]

@njit
def numpy_index(array, item):
    """ Returns index of the first occurence of the item in the array
    If numba is installed, the function is accelerated.
    """
    for idx, val in np.ndenumerate(array):
        if val == item:
            return idx
    # If no item was found return None, other return types might be a problem due to
    # numbas type inference.

if hasattr(functools, 'cached_property'):
    cached_property = functools.cached_property
else:
    #https://github.com/pydanny/cached-property/blob/master/cached_property.py
    class cached_property(object):
        """
        A property that is only computed once per instance and then replaces itself
        with an ordinary attribute. Deleting the attribute resets the property.
        Source: https://github.com/bottlepy/bottle/commit/fa7733e075da0d790d809aa3d2f53071897e6f76
        """  # noqa

        def __init__(self, func):
            self.__doc__ = getattr(func, "__doc__")
            self.func = func

        def __get__(self, obj, cls):
            if obj is None:
                return self

            if asyncio and asyncio.iscoroutinefunction(self.func):
                return self._wrap_in_coroutine(obj)

            value = obj.__dict__[self.func.__name__] = self.func(obj)
            return value

        def _wrap_in_coroutine(self, obj):
            @wraps(obj)
            @asyncio.coroutine
            def wrapper():
                future = asyncio.ensure_future(self.func(obj))
                obj.__dict__[self.func.__name__] = future
                return future

            return wrapper()

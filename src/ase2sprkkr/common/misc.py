import functools
import inspect
import heapq
import copy
import asyncio
import numpy as np
_x = object()
from collections import OrderedDict as _OrderedDict

try:
  from numba import njit
except ImportError:
  def njit(fce):
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

def add_to_signature(func, add=set()):
  """
  Add the arguments in the <func> function to the list of arguments
  of the resulting function (as keyword_only arguments)
  The modified function has to have its arguments defined as
  in the following example:

  def parent(.....):
      ....

  @add_to_signature(parent)
  def child(*args, new_param, **kwargs):
      print(f'New param is {new_param})
      parent(*args, **kwargs)
  """

  signature = inspect.signature(func)
  pars = list( signature.parameters.values())
  P = inspect.Parameter
  names = set((i.name for i in pars))

  def modify(mod_func):
      mod_sig = inspect.signature(mod_func)
      pars2 = mod_sig.parameters.values()
      forbidden = set((P.VAR_KEYWORD, P.VAR_POSITIONAL)) - add
      if forbidden:
         pars2 = (i for i in pars2 if i.kind not in forbidden and i.name not in names)
      result = heapq.merge(pars, pars2, key = lambda x: x.kind)

      @functools.wraps(mod_func)
      def wrapper(*args, **kwargs):
          return mod_func(*args, **kwargs)
      wrapper.__signature__ = mod_sig.replace(parameters = list(result))

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
    Cls:
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

import functools
import inspect
import heapq
import copy
import numpy as np
_x = object()
from collections import OrderedDict as _OrderedDict

try:
  from numba import njit
except ImportError:
  def njit(fce):
      return fce

def lazy_value(fce):
    """ The decorator for only once computed value """
    x = []

    """ Hack for staticmethod decorator, which is in fact binded by the descriptor protocol """
    if isinstance(fce, staticmethod):
        fce = fce.__func__

    @functools.wraps(fce)
    def cached_fce():
        if not x:
           x.append(fce())
        return x[0]
    return cached_fce

def add_to_signature(func, add=set()):
  """ Add the arguments in the <func> function to the list of arguments
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
    """ Returns index of the first occurence of item in the list
    If numba is installed, the function is accelerated.
    """
    for idx, val in np.ndenumerate(array):
        if val == item:
            return idx
    # If no item was found return None, other return types might be a problem due to
    # numbas type inference.

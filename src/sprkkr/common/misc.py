import functools
import copy
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

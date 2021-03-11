import functools
_x = object()
from collections import OrderedDict

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


class classproperty:
    """
    Decorator that converts a method with a single cls argument into a property
    that can be accessed directly from the class.
    """
    def __init__(self, method=None):
        self.fget = method

    def __get__(self, instance, cls=None):
        return self.fget(cls)

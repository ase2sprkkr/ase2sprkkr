""" Various decorators, mainly for class methods """
import functools
import itertools
import inspect
import heapq
import asyncio

if hasattr(functools,'cache'):
   cache = functools.cache
   """ Functools.cache. Python 3.8 and earlier does not have this method, so it is mocked for this version of python. """
else:
   cache = functools.lru_cache(maxsize=None)
   """ Functools.cache. Python 3.8 and earlier does not have this method, so it is mocked for this version of python. """


class cached_class_property:
    """
    Decorator that converts a method with a single cls argument into a property
    that can be accessed directly from the class. The value is computed only
    once.

    Example
    -------

    .. code-block::

       class Cls:
          @class_property
          def cls_property():
              return some_cached_value

          @class_property
          def another_cls_property(cls):
              return another_cached_value

          x = Cls.cls_property
    """
    def __init__(self, method=None):
        self.fget = method
        functools.update_wrapper(self, method)

    def __set_name__(self, owner, name):
        self.owner = owner

    def __get__(self, instance, cls=None):
        if self.fget.__code__.co_argcount:
           out = self.fget(cls)
           owner = cls
        else:
           out = self.fget()
           owner = self.owner
        setattr(owner, self.fget.__name__, out)
        return out


if hasattr(functools, 'cached_property'):
    cached_property = functools.cached_property
    """ Functools.cached_property decorator - the value is computed only once and then stored as an (same-name) instance attribute.
    You can delete the attribute to invalidate the cache.

    Earlier versions of python do not have this decorator, so for these versions is implemented below, otherwise it is taken from functools.
    """
else:
    #https://github.com/pydanny/cached-property/blob/master/cached_property.py
    class cached_property(object):
        """
        A property that is only computed once per instance and then replaces itself
        with an ordinary attribute. Deleting the attribute resets the property.
        Source: `https://github.com/bottlepy/bottle/commit/fa7733e075da0d790d809aa3d2f53071897e6f76`
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

class class_property:
    """
    Decorator that converts a method with a single cls argument into a property
    that can be accessed directly from the class. The value is not cached, see
    :func:`cached_class_property`, if caching is desirable.

    Example
    -------

    .. code-block::

       class Cls:
          @class_property
          def cls_property():
              return some_value

          x = Cls.cls_property
    """
    def __init__(self, method=None):
        self.fget = method

    def __get__(self, instance, cls=None):
        return self.fget(cls)


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

def add_called_class_as_argument(decorator):
    """ If a decorator is used on a method, the information about the defining class is lost.
    As a consequence, the execution of ``super().... failed.

    This decorator can be applied to another decorator: it will add the defining class as a first
    argument of the resulting class.

    The implementation is very tricky - function-style decorator is not able to handle this situation,
    the __set_name__ special method (in class-style decorator) has to be utilized. Moreover, class
    style decorator lose 'self' from decorated method, thus, the self has to be recovered using
    descriptor protocol.
    To speed up the whole thing, the first access replaces the class-style decorator class with
    generated method (that already can know - and knows - everything needed)

    Example

    >>> @add_called_class_as_argument
    ... def decorator(func):
    ...
    ...       def wrapped(cls, self):
    ...            super(cls, self).call()
    ...            func(self)
    ...       return wrapped
    ...
    >>> class A:
    ...    def call(self):
    ...         print('A', end='')
    ...
    >>> class B(A):
    ...      @decorator
    ...      def call(self):
    ...           print('B', end='')
    ...
    >>> B().call()
    AB
    """

    class AddCalledClassAsArgument:

        def __init__(self, function):
            self.func = function

        def __set_name__(self, owner, name):
            """ Catch the defining class """
            self.cls=owner
            self.name=name

        def __call__(self, instance, *args, **kwargs):
             func = self.func
             cls = self.cls

             @functools.wraps(func)
             def wrapped(self, *args, **kwargs):
                 return func(cls, self, *args, **kwargs)

             setattr(cls, self.name, wrapped)
             return wrapped(instance, *args, **kwargs)

        def __get__(self, instance, instancetype):
           """ Descriptor method - via this method the actual calling instance
           (self of the very resulting method) is obtained. """
           out = functools.partial(self.__call__, instance)
           return functools.update_wrapper(out, self.func)

    #@functools.wraps(function)
    def wrapper(function):
        function = decorator(function)
        return AddCalledClassAsArgument(function)

    return wrapper
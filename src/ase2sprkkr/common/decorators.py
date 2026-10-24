"""Various decorators, mainly for class methods"""

import functools
import itertools
import inspect
import heapq
import warnings
from . import backward_compatibility  # NOQA
from functools import cache  # NOQA


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


std_cached_property = functools.cached_property
"""``functools.cached_property`` used as the base implementation."""


class cached_property(std_cached_property):
    """Cached property that allows a setter and deleter."""

    def __init__(self, func, fset=None, fdel=None):
        super().__init__(func)
        self.fset = fset
        self.fdel = fdel

    def __set__(self, obj, value):
        if self.fset:
            self.fset(obj, value)
        else:
            obj.__dict__[self.attrname] = value

    def __delete__(self, obj):
        if self.fdel:
            self.fdel(obj)
        else:
            del obj.__dict__[self.attrname]

    def setter(self, fset):
        return type(self)(self.func, fset, self.fdel)

    def deleter(self, fdel):
        return type(self)(self.func, self.fset, fdel)


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


def add_to_signature(func, prepend=False, old_as_keyword=False, excluding=None, kwargs=False, apply_defaults=False):
    """
    Add the arguments in the ``func`` function to the list of arguments
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
      If True, the new positional arguments will be the first.
      If False, they will be after the parents positional arguments.

      If prepend=True, the order of the same-named parameters is used from the new function.
      Otherwise, it the different order of parameters does not overwrite the order in the
      "parent function"

    old_as_keyword
      If True, the original POSITIONAL_OR_KEYOWRD will be changed to keyword

    excluding
      Do not add arguments of this name(s)

    kwargs
      The resulting function will have kwargs argument

    apply_defaults
      If true, the defaults values from the "parent signature" are applied before
      the function is called. Otherwise, they could not applied (for the sake of
      the simplicity of call-wrapper) and it is presumed,
      that the default will be aplied by calling the func.
    """
    if not excluding:
        excluding = []
    elif isinstance(excluding, str):
        excluding = [excluding]

    P = inspect.Parameter

    def modify(mod_func):
        new_sig = inspect.signature(mod_func)
        new_pars = new_sig.parameters.values()
        new_pars = {str(i.name): i for i in new_pars if (i.kind != P.VAR_KEYWORD or kwargs)}

        signature = inspect.signature(func)
        pars = signature.parameters.values()

        if old_as_keyword:
            pars = [
                i.replace(kind=P.KEYWORD_ONLY)
                if i.kind == P.POSITIONAL_OR_KEYWORD else i
                for i in pars
                if i.name not in excluding
            ]
            pars = { i.name:i for i in pars }
        else:
            pars = {
                i.name: i for i in pars if i.name not in excluding
            }


        if prepend:
            old_pars = new_pars
            add_pars = [ i for name,i in pars.items() if name not in new_pars ]
        else:
            old_pars = { name : (new_pars[name] if name in new_pars else i) for name, i in pars.items() }
            add_pars = [ i for i in new_pars.values() if i.name not in old_pars ]

        result = [*old_pars.values(), *add_pars]
        result.sort(key = lambda p: p.kind.value)
        result = new_sig.replace(parameters=result)
        mod_func.__signature__ = result

        same = P.KEYWORD_ONLY if apply_defaults else True
        for a,b in itertools.zip_longest(result.parameters.values(), new_sig.parameters.values()):
            if a.kind in (P.KEYWORD_ONLY, P.VAR_KEYWORD):
                if b is not None and b.kind not in (P.KEYWORD_ONLY, P.VAR_KEYWORD):
                    same = b.kind
                break
            if b is None:
                same = a.kind
                break
            if a.kind == P.VAR_POSITIONAL:
                if b.kind != P.VAR_POSITIONAL:
                    same = min(a.kind, b.kind)
                break
            if a.name != b.name or a.kind == P.POSITIONAL_OR_KEYWORD and b.kind != P.POSITIONAL_OR_KEYWORD:
                same = min(a.kind, b.kind)
                break

        if same is True:
            return mod_func

        if same > P.POSITIONAL_OR_KEYWORD:
            if apply_defaults:
                def wrapper(*args, **kwargs):
                    ba = result.bind(*args, **kwargs)
                    ba.apply_defaults()
                    return mod_func(*ba.args, **ba.kwargs)
            else:
                def wrapper(*args, **kwargs):
                    ba = result.bind(*args, **kwargs)
                    return mod_func(*ba.args, **ba.kwargs)
            return functools.wraps(mod_func)(wrapper)

        if same == P.POSITIONAL_OR_KEYWORD and not any(i.kind in (P.POSITIONAL_ONLY, P.VAR_POSITIONAL) for i in new_pars.values()):
            if apply_defaults:
                def wrapper(*args, **kwargs):
                    ba = result.bind(*args, **kwargs)
                    ba.apply_defaults()
                    return mod_func(**ba.arguments)
            else:
                def wrapper(*args, **kwargs):
                    ba = result.bind(*args, **kwargs)
                    return mod_func(**ba.arguments)
            return functools.wraps(mod_func)(wrapper)

        @functools.wraps(mod_func)
        def wrapper(*args, **kwargs):
            x=new_pars
            y=pars
            z=old_pars
            gg=add_pars


            ba = result.bind(*args, **kwargs)
            ba.apply_defaults()
            args = ba.arguments
            nargs, nkwargs = [], {}

            for k, p in new_sig.parameters.items():
                if p.kind in (P.POSITIONAL_ONLY, P.POSITIONAL_OR_KEYWORD):
                    nargs.append(args.pop(k))
                elif p.kind == P.VAR_POSITIONAL:
                    np = result.parameters
                    nargs.extend((args.pop(n) for n in list(args.keys()) if np[n].kind in (P.POSITIONAL_ONLY, P.POSITIONAL_OR_KEYWORD) ))
                    nargs.extend(args.pop(k, []))
                elif p.kind == P.VAR_KEYWORD:
                    nkwargs.update(args.pop(k, {}))
                else:
                    nkwargs[k] = args.pop(k)
            nkwargs.update(ba.arguments)
            return mod_func(*nargs, **nkwargs)

        wrapper.__signature__ = result
        return wrapper

    return modify


def add_called_class_as_argument(decorator):
    """If a decorator is used on a method, the information about the defining class is lost.
    As a consequence, the execution of ``super()....`` failed.

    This decorator can be applied to another decorator: it will add the defining class as a first
    argument of the resulting function.

    The implementation is very tricky - function-style decorator is not able to handle this situation,
    the __set_name__ special method (in class-style decorator) has to be utilized. Moreover, class
    style decorator lose 'self' from decorated method, thus, the self has to be recovered using
    descriptor protocol.
    To speed up the whole thing, the first access replaces the class-style decorator class with
    generated method (that can already know - and knows - everything needed)

    Example

    >>> @add_called_class_as_argument
    ... def decorator(func):
    ...
    ...     def wrapped(cls, self):
    ...         super(cls, self).call()
    ...         func(self)
    ...
    ...     return wrapped
    >>> class A:
    ...     def call(self):
    ...         print("A", end="")
    >>> class B(A):
    ...     @decorator
    ...     def call(self):
    ...         print("B", end="")
    >>> B().call()
    AB
    """

    class AddCalledClassAsArgument:
        def __init__(self, function):
            self.func = function

        def __set_name__(self, owner, name):
            """Catch the defining class"""
            self.cls = owner
            self.name = name

        def __call__(self, instance, *args, **kwargs):
            func = self.func
            cls = self.cls

            @functools.wraps(func)
            def wrapped(self, *args, **kwargs):
                return func(cls, self, *args, **kwargs)

            setattr(cls, self.name, wrapped)
            return wrapped(instance, *args, **kwargs)

        def __get__(self, instance, instancetype):
            """Descriptor method - via this method the actual calling instance
            (self of the very resulting method) is obtained."""
            out = functools.partial(self.__call__, instance)
            return functools.update_wrapper(out, self.func)

    def wrapper(function):
        function = decorator(function)
        return AddCalledClassAsArgument(function)

    return wrapper


def warnings_from_here(stacklevel=1):
    stacklevel += 1

    def wrapper(func):

        @functools.wraps(func)
        def wrapped(*args, **kwargs):
            with warnings.catch_warnings(record=True) as warning_list:
                result = func(*args, **kwargs)
            for warning in warning_list:
                warnings.warn(warning.message, warning.category, stacklevel=stacklevel)
            return result

        return wrapped

    return wrapper


class maybeclassmethod(object):
    """This decorator creates a method, that can behawes both as
    classmethod and normal method - it provides both self and cls,
    when self can be None"""

    def __init__(self, method):
        self.method = method
        self.wrapper = functools.wraps(method)

    def __get__(self, instance, cls):
        return self.wrapper(lambda *args, **kw: self.method(instance, cls, *args, **kw))

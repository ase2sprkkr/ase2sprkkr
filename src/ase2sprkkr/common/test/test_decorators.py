import inspect
import pytest

from ..decorators import add_to_signature


def test_prepend_respects_excluding():
    def parent(a=1, b=2):
        return a, b

    @add_to_signature(parent, prepend=True, excluding="b")
    def child(c=3, *args, **kwargs):
        return c, parent(b=5, *args, **kwargs)

    signature = inspect.signature(child)
    assert "a" in signature.parameters
    assert "b" not in signature.parameters
    assert child("c", "a") == ("c", ("a", 5))


def test_override_parent_kind():
    def parent(a: int = 1, b=2):
        return a, b

    @add_to_signature(parent, prepend=True)
    def child(c, *args, a: str = "new", **kwargs):
        return c, parent(a, *args, **kwargs)

    parameter = inspect.signature(child).parameters["a"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default == "new"
    assert parameter.annotation is str
    assert child("c", "b") == ("c", ("new", "b") )

    with pytest.raises(ValueError):
        @add_to_signature(parent)
        def child(c, *args, a: str = "new", **kwargs):
            return c, parent(a, *args, **kwargs)

    @add_to_signature(parent)
    def child(c=3, *args, a: str = "new", **kwargs):
        return c, parent(a, *args, **kwargs)

    parameter = inspect.signature(child).parameters["a"]
    assert parameter.kind is inspect.Parameter.KEYWORD_ONLY
    assert parameter.default == "new"
    assert parameter.annotation is str
    assert child("b", "c") == ("c", ("new", "b") )
    assert child("b", "c", a='a') == ("c", ("a", "b") )


def test_add_to_signature():
    def f(a, b=2, **kwargs):
        return a, b, kwargs

    @add_to_signature(f)
    def ff(c=2, *args, **kwargs):
        return c, f(*args, **kwargs)

    assert ff("a", c="c", d="d", e="e") == ("c", ("a", 2, {"d": "d", "e": "e"}))
    assert ff("a", "b", "c", d="d", e="e") == ("c", ("a", "b", {"d": "d", "e": "e"}))
    assert "a" in inspect.signature(ff).parameters
    assert list(inspect.signature(ff).parameters.values())[0].name == "a"

    @add_to_signature(f, prepend=True)
    def ff(c, *args, **kwargs):
        return c, f(*args, **kwargs)

    assert ff("c", "a", d="d", e="e") == ("c", ("a", 2, {"d": "d", "e": "e"}))

    def f(c=1, d=5):
        return c, d

    @add_to_signature(f, prepend=True)
    def ff(a=2, d="d", *args, **kwargs):
        return a, f(d=d, *args, **kwargs)

    assert ff("a", "c") == ("a", (1, "c"))

    @add_to_signature(f, prepend=True)
    def ff(a=2, *args, d="d", **kwargs):
        return a, f(d=d, *args, **kwargs)

    assert ff("a", "c") == ("a", ("c", "d"))

    def f(*args):
        return args

    @add_to_signature(f)
    def ff(a, *args):
        return a, f(*args)

    assert ff("a", "c", "d") == ("a", ("c", "d"))

    @add_to_signature(lambda: None)
    def ff(*args, c, **kwargs):
        return c

    assert ff(c=1) == 1

    @add_to_signature(lambda x: x)
    def ff(*args, c, **kwargs):
        return c, args

    assert ff(2, c=1) == (1, (2,))

    @add_to_signature(lambda x,y='y': x)
    def ff(*args, c, **kwargs):
        return c, args, kwargs

    assert ff(x=2, c=1) == (1, (2, 'y'), {})

    @add_to_signature(lambda: None)
    def ff(a, *args, c, **kwargs):
        return a, c

    assert ff(2, c=1) == (2, 1)

    def f(a=1, b=2):
        return a, b

    @add_to_signature(f, excluding="b")
    def ff(c=1, *args, **kwargs):
        return c, f(b=5, *args, **kwargs)

    assert ff(3, 7) == (7, (3, 5))
    assert ff(3, c=7) == (7, (3, 5))
    assert ff(a=3, c=7) == (7, (3, 5))




def test_order():
    def parent(a, b):
        pass

    @add_to_signature(parent)
    def child(a, c, b):
        return a,b,c

    assert (1,2,3) == child(1,2,3)

    @add_to_signature(parent, prepend=True)
    def child(a, c, b):
        return a,b,c

    assert (1,3,2) == child(1,2,3)

    @add_to_signature(parent)
    def child(a, c, b, *, d=1):
        return a,b,c,d

    assert (1,2,3,4) == child(1,2,3,d=4)

    @add_to_signature(parent, prepend=True)
    def child(a, c, b, *, d=1):
        return a,b,c,d

    assert (1,3,2,4) == child(1,2,3,d=4)


    def parent(a, b, *, e=3):
        pass

    @add_to_signature(parent)
    def child(a, c, b, *, d=1, **kwargs):
        return a,b,c,d, kwargs

    assert (1,2,3,4,{}) == child(1,2,3,d=4)

    @add_to_signature(parent, apply_defaults=True)
    def child(a, c, b, *, d=1, **kwargs):
        return a,b,c,d, kwargs

    assert (1,2,3,4, {'e':3}) == child(1,2,3,d=4)


    def parent(a, b, *, e=3):
        pass

    @add_to_signature(parent, prepend=True)
    def child(a, c, b, *, d=1):
        return a,b,c,d

    assert (1,2,3,4) == child(1,3,2,d=4)

    def parent(a, /, b=2, *, e=3):
        pass

    @add_to_signature(parent, prepend=True)
    def child(a=1, /, **kwargs) :
        return a, kwargs

    assert child() == (1, { 'b':2, 'e':3})
    assert child(b=4) == (1, { 'b':4, 'e':3})
    assert child(e=4) == (1, { 'b':2, 'e':4})

    @add_to_signature(parent, prepend=True)
    def child(a=1, n=2, /, **kwargs) :
        return a, n, kwargs

    assert child() == (1,2, {'b': 2, 'e': 3})
    assert child(b=4) == (1,2, {'b':4, 'e': 3})
    assert child(e=4) == (1,2, {'b': 2, 'e':4})
    assert child(5,6,e=4) == (5,6, {'b': 2, 'e':4})
    assert child(5,e=4) == (5,2, {'b': 2, 'e':4})

    @add_to_signature(parent)
    def child(n=2, /, **kwargs) :
        return n, kwargs

    assert child(1) == (2, {'a':1, 'b':2, 'e':3})
    assert child(1, b=4) == (2, {'a':1, 'b':4, 'e':3})
    assert child(1, e=4) == (2, {'a':1, 'b':2, 'e':4})
    assert child(5,6,e=4) == (6, {'a':5, 'b':2, 'e':4})
    assert child(5,e=4) == (2, {'a':5, 'b':2, 'e':4})

    @add_to_signature(parent, apply_defaults=True)
    def child(n=2, /, **kwargs) :
        return n, kwargs

    assert child(1) == (2, {'a':1, 'b':2, 'e':3})
    assert child(1, b=4) == (2, {'a':1, 'b':4, 'e':3})
    assert child(1, e=4) == (2, {'a':1, 'b':2, 'e':4})
    assert child(5,6,e=4) == (6, {'a':5, 'b':2, 'e':4})
    assert child(5,e=4) == (2, {'a':5, 'b':2, 'e':4})

    @add_to_signature(parent, prepend=True)
    def child(n, /, **kwargs) :
        return n, kwargs

    with pytest.raises(TypeError):
        child(1)
    with pytest.raises(TypeError):
        child(1, b=4, e=4)
    assert child(5,6,e=4) == (5, {'a':6, 'b':2, 'e':4})

    def parent(a=9, /, b=2, *, e=3):
        pass

    @add_to_signature(parent, prepend=True)
    def child(n, /, **kwargs) :
        return n, kwargs

    assert child(1) == (1, {'a':9, 'b':2, 'e':3})
    assert child(1, b=4) == (1, {'a':9, 'b':4, 'e':3})
    assert child(1, e=4) == (1, {'a':9, 'b':2, 'e':4})
    assert child(5,6,e=4) == (5, {'a':6, 'b':2, 'e':4})
    assert child(5,e=4) == (5, {'a':9, 'b':2, 'e':4})

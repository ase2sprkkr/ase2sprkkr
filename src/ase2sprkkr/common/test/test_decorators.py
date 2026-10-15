import inspect

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


def test_prepend_keyword_only_override_keeps_parent_kind():
    def parent(a: int = 1, b=2):
        return a, b

    @add_to_signature(parent, prepend=True)
    def child(c, *args, a: str = "new", **kwargs):
        return c, a, kwargs

    parameter = inspect.signature(child).parameters["a"]
    assert parameter.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
    assert parameter.default == "new"
    assert parameter.annotation is str
    assert child("c", "value") == ("c", "value", {"b": 2})

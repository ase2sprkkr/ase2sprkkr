import inspect

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ..decorators import add_to_signature  # NOQA: E402


class TestDecorators(TestCase):
    def test_prepend_respects_excluding(self):
        def parent(a=1, b=2):
            return a, b

        @add_to_signature(parent, prepend=True, excluding="b")
        def child(c=3, *args, **kwargs):
            return c, parent(b=5, *args, **kwargs)

        signature = inspect.signature(child)
        self.assertTrue("a" in signature.parameters)
        self.assertFalse("b" in signature.parameters)
        self.assertEqual(("c", ("a", 5)), child("c", "a"))

    def test_prepend_keyword_only_override_keeps_parent_kind(self):
        def parent(a: int = 1, b=2):
            return a, b

        @add_to_signature(parent, prepend=True)
        def child(c, *args, a: str = "new", **kwargs):
            return c, a, kwargs

        parameter = inspect.signature(child).parameters["a"]
        self.assertEqual(inspect.Parameter.POSITIONAL_OR_KEYWORD, parameter.kind)
        self.assertEqual("new", parameter.default)
        self.assertEqual(str, parameter.annotation)
        self.assertEqual(("c", "value", {"b": 2}), child("c", "value"))

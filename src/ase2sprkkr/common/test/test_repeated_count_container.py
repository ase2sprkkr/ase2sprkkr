import pyparsing as pp

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

from ...common.grammar import generate_grammar
from ...input_parameters import input_parameters_definitions as cd  # NOQA: E402
from ...common.repetition import Repeated

V = cd.InputValueDefinition

class S(cd.InputSectionDefinition):
    custom_class = None

class IP(cd.InputParametersDefinition):
    custom_class = None


class TestRepeatedCountParsing(TestCase):
    def test_repeated_count_fixed_int(self):
        item = S("ITEM", members=[V("VAL", int, name_in_grammar=False)])
        parent = S("PARENT", members=[V("N", int), item])
        root = IP([parent])

        # First, test not-repeated at all
        data = "PARENT\n\tN=3\n\tITEM\n\t\t1\n\t"
        parsed = root.read_from_string(data)
        p = parsed.PARENT
        self.assertEqual(p.ITEM.VAL(), 1)

        data = "PARENT\n\tITEM\n\t\t10\n\tITEM\n\t\t20\n"
        with self.assertRaises(Exception):
            parent.parse(data)

        # Then, test repeated, any number of repetition
        item = S("ITEM", members=[V("VAL", int, name_in_grammar=False)], is_repeated=Repeated.REPEATED, repeated_with_name=True)
        parent = S("PARENT", members=[V("N", int), item])
        root = IP([parent])

        # repeated-header style: two ITEM blocks each with one unnamed integer
        data = "PARENT\n\tITEM\n\t\t10\n\tITEM\n\t\t20\n"
        parsed = root.read_from_string(data)
        data = "PARENT\n\tITEM\n\t\t10\n\tITEM\n\t\t20\n\tITEM\n\t\t30\n"
        parsed = root.read_from_string(data)

        p = parsed.PARENT
        items = p.ITEM
        self.assertEqual(len(items), 3)
        self.assertEqual(items[0]['VAL'](), 10)
        self.assertEqual(items[1]['VAL'](), 20)
        self.assertEqual(items[2]['VAL'](), 30)

        # End-to-end parse: ITEM is repeated exactly twice
        # Define VAL as a nameless line inside ITEM to avoid duplicate-name collisions
        item = S("ITEM", members=[V("VAL", int, name_in_grammar=False)], is_repeated=Repeated.REPEATED, repeated_count=2, repeated_with_name=True)
        parent = S("PARENT", members=[V("N", int), item])
        root = IP([parent])

        # repeated-header style: two ITEM blocks each with one unnamed integer
        data = "PARENT\n\tITEM\n\t\t10\n\tITEM\n\t\t20\n"
        parsed = root.read_from_string(data)
        p = parsed.PARENT
        items = p.ITEM
        self.assertEqual(len(items), 2)
        self.assertEqual(items[0]['VAL'](), 10)
        self.assertEqual(items[1]['VAL'](), 20)

        # too few items -> expect parse error
        data_few = "PARENT\n\tITEM\n\t\t10\n"
        with self.assertRaises(Exception):
            parent.parse(data_few)

        # too many items -> expect parse error
        data_many = "PARENT\n\tITEM\n\t\t10\n\tITEM\n\t\t20\n\tITEM\n\t\t30\n"
        with self.assertRaises(Exception):
            parent.parse(data_many)

    def test_repeated_count_from_dependency(self):
        # End-to-end parse: ITEM repetition count depends on N
        item = S("ITEM", members=[V("VAL", int, name_in_grammar=False)], is_repeated=Repeated.REPEATED, repeated_count="N", repeated_with_name=True)
        parent = S("PARENT", members=[V("N", int), item])
        root = IP([parent])

        # repeated-header style: N before, then N ITEM blocks each with one unnamed integer
        data = "PARENT\n\tN=3\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n\tITEM\n\t\t3\n"
        parsed = root.read_from_string(data)
        p = parsed.PARENT
        items = p.ITEM
        self.assertEqual(len(items), 3)
        self.assertEqual([it['VAL']() for it in items.values()], [1, 2, 3])

        # too few
        data_few = "PARENT\n\tN=3\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n"
        with self.assertRaises(Exception):
            parent.parse(data_few)

        # too many
        data_many = "PARENT\n\tN=3\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n\tITEM\n\t\t3\n\tITEM\n\t\t4\n"
        with self.assertRaises(Exception):
            parent.parse(data_many)

    def test_repeated_count_callable(self):
        # End-to-end parse: callable doubles the given N -> 2 -> 4 items
        def double_n(N):
            return int(N) * 2

        item = S("ITEM", members=[V("VAL", int, name_in_grammar=False)], is_repeated=Repeated.REPEATED, repeated_count=double_n, repeated_with_name=True)
        parent = S("PARENT", members=[V("N", int), item])
        root = IP([parent])

        # repeated-header style: N=2 then 4 ITEM blocks in total because callable doubles N
        data = "PARENT\n\tN=2\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n\tITEM\n\t\t3\n\tITEM\n\t\t4\n"
        parsed = root.read_from_string(data)
        p = parsed.PARENT
        items = p.ITEM
        self.assertEqual(len(items), 4)
        self.assertEqual([it['VAL']() for it in items.values()], [1, 2, 3, 4])

        # too few
        data_few = "PARENT\n\tN=2\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n"
        with self.assertRaises(Exception):
            parent.parse(data_few)

        # too many
        data_many = "PARENT\n\tN=2\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n\tITEM\n\t\t3\n\tITEM\n\t\t4\n\tITEM\n\t\t5\n"
        with self.assertRaises(Exception):
            parent.parse(data_many)

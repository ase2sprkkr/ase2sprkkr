import pytest

from ...common.repetition import Repeated
from ...input_parameters import input_parameters_definitions as cd


V = cd.InputValueDefinition


class S(cd.InputSectionDefinition):
    custom_class = None


class IP(cd.InputParametersDefinition):
    custom_class = None


class TestRepeatedCountParsing:
    def test_container_callbacks_are_deferred_until_first_use(self):
        class CountingValue(V):
            def __init__(self, *args, **kwargs):
                self.added_count = 0
                super().__init__(*args, **kwargs)

            def added_to_container(self, container):
                self.added_count += 1
                super().added_to_container(container)

        value = CountingValue("VALUE", int)
        section = S("SECTION", members=[value])

        assert value.container is None
        assert value.added_count == 0

        root = IP([section])
        assert value.container is section
        assert value.added_count == 1

        root.create_object()
        assert value.added_count == 1

    def test_repeated_count_fixed_int(self):
        item = S("ITEM", members=[V("VAL", int, name_in_grammar=False)])
        parent = S("PARENT", members=[V("N", int), item])
        root = IP([parent])

        data = "PARENT\n\tN=3\n\tITEM\n\t\t1\n\t"
        parsed = root.read_from_string(data)
        assert parsed.PARENT.ITEM.VAL() == 1

        data = "PARENT\n\tITEM\n\t\t10\n\tITEM\n\t\t20\n"
        with pytest.raises(Exception):
            parent.parse(data)

        item = S(
            "ITEM",
            members=[V("VAL", int, name_in_grammar=False)],
            is_repeated=Repeated.REPEATED,
            repeated_with_name=True,
        )
        parent = S("PARENT", members=[V("N", int), item])
        root = IP([parent])

        data = "PARENT\n\tN=2\n\tITEM\n\t\t10\n\tITEM\n\t\t20\n\tITEM\n\t\t30\n"
        parsed = root.read_from_string(data)
        items = parsed.PARENT.ITEM
        assert len(items) == 3
        assert items[0]["VAL"]() == 10
        assert items[1]["VAL"]() == 20
        assert items[2]["VAL"]() == 30

        item = S(
            "ITEM",
            members=[V("VAL", int, name_in_grammar=False)],
            is_repeated=Repeated.REPEATED,
            repeated_count=2,
            repeated_with_name=True,
        )
        parent = S("PARENT", members=[V("N", int), item])
        root = IP([parent])

        data = "PARENT\n\tN=2\n\tITEM\n\t\t10\n\tITEM\n\t\t20\n"
        parsed = root.read_from_string(data)
        items = parsed.PARENT.ITEM
        assert len(items) == 2
        assert items[0]["VAL"]() == 10
        assert items[1]["VAL"]() == 20

        data_few = "PARENT\n\tITEM\n\t\t10\n"
        with pytest.raises(Exception):
            parent.parse(data_few)

        data_many = "PARENT\n\tN=2\n\tITEM\n\t\t10\n\tITEM\n\t\t20\n\tITEM\n\t\t30\n"
        with pytest.raises(Exception):
            parent.parse(data_many)

    def test_repeated_count_from_dependency(self):
        item = S(
            "ITEM",
            members=[V("VAL", int, name_in_grammar=False)],
            is_repeated=Repeated.REPEATED,
            repeated_count="N",
            repeated_with_name=True,
        )
        parent = S("PARENT", members=[V("N", int), item])
        root = IP([parent])

        data = "PARENT\n\tN=3\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n\tITEM\n\t\t3\n"
        parsed = root.read_from_string(data)
        items = parsed.PARENT.ITEM
        assert len(items) == 3
        assert [it["VAL"]() for it in items.values()] == [1, 2, 3]

        data_few = "PARENT\n\tN=3\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n"
        with pytest.raises(Exception):
            parent.parse(data_few)

        data_many = "PARENT\n\tN=3\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n\tITEM\n\t\t3\n\tITEM\n\t\t4\n"
        with pytest.raises(Exception):
            parent.parse(data_many)

    def test_repeated_count_one_can_be_parsed_repeatedly(self):
        item = S(
            "ITEM",
            members=[V("VAL", int, name_in_grammar=False)],
            is_repeated=Repeated.REPEATED,
            repeated_count="N",
            repeated_with_name=True,
        )
        root = IP([S("PARENT", members=[V("N", int), item])])
        data = "PARENT\n\tN=1\n\tITEM\n\t\t10\n"

        first = root.read_from_string(data)
        second = root.read_from_string(data)

        assert first.PARENT.ITEM[0].VAL() == 10
        assert second.PARENT.ITEM[0].VAL() == 10

    def test_repeated_count_callable(self):
        def double_n(N):
            return int(N) * 2

        item = S(
            "ITEM",
            members=[V("VAL", int, name_in_grammar=False)],
            is_repeated=Repeated.REPEATED,
            repeated_count=double_n,
            repeated_with_name=True,
        )
        parent = S("PARENT", members=[V("N", int), item])
        root = IP([parent])

        data = "PARENT\n\tN=2\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n\tITEM\n\t\t3\n\tITEM\n\t\t4\n"
        parsed = root.read_from_string(data)
        items = parsed.PARENT.ITEM
        assert len(items) == 4
        assert [it["VAL"]() for it in items.values()] == [1, 2, 3, 4]

        data_few = "PARENT\n\tN=2\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n"
        with pytest.raises(Exception):
            parent.parse(data_few)

        data_many = "PARENT\n\tN=2\n\tITEM\n\t\t1\n\tITEM\n\t\t2\n\tITEM\n\t\t3\n\tITEM\n\t\t4\n\tITEM\n\t\t5\n"
        with pytest.raises(Exception):
            parent.parse(data_many)

    def test_repeated_count_callable_with_multiple_dependencies(self):
        def total(N, M):
            return int(N) + int(M)

        item = S(
            "ITEM",
            members=[V("VAL", int, name_in_grammar=False)],
            is_repeated=Repeated.REPEATED,
            repeated_count=total,
            repeated_with_name=True,
        )
        parent = S("PARENT", members=[V("N", int), V("M", int), item])
        root = IP([parent])

        parsed = root.read_from_string(
            "PARENT\n\tN=2\n\tM=1\n"
            "\tITEM\n\t\t10\n"
            "\tITEM\n\t\t20\n"
            "\tITEM\n\t\t30\n"
        )
        assert [instance["VAL"]() for instance in parsed.PARENT.ITEM.values()] == [10, 20, 30]

    def test_repeated_count_dependency_from_parent_section(self):
        item = S(
            "ITEM",
            members=[V("VAL", int, name_in_grammar=False)],
            is_repeated=Repeated.REPEATED,
            repeated_count="..N",
            repeated_with_name=True,
        )
        child = S("CHILD", members=[V("N", int), item])
        parent = S("PARENT", members=[V("N", int), child])
        root = IP([parent])

        parsed = root.read_from_string(
            "PARENT\n\tN=2\n\tCHILD\n"
            "\t\tN=4\n"
            "\t\tITEM\n\t\t\t10\n"
            "\t\tITEM\n\t\t\t20\n"
        )

        assert child.get_member("N") is child["N"]
        assert child.get_member("..N") is parent["N"]
        assert len(parent["N"].grammar_hooks) == 1
        assert child["N"].grammar_hooks == []
        assert parsed.PARENT.CHILD.get_member("N") is parsed.PARENT.CHILD["N"]
        assert parsed.PARENT.CHILD.get_member("..N") is parsed.PARENT["N"]
        assert parsed.PARENT.CHILD.ITEM.get_member("..N") is parsed.PARENT.CHILD["N"]
        assert parsed.PARENT.CHILD.ITEM.get_member("....N") is parsed.PARENT["N"]
        assert parsed.PARENT.CHILD.ITEM[0].get_member("..N") is parsed.PARENT.CHILD["N"]
        assert parsed.PARENT.CHILD.ITEM.get_member(0) is parsed.PARENT.CHILD.ITEM[0]
        assert parsed.PARENT.CHILD.ITEM.get_member("0.VAL") is parsed.PARENT.CHILD.ITEM[0]["VAL"]
        assert [instance["VAL"]() for instance in parsed.PARENT.CHILD.ITEM.values()] == [10, 20]

        rows = S(
            "ROWS",
            members=[V("VAL", int)],
            is_repeated=Repeated.DICT_SECTION,
        )
        dictionary = IP([rows]).create_object()
        dictionary.ROWS.set({1: {"VAL": 10}})
        assert dictionary.ROWS.get_member(1) is dictionary.ROWS[1]
        assert dictionary.ROWS.get_member("1.VAL") is dictionary.ROWS[1]["VAL"]

import pytest

from ..debug import backtrace_pyparsing, check_whitespaces


class GrammarNode:
    def __init__(self, *, expr=None, exprs=None, white_chars=None):
        if expr is not None:
            self.expr = expr
        if exprs is not None:
            self.exprs = exprs
        self.whiteChars = frozenset({"\t", " ", "\r"}) if white_chars is None else white_chars


def test_backtrace_pyparsing_keeps_callback_for_nested_terms():
    leaf = GrammarNode()
    branch = GrammarNode(expr=leaf)
    root = GrammarNode(exprs=[branch])

    visited = []
    backtrace_pyparsing(root, visited.append)

    assert visited == [root, branch, leaf]


def test_backtrace_pyparsing_handles_recursive_grammar():
    root = GrammarNode()
    root.expr = root

    visited = []
    backtrace_pyparsing(root, visited.append)

    assert visited == [root]


def test_check_whitespaces_checks_nested_terms():
    bad = GrammarNode(white_chars={" "})
    root = GrammarNode(exprs=[GrammarNode(expr=bad)])

    with pytest.raises(ValueError, match="Bad white chars"):
        check_whitespaces(root)


def test_check_whitespaces_accepts_matching_tree():
    root = GrammarNode(exprs=[GrammarNode(expr=GrammarNode())])

    check_whitespaces(root)

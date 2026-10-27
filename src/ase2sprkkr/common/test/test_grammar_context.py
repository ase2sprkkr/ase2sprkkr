from concurrent.futures import ThreadPoolExecutor
from threading import Event

import pyparsing as pp
import pytest

from ..grammar import generate_grammar


def parser_defaults():
    return pp.ParserElement.DEFAULT_WHITE_CHARS, pp.Keyword.DEFAULT_KEYWORD_CHARS


def test_nested_grammar_context_restores_defaults_after_exception():
    original = parser_defaults()
    with pytest.raises(ValueError, match="grammar failed"), generate_grammar():
        active = parser_defaults()
        assert active[0] == " \t\r"
        assert "-" in active[1]
        with generate_grammar():
            assert parser_defaults() == active
        assert parser_defaults() == active
        raise ValueError("grammar failed")
    assert parser_defaults() == original


def test_overlapping_grammar_contexts_restore_only_after_last_exit():
    original = parser_defaults()
    first_entered, second_entered, first_exited = Event(), Event(), Event()

    def first():
        with generate_grammar():
            first_entered.set()
            assert second_entered.wait(5)
        first_exited.set()

    def second():
        assert first_entered.wait(5)
        with generate_grammar():
            second_entered.set()
            assert first_exited.wait(5)
            # Construction still sees SPRKKR defaults after the first exit.
            assert "\n" not in pp.Word("x").whiteChars
            assert "-" in pp.Keyword("x").identChars

    with ThreadPoolExecutor(max_workers=2) as workers:
        futures = [workers.submit(first), workers.submit(second)]
        for future in futures:
            future.result(timeout=15)
    assert parser_defaults() == original

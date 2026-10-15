import pickle

from ..grammar_types.basic import Keyword


def test_keyword_pickles():
    pickle.dumps(Keyword("A", "B"))

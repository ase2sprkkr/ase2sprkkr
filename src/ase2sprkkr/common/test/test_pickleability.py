import pickle

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from ..grammar_types.basic import Keyword


class TestGrammarTypesPickle(TestCase):

    def test_keyword_pickles(self):
        pickle.dumps(Keyword('A', 'B'))

import inspect
import pickle

if __package__:
    from .init_tests import TestCase, patch_package
else:
    from init_tests import TestCase, patch_package
__package__, __name__ = patch_package(__package__, __name__)

if True:
    from ..definitions.jxc import input_parameters as jxc_input_parameters
    from ..definitions.scf import input_parameters as scf_input_parameters
    from ..definitions.sections import ENERGY, TAU
    from ..definitions.torque import input_parameters as torque_input_parameters
    from ..input_parameters import InputParameters
    from ...potentials.definitions.potential import potential_definition


DEFINITION_HOOK_ATTRIBUTES = (
    '_default_value',
    'is_optional',
    'is_required',
    'write_condition',
    'condition',
    'warning_condition',
)

TYPE_HOOK_ATTRIBUTES = (
    'condition',
    'after_convert',
    'free_header',
)


def _iter_definition_nodes(root):
    seen = set()
    stack = [root]

    while stack:
        node = stack.pop()
        if id(node) in seen:
            continue
        seen.add(id(node))
        yield node

        members = getattr(node, 'members', None)
        if callable(members):
            stack.extend(members())

        for attr in ('type', 'grammar_type'):
            child = getattr(node, attr, None)
            if child is not None and not inspect.isclass(child):
                stack.append(child)


def _iter_callback_hooks(root):
    for node in _iter_definition_nodes(root):
        node_name = getattr(node, 'name', node.__class__.__name__)
        for attr in DEFINITION_HOOK_ATTRIBUTES:
            hook = getattr(node, attr, None)
            if callable(hook):
                yield f'{node_name}.{attr}', hook

        node_type = getattr(node, 'type', None)
        if node_type is not None:
            for attr in TYPE_HOOK_ATTRIBUTES:
                hook = getattr(node_type, attr, None)
                if callable(hook):
                    yield f'{node_name}.type.{attr}', hook


def _assert_importable_hook(path, hook):
    assert inspect.isfunction(hook), \
        f'{path} should reference a module-level function, got {type(hook)!r}'
    assert hook.__name__ != '<lambda>', \
        f'{path} still references a lambda'
    assert '<locals>' not in hook.__qualname__, \
        f'{path} still references a local function: {hook.__qualname__}'
    pickle.dumps(hook)


class TestDefinitionCallbacks(TestCase):

    def test_callbacks_use_importable_functions(self):
        roots = [
            potential_definition,
            TAU,
            ENERGY(emin=(0.7, 'the energy to compute the BSF', None), emax='emin'),
            ENERGY(
                emin=(None, 'Minimum of the energy window in eV with respect to the Fermi level', -8.0),
                emax=(None, 'Maximum of the energy window in eV with respect to the Fermi level', 5.0),
            ),
            torque_input_parameters(),
            jxc_input_parameters(),
        ]

        hooks = []
        for root in roots:
            hooks.extend(_iter_callback_hooks(root))

        assert hooks, 'Expected to find callback hooks in the definition trees'
        for path, hook in hooks:
            _assert_importable_hook(path, hook)


class TestInputParametersPickle(TestCase):

    def test_scf_input_parameters_pickle(self):
        pickle.dumps(scf_input_parameters())

    def test_uninitialized_container_dunder_lookup_does_not_recurse(self):
        obj = object.__new__(InputParameters)

        try:
            getattr(obj, '__setstate__')
        except AttributeError:
            pass
        else:
            raise AssertionError('Expected missing __setstate__ on an uninitialized container')

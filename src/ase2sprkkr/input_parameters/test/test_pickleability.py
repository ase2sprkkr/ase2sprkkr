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
    from ..input_parameters_definitions import (
        InputParametersDefinition,
        InputSectionDefinition,
        InputValueDefinition,
    )
    from ...common.generated_configuration_definitions import Length
    from ...common.grammar_types import Array
    from ...common.warnings import DataValidityError
    from ...potentials.definitions.potential import potential_definition


def _copy_validator(_item, _values, _why):
    return None


DEFINITION_HOOK_ATTRIBUTES = (
    '_default_value',
    'is_optional',
    'is_required',
    'write_condition',
    'condition',
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

        for index, hook in enumerate(getattr(node, 'validators', ())):
            yield f'{node_name}.validators[{index}]', hook

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

    def test_root_validators_survive_copy(self):
        definition = InputParameters.create("SCF")._definition.copy(
            validators=(_copy_validator,)
        )
        copied = definition.copy()

        self.assertEqual(copied.validators, (_copy_validator,))

    def test_length_definition_pickle(self):
        definition = InputValueDefinition('NVALUES', Length('VALUES'))

        restored = pickle.loads(pickle.dumps(definition))

        self.assertEqual(restored.name, 'NVALUES')
        self.assertEqual(restored._length_of, ('VALUES',))
        self.assertEqual(len(restored._modifiers), 1)
        assert isinstance(restored._modifiers[0], Length)

    def test_length_behavior_survives_pickle(self):
        definition = InputParametersDefinition(
            [
                InputSectionDefinition(
                    'CONTROL',
                    [
                        InputValueDefinition('NVALUES', Length('VALUES')),
                        InputValueDefinition('VALUES', Array(int)),
                    ],
                )
            ]
        )

        restored = pickle.loads(pickle.dumps(definition)).create_object()
        restored.CONTROL.VALUES.set([1, 2])

        self.assertEqual(restored.CONTROL.NVALUES(), 2)
        with self.assertRaises(DataValidityError):
            restored.CONTROL.NVALUES.set(3)

    def test_length_modifier_survives_copy_and_pickle(self):
        definition = InputValueDefinition(
            'NVALUES', Length('VALUES', default_values=0)
        )

        copied = definition.copy()
        restored = pickle.loads(pickle.dumps(copied))

        assert restored._base_classes == (Length, InputValueDefinition)
        assert restored.__class__.__bases__ == restored._base_classes
        assert len(restored._modifiers) == 1
        assert isinstance(restored._modifiers[0], Length)
        assert restored._length_has_default_values
        self.assertEqual(restored._length_defaults, (0,))

    def test_length_copy_with_validators_keeps_length_validation(self):
        length = InputValueDefinition(
            'NVALUES', Length('FIRST', 'SECOND')
        ).copy(validators=(_copy_validator,))
        definition = InputParametersDefinition(
            [
                InputSectionDefinition(
                    'CONTROL',
                    [
                        InputValueDefinition('FIRST', Array(int)),
                        InputValueDefinition('SECOND', Array(int)),
                        length,
                    ],
                )
            ]
        )
        parameters = definition.create_object()
        parameters.CONTROL.set(FIRST=[1], SECOND=[1])

        with self.assertRaises(DataValidityError):
            parameters.CONTROL.SECOND.set([1, 2])

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

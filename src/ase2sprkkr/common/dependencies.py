"""Values derived from other configuration definitions."""

import copy
import inspect


class DependentValue:
    """Resolve and evaluate a value derived from other configuration items.

    A plain string denotes one exact definition path and can have an optional
    transform. For a callable, its parameter names denote exact local paths.
    """

    _MISSING = object()

    def __init__(self, value, transform=None):
        if isinstance(value, str):
            if transform is not None and not callable(transform):
                raise TypeError("A dependent-value transform has to be callable")
            self.paths = (value,)
            self.function = transform
        elif callable(value) and transform is None:
            self.function = value
            self.paths = tuple(inspect.signature(value).parameters)
        else:
            raise TypeError(
                "A dependent value has to be a path or a callable"
            )

        if not self.paths:
            raise ValueError("A dependent value needs at least one source")

        self.transforms = ()
        self.hooks = []

    @classmethod
    def create(cls, value):
        """Return an independent dependency description for ``value``."""

        if isinstance(value, cls):
            return value.copy()
        return cls(value)

    def copy(self):
        """Copy the description without copying installed grammar hooks."""

        out = copy.copy(self)
        out.hooks = []
        return out

    def mapped(self, transform):
        """Return a dependency applying ``transform`` to this result."""

        out = self.copy()
        out.transforms = (*self.transforms, transform)
        return out

    def __call__(self, *values):
        value = self.function(*values) if self.function else values[0]
        for transform in self.transforms:
            value = transform(value)
        return value

    @staticmethod
    def _parsed_value(tokens):
        """Extract the value from a configuration-definition parse result."""

        return tokens[0][1]

    def unbind(self):
        """Remove grammar hooks previously installed by :meth:`bind`."""

        for source, hook in self.hooks:
            source.remove_grammar_hook(hook)
        self.hooks = []

    def bind(self, container, callback):
        """Call ``callback`` with the derived value while parsing sources."""

        self.unbind()
        if container is None:
            return False

        try:
            definitions = tuple(
                container.get_member(path) for path in self.paths
            )
        except KeyError:
            return False

        values = [self._MISSING] * len(definitions)

        for index, source in enumerate(definitions):

            def parse_action(string, location, tokens, index=index):
                values[index] = self._parsed_value(tokens)
                if index == len(values) - 1:
                    if any(value is self._MISSING for value in values):
                        missing = ", ".join(
                            path
                            for path, value in zip(self.paths, values)
                            if value is self._MISSING
                        )
                        raise KeyError(
                            f"Dependent values have not been parsed: {missing}"
                        )
                    callback(self(*values))
                return tokens

            def grammar_hook(grammar, parse_action=parse_action):
                grammar.add_parse_action(parse_action)

            source.add_grammar_hook(grammar_hook)
            self.hooks.append((source, grammar_hook))
        return True

    def runtime_value(self, item):
        """Evaluate the dependency for a runtime option or container."""

        return self(
            *(
                item._container.get_member(path, unknown="fail")()
                for path in self.paths
            )
        )

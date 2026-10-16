"""Declarative metadata for SPR-KKR executables."""

from collections.abc import Mapping
from dataclasses import dataclass

from ..common.module_registry import RegistryDefinitionError, RegistryEntryNotFound


def _nonempty_string(value, field):
    if not isinstance(value, str):
        raise TypeError(f"{field} must be a string")
    value = value.strip()
    if not value:
        raise ValueError(f"{field} must not be empty")
    return value


@dataclass(frozen=True)
class ExecutableDefinition:
    """Metadata shared by tasks executed by the same SPR-KKR binary."""

    name: str
    default_output_reader: str = "DEFAULT"
    description: str = ""

    def __post_init__(self):
        name = _nonempty_string(self.name, "name")
        reader = _nonempty_string(
            self.default_output_reader, "default_output_reader"
        ).upper()
        if not isinstance(self.description, str):
            raise TypeError("description must be a string")

        object.__setattr__(self, "name", name)
        object.__setattr__(self, "default_output_reader", reader)


class ExecutableRegistry(Mapping):
    """Registry of known SPR-KKR executable metadata.

    Task definitions are deliberately not validated against this registry.
    A custom task may name any executable; the registry is consulted only when
    output-reader selection needs that executable's default.
    """

    def __init__(self, definitions=()):
        self._definitions = {}
        for definition in definitions:
            self.register(definition)

    def register(self, definition):
        if not isinstance(definition, ExecutableDefinition):
            raise TypeError("definition must be an ExecutableDefinition")
        if definition.name in self._definitions:
            raise RegistryDefinitionError(
                f"Duplicate executable definition {definition.name!r}"
            )
        self._definitions[definition.name] = definition
        return definition

    def __getitem__(self, name):
        try:
            return self._definitions[name]
        except KeyError:
            raise RegistryEntryNotFound(name) from None

    def __iter__(self):
        return iter(self._definitions)

    def __len__(self):
        return len(self._definitions)


executable_registry = ExecutableRegistry(
    (
        ExecutableDefinition(
            "kkrscf",
            description="SPR-KKR potential and self-consistency executable.",
        ),
        ExecutableDefinition(
            "kkrgen",
            description="General-purpose SPR-KKR property executable.",
        ),
        ExecutableDefinition(
            "kkrspec",
            default_output_reader="SPEC",
            description="SPR-KKR spectroscopy executable.",
        ),
        ExecutableDefinition(
            "kkrchi",
            description="SPR-KKR response-function executable.",
        ),
    )
)

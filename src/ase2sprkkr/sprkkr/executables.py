"""Declarative metadata for SPR-KKR executables."""

from collections.abc import Mapping
from dataclasses import dataclass

from ..common.module_registry import RegistryDefinitionError, RegistryEntryNotFound


__all__ = ("ExecutableDefinition", "executables")


def _nonblank_string(value, field):
    if not isinstance(value, str):
        raise TypeError(f"{field} must be a string")
    if not value.strip():
        raise ValueError(f"{field} must not be empty")
    return value


def _normalize_reader_name(value):
    return _nonblank_string(value, "default_output_reader").strip().upper()


@dataclass(frozen=True)
class ExecutableDefinition:
    """Metadata shared by tasks executed by the same SPR-KKR binary."""

    name: str
    default_output_reader: str = "DEFAULT"
    description: str = ""

    def __post_init__(self):
        name = _nonblank_string(self.name, "name")
        reader = _normalize_reader_name(self.default_output_reader)
        if not isinstance(self.description, str):
            raise TypeError("description must be a string")

        object.__setattr__(self, "name", name)
        object.__setattr__(self, "default_output_reader", reader)


class _ExecutableRegistry(Mapping):
    """Read-only mapping of known SPR-KKR executable metadata."""

    def __init__(self, definitions=()):
        items = {}
        for definition in definitions:
            if not isinstance(definition, ExecutableDefinition):
                raise TypeError("definition must be an ExecutableDefinition")
            if definition.name in items:
                raise RegistryDefinitionError(
                    f"Duplicate executable definition {definition.name!r}"
                )
            items[definition.name] = definition
        self._definitions = items

    def __getitem__(self, name):
        try:
            return self._definitions[name]
        except KeyError:
            raise RegistryEntryNotFound(name) from None

    def __iter__(self):
        return iter(self._definitions)

    def __len__(self):
        return len(self._definitions)


executables = _ExecutableRegistry(
    (
        ExecutableDefinition(
            "kkrscf",
            description="Potential generation and self-consistent calculations.",
        ),
        ExecutableDefinition(
            "kkrgen",
            description="Electronic-structure and general property calculations.",
        ),
        ExecutableDefinition(
            "kkrspec",
            default_output_reader="SPEC",
            description="Spectroscopy and photoemission calculations.",
        ),
        ExecutableDefinition(
            "kkrchi",
            description="Response-function calculations.",
        ),
    )
)

"""YAML persistence helpers for the ASE2SPRKKR configuration.

This module deliberately only translates between YAML data and the existing
configuration model. Selecting the user configuration file and migrating the
legacy Python file are handled separately.
"""

from collections.abc import Mapping
from pathlib import Path

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap


def _yaml():
    yaml = YAML(typ="rt")
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)
    return yaml


def load_yaml_document(file):
    """Load a YAML mapping while retaining comments and formatting metadata."""
    path = Path(file)
    with path.open("r", encoding="utf8") as stream:
        data = _yaml().load(stream)
    if data is None:
        return CommentedMap()
    if not isinstance(data, Mapping):
        raise ValueError("ASE2SPRKKR configuration YAML must contain a mapping")
    return data


def dump_yaml_document(file, data):
    """Write a YAML document using ruamel's round-trip representation."""
    path = Path(file)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf8") as stream:
        _yaml().dump(data, stream)


def explicit_configuration(config):
    """Return only values explicitly stored in ASE2SPRKKR configuration.

    Generated compatibility values such as the virtual ``legacy`` profile are
    intentionally omitted. Virtual profiles are likewise not part of the
    persistent repeated section.
    """
    return config.as_dict(only_changed="explicit", generated=False) or {}


def apply_yaml_document(config, data):
    """Apply a YAML mapping to a configuration object."""
    if not isinstance(data, Mapping):
        raise ValueError("ASE2SPRKKR configuration YAML must contain a mapping")
    config.set(data, unknown="fail")


def load_configuration_yaml(config, file):
    """Load ``file`` and apply it to ``config``; return the round-trip document."""
    data = load_yaml_document(file)
    apply_yaml_document(config, data)
    return data


def save_configuration_yaml(config, file, document=None):
    """Persist explicit config values to YAML and preserve existing comments.

    When a round-trip ``document`` is supplied, mappings are updated in place so
    comments attached to still-existing keys survive the save.
    """
    values = explicit_configuration(config)
    if document is None:
        document = CommentedMap()
    _update_mapping(document, values)
    dump_yaml_document(file, document)
    return document


def _update_mapping(target, values):
    """Synchronize a CommentedMap recursively without replacing its nodes."""
    for key in tuple(target):
        if key not in values:
            del target[key]

    for key, value in values.items():
        current = target.get(key)
        if isinstance(value, Mapping):
            if not isinstance(current, Mapping):
                current = CommentedMap()
                target[key] = current
            _update_mapping(current, value)
        else:
            target[key] = value

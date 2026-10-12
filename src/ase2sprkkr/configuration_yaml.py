"""YAML persistence helpers for the ASE2SPRKKR configuration.

The YAML document is kept separate from the typed configuration. Loading first
tries the whole document atomically; if a user-data error is found, valid
fragments are recovered while invalid/unknown raw nodes remain available for
round-trip persistence.
"""

import copy
from collections.abc import Mapping
from pathlib import Path
import shlex
import warnings

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

from .common.warnings import DataValidityError


class YamlConfigurationWarning(UserWarning):
    """Problems found while recovering a syntactically valid YAML config."""


_EXPECTED_USER_ERRORS = (KeyError, TypeError, ValueError, DataValidityError)
_MISSING = object()


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


def _shell_string(value):
    """Serialize an argv-like value as one shell-escaped command string."""
    if isinstance(value, (list, tuple)):
        return shlex.join([str(item) for item in value])
    return value


def _storage_configuration(values):
    """Convert command argv values to the canonical YAML representation."""
    profiles = values.get("profiles") if isinstance(values, Mapping) else None
    if not isinstance(profiles, Mapping):
        return values

    for profile in profiles.values():
        if not isinstance(profile, Mapping):
            continue
        if "mpi" in profile:
            profile["mpi"] = _shell_string(profile["mpi"])
        if "ssh_options" in profile:
            profile["ssh_options"] = _shell_string(profile["ssh_options"])
    return values


def _typed_configuration(data):
    """Convert shell-string YAML values to the typed configuration shape."""
    out = copy.deepcopy(data)
    profiles = out.get("profiles") if isinstance(out, Mapping) else None
    if not isinstance(profiles, Mapping):
        return out

    for profile in profiles.values():
        if not isinstance(profile, Mapping):
            continue
        ssh_options = profile.get("ssh_options")
        if isinstance(ssh_options, str):
            profile["ssh_options"] = shlex.split(ssh_options)
    return out


def explicit_configuration(config):
    """Return explicitly stored values in their canonical YAML representation."""
    values = config.as_dict(only_changed="explicit", generated=False) or {}
    return _storage_configuration(values)


def _validate_yaml_shape(data):
    if not isinstance(data, Mapping):
        raise ValueError("ASE2SPRKKR configuration YAML must contain a mapping")


def _stage_persisted_configuration(config, transaction, data):
    """Stage persisted data, including the private legacy profile."""
    profiles = data.get("profiles", _MISSING)
    values = {name: value for name, value in data.items() if name != "profiles"}
    config.stage(transaction, values, unknown="fail")
    if profiles is not _MISSING:
        config.profiles._stage_persisted(
            transaction,
            profiles,
            unknown="fail",
        )


def _set_persisted_configuration(config, data):
    """Apply persisted data atomically through the privileged profile path."""
    with config._mutation("set") as (transaction, _policy):
        _stage_persisted_configuration(config, transaction, data)


def apply_yaml_document(config, data):
    """Apply a complete YAML mapping atomically to a configuration object."""
    _validate_yaml_shape(data)
    _set_persisted_configuration(config, _typed_configuration(data))


def _recover_profiles(config, value, accepted, issues):
    if not isinstance(value, Mapping):
        issues.append(("profiles", ValueError("profiles must be a mapping")))
        return

    accepted_profiles = CommentedMap()
    for name, profile_data in value.items():
        try:
            with config.profiles._mutation("set") as (transaction, _policy):
                config.profiles._stage_persisted(
                    transaction,
                    {name: profile_data},
                    unknown="fail",
                    merge=True,
                )
        except _EXPECTED_USER_ERRORS as error:
            issues.append((f"profile {name!r}", error))
        else:
            accepted_profiles[name] = profile_data
    accepted["profiles"] = accepted_profiles


def _recover_section(member, path, value, accepted, issues):
    """Recover a normal section, then fall back to individual options."""
    try:
        member.set(value, unknown="fail")
    except _EXPECTED_USER_ERRORS:
        if not isinstance(value, Mapping):
            issues.append((path, ValueError(f"{path} must be a mapping")))
            return
    else:
        accepted[path] = value
        return

    accepted_section = CommentedMap()
    for name, child_value in value.items():
        child_path = f"{path}.{name}"
        try:
            child = member.get_member(
                name, unknown="fail", is_option=None, lower_case=False
            )
            child.set(child_value, unknown="fail")
        except _EXPECTED_USER_ERRORS as error:
            issues.append((child_path, error))
        else:
            accepted_section[name] = child_value
    accepted[path] = accepted_section


def _valid_default_profile(config):
    """Return an error if the configured implicit profile cannot be resolved."""
    name = config.default_profile()
    if name is None or name in config.profiles:
        return None
    return ValueError(f"unknown execution profile {name!r}")


def _discard_invalid_default_profile(config, accepted, issues):
    error = _valid_default_profile(config)
    if error is None:
        return
    config.default_profile.set(None)
    accepted.pop("default_profile", None)
    issues.append(("default_profile", error))


def _recover_configuration(config, data):
    """Recover valid logical units after an atomic load has failed."""
    accepted = CommentedMap()
    issues = []

    if "profiles" in data:
        _recover_profiles(config, data["profiles"], accepted, issues)

    keys = [key for key in data if key not in ("profiles", "default_profile")]
    if "default_profile" in data:
        keys.append("default_profile")

    for name in keys:
        value = data[name]
        try:
            member = config.get_member(
                name, unknown="fail", is_option=None, lower_case=False
            )
        except KeyError as error:
            issues.append((str(name), error))
            continue

        if member._definition.is_option:
            try:
                member.set(value, unknown="fail")
            except _EXPECTED_USER_ERRORS as error:
                issues.append((str(name), error))
            else:
                accepted[name] = value
        else:
            _recover_section(member, str(name), value, accepted, issues)

    _discard_invalid_default_profile(config, accepted, issues)
    return accepted, issues


def _format_issues(file, issues):
    lines = [
        f"Problems found in ASE2SPRKKR configuration '{file}':",
        "",
    ]
    for path, error in issues:
        lines.append(f"- {path}: {error}")
    lines.extend(("", "Other valid configuration values were loaded."))
    return "\n".join(lines)


def _apply_with_recovery(config, data, *, file=None, warn=True):
    _validate_yaml_shape(data)
    typed_data = _typed_configuration(data)
    try:
        _set_persisted_configuration(config, typed_data)
    except _EXPECTED_USER_ERRORS:
        accepted, issues = _recover_configuration(config, typed_data)
    else:
        accepted = CommentedMap(data)
        issues = []
        _discard_invalid_default_profile(config, accepted, issues)

    if warn and issues:
        warnings.warn(
            _format_issues(file or "<YAML>", issues),
            YamlConfigurationWarning,
            stacklevel=3,
        )
    return accepted, issues


def load_configuration_yaml(config, file):
    """Load YAML strictly first, then recover independent valid fragments."""
    data = load_yaml_document(file)
    _apply_with_recovery(config, data, file=file)
    return data


def _classify_document(config, document):
    """Return the subset of ``document`` which the current schema accepts."""
    fresh = config._definition.create_object()
    accepted, _issues = _apply_with_recovery(fresh, document, warn=False)
    return accepted


def save_configuration_yaml(config, file, document=None):
    """Persist explicit values without destroying unknown or invalid YAML.

    Existing YAML is classified through the same loader. Valid nodes are
    synchronized with ``config``; unknown or invalid nodes are preserved as raw
    round-trip YAML unless the caller explicitly replaces the same logical unit.
    """
    values = explicit_configuration(config)
    if document is None:
        document = CommentedMap()
        accepted = CommentedMap()
    else:
        accepted = _classify_document(config, document)
    _update_mapping(document, values, accepted)
    dump_yaml_document(file, document)
    return document


def _update_mapping(target, values, accepted):
    """Synchronize known values while leaving unaccepted raw nodes untouched."""
    for key in tuple(accepted):
        if key not in values and key in target:
            del target[key]

    for key, value in values.items():
        current = target.get(key, _MISSING)
        known_current = accepted.get(key, _MISSING)

        if isinstance(value, Mapping):
            if (
                current is _MISSING
                or not isinstance(current, Mapping)
                or known_current is _MISSING
            ):
                current = CommentedMap()
                target[key] = current
            known_children = (
                known_current
                if isinstance(known_current, Mapping)
                else CommentedMap()
            )
            _update_mapping(current, value, known_children)
        else:
            target[key] = value

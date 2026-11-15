"""Migrate legacy ASE2SPRKKR Python configuration to YAML."""

from . import config as config_command


help = "Migrate legacy ASE2SPRKKR configuration to config.yaml."
description = (
    "Convert declarative legacy Python configuration to YAML, migrate "
    "execution settings into profiles, and retain unrelated Python code."
)


def parser(parser):
    config_command.add_migration_cleanup_arguments(parser)


def run(args, global_args):
    from ase2sprkkr.configuration import (
        legacy_user_preferences_file,
        user_preferences_file,
    )
    config_command._migrate(
        args,
        legacy_user_preferences_file(),
        user_preferences_file(),
    )

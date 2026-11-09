"""Delegate the user-facing remote CLI to KKRunner, only when requested."""

help = "Inspect, follow or cancel an execution handle on this host."
description = "User commands delegated to the KKRunner remote service."

# Discovery and top-level --help must not load the runner.
lazy_parser = True


def parser(argument_parser):
    from kkrunner_cli.main import user_parser
    return user_parser(argument_parser)


def run(args, global_args):
    from kkrunner_cli.main import run as runner_run
    return runner_run(args, global_args)

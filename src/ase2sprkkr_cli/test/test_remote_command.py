"""The end-user CLI delegates to the runner's single remote implementation."""
import argparse
import json

import pytest
from kkrunner_cli import main as cli
from kkrunner.remote import service as remote_service
from ase2sprkkr_cli.commands import remote


@pytest.mark.parametrize('action', ['status', 'cancel'])
def test_user_command_delegates_without_protocol_arguments(monkeypatch, capsys, action):
    calls = []
    monkeypatch.setattr(remote_service, action + '_remote',
                        lambda handle: calls.append(handle) or {'state': 'running'})
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    remote.run(parser.parse_args([action, 'opaque-handle']), None)
    assert calls == ['opaque-handle']
    assert json.loads(capsys.readouterr().out) == {'state': 'running'}


def test_follow_exposes_user_stream_options():
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    args = parser.parse_args(['follow', 'opaque', '--live'])
    assert args.live
    assert args.remote_id == 'opaque'


@pytest.mark.parametrize('action', ['submit', 'recover', 'ack', 'ack-submit', 'resolve-mpi'])
def test_internal_protocol_commands_belong_to_runner(action):
    parser = argparse.ArgumentParser()
    remote.parser(parser)
    with pytest.raises(SystemExit):
        parser.parse_args([action])


def test_parser_loads_without_calculation_or_execution_library():
    import subprocess
    import sys

    code = """\
import sys
import argparse
from ase2sprkkr_cli.commands import remote
p = argparse.ArgumentParser()
remote.parser(p)
p.parse_args(["status", "opaque"])
assert "ase2sprkkr" not in sys.modules
assert "kkrunner.remote.service" not in sys.modules
assert "kkrunner.runs.manager" not in sys.modules
"""
    subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)

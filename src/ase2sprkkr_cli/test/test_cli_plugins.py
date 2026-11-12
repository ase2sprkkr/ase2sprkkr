"""Third-party CLI command packages use the same interface as native commands."""

import importlib

import pytest

from ase2sprkkr_cli.main import run


def _install_plugin(monkeypatch, tmp_path, package, modules, *, target=None):
    root = tmp_path / package
    commands = root / "commands"
    commands.mkdir(parents=True)
    (root / "__init__.py").write_text("", encoding="utf8")
    (commands / "__init__.py").write_text("", encoding="utf8")
    for filename, source in modules.items():
        (commands / filename).write_text(source, encoding="utf8")
    dist = tmp_path / f"{package}-1.0.dist-info"
    dist.mkdir()
    (dist / "METADATA").write_text(
        f"Metadata-Version: 2.1\nName: {package}\nVersion: 1.0\n",
        encoding="utf8",
    )
    (dist / "entry_points.txt").write_text(
        f"[ase2sprkkr]\ncli = {target or package + '.commands'}\n",
        encoding="utf8",
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    importlib.invalidate_caches()


_COMMAND = """
help = "Example plugin command."
description = "Provided by a separate distribution."
def parser(parser):
    parser.add_argument("--message", default="hello")
def run(args, global_args):
    print(args.message)
"""


def test_installed_package_exposes_commands_and_aliases(monkeypatch, tmp_path, capsys):
    _install_plugin(monkeypatch, tmp_path, "ase2sprkkr_plugin_echo", {
        "hello_world.py": """
help = "Hello from a plugin."
description = "Plugin detail."
aliases = ["greet"]
def extra_help():
    return "PLUGIN EXTENDED HELP"
def parser(parser):
    parser.add_argument("--message", required=True)
def run(args, global_args):
    print(f"plugin: {args.message}; debug={global_args['debug']}")
""",
        "support.py": "VALUE = 42\n",  # Not a command, therefore not registered.
    })

    with pytest.raises(SystemExit) as exc:
        run(["--help"])
    assert exc.value.code == 0
    output = capsys.readouterr().out
    assert "hello-world" in output
    assert "greet" in output
    assert "support" not in output

    run(["greet", "--message", "from alias"])
    assert "plugin: from alias; debug=False" in capsys.readouterr().out

    with pytest.raises(SystemExit) as exc:
        run(["hello-world", "--help"])
    assert exc.value.code == 0
    assert "PLUGIN EXTENDED HELP" in capsys.readouterr().out


def test_plugin_parsers_are_only_initialized_on_selection(monkeypatch, tmp_path, capsys):
    _install_plugin(monkeypatch, tmp_path, "ase2sprkkr_plugin_lazy", {
        "lazy_command.py": """
help = "Lazy plugin."
description = "Do not configure unless selected."
def parser(parser):
    raise RuntimeError("plugin parser was selected")
def run(args, global_args):
    pass
""",
    })
    with pytest.raises(SystemExit) as exc:
        run(["--help"])
    assert exc.value.code == 0
    assert "lazy-command" in capsys.readouterr().out

    with pytest.raises(RuntimeError, match="plugin parser was selected"):
        run(["lazy-command"])


def test_native_names_and_aliases_cannot_be_overridden(monkeypatch, tmp_path, capsys):
    _install_plugin(monkeypatch, tmp_path, "ase2sprkkr_plugin_conflicts", {
        "run.py": _COMMAND,
        "custom.py": 'aliases = ["view"]\n' + _COMMAND,
        "safe.py": _COMMAND,
    })

    run(["safe", "--message", "available"])
    captured = capsys.readouterr()
    assert "available" in captured.out
    assert "run" in captured.err
    assert "view" in captured.err
    assert "conflicts with an existing command" in captured.err


def test_broken_plugin_modules_do_not_hide_working_commands(monkeypatch, tmp_path, capsys):
    _install_plugin(monkeypatch, tmp_path, "ase2sprkkr_plugin_broken_module", {
        "broken.py": 'raise ImportError("bad optional dependency")\n',
        "incomplete.py": 'help = "Bad"\ndescription = "No run"\ndef parser(p): pass\n',
        "working.py": _COMMAND,
    })
    # The entry point itself can be invalid independently of its modules.
    invalid = tmp_path / "ase2sprkkr_plugin_invalid_entry-1.0.dist-info"
    invalid.mkdir()
    (invalid / "METADATA").write_text(
        "Metadata-Version: 2.1\nName: invalid-entry\nVersion: 1.0\n", encoding="utf8"
    )
    (invalid / "entry_points.txt").write_text(
        "[ase2sprkkr]\ncli = math\n", encoding="utf8"
    )

    run(["working", "--message", "still works"])
    captured = capsys.readouterr()
    assert "still works" in captured.out
    assert "bad optional dependency" in captured.err
    assert "incomplete" in captured.err
    assert "must expose a Python package" in captured.err


def test_plugin_unknown_argument_passthrough(monkeypatch, tmp_path, capsys):
    _install_plugin(monkeypatch, tmp_path, "ase2sprkkr_plugin_pass", {
        "pass_through.py": """
help = "Passthrough."
description = "Accept additional flags."
unknowns = "remaining"
def parser(parser):
    parser.set_defaults(remaining=[])
def run(args, global_args):
    print(repr(args.remaining))
""",
    })
    run(["pass-through", "--custom-flag", "value"])
    assert "['--custom-flag', 'value']" in capsys.readouterr().out


def test_python39_entry_points_mapping_is_supported(monkeypatch, tmp_path, capsys):
    _install_plugin(monkeypatch, tmp_path, "ase2sprkkr_plugin_legacy", {
        "legacy.py": _COMMAND,
    })

    class Entry:
        name = "cli"
        value = "ase2sprkkr_plugin_legacy.commands"

        def load(self):
            return importlib.import_module(self.value)

    class Unrelated(Entry):
        name = "other"

    import importlib.metadata
    monkeypatch.setattr(
        importlib.metadata, "entry_points",
        lambda: {"ase2sprkkr": [Entry(), Unrelated()]},
    )

    run(["legacy", "--message", "mapping works"])
    assert "mapping works" in capsys.readouterr().out

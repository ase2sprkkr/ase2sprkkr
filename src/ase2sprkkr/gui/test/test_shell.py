"""Tests for optional interactive-shell dependencies."""

import builtins
import pytest

from .shell import JUPYTER_INSTALL_COMMAND, JupyterDependencyError, load_jupyter_dependencies


def test_missing_jupyter_dependency_has_install_command(monkeypatch):
    real_import = builtins.__import__

    def without_jupyter(name, *args, **kwargs):
        if name in {"nbformat", "nbconvert.preprocessors", "jupyterlab"}:
            raise ModuleNotFoundError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_jupyter)

    with pytest.raises(JupyterDependencyError) as exc_info:
        load_jupyter_dependencies()

    message = str(exc_info.value)
    assert "Jupyter support is not installed" in message
    assert JUPYTER_INSTALL_COMMAND in message

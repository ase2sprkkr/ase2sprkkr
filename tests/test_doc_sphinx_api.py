"""Regression checks for generated Sphinx API pages.

These tests inspect the generation templates and the skip callback without
importing SPR-KKR or Matplotlib. A full Sphinx build remains the integration
check for actual cross-reference resolution.
"""

import ast
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "doc_src"


def _skip_member():
    tree = ast.parse((DOCS / "conf.py").read_text(encoding="utf-8"))
    node = next(
        node for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "skip_member"
    )
    namespace = {"InputParameters": object(), "Potential": object()}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(DOCS / "conf.py"), "exec"), namespace)
    return namespace["skip_member"]


def test_skip_non_callable_cached_descriptor():
    from types import SimpleNamespace

    skip = _skip_member()
    descriptor = SimpleNamespace(__module__="ase2sprkkr.common.decorators")
    assert skip(None, "method", "_validate_grammar", descriptor, False, {}) is True
    assert skip(None, "attribute", "BaseString._validate_grammar", descriptor, False, {}) is True
    assert skip(None, "attribute", "CubicSpline", descriptor, False, {}) is True


def test_skip_external_inherited_library_methods():
    from types import SimpleNamespace

    skip = _skip_member()
    matplotlib_set = SimpleNamespace(__module__="matplotlib.artist")
    scipy_method = SimpleNamespace(__module__="scipy.interpolate")
    own_method = SimpleNamespace(__module__="ase2sprkkr.gui.k_path")
    assert skip(None, "method", "set", matplotlib_set, False, {}) is True
    assert skip(None, "method", "spline", scipy_method, False, {}) is True
    assert skip(None, "method", "draw", own_method, False, {}) is False


@pytest.mark.parametrize("class_name", ("Arrow3D", "BaseString"))
def test_autosummary_class_template_excludes_problem_members(class_name):
    jinja2 = pytest.importorskip("jinja2")

    template = jinja2.Environment(
        loader=jinja2.FileSystemLoader(DOCS / "_templates"),
    )
    template.filters["underline"] = lambda name: str(name) + "\n" + "=" * len(name)
    doc = template.get_template("custom-class-template.rst").render(
        name=class_name,
        objname=class_name,
        module=(
            "ase2sprkkr.gui.k_path"
            if class_name == "Arrow3D"
            else "ase2sprkkr.common.grammar_types.basic"
        ),
        fullname=class_name,
        _=lambda text: text,
        own_methods=["draw", "set", "_validate_grammar"],
        inherited_methods=["set", "_validate_grammar"],
        own_attributes=[],
        inherited_attributes=["axes"],
    )
    assert ":exclude-members: _validate_grammar, CubicSpline" in doc
    assert "~" + class_name + "._validate_grammar" not in doc
    if class_name == "Arrow3D":
        assert ":exclude-members: _validate_grammar, CubicSpline, set" in doc
        assert ":inherited-members:" not in doc
        assert "Inherited methods" not in doc
        assert "~Arrow3D.set" not in doc
    else:
        # Inherited members can be linked, not registered as new objects.
        assert "Inherited methods" in doc
        assert "~BaseString.draw" in doc
    assert not any(
        line.strip().startswith(":inherited-members:")
        for line in doc.splitlines()
    )


def test_inherited_members_not_registered_again():
    """Document inherited APIs on defining classes, without duplicates."""
    conf = (DOCS / "conf.py").read_text(encoding="utf-8")
    parsed = ast.parse(conf)
    assignment = next(
        node for node in parsed.body
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name) and target.id == "autodoc_default_options"
            for target in node.targets
        )
    )
    options = ast.literal_eval(assignment.value)
    assert "inherited-members" not in options

    template = (DOCS / "_templates" / "custom-class-template.rst").read_text(
        encoding="utf-8"
    )
    # The class hierarchy and optional *links* to inherited members are
    # useful, but an :inherited-members: directive duplicates targets
    # throughout the generated API.
    assert ":inherited-members:" not in template
    assert "inheritance-diagram" in template


def test_type_annotations_are_not_ambiguous_short_names():
    sources = {
        "repetition": ROOT / "src/ase2sprkkr/common/repetition.py",
        "transaction": ROOT / "src/ase2sprkkr/common/configuration_transaction.py",
        "container": ROOT / "src/ase2sprkkr/common/configuration_containers.py",
    }
    text = {key: file.read_text(encoding="utf-8") for key, file in sources.items()}
    assert '-> "Repeated.Numbering"' not in text["repetition"]
    assert 'configuration: "Configuration"' not in text["transaction"]
    assert "Configuration as BaseConfiguration" in text["transaction"]
    assert "Configuration as BaseConfiguration" in text["container"]

def test_occupation_annotations_resolve_physical_atomic_type():
    """Do not emit the ambiguous bare AtomicType class in API type hints."""
    source = (
        ROOT / "src/ase2sprkkr/sprkkr/occupations.py"
    ).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert "from . import atomic_types" in source

    annotations = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        parameters = (
            node.args.posonlyargs + node.args.args
            + node.args.kwonlyargs
        )
        annotations.extend(arg.annotation for arg in parameters if arg.annotation)
        if node.returns:
            annotations.append(node.returns)
    assert annotations
    assert any(
        isinstance(member, ast.Attribute)
        and member.attr == "AtomicType"
        and isinstance(member.value, ast.Name)
        and member.value.id == "atomic_types"
        for annotation in annotations for member in ast.walk(annotation)
    )
    assert not any(
        isinstance(member, ast.Name) and member.id == "AtomicType"
        for annotation in annotations for member in ast.walk(annotation)
    )
    assert "except_from: ase2sprkkr.sprkkr.atomic_types.AtomicType" in source


def test_raw_data_constructor_preserves_variadic_keywords():
    """Docstring edits must not turn **kwargs into a required argument."""
    source = ROOT / "src/ase2sprkkr/common/grammar_types/data.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))
    raw_data = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "RawData")
    constructor = next(node for node in raw_data.body if isinstance(node, ast.FunctionDef) and node.name == "__init__")
    assert constructor.args.kwarg is not None
    assert constructor.args.kwarg.arg == "kwargs"

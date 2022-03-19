# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../src/'))


# -- Project information -----------------------------------------------------

project = 'ASE2SPRKKR'
copyright = '2021, Matyáš Novák & Jano Minár'
author = 'Matyáš Novák & Jano Minár'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx_toolbox.sidebar_links',
    'sphinx.ext.viewcode',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.inheritance_diagram'
]

autodoc_default_options = {
    'member-order': 'bysource',
    'private-members': True,
    'undoc-members': True,
    'inherited-members' : False
}

napoleon_attr_annotations = True
autodoc_member_order = 'bysource'
autosummary_generate = True

#autoclass_content = "both"  # Add __init__ doc (ie. params) to class summaries

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
autodoc_typehints = 'description'
github_url =  'https://github.com/ase2sprkkr/ase2sprkkr/'

html_css_files = [
    'mods.css',
]

def skip_member(app, what, name, obj, skip, options):
    if skip:
       return True
    if what=='module' or what=='package':
       if name == 'test':
          return True
       if name.endswith('.test'):
          return True
    if what=='method':
       if getattr(obj, '__objclass__', None) is object:
           return True
       if getattr(obj, '__module__', None) is None:
           return True
    if name == '__subclasshook__':
       breakpoint()
    if what=='attribute' and name in ['__module__','__weakref__', '__doc__', '__dict__']:
       return True
    return False

def setup(app):
    app.connect('autodoc-skip-member', skip_member)


exclude_patterns = ['*/test/*', '*/test']

inheritance_graph_attrs = dict(size='"10.0, 10.0"',
    fontsize=14, ranksep=0.3 )
inheritance_edge_attrs = { 'arrowsize': 0.7 }
inheritance_parent_node_attrs = { 'color' : 'gray28', 'fontcolor' : 'gray28', 'style' : 'solid' }
inheritance_node_attrs = { 'style' : '"filled"', 'fillcolor' : 'lightgray' }

def inheritance_style_callback(class_names, current_module, styles):
    styles['graph']['rankdir'] = 'TB' if len(class_names) <= 2 else 'LR'

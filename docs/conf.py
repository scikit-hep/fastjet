# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# Warning: do not change the path here. To use autodoc, you need to install the
# package first.

import os
import sys
from typing import List

# -- Project information -----------------------------------------------------

project = "fastjet"
copyright = "2021, The Scikit-HEP admins"
author = "Aryan Roy"


curdir = os.path.dirname(__file__)
finalpath = os.path.join(curdir, "..", "src")
sys.path.insert(0, finalpath)

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
]

autodoc_mock_imports = ["fastjet._ext", "fastjet._pyjet", "fastjet._swig"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "**.ipynb_checkpoints", "Thumbs.db", ".DS_Store", ".env"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

html_show_sourcelink = False
html_logo = "../docs/logo.svg"
html_theme_options = {"logo_only": True, "style_nav_header_background": "#fcfcfc"}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path: list[str] = []

master_doc = "index"

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

import jupytext

sys.path.insert(0, os.path.abspath("../.."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "xenopict"
copyright = "2024, S. Joshua Swamidass"
author = "S. Joshua Swamidass"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.githubpages",
    "sphinx_autodoc_typehints",
    "nbsphinx",
]

apidoc_separate_modules = True

nbsphinx_custom_formats = {
    ".py": lambda s: jupytext.reads(s, fmt="py:percent"),
}

# source_suffix = [".rst", ".md"]
source_suffix = [".rst", ".md"]

templates_path = ["_templates"]
exclude_patterns = ["build"]

nbsphinx_execute = "always"
nbsphinx_allow_errors = True

autosummary_generate = True

html_sourcelink_suffix = ""

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

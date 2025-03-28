# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


import sys
import os


# source_suffix = [".rst", ".md"]
source_suffix = [".rst", ".md"]

sys.path.insert(0, os.path.abspath("../.."))


project = "XenoPict"
copyright = "2022, S. Joshua Swamidass"
author = "S. Joshua Swamidass"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx_rtd_theme",
    "sphinxcontrib.apidoc",
    "nbsphinx",
]

apidoc_module_dir = "../../xenopict"
apidoc_output_dir = "api"
apidoc_toc_file = False
apidoc_excluded_paths = ["tests", "monkey", "magic"]
apidoc_separate_modules = True

import jupytext

nbsphinx_custom_formats = {
    ".md": lambda s: jupytext.reads(s, ".md"),
}


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

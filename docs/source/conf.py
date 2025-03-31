# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

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
    "myst_nb",
]

apidoc_separate_modules = True

# MyST configuration
myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
]

# MyST-NB settings
nb_execution_mode = "force"
nb_execution_timeout = 300
nb_execution_allow_errors = False
nb_execution_raise_on_error = True

templates_path = ["_templates"]
exclude_patterns = ["build", "conf.py"]

# Configure image handling
nb_output_stderr = "remove"

# Add nbsphinx kernel name for Python notebooks
nbsphinx_kernel_name = "python3"

# Configure image handling
nbsphinx_execute_arguments = [
    "--InlineBackend.figure_formats={'svg', 'pdf'}",
    "--InlineBackend.rc={'figure.dpi': 96}",
]

autosummary_generate = True

html_sourcelink_suffix = ""

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = []

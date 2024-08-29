"""Sphinx configuration."""

from datetime import datetime

project = "python-Flint"
author = "Fredrik Johansson"
copyright = f"{datetime.now().year}, {author}"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
]
#autodoc_typehints = "description"
html_theme = "sphinx_rtd_theme"

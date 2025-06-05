# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'eikonalfm'
copyright = '2025, Kevin Ganster'
author = 'Kevin Ganster'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # "myst_nb",              # markdown/notebook support
    "myst_parser",          # markdown support (requires Sphinx >= 2.1)
    "sphinx_rtd_theme",     # 'read the docs' theme
    "autoapi.extension",    # parse source code and docstrings to create API reference sheet
    "sphinx.ext.napoleon",  # parse numpydoc style docstrings
    "sphinx.ext.viewcode",  # adds a helpful link to the source code of each object in the API reference sheet
    # "sphinx.ext.mathjax"    # use JS to render math
]

autoapi_dirs = ["../src"]

myst_enable_extensions = ["dollarmath"]


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

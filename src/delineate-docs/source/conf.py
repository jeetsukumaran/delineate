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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'DELINEATE'
copyright = '2020, Jeet Sukumaran and Mark T. Holder'
author = 'Jeet Sukumaran and Mark T. Holder'

# The full version, including alpha/beta/rc tags
release = '1.2.3'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.githubpages",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

def setup(app):
  app.add_stylesheet("css/docs.css")

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

rst_epilog = """
.. |Anaconda| replace:: Anaconda
.. _Anaconda: https://www.anaconda.com/
.. |BPP| replace:: BP&P
.. _BPP: https://github.com/bpp/bpp
.. |delineate| replace:: DELINEATE
.. |pip| replace:: pip
.. _pip: https://docs.python.org/3/installing/index.html
.. |Python| replace:: Python
.. _Python: https://www.python.org/
.. |StarBeast2| replace:: StarBeast2
.. _StarBeast2: https://taming-the-beast.org/tutorials/starbeast2-tutorial/
"""

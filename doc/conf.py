#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# RainbowLasso documentation build configuration file.
#
# This file is execfile()d with the current directory set to its
# containing dir.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- General configuration ---------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    'sphinx.ext.autosectionlabel',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['templates']

# The suffix(es) of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'RainbowLasso'
copyright = u"2024-2026, Johannes Buchner"
author = u"Johannes Buchner"

# The version info for the project you're documenting.
version = u''
release = u''

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'build', 'Thumbs.db', '.DS_Store',
    '**.ipynb_checkpoints']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

autosectionlabel_prefix_document = True

# -- Options for HTML output -------------------------------------------

# The theme to use for HTML and HTML Help pages.
html_theme = "sphinx_book_theme"

html_baseurl = 'https://johannesbuchner.github.io/RainbowLasso/'
html_logo = "../logo.png"
html_show_sourcelink = False
html_favicon = None

html_show_sphinx = False

# Add any paths that contain custom static files here.
html_static_path = []

# -- Options for HTMLHelp output ---------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'RainbowLassodoc'


# -- Options for LaTeX output ------------------------------------------

latex_elements = {}

# Grouping the document tree into LaTeX files.
latex_documents = [
    (master_doc, 'RainbowLasso.tex',
     u'RainbowLasso Documentation',
     u'Johannes Buchner', 'manual'),
]


# -- Options for manual page output ------------------------------------

man_pages = [
    (master_doc, 'rainbowlasso',
     u'RainbowLasso Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------

texinfo_documents = [
    (master_doc, 'RainbowLasso',
     u'RainbowLasso Documentation',
     author,
     'RainbowLasso',
     'Compiles reliable ultraviolet to infrared photometry for '
     'modelling stars, galaxies and AGN.',
     'Miscellaneous'),
]

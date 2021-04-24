"""
Mpmath documentation build configuration file.

This file is execfile()d with the current directory set to its
containing dir.

The contents of this file are pickled, so don't put values in the
namespace that aren't pickleable (module imports are okay, they're
removed automatically).
"""

import mpmath


# Add any Sphinx extension module names here, as strings.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.mathjax',
              'sphinx.ext.intersphinx']

# The master toctree document.
master_doc = 'index'

# Project information.
project = 'mpmath'
copyright = '2007-2021, Fredrik Johansson and mpmath developers'
version = mpmath.__version__
release = version

# Define how the current time is formatted using time.strftime().
today_fmt = '%B %d, %Y'

# The "theme" that the HTML output should use.
html_theme = 'classic'

# Output file base name for HTML help builder.
htmlhelp_basename = 'mpmathdoc'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = [(master_doc, 'mpmath.tex', 'mpmath documentation',
                    r'Fredrik Johansson \and mpmath contributors', 'manual')]

# Additional stuff for the LaTeX preamble.
latex_preamble = r'\usepackage{amsfonts}'

# The name of default reST role, that is, for text marked up `like this`.
default_role = 'math'

# Contains mapping the locations and names of other projects that
# should be linked to in this documentation.
intersphinx_mapping = {
    'python3': ('https://docs.python.org/3/', None),
}

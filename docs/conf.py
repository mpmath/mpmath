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
              'sphinx.ext.intersphinx', 'sphinxcontrib.autoprogram',
              'matplotlib.sphinxext.plot_directive']

# Sphinx will warn about all references where the target cannot be found.
nitpicky = True

# Project information.
project = mpmath.__name__
copyright = '2007-2025, Fredrik Johansson and mpmath developers'
release = version = mpmath.__version__

# Define how the current time is formatted using time.strftime().
today_fmt = '%B %d, %Y'

# The "theme" that the HTML output should use.
html_theme = 'classic'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = [('index', 'mpmath.tex', 'mpmath documentation',
                    r'Fredrik Johansson \and mpmath contributors', 'manual')]

# The name of default reST role, that is, for text marked up `like this`.
default_role = 'math'

# Contains mapping the locations and names of other projects that
# should be linked to in this documentation.
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sympy': ('https://docs.sympy.org/latest/', None),
}

plot_include_source = True
plot_formats = [('png', 96), 'pdf']
plot_html_show_formats = False
plot_html_show_source_link = False

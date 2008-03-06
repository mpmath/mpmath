from docutils.core import publish_file

publish_file(
  source_path='manual.rst',
  destination_path='manual.html',
  writer_name='html')

"""
publish_file(
  source_path='manual.rst',
  destination_path='manual.tex',
  writer_name='latex')

import os
s = os.system('pdflatex manual.tex')
if not s:
    for e in ['aux','log','out','tex']:
        os.system('rm manual.%s' % e)
"""
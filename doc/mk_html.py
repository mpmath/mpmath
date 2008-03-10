from docutils.core import publish_file

publish_file(
  source_path='manual.rst',
  destination_path='manual.html',
  writer_name='html')

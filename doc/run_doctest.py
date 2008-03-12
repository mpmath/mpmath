import os
import sys

root_dir = os.path.abspath(os.path.dirname(__file__))
path = lambda *paths: os.path.abspath(os.path.join(*((root_dir,)+paths)))

sys.path.insert(0, path('..'))

import doctest
doctest.testfile(path("manual.rst"), module_relative=False)

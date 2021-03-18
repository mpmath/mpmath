#!/usr/bin/env python

import os
import setuptools

# where version is stored relative to setup.py
_version_path = os.path.join(
    os.path.dirname(__file__),
    'mpmath/__init__.py')
# get version without importing anything
# note that this is relying on __version__
# being the *first* line of the __init__ file
with open(_version_path, 'r') as f:
    # use eval to get a clean string of version from file
    __version__ = eval(f.readline().strip().split('=')[-1])

setuptools.setup(version=__version__)

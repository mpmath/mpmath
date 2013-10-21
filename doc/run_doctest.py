#!/usr/bin/env python

import os
import os.path
import sys

sys.path.insert(0, '..')

path = "source"

import doctest
for f in os.listdir(path):
    if f.endswith(".txt"):
        print f
        doctest.testfile(os.path.join(path, f), module_relative=False)

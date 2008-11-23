#!/usr/bin/env python

import os
if not os.path.exists("build"):
    os.mkdir("build")
os.system("sphinx-build -E source build")

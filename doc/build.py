#!/usr/bin/env python

import os
if not os.path.exists("build"):
    os.mkdir("build")
os.system("sphinx-build -E source build")

# for formulas (XXX)
imgstyle = """
/* inline formulas */
p img { vertical-align: middle; }

/* center display formulas */
div.math img {
    display: block;
    margin-left: auto;
    margin-right: auto;
}

"""
f = open("build/_static/default.css", "r+w")
f.seek(0)
if imgstyle not in f.read():
    f.seek(0, 2)
    f.write(imgstyle)

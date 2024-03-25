"""Airy function Bi(z) in the complex plane"""
from mpmath import *
cplot(airybi, [-8,8], [-8,8], points=50000)

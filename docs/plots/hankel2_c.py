"""
Hankel function `H2_n(z)` in the complex plane.
----------------------------------------------------
"""
from mpmath import *
cplot(lambda z: hankel2(1,z), [-8,8], [-8,8], points=50000)

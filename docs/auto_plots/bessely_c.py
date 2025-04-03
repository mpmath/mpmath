"""
Bessel function of 2nd kind `Y_n(z)` in the complex plane.
---------------------------------------------------------------
"""
from mpmath import *
cplot(lambda z: bessely(1,z), [-8,8], [-8,8], points=50000)

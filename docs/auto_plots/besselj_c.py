"""
Bessel function `J_n(z)` in the complex plane.
---------------------------------------------------
"""
from mpmath import *
cplot(lambda z: besselj(1,z), [-8,8], [-8,8], points=50000)

"""
Irregular Coulomb wave function in the complex plane.
------------------------------------------------------------
"""
from mpmath import *
cplot(lambda z: coulombg(1,1,z), points=50000)

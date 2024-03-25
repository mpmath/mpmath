"""Regular Coulomb wave function in the complex plane"""
from mpmath import *
cplot(lambda z: coulombf(1,1,z), points=50000)

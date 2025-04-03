"""
Principal branch of the Lambert W function `W(z)`.
--------------------------------------------------
"""
from mpmath import *
cplot(lambertw, [-1,1], [-1,1], points=50000)

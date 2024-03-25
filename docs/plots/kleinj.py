"""Klein J-function as function of the number-theoretic nome"""
from mpmath import *
fp.cplot(lambda q: fp.kleinj(qbar=q), [-1,1], [-1,1], points=50000)

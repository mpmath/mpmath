"""Scorer function Hi(x) and Hi'(x) on the real line"""
from mpmath import *
plot([scorerhi, diffun(scorerhi)], [-10,2], [0,2])

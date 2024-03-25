"""Scorer function Gi(x) and Gi'(x) on the real line"""
from mpmath import *
plot([scorergi, diffun(scorergi)], [-10,10])

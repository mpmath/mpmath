"""
Mathematical constants
"""
import math

from util import *
from floatop import *
from squareroot import *

def constant_memo(f):
    """Cache computer values of mathematical constants"""
    f.memo_prec = -1
    f.memo_val = None
    def g(prec):
        if prec == f.memo_prec:
            return f.memo_val
        if prec < f.memo_prec:
            return f.memo_val >> (f.memo_prec-prec)
        f.memo_val = f(prec)
        f.memo_prec = prec
        return f.memo_val
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g

def acot(n, prec, hyperbolic):
    """Compute acot of an integer using fixed-point arithmetic. With
    hyperbolic=True, compute acoth. The standard Taylor series
    is used."""
    s = t = (1 << prec) // n  # 1 / n
    k = 3
    while 1:
        # Repeatedly divide by k * n**2, and add
        t //= (n*n)
        term = t // k
        if not term:
            break
        # Alternate signs
        if hyperbolic or not k & 2:
            s += term
        else:
            s -= term
        k += 2
    return s

def machin(coefs, prec, hyperbolic=False):
    """
    Evaluate a Machin-like formula, i.e., a linear combination of
    acot(n) or acoth(n) for specific integer values of n, using fixed-
    point arithmetic.

    The input should be a list [(c, n), ...], giving c*acot[h](n) + ...
    """
    extraprec = 10
    s = 0
    for a, b in coefs:
        s += a * acot(b, prec+extraprec, hyperbolic)
    return (s >> extraprec)

@constant_memo
def pi_fixed(prec):
    """
    Compute floor(pi * 2**prec) as a big integer.

    For low precisions, Machin's formula pi = 16*acot(5)-4*acot(239)
    is used. For high precisions, the more efficient arithmetic-
    geometric mean iteration is used.
    """
    return machin([(16, 5), (-4, 239)], prec)

def fpi(prec, rounding):
    """Compute a floating-point approximation of pi"""
    return from_man_exp(pi_fixed(prec+5), -prec-5, prec, rounding)

@constant_memo
def log2_fixed(prec):
    return machin([(18, 26), (-2, 4801), (8, 8749)], prec, True)

def flog2(prec, rounding):
    return from_man_exp(log2_fixed(prec+5), -prec-5, prec, rounding)

@constant_memo
def log10_fixed(prec):
    return machin([(46, 31), (34, 49), (20, 161)], prec, True)

def flog10(prec, rounding):
    return from_man_exp(log10_fixed(prec+5), -prec-5, prec, rounding)


"""
Euler's constant (gamma) is computed using the Brent-McMillan formula,
gamma ~= A(n)/B(n) - log(n), where

  A(n) = sum_{k=0,1,2,...} (n**k / k!)**2 * H(k)
  B(n) = sum_{k=0,1,2,...} (n**k / k!)**2
  H(k) = 1 + 1/2 + 1/3 + ... + 1/k

The error is bounded by O(exp(-4n)). Choosing n to be a power
of two, 2**p, the logarithm becomes particularly easy to calculate.

Reference:
Xavier Gourdon & Pascal Sebah, The Euler constant: gamma
http://numbers.computation.free.fr/Constants/Gamma/gamma.pdf
"""

@constant_memo
def gamma_fixed(prec):
    prec += 30
    # choose p such that exp(-4*(2**p)) < 2**-n
    p = int(math.log((prec/4) * math.log(2), 2)) + 1
    n = 1<<p
    r = one = 1<<prec
    H, A, B, npow, k, d = 0, 0, 0, 1, 1, 1
    while r:
        A += (r * H) >> prec
        B += r
        r = r * (n*n) // (k*k)
        H += one // k
        k += 1
    S = ((A<<prec) // B) - p*log2_fixed(prec)
    return S >> 30

def fgamma(prec, rounding):
    return from_man_exp(gamma_fixed(prec+5), -prec-5, prec, rounding)

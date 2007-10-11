from util import *
from floatop import *


import math
import decimal

LOG2_10 = math.log(10,2)  # 3.3219...


# TODO: only binary_to_decimal and decimal_to_binary are used currently.
# Things could be sped up by using the other functions below, currently used
# only by the pidigits.py demo

getctx = decimal.getcontext
Dec = decimal.Decimal

def binary_to_decimal(s, n):
    """Represent as a decimal string with at most n digits"""
    man, exp, bc = s
    prec_ = getctx().prec
    getctx().prec = n + 10
    d = Dec(man) * Dec(2)**exp
    getctx().prec = n
    a = str(+d)
    getctx().prec = prec_
    return a

def decimal_to_binary(x, prec, rounding):
    dps = int(prec*LOG2_10) + 10
    prec_ = getctx().prec
    getctx().prec = dps
    d = Dec(x).normalize()
    sgn, digits, dexp = d.as_tuple()
    d = d * Dec(10)**(-dexp)
    power = fpow(ften, -dexp, prec+10, ROUND_FLOOR)
    y = fdiv(float_from_int(int(d), prec+10, ROUND_FLOOR), power, prec+10, ROUND_HALF_EVEN)
    getctx().prec = prec_
    y = normalize(y[0], y[1], prec, rounding)
    return y


def float_from_int(n, prec, rounding):
    return normalize(n, 0, prec, rounding)

def float_from_rational(p, q, prec, rounding):
    """Create floating-point number from a rational number p/q"""
    n = prec + bitcount(q) + 2
    return normalize((p<<n)//q, -n, prec, rounding)

def float_from_pyfloat(x, prec, rounding):
    # We assume that a float mantissa has 53 bits
    m, e = math.frexp(x)
    return normalize(int(m*(1<<53)), e-53, prec, rounding)

def float_to_int(s):
    man, exp, bc = s
    return rshift(man, -exp, ROUND_DOWN)

def float_to_pyfloat(s):
    """Convert to a Python float. May raise OverflowError."""
    man, exp, bc = s
    try:
        return math.ldexp(man, exp)
    except OverflowError:
        # Try resizing the mantissa. Overflow may still happen here.
        n = bc - 53
        m = man >> n
        return math.ldexp(m, exp + n)

def float_to_rational(s):
    """Return (p, q) such that s = p/q exactly. p and q are not reduced
    to lowest terms."""
    man, exp, bc = s
    if exp > 0:
        return man * 2**exp, 1
    else:
        return man, 2**-exp


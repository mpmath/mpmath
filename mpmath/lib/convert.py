from util import *
from floatop import *
from constants import flog2, flog10
import math

LOG2_10 = math.log(10,2)  # 3.3219...


def f2s(s, dps):
    """
    Helper function for representing the floating-point number s as a
    decimal with dps places. Returns (string, exponent) containing
    the decimal digits as a string and an integer for the decimal
    exponent.

    If inexact, the decimal representation is rounded to the floor.
    """
    RF = ROUND_FLOOR
    bitprec = int(dps * LOG2_10) + 10
    man, exp, bc = s

    # Cut down to size
    if abs(exp) > 500:
        # Set b = int(exp * log(2)/log(10))
        # If exp is huge, we must use high-precision arithmetic to
        # find the nearest power of ten
        expprec = bitcount(exp) + 5
        tmp = float_from_int(exp, expprec, RF)
        tmp = fmul(tmp, flog2(expprec, RF), expprec, RF)
        tmp = fdiv(tmp, flog10(expprec, RF), expprec, RF)
        b = float_to_int(tmp)
        s = fdiv(s, fpow(ften, b, bitprec, ROUND_FLOOR), bitprec, ROUND_FLOOR)
        man, exp, bc = s
        exponent = b
    else:
        exponent = 0

    # First, calculate mantissa digits by converting to a binary
    # fixed-point number and then converting that number to
    # a decimal fixed-point number.
    fixprec = max(bitprec - exp, 0)
    fixdps = int(fixprec / LOG2_10 + 0.5)
    sf = make_fixed(s, fixprec)
    sd = bin_to_radix(sf, fixprec, 10, fixdps)
    digits = numeral(sd, base=10, size=dps)

    # Next, figure out exponent and placement of decimal point
    exponent += len(digits) - fixdps - 1
    return digits, exponent


def binary_to_decimal(s, dps):
    """
    Represent as a decimal floating-point string that can be parsed by
    float.__init__ or Decimal.__init__ (or for that matter, by a human).

    The dps option specifies how many digits to include.
    """
    digits, exponent = f2s(s, dps)

    # Prettify numbers close to unit magnitude
    if -5 < exponent < 17:
        if exponent < 0:
            digits = ("0"*(-exponent)) + digits
            split = 1
        else:
            split = exponent + 1
        exponent = 0
    else:
        split = 1
    digits = (digits[:split] + "." + digits[split:dps])

    # Clean up trailing zeros
    digits = digits.rstrip('0')
    if digits[-1] == ".":
        digits += "0"

    if exponent == 0:
        return digits
    else:
        return digits + "e" + str(exponent)


def decimal_to_binary(x, prec, rounding):
    """
    Create a floating-point number from a decimal float literal.
    """
    # Verify that the input is a valid float literal
    float(x)
    # Split into mantissa, exponent
    x = x.lower()
    parts = x.split('e')
    if len(parts) == 1:
        exp = 0
    else: # == 2
        x = parts[0]
        exp = int(parts[1])
    # Look for radix point in mantissa
    parts = x.split('.')
    if len(parts) == 2:
        a, b = parts[0], parts[1].rstrip('0')
        exp -= len(b)
        x = a + b
    x = int(x)
    s = float_from_int(x, prec+10, ROUND_FLOOR)
    s = fmul(s, fpow(ften, exp, prec+10, ROUND_FLOOR), prec, rounding)
    return s

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


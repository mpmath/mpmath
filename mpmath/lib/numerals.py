from floatop import *
from constants import *

#----------------------------------------------------------------------
# Strings
#


def make_fixed(s, prec):
    """Convert a floating-point number to a fixed-point big integer"""
    man, exp, bc = s
    offset = exp + prec
    if offset >= 0:
        return man << offset
    else:
        return man >> (-offset)

# TODO: speed up for bases 2, 4, 8, 16, ...

def bin_to_radix(x, xbits, base, bdigits):
    """
    Radix conversion for fixed-point numbers. That is, convert
    x * 2**xbits to floor(x * 10**bdigits).
    """
    return x * (base**bdigits) >> xbits

stddigits = '0123456789abcdefghijklmnopqrstuvwxyz'

def small_numeral(n, base=10, digits=stddigits):
    """
    Return the string numeral of a positive integer in an arbitrary
    base. Most efficient for small input.
    """
    if base == 10:
        return str(n)
    digs = []
    while n:
        n, digit = divmod(n, base)
        digs.append(digits[digit])
    return "".join(digs[::-1])

def numeral(n, base=10, size=0, digits=stddigits):
    """
    Represent the integer n as a string of digits in the given base.
    Recursive division is used to make this function about 3x faster
    than Python's str() for converting integers to decimal strings.

    The 'size' parameters specifies the number of digits in n; this
    number is only used to determine splitting points and need not
    be exact.
    """

    if n < 0:
        return "-" + numeral(-n, base, size, digits)

    # Fast enough to do directly
    if size < 250:
        return small_numeral(n, base, digits)

    # Divide in half
    half = (size // 2) + (size & 1)
    A, B = divmod(n, base**half)
    ad = numeral(A, base, half, digits)
    bd = numeral(B, base, half, digits).rjust(half, "0")
    return ad + bd



def to_digits_exp(s, dps):
    """
    Helper function for representing the floating-point number s as a
    decimal with dps places. Returns (sign, string, exponent)
    containing '' or '-', the decimal digits as a string, and an
    integer for the decimal exponent.

    If inexact, the decimal representation is rounded toward zero.
    """

    # Extract sign so it doesn't mess up the string digit count
    sign = ''
    if s[0] < 0:
        sign = '-'
    s = fabs(s)
    man, exp, bc = s

    if not man:
        return '', '0', 0

    bitprec = int(dps * math.log(10,2)) + 10

    # Cut down to size
    # TODO: account for precision when doing this
    exp_from_1 = exp + bc
    if abs(exp) > 2500:
        # Set b = int(exp * log(2)/log(10))
        # If exp is huge, we must use high-precision arithmetic to
        # find the nearest power of ten
        expprec = bitcount(exp) + 5
        RF = round_floor
        tmp = from_int(exp, expprec, RF)
        tmp = fmul(tmp, flog2(expprec, RF), expprec, RF)
        tmp = fdiv(tmp, flog10(expprec, RF), expprec, RF)
        b = to_int(tmp)
        s = fdiv(s, fpow(ften, b, bitprec, RF), bitprec, RF)
        man, exp, bc = s
        exponent = b
    else:
        exponent = 0

    # First, calculate mantissa digits by converting to a binary
    # fixed-point number and then converting that number to
    # a decimal fixed-point number.
    fixprec = max(bitprec - exp, 0)
    fixdps = int(fixprec / math.log(10,2) + 0.5)
    sf = make_fixed(s, fixprec)
    sd = bin_to_radix(sf, fixprec, 10, fixdps)
    digits = numeral(sd, base=10, size=dps)

    exponent += len(digits) - fixdps - 1
    return sign, digits, exponent


def to_str(s, dps):
    """Convert a raw mpf to a decimal floating-point literal with at
    most `dps` decimal digits in the mantissa (not counting extra zeros
    that may be inserted for visual purposes).

    The literal is formatted so that it can be parsed back to a number
    by to_str, float() or Decimal()."""

    # Special numbers
    if s[2] == -1:
        return s[0]

    # to_digits_exp rounds to floor.
    # This sometimes kills some instances of "...00001"
    sign, digits, exponent = to_digits_exp(s, dps+3)

    # Rounding up kills some instances of "...99999"
    if len(digits) > dps and digits[dps] in '56789':
        digits2 = str(int(digits[:dps]) + 1)
        if len(digits2) > dps:
            digits2 = digits2[:dps]
            exponent += 1
        digits = digits2
    else:
        digits = digits[:dps]

    # Prettify numbers close to unit magnitude
    if -(dps//3) < exponent < dps:
        if exponent < 0:
            digits = ("0"*(-exponent)) + digits
            split = 1
        else:
            split = exponent + 1
        exponent = 0
    else:
        split = 1

    digits = (digits[:split] + "." + digits[split:])

    # Clean up trailing zeros
    digits = digits.rstrip('0')
    if digits[-1] == ".":
        digits += "0"

    if exponent == 0: return sign + digits
    if exponent > 0: return sign + digits + "e+" + str(exponent)
    if exponent < 0: return sign + digits + "e" + str(exponent)


def str_to_man_exp(x, base=10):
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
    x = int(x, base)
    return x, exp

special_str = {'inf':finf, '+inf':finf, '-inf':fninf, 'nan':fnan}

def from_str(x, prec, rounding):
    """Create a raw mpf from a decimal literal, rounding in the
    specified direction if the input number cannot be represented
    exactly as a binary floating-point number with the given number of
    bits. The literal syntax accepted is the same as for Python
    floats.

    TODO: the rounding does not work properly for large exponents.
    """
    if x in special_str:
        return special_str[x]

    man, exp = str_to_man_exp(x, base=10)

    # XXX: appropriate cutoffs & track direction
    # note no factors of 5
    if abs(exp) > 400:
        s = from_int(man, prec+10, round_floor)
        s = fmul(s, fpow(ften, exp, prec+10, round_floor), prec, rounding)
    else:
        if exp >= 0:
            s = from_int(man * 10**exp, prec, rounding)
        else:
            s = from_rational(man, 10**-exp, prec, rounding)
    return s

# Binary string conversion. These are currently mainly used for debugging
# and could use some improvement in the future

def from_bstr(x):
    man, exp = str_to_man_exp(x, base=2)
    bc = bitcount(man)
    return normalize(man, exp, bc, bc, round_floor)

def to_bstr(x):
    man, exp, bc = x
    return numeral(man, size=bitcount(man), base=2) + ("e%i" % exp)

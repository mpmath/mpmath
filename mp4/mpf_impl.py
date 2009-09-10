from lib import MPZ, MPZ_0, MPZ_1, bitcount2, bitcount
from lib import (i_trim, i_add, i_div, i_sqrt, i_exp, i_ln, i_atan2, i_pi, i_ln2,
    i_pow, i_pow_n, i_cos_sin, i_ln_hypot)

from math import ldexp as math_ldexp
from math import frexp as math_frexp

inttypes = tuple(set((int, long, type(MPZ(1)))))

# Note: the literal 1e1000 is broken in Python 2.4; this seems to work
py_inf = 1e300 * 1e300
py_ninf = -py_inf
py_nan = py_inf - py_inf

# Flags for special values
S_NORMAL = 0

S_HAS_INF = 1
S_HAS_NAN = 2
S_HAS_RE = 4
S_HAS_IM = 8

S_REAL_INF = S_HAS_INF | S_HAS_RE | 16
S_IMAG_INF = S_HAS_INF | S_HAS_IM | 32
S_UNSIGNED_INF = S_HAS_INF | S_HAS_RE | S_HAS_IM | 64

S_NAN = S_HAS_NAN | S_HAS_RE | S_HAS_IM | 128

# Constants
V_0 = (MPZ_0, 0, MPZ_0, 0, S_NORMAL)
V_1 = (MPZ_1, 0, MPZ_0, 0, S_NORMAL)
V_2 = (MPZ_1, 1, MPZ_0, 0, S_NORMAL)
V_n1 = (-MPZ_1, 0, MPZ_0, 0, S_NORMAL)
V_half = (MPZ_1, -1, MPZ_0, 0, S_NORMAL)
V_j = (MPZ_0, 0, MPZ_1, 0, S_NORMAL)
V_nj = (MPZ_0, 0, -MPZ_1, 0, S_NORMAL)
V_inf = (MPZ_1, 0, MPZ_0, 0, S_REAL_INF)
V_ninf = (-MPZ_1, 0, MPZ_0, 0, S_REAL_INF)
V_infj = (MPZ_0, 0, MPZ_1, 0, S_IMAG_INF)
V_ninfj = (MPZ_0, 0, -MPZ_1, 0, S_IMAG_INF)
V_uinf = (MPZ_0, 0, MPZ_0, 0, S_UNSIGNED_INF)
V_nan = (MPZ_0, 0, MPZ_0, 0, S_NAN)

# TODO: implement here
from libmpf import from_man_exp, to_str, finf, fninf, fnan

def V_to_str(x, dps, **kwargs):
    am, ae, bm, be, special = x
    if special:
        if x == V_inf: return "+inf"
        if x == V_ninf: return "-inf"
        if x == V_infj: return "+inf*j"
        if x == V_ninfj: return "-inf*j"
        if x == V_uinf: return "uinf"
        if x == V_nan: return "nan"
        raise NotImplementedError
    if am:
        ar = from_man_exp(am, ae)
        astr = to_str(ar, dps, **kwargs)
        if bm:
            if bm < 0:
                br = from_man_exp(-bm, be)
                return "(" + astr + " - " + to_str(br, dps, **kwargs) + "j)"
            else:
                br = from_man_exp(bm, be)
                return "(" + astr + " + " + to_str(br, dps, **kwargs) + "j)"
        return to_str(ar, dps, **kwargs)
    if bm:
        br = from_man_exp(bm, be)
        bstr = to_str(br, dps, **kwargs)
        return bstr + "j"
    return "0.0"

def V_to_bstr(x, base=2, prefix=True, exp="e", imag="j"):
    r"""
    Represent verbatim as a literal with a binary exponent.
    The base for the mantissa can be 2, 10 or 16. The exponent is
    always represented in decimal. With *prefix=True*, a binary or hex
    mantissa is prefixed with '0b' or '0x'.

    >>> V_to_bstr(V_from(0.48291015625))
    '0b1111011101e-11'
    >>> V_to_bstr(V_from(0.48291015625), base=10)
    '989e-11'
    >>> V_to_bstr(V_0)
    '0b0e0'
    >>> V_to_bstr(V_from(5j))
    '0b101e0j'
    >>> V_to_bstr(V_from(-5j))
    '-0b101e0j'
    >>> V_to_bstr(V_from(3+4j))
    '0b11e0 + 0b1e0j'
    >>> V_to_bstr(V_from(-9-7j))
    '-0b1001e0 - 0b111e0j'
    >>> V_to_bstr(V_from(1.25))
    '0b101e-2'
    >>> V_to_bstr(V_from(1.25+3.5j), prefix=False, exp="*2**", imag="*I")
    '101*2**-2 + 111*2**-2*I'
    >>> V_to_bstr(V_from(0), prefix="", exp="*2**")
    '0*2**0'
    >>> V_to_bstr(V_from(1.25+3.5j), base=10, prefix=False, exp="*2**", imag="*I")
    '5*2**-2 + 7*2**-2*I'

    """
    a, b, c, d, special = x
    if prefix:
        if   base == 2: prefix = "0b"
        elif base == 16: prefix = "0x"
        elif base == 10: prefix = ""
        else: raise ValueError("expected base 2, 10, or 16")
    else:
        prefix = ""
    if special:
        raise NotImplementedError
    if a or not c:
        s = "%s%s%s%s%s" % \
            (['','-'][a<0], prefix, numeral(abs(a), base), exp, b)
        if c:
            s += " %s %s%s%s%s%s" % \
                (['+','-'][c<0], prefix, numeral(abs(c), base), exp, b, imag)
    else:
        s = "%s%s%s%s%s%s" % \
            (['','-'][c<0], prefix, numeral(abs(c), base), exp, d, imag)
    return s

_illegal_return_values = set([0.0, py_inf, py_ninf])

def V_to_float(x, catch_overflow=False, catch_underflow=False, catch_complex=True):
    r"""
    Convert value to a Python float. The result is exact if the mantissa
    bitcount is `\le 53` and no underflow/overflow occurs. Rounding,
    if necessary, happens in an arbitrary direction that depends on
    the Python implementation's :func:`math.ldexp` (to obtain a desired
    rounding, pre-round the number to 53 bits or less).

    If the number is too large or too small to represent as a regular
    float, it will be converted to an infinity or 0.0. An *OverflowError*
    can be forced instead by setting the appropriate keyword argument.
    Complex-valued numbers raise *ValueError* by default; otherwise
    a float NaN is returned.

    >>> V_to_float(V_from(2.23))
    2.23
    >>> V_to_float(V_from(0.0))
    0.0
    >>> V_to_float(V_from(-1))
    -1.0
    >>> V_to_float(V_from(100000))
    100000.0
    >>> V_to_float(V_inf) == py_inf
    True
    >>> V_to_float(V_ninf) == py_ninf
    True
    >>> n = V_to_float(V_nan)
    >>> n != n
    True

    Overflow handling:

    >>> big = (MPZ(3), 2000, MPZ_0, 0, S_NORMAL)
    >>> negbig = (MPZ(-3), 2000, MPZ_0, 0, S_NORMAL)
    >>> big2 = (MPZ(3)**1000, 50, MPZ_0, 0, S_NORMAL)
    >>> negbig2 = (-MPZ(3)**1000, 50, MPZ_0, 0, S_NORMAL)
    >>> small = (MPZ(3), -2000, MPZ_0, 0, S_NORMAL)
    >>> wide = (MPZ(3)**1000, -1585, MPZ_0, 0, S_NORMAL)
    >>> wide_large = (MPZ(3)**1000, 10000, MPZ_0, 0, S_NORMAL)
    >>> wide_small = (MPZ(3)**1000, -10000, MPZ_0, 0, S_NORMAL)
    >>> cplx = V_from(2+3j)
    >>> V_to_float(big) == py_inf
    True
    >>> V_to_float(negbig) == py_ninf
    True
    >>> V_to_float(small)
    0.0
    >>> V_to_float(big, catch_overflow=True)
    Traceback (most recent call last):
      ...
    OverflowError: value is too large to convert to a float
    >>> V_to_float(negbig, catch_overflow=True)
    Traceback (most recent call last):
      ...
    OverflowError: value is too large to convert to a float
    >>> V_to_float(small)
    0.0
    >>> V_to_float(small, catch_underflow=True)
    Traceback (most recent call last):
      ...
    OverflowError: value is too small to convert to a float
    >>> V_to_float(cplx)
    Traceback (most recent call last):
      ...
    ValueError: float() requires a real number
    >>> n = V_to_float(cplx, catch_complex=False)
    >>> n != n
    True
    >>> V_to_float(wide)
    0.97434237824355119
    >>> V_to_float(wide_large) == py_inf
    True
    >>> V_to_float(wide_small)
    0.0
    >>> V_to_float(big2) == py_inf
    True
    >>> V_to_float(negbig2) == py_ninf
    True
    >>> V_to_float(big2, catch_overflow=True)
    Traceback (most recent call last):
      ...
    OverflowError: value is too large to convert to a float
    >>> V_to_float(negbig2, catch_overflow=True)
    Traceback (most recent call last):
      ...
    OverflowError: value is too large to convert to a float

    """
    a, b, c, d, special = x
    if c or special:
        if c:
            if catch_complex:
                raise ValueError("float() requires a real number")
            else:
                return py_nan
        if special == S_POS_INF: return py_inf
        if special == S_NEG_INF: return py_ninf
        if special == S_NAN: return py_nan
    if not a:
        return 0.0
    try:
        # There is a high probability that the following will succeed
        if abs(b) < 1000:
            v = math_ldexp(a, b)
            if v not in _illegal_return_values:
                return v
    except OverflowError:
        pass
    bc = bitcount2(a)
    try:
        # Try resizing the mantissa. Overflow may still happen here.
        if bc > 100:
            n = bc - 53
            a >>= n
            b += n
        v = math_ldexp(a, b)
        if v not in _illegal_return_values:
            return v
    except OverflowError:
        pass
    # Overflow to infinity
    if b + bc > 0:
        if catch_overflow:
            raise OverflowError("value is too large to convert to a float")
        if a < 0:
            return py_ninf
        else:
            return py_inf
    else:
        # XXX: never reached!
        if catch_underflow:
            raise OverflowError("value is too small to convert to a float")
        # Underflow to zero
        return 0.0

def V_to_complex(x, **kwargs):
    a, b, c, d, special = x
    if special:
        raise NotImplementedError
    re = V_to_float((a, b, MPZ_0, 0, special), **kwargs)
    im = V_to_float((c, d, MPZ_0, 0, special), **kwargs)
    return complex(re, im)

def V_from(x):
    """
    Create exact value from a given number.

    Recognized types are ``int``, ``long``, ``MPZ``, ``float``, ``complex``.
    Any type that supports exact conversion to a builtin Python number type and
    back (e.g. NumPy floats; presumably any integer type) will also work.
    For unknown types, conversion in both directions is attempted in order
    not to accidentally perform a lossy conversion. Strings are rejected.

    >>> V_from(py_inf) == V_inf
    True
    >>> V_from(py_ninf) == V_ninf
    True
    >>> for typ in (int, long, MPZ, float, complex):
    ...     assert V_from(typ(3)) == (3, 0, 0, 0, 0)
    ...     assert V_from(typ(10)) == (5, 1, 0, 0, 0)
    ...     assert V_from(typ(-2000)) == (-125, 4, 0, 0, 0)
    ...
    >>> V_from(1.125 + 25.5j) == (9, -3, 51, -1, 0)
    True
    >>> try:
    ...     from numpy import float64, complex128
    ...     assert V_from(float64(10.5)) == V_from(10.5)
    ...     assert V_from(complex128(10.5+7.5j)) == V_from(10.5+7.5j)
    ... except ImportError:
    ...     pass
    ...
    >>> V_from('0.25')
    Traceback (most recent call last):
      ...
    TypeError: unable to create number from '0.25'


    """
    if x in V_cache:
        return V_cache[x]
    tp = type(x)
    if tp in inttypes:
        a = MPZ(x)
        if a & 1:
            b = 0
        else:
            a, b = i_trim(a, 0)
        return a, b, MPZ_0, 0, S_NORMAL
    elif tp is float:
        if x != x:
            return V_nan
        # Note: inf, -inf already caught in cache
        a, b = math_frexp(x)
        a = int(a*(2**53))
        b -= 53
        if not a & 1:
            a, b = i_trim(a, b)
        return MPZ(a), b, MPZ_0, 0, S_NORMAL
    elif tp is complex:
        a, b, _, _, s1 = V_from(x.real)
        i = x.imag
        if i:
            c, d, _, _, s2 = V_from(i)
        else:
            return a, b, MPZ_0, 0, S_NORMAL
        if s1 or s2:
            raise NotImplementedError
        return a, b, c, d, S_NORMAL
    elif tp is tuple:
        if len(x) == 2:
            a, b = i_trim(*x)
            return a, b, MPZ_0, 0, S_NORMAL
        elif len(x) == 4:
            return V_from_mpf(x)
    for builtin_type in [int, float, complex]:
        try:
            y = builtin_type(x)
            if tp(y) == x and x == y:
                return V_from(y)
        except (TypeError, ValueError, OverflowError):
            pass
    raise TypeError("unable to create number from %s" % repr(x))

# Precache common int, float & complex values
# This is particularly intended to speed up operations with
# a Python number literal, e.g. 0.5*x
V_cache = {0:V_0, py_inf:V_inf, py_ninf:V_ninf}

for _x in (range(-256, 257) + [
    0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3,
    -0.1, -0.2, -0.3, -0.4, -0.6, -0.7, -0.8, -0.9, -1.1, -1.2, -1.3,
    -2.5, -1.5, -0.5, 0.5, 1.5, 2.5,
    -1.25, -0.75, -0.25, 0.25, 0.75, 1.25,
     2+2j,  2+1j,  2+0.5j,  2-0.5j,  2-1j,  2-2j,
     1+2j,  1+1j,  1+0.5j,  1-0.5j,  1-1j,  1-2j,
      +2j,   +1j,   +0.5j,   -0.5j,   -1j,   -2j,
    -1+2j, -1+1j, -1-0.5j, -1+0.5j, -1-1j, -1-2j,
    -2+2j, -2+1j, -2-0.5j, -2+0.5j, -2-1j, -2-2j]):
    V_cache[_x] = V_from(_x)

def V_to_old(x):
    """
    Convert to old format for an mpmath number.
    """
    pass

def V_from_old(x):
    """
    Convert from old format for an mpmath number.
    """
    pass

def V_from_mpf(x):
    sign, man, exp, bc = x
    if man:
        if sign:
            man = -man
        return man, exp, MPZ_0, 0, S_NORMAL
    elif exp:
        if x == finf: return V_inf
        if x == fninf: return V_ninf
        if x == fnan: return V_nan
    else:
        return V_0

def V_from_mpc(x):
    re, im = x
    rsign, rman, rexp, rbc = re
    isign, iman, iexp, ibc = im
    if rman:
        if rsign:
            rman = -rman
    elif rexp:
        raise NotImplementedError
    if iman:
        if isign:
            iman = -iman
    elif iexp:
        raise NotImplementedError
    return rman, rexp, iman, iexp, S_NORMAL

from libmpf import finf, fninf, fnan

def V_to_mpf(x):
    a, b, c, d, special = x
    if special:
        if x == V_inf: return finf
        if x == V_ninf: return fninf
        if x == V_nan: return fnan
        raise NotImplementedError
    if c:
        raise ValueError
    sign = int(a < 0)
    if sign:
        a = -a
    return sign, a, b, bitcount(a)

def V_to_mpc(x):
    a, b, c, d, special = x
    if special:
        raise NotImplementedError
    asign = int(a < 0)
    if asign:
        a = -a
    csign = int(c < 0)
    if csign:
        c = -c
    return (asign, a, b, bitcount(a)), (csign, c, d, bitcount(c))

def V_hash(x):
    a, b, c, d, special = x
    if c or special:
        if special:
            # XXX
            if x == V_inf: return hash(py_inf)
            if x == V_ninf: return hash(py_ninf)
        return hash(V_to_complex(x))
    if 0 <= b < 1000:
        return hash(a<<b)
    try:
        return hash(V_to_float(x))
    except:
        return hash(x)

import re
get_complex = re.compile(r'^\(?(?P<re>[\+\-]?\d*\.?\d*(e[\+\-]?\d+)?)??'
                         r'(?P<im>[\+\-]?\d*\.?\d*(e[\+\-]?\d+)?j)?\)?$')

from libmpf import from_str

def V_from_str(s, prec):
    if '/' in s:
        fract = s.split('/')
        p, q = fract
        return V_from_rational(MPZ(p), q, prec)
    if 'j' in s.lower():
        s = s.lower().replace(' ', '')
        match = get_complex.match(s)
        re = match.group('re')
        if not re:
            re = '0'
        im = match.group('im').rstrip('j')
        re = V_from_str(re, prec)
        im = V_from_str(im, prec)
        return V_add(re, V_mul(im, V_j), prec)
    return V_from_mpf(from_str(s, prec, 'n'))

#------------------------------------------------------------------------------#
# Test and comparison functions                                                #
#------------------------------------------------------------------------------#

def V_cmp(x, y):
    if x == y:
        return 0
    am, ae, bm, be, xspecial = x
    cm, ce, dm, de, yspecial = y
    if bm or dm:
        raise ValueError("cannot compare complex numbers")
    if xspecial or yspecial:
        if y == V_inf:
            #if x == V_inf:
            #    return 0
            return -1
        if x == V_inf:
            #if y == V_inf:
            #    return 0
            return 1
        if x == V_ninf:
            #if y == V_ninf:
            #    return 0
            return -1
        if y == V_ninf:
            #if x == V_ninf:
            #    return 0
            return 1
        return -1
        #raise NotImplementedError
    # First of all compare signs
    if am > 0:
        if cm <= 0:
            return 1
    elif am < 0:
        if cm >= 0:
            return -1
    else:
        if cm < 0: return 1
        if cm > 0: return -1
    # If not too huge shifts involved, shift and subtract
    if ae > ce:
        offset = ae - ce
        if offset < 1000:
            return cmp((am<<offset) - cm, 0)
    elif ae < ce:
        offset = ce - ae
        if offset < 1000:
            return cmp(am - (cm<<offset), 0)
    else:
        return cmp(am, cm)
    # Fall back to addition
    m, e = i_add(am, ae, -cm, ce, 1, 'd')
    return cmp(m, 0)

def V_is_int(x, gaussian=False):
    """
    Determine whether *x* is integer-valued. Gaussian (complex) integers
    are recognized with *gaussian=True*.

    >>> V_is_int(V_from(0))
    True
    >>> V_is_int(V_from(10))
    True
    >>> V_is_int(V_from(-5))
    True
    >>> V_is_int(V_from(0.25))
    False
    >>> V_is_int(V_from(1+1j))
    False
    >>> V_is_int(V_from(1+1j), gaussian=True)
    True
    >>> V_is_int(V_inf); V_is_int(V_ninf); V_is_int(V_nan)
    False
    False
    False

    """
    a, b, c, d, special = x
    if special:
        return False
    if c:
        return gaussian and b >= 0 and d >= 0
    return b >= 0

def V_is_even(x):
    """
    Determine whether *x* is an even integer.

    >>> V_is_even(V_from(0))
    True
    >>> V_is_even(V_from(1))
    False
    >>> V_is_even(V_from(-27))
    False
    >>> V_is_even(V_from(64))
    True
    >>> V_is_even(V_from(0.5))
    False
    >>> V_is_even(V_from(2+4j))
    False
    >>> V_is_even(V_inf)
    False
    """
    a, b, c, d, special = x
    return bool((b > 1 or not a) and not (c or special))

def V_is_odd(x):
    """
    Determine whether *x* is an odd integer.

    >>> V_is_odd(V_from(0))
    False
    >>> V_is_odd(V_from(1))
    True
    >>> V_is_odd(V_from(-27))
    True
    >>> V_is_odd(V_from(64))
    False
    >>> V_is_odd(V_from(0.5))
    False
    >>> V_is_odd(V_from(3+5j))
    False
    >>> V_is_odd(V_inf)
    False
    """
    a, b, c, d, special = x
    return bool(b == 0 and a and not (c or special))

def V_is_halfint(x, strict=False):
    """
    Determine whether `2x` is an integer, i.e. whether `x` is of the
    form `n/2` for an integer `n`. With *strict=True*, returns True
    only if `x = n/2` with `n` odd.

    >>> V_is_halfint(V_from(-3.5))
    True
    >>> V_is_halfint(V_from(0.25))
    False
    >>> V_is_halfint(V_from(0))
    True
    >>> V_is_halfint(V_from(0), strict=True)
    False
    >>> V_is_halfint(V_from(4))
    True
    >>> V_is_halfint(V_from(4), strict=True)
    False
    >>> V_is_halfint(V_from(1+1j))
    False
    >>> V_is_halfint(V_inf)
    False

    """
    a, b, c, d, special = x
    if c or special:
        return False
    if strict:
        return b == -1
    return b >= -1

def V_is_inf(x):
    raise NotImplementedError

def V_is_nan(x):
    return x[-1] == S_NAN

def V_mag(x):
    r"""
    Quickly returns an integer providing an upper bound for `\log_2(|x|)`.
    If *x* is pure real or pure imaginary, this returns exactly
    `\lfloor \log_2(|x|) \rfloor + 1`. For general complex *x*, it returns a
    number between `\lceil \log_2(|x|) \rceil` and
    `\lceil \log_2(|x|) \rceil + 2` (a tighter bound could be given).

    If *x* is zero or infinite, an infinite-valued Python float is returned.

    >>> from math import log, ceil
    >>> for x in [0.25, 125, 0.001, -0.001j, 3+4j, 1000+500j]:
    ...     print V_mag(V_from(x)), int(ceil(log(abs(x), 2)))
    ...
    -1 -2
    7 7
    -9 -9
    -9 -9
    4 3
    11 11
    >>> V_mag(V_inf) == py_inf
    True
    >>> V_mag(V_ninf) == py_inf
    True
    >>> V_mag(V_0) == py_ninf
    True
    >>> V_mag(V_nan)
    Traceback (most recent call last):
      ...
    ValueError

    """
    a, b, c, d, special = x
    if special:
        if special & S_HAS_NAN:
            raise ValueError
        return py_inf
    if a:
        m = bitcount2(a) + b
        if c:
            m = max(m, bitcount2(c) + d) + 1
    elif c:
        m = bitcount2(c) + d
    else:
        return py_ninf
    return m

def V_nint_distance(x):
    a, b, c, d, special = x
    if special:
        raise ValueError

    if c:
        im_dist = d + bitcount2(c)
    else:
        im_dist = py_ninf

    if a:
        bc = bitcount2(a)
        shift = b + bc
        if shift < -1:
            n = 0
            re_dist = shift
        else:
            if b >= 0:
                n = a << b
                re_dist = py_ninf
            else:
                if shift >= 0:
                    xfixed = a << shift
                else:
                    xfixed = a >> (-shift)
                n1 = xfixed >> bc
                n2 = -((-xfixed) >> bc)
                dist1 = abs(xfixed - (n1<<bc))
                dist2 = abs(xfixed - (n2<<bc))
                if dist1 < dist2:
                    re_dist = dist1
                    n = n1
                else:
                    re_dist = dist2
                    n = n2
                if re_dist:
                    re_dist = bitcount2(re_dist) - bc
                else:
                    re_dist = py_ninf
    else:
        re_dist = py_ninf
        n = 0

    return n, max(re_dist, im_dist)

def V_sign(x, prec=0):
    a, b, c, d, special = x
    if special:
        if x == V_inf: return V_1
        if x == V_ninf: return V_n1
        if x == V_infj: return V_j
        if x == V_ninfj: return V_nj
        return V_nan
    if a:
        if c:
            if not prec:
                raise ValueError
            # sign(x) = x / abs(x)
            wp = prec+20
            tm, te = i_add(a*a, b+b, c*c, d+d, wp)
            tm, te = i_sqrt(tm, te, wp)
            rm, re = i_div(a, b, tm, te, prec, 'n')
            im, ie = i_div(c, d, tm, te, prec, 'n')
            return rm, re, im, ie, S_NORMAL
        if a > 0:
            return V_1
        return V_n1
    if c:
        if c > 0:
            return V_j
        return V_nj
    return V_0

"""
def V_sign(x):
    a, b, c, d, special = x
    if c or special:
        raise ValueError
    return cmp(a, 0)
"""

def V_sign2(x):
    """
    Returns (sign(re(x)), sign(im(x))).
    """
    a, b, c, d, special = x
    if special:
        raise ValueError
    return (cmp(a, 0), cmp(c, 0))

def V_isign(x):
    a, b, c, d, special = x
    if special:
        if special == S_REAL_INF:
            if a > 0: return 0
            if a < 0: return 2
        if special == S_IMAG_INF:
            if c > 0: return 1
            if c < 0: return 3
    elif (a or c) and (not (a and c)):
        if a > 0: return 0
        if a < 0: return 2
        if c > 0: return 1
        if c < 0: return 3
    return None


def V_re(x, prec=0):
    a, b, c, d, special = x
    if special:
        if x == V_inf or x == V_ninf:
            return x
        raise NotImplementedError
    return a, b, MPZ_0, 0, S_NORMAL

def V_im(x, prec=0):
    a, b, c, d, special = x
    if special:
        if x == V_infj: return V_inf
        if x == V_ninfj: return V_ninf
        if x == V_nan: return V_nan
        return V_0 # XXX
    @property
    def imag(x):
        am, ae, bm, be, special = v = x._v
        if special:
            if v == V_infj: return x._ctx.infj
            if v == V_ninfj: return x._ctx.ninfj
            if v == V_nan: return x._ctx.nan

        raise NotImplementedError
    return c, d, MPZ_0, 0, S_NORMAL

#------------------------------------------------------------------------------#
# Arithmetic functions                                                         #
#------------------------------------------------------------------------------#



def V_add(xval, yval, prec):
    """

    """
    am, ae, bm, be, xspecial = xval
    cm, ce, dm, de, yspecial = yval
    if xspecial or yspecial:
        if (xspecial and yspecial) and xval != yval:
            return V_nan
        elif xspecial:
            return xval
        else:
            return yval
    if cm:
        if am:
            am, ae = i_add(am, ae, cm, ce, prec, 'n')
        else:
            am = cm
            ae = ce
    if dm:
        if bm:
            bm, be = i_add(bm, be, dm, de, prec, 'n')
        else:
            bm = dm
            be = de
    return am, ae, bm, be, xspecial

def V_sub(xval, yval, prec):
    """

    """
    am, ae, bm, be, xspecial = xval
    cm, ce, dm, de, yspecial = yval
    if xspecial or yspecial:
        return V_add(xval, (-cm, ce, -dm, de, yspecial), prec)
    if cm:
        if am:
            am, ae = i_add(am, ae, -cm, ce, prec, 'n')
        else:
            am = -cm
            ae = ce
    if dm:
        if bm:
            bm, be = i_add(bm, be, -dm, de, prec, 'n')
        else:
            bm = -dm
            be = de
    return am, ae, bm, be, xspecial

# XXX
from gmpy import _mpmath_mult as i_mul

def V_mul(x, y, prec=0):
    """
    Multiply *x* and *y*, rounding the mantissas to a precision of *prec*
    bits in the given directions. Correct rounding is used in all cases.

    With *prec=0*, the multiplication is performed exactly. Note that
    this can generate very wide mantissas and possibly overflow if multiplying
    complex numbers of vastly different magnitude. Exact multiplication will
    also generate a large mantissa (and thus be slow) if used for several
    successive multiplications. A single exact real multiplication is, however,
    significantly faster than a rounded multiplication at low precision.

    Multiplication of an infinity by a pure real or pure imaginary number
    is done with the convention that the result is a real or imaginary
    infinity and that the orthogonal component is zero (rather than NaN,
    which would be the result using the general rectangular multiplication
    formula for complex numbers).

    >>> cases = [0, 1.5, 7j, 2+3.25j]
    >>> for x in cases:
    ...     for y in cases:
    ...         assert V_mul(V_from(x), V_from(y), 53) == V_from(x*y)
    ...
    >>>


    """
    a, b, c, d, xspecial = x
    e, f, g, h, yspecial = y
    if xspecial or yspecial:
        if (xspecial & S_HAS_NAN) or (yspecial & S_HAS_NAN):
            return V_nan
        if x == V_0 or y == V_0:
            return V_nan
        s1 = V_isign(x)
        s2 = V_isign(y)
        if s1 is None or s2 is None:
            raise ValueError
        sign = (s1+s2) % 4
        if sign == 0: return V_inf
        if sign == 1: return V_infj
        if sign == 2: return V_ninf
        if sign == 3: return V_ninfj
        raise NotImplementedError
    # Have any imaginary parts?
    if c or g:
        # x is a full complex number
        if a and c:
            # full complex-complex multiply
            if e and g:
                sm, se = i_add(a*e, b+f, -c*g, d+h, prec, 'n')
                tm, te = i_add(a*g, b+h, c*e, d+f, prec, 'n')
                return sm, se, tm, te, S_NORMAL
            elif e: sm, se, tm, te = a*e, b+f, c*e, d+f
            elif g: sm, se, tm, te = -c*g, d+h, a*g, b+h
            else:   return V_0
        # y is a full complex number
        elif e and g:
            if a:   sm, se, tm, te = a*e, b+f, a*g, b+h
            elif c: sm, se, tm, te = -c*g, d+h, c*e, d+f
            else:   return V_0
        elif c:
            if e:   sm, se, tm, te = MPZ_0, 0, c*e, d+f
            elif g: sm, se, tm, te = -c*g, d+h, MPZ_0, 0
            else:   return V_0
        elif g and a: sm, se, tm, te = MPZ_0, 0, a*g, b+h
        else:         return V_0
    elif a and e:
        if prec:
            sm, se = i_trim(a*e, b+f, prec, 'n')
            #sm, se = i_mul(a, b, e, f, prec, 'n')
        else:
            sm = a*e
            se = b+f
        return sm, se, MPZ_0, 0, S_NORMAL
        #sm, se, tm, te = a*e, b+f, MPZ_0, 0
    else:
        return V_0
    if sm: sm, se = i_trim(sm, se, prec, 'n')
    if tm: tm, te = i_trim(tm, te, prec, 'n')
    return sm, se, tm, te, S_NORMAL

def V_sum(terms, prec=0, absolute=False):
    #s = V_0
    #for term in terms:
    #    s = V_add(s, term, prec+20)
    #return s

    if prec:
        wp = prec + 20
    else:
        wp = 0
    rem = MPZ_0
    ree = 0
    imm = MPZ_0
    ime = 0
    maxshift = max(200, 2*prec)
    for am, ae, bm, be, special in terms:
        if special:
            raise NotImplementedError
        if absolute:
            if am and bm:
                raise NotImplementedError
            if am < 0: am = -am
            if bm < 0: bm = -bm
        if am:
            offset = ree - ae
            if offset >= 0:
                if offset > maxshift:
                    rem, ree = i_add(rem, ree, am, ae, wp)
                else:
                    rem = (rem<<offset) + am
                    ree = ae
            else:
                if offset < -maxshift:
                    rem, ree = i_add(rem, ree, am, ae, wp)
                else:
                    rem += (am<<(-offset))

        if bm:
            offset = ime - be
            if offset >= 0:
                if offset > maxshift:
                    imm, ime = i_add(imm, ime, bm, be, wp)
                else:
                    imm = (imm<<offset) + bm
                    ime = be
            else:
                if offset < -maxshift:
                    imm, ime = i_add(imm, ime, bm, be, wp)
                else:
                    imm += (bm<<(-offset))

    rem, ree = i_trim(rem, ree, prec, 'n')
    imm, ime = i_trim(imm, ime, prec, 'n')

    return rem, ree, imm, ime, S_NORMAL

def V_dot(pairs, prec=0):
    def iterpairs():
        for x, y in pairs:
            a, b, c, d, s1 = x
            e, f, g, h, s2 = y
            if s1 or s2:
                yield V_mul(x, y)
            if c or g:
                if a and c:
                    if e and g:
                        yield a*e, b+f, a*g, b+h, S_NORMAL
                        yield -c*g, d+h, c*e, d+f, S_NORMAL
                    elif e: yield a*e, b+f, c*e, d+f, S_NORMAL
                    elif g: yield -c*g, d+h, a*g, b+h, S_NORMAL
                    else:   continue
                elif e and g:
                    if   a: yield a*e, b+f, a*g, b+h, S_NORMAL
                    elif c: yield -c*g, d+h, c*e, d+f, S_NORMAL
                    else:   continue
                elif c:
                    if   e: yield MPZ_0, 0, c*e, d+f, S_NORMAL
                    elif g: yield -c*g, d+h, MPZ_0, 0, S_NORMAL
                    else:   continue
                elif g and a:
                    yield MPZ_0, 0, a*g, b+h, S_NORMAL
                else:
                    continue
            elif a and e:
                yield a*e, b+f, MPZ_0, 0, S_NORMAL
            else:
                continue
    return V_sum(iterpairs(), prec)



def V_sqr(x, prec=0):
    """
    Square *x*. This is slightly faster than but otherwise fully equivalent
    to *V_mul(x, x, ...)*.

    >>> cases = [0, 1.5, -3.25, 2j, -3j, 2.25+1j]
    >>> for x in cases:
    ...     assert V_sqr(V_from(x), 53) == V_from(x**2)
    ...
    >>>

    """
    a, b, c, d, special = x
    if special:
        return V_mul(x, x, prec)
    # nonreal
    if c:
        if a:
            # (a+c)^2 = (a^2 - c^2) + 2aci
            sm, se = i_add(a*a, b+b, -c*c, d+d, prec, 'n')
            if prec:
                tm, te = i_trim(a*c, b+d+1, prec, 'n')
            else:
                tm, te = a*c, b+d+1
            return sm, se, tm, te, special
        else:
            if prec:
                sm, se = i_trim(-c*c, d+d, prec, 'n')
                return sm, se, a, b, special
            else:
                return -c*c, d+d, a, b, special
    elif a:
        if prec:
            sm, se = i_trim(a*a, b+b, prec, 'n')
            return sm, se, c, d, special
        else:
            return a*a, b+b, c, d, special
    # zero
    else:
        return x

def V_cub(x, prec=0):
    a, b, c, d, special = x
    # TODO: fast complex
    if special or c:
        return V_mul(x, V_mul(x, x, prec+20), prec)
    elif a:
        if prec:
            sm, se = i_trim(a**3, b*3, prec, 'n')
            return sm, se, c, d, special
        else:
            return a**3, b*3, c, d, special
    # zero
    else:
        return x


def V_recip(x, prec, raises=False):
    """
    Multiplicative inverse, 1/x.

    >>> print V_to_complex(V_recip(V_from(5), 53, 'n'))
    (0.2+0j)
    >>> print V_to_complex(V_recip(V_from(0.2), 53, 'n'))
    (5+0j)
    >>> print V_to_complex(V_recip(V_from(5j), 53, 'n'))
    -0.2j
    >>> print V_to_complex(V_recip(V_from(-5j), 53, 'n'))
    0.2j
    >>> print V_to_complex(V_recip(V_from(2+5j), 53, 'n'))
    (0.0689655172414-0.172413793103j)
    >>> V_recip(V_0, 53, raises=True)
    Traceback (most recent call last):
      ...
    ZeroDivisionError: floating-point division by zero

    """
    am, ae, bm, be, xspecial = x
    if xspecial:
        return V_div(V_1, x, prec, raises)
    if bm:
        # 1/(a+bi) = (a-ib) / (a^2 + b^2)
        if am:
            # XXX: rounding direction
            mm, me = i_add(am*am, ae+ae, bm*bm, be+be, prec+20)
            rm, re = i_div(am, ae, mm, me, prec, 'n')
            tm, te = i_div(-bm, be, mm, me, prec, 'n')
            return rm, re, tm, te, S_NORMAL
        # 1/(bi) = -i/b
        sm, se = i_div(-MPZ_1, 0, bm, be, prec, 'n')
        return MPZ_0, 0, sm, se, S_NORMAL
    # 1/a
    if am:
        sm, se = i_div(MPZ_1, 0, am, ae, prec, 'n')
        return sm, se, MPZ_0, 0, S_NORMAL
    if raises:
        raise ZeroDivisionError #("floating-point division by zero")
    return V_uinf

def V_div(x, y, prec, raises=True):
    """
    x/y
    """
    am, ae, bm, be, xspecial = x
    cm, ce, dm, de, yspecial = y
    if xspecial or yspecial:
        if (xspecial & S_HAS_NAN) or (yspecial & S_HAS_NAN):
            return V_nan
        if y == V_0:
            if raises:
                raise ZeroDivisionError
            return V_uinf
        # Now y != 0
        if xspecial & S_HAS_INF:
            if yspecial & S_HAS_INF:
                return V_nan
            s1 = V_isign(x)
            s2 = V_isign(y)
            if s1 is None or s2 is None:
                raise ValueError
            sign = (s1-s2) % 4
            if sign == 0: return V_inf
            if sign == 1: return V_infj
            if sign == 2: return V_ninf
            if sign == 3: return V_ninfj
            # here: sign!
        if yspecial & S_HAS_INF:
            return V_0
        raise NotImplementedError
    if cm:
        if dm:
            # TODO: faster when pure real or pure imag numerator ...
            # (a+bi) / (c+di) = [(ac+bd) + i(bc-ad)] / (c^d + d^2)
            wp = prec + 20
            mm, me = i_add(cm*cm, ce+ce, dm*dm, de+de, wp)
            mm, me = i_div(MPZ_1, 0, mm, me, wp)
            rm, re = i_add(am*cm, ae+ce, bm*dm, be+de, wp)
            rm, re = i_trim(rm*mm, re+me, prec, 'n')
            im, ie = i_add(bm*cm, be+ce, -am*dm, ae+de, wp)
            im, ie = i_trim(im*mm, ie+me, prec, 'n')
            return rm, re, im, ie, S_NORMAL

        # (a+bi) / c = a/c + (b/c)i
        if am: sm, se = i_div(am, ae, cm, ce, prec, 'n')
        else:  sm, se = MPZ_0, 0
        if bm: tm, te = i_div(bm, be, cm, ce, prec, 'n')
        else:  tm, te = MPZ_0, 0
        return sm, se, tm, te, S_NORMAL
    if dm:
        # (a+bi)/(di) = (b-ai)/d
        if bm: sm, se = i_div(bm, be, dm, de, prec, 'n')
        else:  sm, se = MPZ_0, 0
        if am: tm, te = i_div(-am, ae, dm, de, prec, 'n')
        else:  tm, te = MPZ_0, 0
        return sm, se, tm, te, S_NORMAL
    # Division by zero
    if raises:
        raise ZeroDivisionError #("floating-point division by zero")
    if am or bm:
        return V_uinf
    return V_nan

def V_from_rational(p, q, prec):
    a, b = i_div(MPZ(p), 0, MPZ(q), 0, prec, 'n')
    return a, b, MPZ_0, 0, S_NORMAL

from libmpf import mpf_mod

def V_mod(x, y, prec):
    #am, ae, bm, be, xspecial = x
    #cm, ce, dm, de, yspecial = y
    #if bm or dm or xspecial or yspecial:
    #    raise ValueError
    return V_from_mpf(mpf_mod(V_to_mpf(x), V_to_mpf(y), prec, 'n'))


def V_mul_i(x, n=1):
    """
    Multiply *x* exactly by `i^n` for integer `n`, i.e. rotate counterclockwise
    by `n 90^\circ` in the complex plane. In particular, *n = -1* divides
    by `i`.

    """
    n %= 4
    raise NotImplementedError

def V_mul_2exp(x, n):
    """
    Multiply *x* exactly by `2^n` for integer `n`.
    """
    a, b, c, d, s = x
    if s or not n:
        return x
    if a:
        if c:
            return a, b+n, c, d+n, s
        return a, b+n, c, d, s
    if c:
        return a, b, c, d+n, s
    return x

def V_sqrt(x, prec, rounding='n'):
    """
    Calculate the principal square root of *x*. Rounding is correct for
    pure real or pure imaginary square roots.

    >>> V_to_float(V_sqrt(V_from(0), 53, 'n'))
    0.0
    >>> V_to_float(V_sqrt(V_from(1), 53, 'n'))
    1.0
    >>> V_to_float(V_sqrt(V_from(2), 53, 'n'))
    1.4142135623730951
    >>> V_to_complex(V_sqrt(V_from(-1), 53, 'n'))
    1j
    >>> V_to_complex(V_sqrt(V_from(-2), 53, 'n'))
    1.4142135623730951j
    >>> V_to_complex(V_sqrt(V_from(2j), 53, 'n'))
    (1+1j)
    >>> V_to_complex(V_sqrt(V_from(-2j), 53, 'n'))
    (1-1j)
    >>> V_to_complex(V_sqrt(V_from(3j), 53, 'n'))
    (1.2247448713915889+1.2247448713915889j)
    >>> V_to_complex(V_sqrt(V_from(-3j), 53, 'n'))
    (1.2247448713915889-1.2247448713915889j)
    >>> V_to_complex(V_sqrt(V_from(2+3j), 53, 'n'))
    (1.6741492280355401+0.89597747612983814j)
    >>> V_to_complex(V_sqrt(V_from(2-3j), 53, 'n'))
    (1.6741492280355401-0.89597747612983814j)
    >>> V_to_complex(V_sqrt(V_from(-2+3j), 53, 'n'))
    (0.89597747612983814+1.6741492280355401j)
    >>> V_to_complex(V_sqrt(V_from(-2-3j), 53, 'n'))
    (0.89597747612983814-1.6741492280355401j)

    """
    am, ae, bm, be, special = x
    if special:
        if x == V_inf: return x
        if x == V_ninf: return V_infj
        return V_nan
    if am:
        if bm:
            # General complex square root: we use the formula
            # sqrt(a+bi) = u/2 + (b/u)i
            # u = sqrt(2(sqrt(a^2+b^2) + a))
            # Note: u is guaranteed to be real because sqrt(a^2+b^2) > |a|. We
            # round upwards to ensure that this also holds numerically.
            # TODO: check for cancellation when adding sqrt(a^2+b^2) to a
            wp = prec + 30
            um, ue = i_add(am*am, ae+ae, bm*bm, be+be, wp, 'u')
            um, ue = i_sqrt(um, ue, wp, 'u')
            um, ue = i_add(um, ue, am, ae, wp, 'u')
            um, ue = i_sqrt(um, ue+1, wp, 'f')
            cman, cexp = i_trim(um, ue-1, prec, rounding)
            dman, dexp = i_div(bm, be, um, ue, prec, rounding)
            return cman, cexp, dman, dexp, S_NORMAL
        if am < 0:
            bm, be = i_sqrt(-am, ae, prec, rounding)
            return MPZ_0, 0, bm, be, S_NORMAL
        am, ae = i_sqrt(am, ae, prec, rounding)
        return am, ae, bm, be, S_NORMAL
    # Pure imaginary number:
    # sqrt(i*b) = (1+i) * sqrt(b/2)
    if bm:
        if bm > 0:
            am, ae = i_sqrt(bm, be-1, prec, rounding)
            return am, ae, am, ae, S_NORMAL
        else:
            am, ae = i_sqrt(-bm, be-1, prec, rounding)
            return am, ae, -am, ae, S_NORMAL
    return V_0


# FIXME: the following needs huge amounts of cleanup

def V_pow(x, y, prec=0):
    am, ae, bm, be, xspecial = x
    cm, ce, dm, de, yspecial = y
    if xspecial or yspecial:
        if y == V_inf and not (xspecial or bm):
            c = V_cmp(x, V_1)
            if c < 0: return V_0
            if c > 0: return V_inf
        if x == V_inf and not yspecial:
            if cm < 0:
                return V_0
            if cm > 0:
                return V_inf
        if x == V_ninf and not yspecial:
            # ^ integer
            if cm > 0 and ce >= 0 and dm == 0:
                if ce == 0:
                    return V_ninf
                return V_inf
        if xspecial & S_HAS_INF and not yspecial:
            if cm < 0:
                return V_0

        return V_nan

        raise NotImplementedError

    # Pure real exponent
    if not dm:
        #(a+bi)^c
        if cm:
            # y = c is an integer
            if ce >= 0:
                n = cm << ce
                if n < 4:
                    if n == 0: return V_1
                    if n == 1: return x
                    if n == 2: return V_sqr(x, prec)
                    if n == 3: return V_cub(x, prec)
                    if n == -1: return V_recip(x, prec)
                    # if n == -2:
                # a^n or (bi)^n
                if not (am and bm):
                    # a^n
                    if am:
                        am, ae = i_pow_n(am, ae, n, prec, 'n')
                        return am, ae, bm, be, S_NORMAL
                    # (bi)^n = i^n * b^n
                    else:
                        # TODO: if m in (2, 3), change rounding direction
                        m = n & 3
                        bm, be = i_pow_n(bm, be, n, prec, 'n')
                        if m == 0: return bm, be, MPZ_0, 0, S_NORMAL
                        if m == 1: return MPZ_0, 0, bm, be, S_NORMAL
                        if m == 2: return -bm, be, MPZ_0, 0, S_NORMAL
                        if m == 3: return MPZ_0, 0, -bm, be, S_NORMAL

            # Square root
            # TODO: generalize
            elif ce == -1 and cm == 1:
                return V_sqrt(x, prec)

            # (a+bi)^c, general real c
            elif bm:
                pass

            # a^c
            else:
                if am > 0:
                    am, ae = i_pow(am, ae, cm, ce, prec, 'n')
                    return am, ae, MPZ_0, 0, S_NORMAL
                #pass
        else:
            # x^0 = 1
            return V_1

    if not (am or bm):
        if cm <= 0:
            raise ZeroDivisionError
        if cm > 0:
            return V_0

    pure_imag = (am < 0) and (ce == -1) and (not (bm or dm))

    # XXX
    if bm or dm or 1:
        wp = prec + 20
        # r^2
        r2m, r2e = i_add(am*am, ae+ae, bm*bm, be+be, wp)
        # t (angle)
        tm, te = i_atan2(am, ae, bm, be, wp)
        # 2*log(r) = log(r^2)
        lm, le = i_ln(r2m, r2e, wp)

        um, ue = i_add(dm*lm, de+le-1, cm*tm, ce+te, wp)

        cxm, cxe, sxm, sxe = i_cos_sin(um, ue, wp)

        # r^c = (r^2)^(c/2)
        # rcm, rce = i_pow(r2m, r2e, cm, ce-1, wp)
        # exp(c*log(r) - d*t)
        mm, me = i_add(cm*lm, ce+le-1, -dm*tm, de+te, wp)

        mm, me = i_exp(mm, me, wp)

        rem, ree = i_trim(mm*cxm, me+cxe, prec, 'n')
        imm, ime = i_trim(mm*sxm, me+sxe, prec, 'n')

        if pure_imag:
            rem, ree = MPZ_0, 0

        return rem, ree, imm, ime, S_NORMAL


def V_cbrt(x, **kwargs):
    pass

def V_pow_n(x, n, **kwargs):
    """
    Compute the power `x^n` for integer *n*.
    """
    pass

def V_n_pow(n, x, **kwargs):
    """
    Compute the power `n^x` for integer *n*.
    """
    pass


def V_pow_r(x, p, q, **kwargs):
    """
    Compute the power `x^{p/q}` for integers *p*, *q*.
    """
    pass

from libelefun import mpf_nthroot
from libmpc import mpc_nthroot

def V_nthroot(x, n, prec):
    if n < 0: return V_recip(V_nthroot(x, -n, prec+20), prec)
    if n == 0: return V_1
    if n == 1: return x
    if n == 2: return V_sqrt(x, prec)
    a, b, c, d, special = x
    if not c:
        try:
            return V_from_mpf(mpf_nthroot(V_to_mpf(x), n, prec, 'n'))
        except:
            pass
    return V_from_mpc(mpc_nthroot(V_to_mpc(x), n, prec, 'n'))

    """
    # if n == 3: return V_cbrt(x, prec)
    a, b, c, d, special = x
    if b or a < 0:
        raise NotImplementedError
        wp = prec + 20
        return V_exp(V_div(V_ln(x, wp), wp), V_from(n), prec)
    a, b = i_nthroot(x, n, prec, 'n')
    return a, b, MPZ_0, 0, S_NORMAL
    """


def V_exp(x, prec):
    a, b, c, d, special = x
    if special:
        if x == V_inf: return x
        if x == V_ninf: return V_0
        return V_nan
        #raise NotImplementedError
    if c:
        if a:
            wp = prec + 10
            mm, me = i_exp(a, b, wp)
            cm, ce, sm, se = i_cos_sin(c, d, wp)
            a, b = i_trim(mm*cm, me+ce, prec, 'n')
            c, d = i_trim(mm*sm, me+se, prec, 'n')
            return a, b, c, d, S_NORMAL
        else:
            a, b, c, d = i_cos_sin(c, d, prec, 'n')
            return a, b, c, d, S_NORMAL
    if a:
        a, b = i_exp(a, b, prec, 'n')
        return a, b, MPZ_0, 0, S_NORMAL
    return V_1


def V_pi(prec):
    wp = prec+20
    m, e = i_trim(i_pi(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

from libelefun import ln2_fixed, pi_fixed, e_fixed, ln10_fixed, phi_fixed, degree_fixed
from gammazeta import euler_fixed, catalan_fixed, glaisher_fixed, apery_fixed, khinchin_fixed, twinprime_fixed, mertens_fixed

def V_ln2(prec):
    wp = prec+20
    m, e = i_trim(ln2_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_e(prec):
    wp = prec+20
    m, e = i_trim(e_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_ln10(prec):
    wp = prec+20
    m, e = i_trim(ln10_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_phi(prec):
    wp = prec+20
    m, e = i_trim(phi_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_degree(prec):
    wp = prec+20
    m, e = i_trim(degree_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_euler(prec):
    wp = prec+20
    m, e = i_trim(euler_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_catalan(prec):
    wp = prec+20
    m, e = i_trim(catalan_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_glaisher(prec):
    wp = prec+20
    m, e = i_trim(glaisher_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_apery(prec):
    wp = prec+20
    m, e = i_trim(e_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_khinchin(prec):
    wp = prec+20
    m, e = i_trim(khinchin_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_twinprime(prec):
    wp = prec+20
    m, e = i_trim(twinprime_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL

def V_mertens(prec):
    wp = prec+20
    m, e = i_trim(mertens_fixed(wp), -wp, prec, 'n')
    return m, e, MPZ_0, 0, S_NORMAL


def V_exppi(x, **kwargs): pass
def V_cos(x, **kwargs): pass
def V_sin(x, **kwargs): pass
def V_tan(x, **kwargs): pass
def V_cospi(x, **kwargs): pass
def V_sinpi(x, **kwargs): pass
def V_tanpi(x, **kwargs): pass
def V_sinc(x, **kwargs): pass
def V_sincpi(x, **kwargs): pass

def V_ln(x, prec):
    a, b, c, d, special = x
    if special:
        raise NotImplementedError
    rm, re = i_ln_hypot(a, b, c, d, prec, 'n')
    im, ie = i_atan2(a, b, c, d, prec, 'n')
    return rm, re, im, ie, S_NORMAL

def V_arg(x, prec):
    a, b, c, d, special = x
    if special:
        if x == V_nan:
            return x
    rm, re = i_atan2(a, b, c, d, prec, 'n')
    return rm, re, MPZ_0, 0, S_NORMAL

def V_log(x, b=None, **kwargs):
    pass



# Special functions?

def V_gamma(x, prec): pass
def V_fac(x, prec): pass
def V_lgamma(x, prec): pass
def V_rgamma(x, prec): pass
def V_hypsum(a_fp, a_rat, b_fp, b_rat, x, prec): pass

# XXX: real only?
def V_erf(x, prec): pass
def V_erfc(x, prec): pass
def V_erfi(x, prec): pass

# XXX!

from libmpf import mpf_hypot
from libelefun import mpf_atan2
from libhyper import mpf_expint, mpf_ei, mpc_ei, mpf_e1, mpc_e1, mpf_agm, mpc_agm, mpf_besseljn, mpc_besseljn
from gammazeta import mpf_psi, mpc_psi, mpf_bernoulli

def V_besselj(n, x, prec):
    try:
        return V_from_mpf(mpf_besseljn(n, V_to_mpf(x), prec, 'n'))
    except:
        return V_from_mpc(mpc_besseljn(n, V_to_mpc(x), prec, 'n'))

def V_ei(x, prec):
    try:
        return V_from_mpf(mpf_ei(V_to_mpf(x), prec, 'n'))
    except:
        return V_from_mpc(mpc_ei(V_to_mpc(x), prec, 'n'))

def V_e1(x, prec):
    try:
        return V_from_mpf(mpf_e1(V_to_mpf(x), prec, 'n'))
    except:
        return V_from_mpc(mpc_e1(V_to_mpc(x), prec, 'n'))

def V_gamma_upper_int(n, z, prec):
    n = int(n)
    if n == 0:
        return V_e1(z, prec)
    if z[2]:
        raise NotImplementedError
    real, imag = mpf_expint(n, V_to_mpf(z), prec, 'n', gamma=True)
    if imag is None:
        return V_from_mpf(real)
    else:
        return V_from_mpc((real, imag))

def V_expint_int(n, z, prec):
    n = int(n)
    if n == 1:
        return V_e1(z, prec)
    if z[2]:
        raise NotImplementedError
    real, imag = mpf_expint(n, V_to_mpf(z), prec, 'n')
    if imag is None:
        return V_from_mpf(real)
    else:
        return V_from_mpc((real, imag))

def V_agm(a, b, prec):
    try:
        aa = V_to_mpf(a)
        bb = V_to_mpf(b)
        v = mpf_agm(aa, bb, prec, 'n')
        return V_from_mpf(v)
    except:
        aa = V_to_mpc(a)
        bb = V_to_mpc(b)
        v = mpc_agm(aa, bb, prec, 'n')
        return V_from_mpc(v)

def V_hypot(x, y, prec):
    return V_from_mpf(mpf_hypot(V_to_mpf(x), V_to_mpf(y), prec, 'n'))

def V_atan2(x, y, prec):
    return V_from_mpf(mpf_atan2(V_to_mpf(x), V_to_mpf(y), prec, 'n'))

def V_psi(n, x, prec):
    try:
        return V_from_mpf(mpf_psi(n, V_to_mpf(x), prec, 'n'))
    except:
        return V_from_mpc(mpc_psi(n, V_to_mpc(x), prec, 'n'))

def V_bernoulli(n, prec):
    return V_from_mpf(mpf_bernoulli(int(n), prec, 'n'))

if __name__ == "__main__":
    import doctest, trace, sys
    tracer = trace.Trace(ignoredirs=[sys.prefix, sys.exec_prefix],
        trace=0, count=1)
    tracer.run('doctest.testmod()')
    r = tracer.results()
    r.write_results(show_missing=True, summary=True, coverdir=".")

    #doctest.testmod()
    print "Passed tests"

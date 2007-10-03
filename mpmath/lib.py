"""
This module contains "low-level" functions for multiprecision floating-
point arithmetic implemented in pure Python. The code is written in a
functional style for simplicity and speed.

A floating-point number x = man * 2**exp is represented by the tuple
(man, exp, bc) where man is the mantissa, exp is the exponent, and bc
is the number of bits in the mantissa. To simplify equality testing,
the mantissa always gets normalized by removing trailing zero bits.

The bitcount is slightly redundant to store in the number, but may as
well be reused since it always gets computed during normalization,
and slightly speeds up subsequent operations on a number.

"""


#----------------------------------------------------------------------------#
#                                                                            #
#                             General utilities                              #
#                                                                            #
#----------------------------------------------------------------------------#

import math
import decimal

# Same as standard Python float
STANDARD_PREC = 53


# All supported rounding modes. We define them as integer constants for easy
# management, but change __repr__ to give more information on inspection

class RoundingMode(int):
    def __new__(cls, level, name):
        a = int.__new__(cls, level)
        a.name = name
        return a
    def __repr__(self): return self.name

ROUND_DOWN    = RoundingMode(1, 'ROUND_DOWN')
ROUND_UP      = RoundingMode(2, 'ROUND_UP')
ROUND_FLOOR   = RoundingMode(3, 'ROUND_FLOOR')
ROUND_CEILING = RoundingMode(4, 'ROUND_CEILING')
ROUND_HALF_UP = RoundingMode(5, 'ROUND_HALF_UP')
ROUND_HALF_DOWN = RoundingMode(6, 'ROUND_HALF_DOWN')
ROUND_HALF_EVEN = RoundingMode(7, 'ROUND_HALF_EVEN')


def giant_steps(start, target):
    """Generate a list of precisions ranging from 'start' to 'target'
    that doubles with each step. This is used by quadratically
    convergent iterations (that is, Newton iterations), where we want
    to keep the precision at the same level as the accuracy in each
    step to minimize work.

    For example, to find a sequence of precisions to reach 1000 bits
    starting from a 53-bit estimate, giant_steps(53, 1000) gives

        [64, 126, 251, 501, 1000]

    So, if we start Newton's method with a 53-bit accurate initial
    guess, the first iteration should be carried out at 64-bit
    precision, the second at 126-bit precision, and so on.

    Note the conservative rounding (1000 to 501, etc); this is used
    guard against unit errors in the last place."""
    L = [target]
    while L[-1] > start*2:
        L = L + [L[-1]//2 + 1]
    return L[::-1]


#----------------------------------------------------------------------------#
#                                                                            #
#                             Radix conversion                               #
#                                                                            #
#----------------------------------------------------------------------------#

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

def decimal_to_binary(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    dps = int(prec*LOG2_10) + 5
    prec_ = getctx().prec
    getctx().prec = dps
    d = Dec(x).normalize()
    sgn, digits, dexp = d.as_tuple()
    d = d * Dec(10)**(-dexp)
    power = fpow(ften, -dexp, prec+5)
    y = fdiv(float_from_int(int(d), prec+5), power, prec, rounding)
    getctx().prec = prec_
    return y

def bin_to_radix(x, xbits, base, bdigits):
    return x * (base**bdigits) >> xbits

def small_numeral(n, base=10, digits='0123456789abcdefghijklmnopqrstuvwxyz'):
    # Calculate numeral of n*(base**digits) in the given base
    if base == 10:
        return str(n)
    digs = []
    while n:
        n, digit = divmod(n, base)
        digs.append(digits[digit])
    return "".join(digs[::-1])

# TODO: speed up for bases 2, 4, 8, 16, ...
def fixed_to_str(x, base, digits, verbose=False):
    if digits < 789:
        return small_numeral(x, base)
    half = (digits // 2) + (digits & 1)
    if verbose and half > 50000: print "  dividing..."
    A, B = divmod(x, base**half)
    ad = fixed_to_str(A, base, half)
    bd = fixed_to_str(B, base, half).rjust(half, "0")
    return ad + bd


#----------------------------------------------------------------------------#
#                                                                            #
#                          Bit manipulation, etc                             #
#                                                                            #
#----------------------------------------------------------------------------#

def make_fixed(s, prec):
    """Convert a floating-point number to a fixed-point big integer"""
    man, exp, bc = s
    offset = exp + prec
    if offset >= 0:
        return man << offset
    else:
        return man >> (-offset)

def bitcount(n, log=math.log, table=(0,1,2,2,3,3,3,3,4,4,4,4,4,4,4,4)):
    """Give size of n in bits; i.e. the position of the highest set bit
    in n. If n is negative, the absolute value is used. The bitcount of
    zero is taken to be 0."""

    if not n: return 0
    if n < 0: n = -n

    # math.log gives a good estimate, and never overflows, but
    # is not always exact. Subtract 2 to underestimate, then
    # count remaining bits by table lookup
    bc = int(log(n, 2)) - 2
    if bc < 0:
        bc = 0
    return bc + table[n >> bc]

# from decimal.py -- faster for small precs
def bitcount2(n, correction = {
        '0': 4, '1': 3, '2': 2, '3': 2,
        '4': 1, '5': 1, '6': 1, '7': 1,
        '8': 0, '9': 0, 'a': 0, 'b': 0,
        'c': 0, 'd': 0, 'e': 0, 'f': 0}):
    if n < 0:
        n = -n
    hex_n = "%x" % n
    return 4*len(hex_n) - correction[hex_n[0]]

def trailing_zeros(n):
    """Count trailing zero bits in an integer. If n is negative, it is
    replaced by its absolute value."""
    if n & 1: return 0
    if not n: return 0
    if n < 0: n = -n
    t = 0
    while not n & 0xffffffffffffffff: n >>= 64; t += 64
    while not n & 0xff: n >>= 8; t += 8
    while not n & 1: n >>= 1; t += 1
    return t

def rshift_quick(x, n):
    """For an integer x, calculate x >> n with the fastest (floor)
    rounding. Unlike the plain Python expression (x >> n), n is
    allowed to be negative, in which case a left shift is performed."""
    if n >= 0: return x >> n
    else:      return x << (-n)

def lshift_quick(x, n):
    """For an integer x, calculate x << n. Unlike the plain Python
    expression (x << n), n is allowed to be negative, in which case a
    right shift with default (floor) rounding is performed."""
    if n >= 0: return x << n
    else:      return x >> (-n)

def rshift(x, n, rounding):
    """Shift x (a plain Python integer) n bits to the right (i.e.,
    calculate x/(2**n)), and round to the nearest integer in accordance
    with the specified rounding mode. The exponent n may be negative,
    in which case x is shifted to the left (and no rounding is
    necessary)."""

    if not n or not x:
        return x
    # Support left-shifting (no rounding needed)
    if n < 0:
        return x << -n

    # To get away easily, we exploit the fact that Python rounds positive
    # integers toward zero and negative integers away from zero when dividing
    # or shifting. The simplest roundings can be handled entirely
    # through shifts:
    if rounding == ROUND_FLOOR:
        return x >> n
    elif rounding < ROUND_HALF_UP:
        if rounding == ROUND_DOWN:
            if x > 0: return x >> n
            else:     return -((-x) >> n)
        if rounding == ROUND_UP:
            if x > 0: return -((-x) >> n)
            else:     return x >> n
        if rounding == ROUND_CEILING:
            return -((-x) >> n)

    # Here we need to inspect the bits around the cutoff point
    if x > 0: t = x >> (n-1)
    else:     t = (-x) >> (n-1)
    if t & 1:
        if rounding == ROUND_HALF_UP or \
           (rounding == ROUND_HALF_DOWN and x & ((1<<(n-1))-1)) or \
           (rounding == ROUND_HALF_EVEN and (t&2 or x & ((1<<(n-1))-1))):
            if x > 0:  return (t>>1)+1
            else:      return -((t>>1)+1)
    if x > 0: return t>>1
    else:     return -(t>>1)

def normalize(man, exp, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """Normalize the binary floating-point number represented by
    man * 2**exp to the specified precision level, rounding if the
    number of bits in the mantissa exceeds prec. The mantissa is also
    stripped of trailing zero bits, and its bits are counted. The
    returned value is a tuple (man, exp, bc)."""

    if not man:
        return 0, 0, 0
    if prec < 100:
        bc = bitcount2(man)
    else:
        bc = bitcount(man)
    if bc > prec:
        # Right shifting by bc-prec nearly always guarantees that
        # the result has at most prec bits. There is one exceptional
        # case: if abs(man) is 1 less than a power of two and rounding
        # is done away from zero, it turns into the higher power of two.
        # This case must be handled separately; otherwise a bc of 0
        # gets returned.

        # TODO: by comparing sign and rounding mode, we could just
        # return (+/- 1, exp+bc, 1) right here
        absman = man
        if absman < 0:
            absman = -absman
        if not ((absman+1) & absman):
            man = rshift(man, bc-prec, rounding)
            exp += (bc - prec)
            bc = bitcount(man)
        # Default case
        else:
            man = rshift(man, bc-prec, rounding)
            exp += (bc - prec)
            bc = prec
    # Strip trailing zeros
    if not man & 1:
        tr = trailing_zeros(man)
        if tr:
            man >>= tr
            exp += tr
            bc -= tr
    #assert bitcount(man) <= prec
    if not man:
        return 0, 0, 0
    return man, exp, bc

#----------------------------------------------------------------------------#
#                                                                            #
#                            Type conversion                                 #
#                                                                            #
#----------------------------------------------------------------------------#

def float_from_int(n, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    return normalize(n, 0, prec, rounding)

def float_from_rational(p, q, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """Create floating-point number from a rational number p/q"""
    n = prec + bitcount(q) + 2
    return normalize((p<<n)//q, -n, prec, rounding)

def float_from_pyfloat(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
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


fzero = float_from_int(0)
fone = float_from_int(1)
ftwo = float_from_int(2)
ften = float_from_int(10)
fhalf = float_from_rational(1, 2)
assert fhalf == float_from_pyfloat(0.5)


#----------------------------------------------------------------------------#
#                                                                            #
#                                Comparison                                  #
#                                                                            #
#----------------------------------------------------------------------------#

def feq(s, t):
    """Floating-point equality testing. The numbers are assumed to
    be normalized, meaning that this simply performs tuple comparison."""
    return s == t

def fcmp(s, t):
    """Floating-point comparison. Return -1 if s < t, 0 if s == t,
    and 1 if s > t."""

    # An inequality between two numbers s and t is determined by looking
    # at the value of s-t. A full floating-point subtraction is relatively
    # slow, so we first try to look at the exponents and signs of s and t.
    sman, sexp, sbc = s
    tman, texp, tbc = t

    # Very easy cases: check for 0's and opposite signs
    if not tman: return cmp(sman, 0)
    if not sman: return cmp(0, tman)
    if sman > 0 and tman < 0: return 1
    if sman < 0 and tman > 0: return -1

    # In this case, the numbers likely have the same magnitude
    if sexp == texp: return cmp(sman, tman)

    # The numbers have the same sign but different exponents. In this
    # case we try to determine if they are of different magnitude by
    # checking the position of the highest set bit in each number.
    a = sbc + sexp
    b = tbc + texp
    if sman > 0:
        if a < b: return -1
        if a > b: return 1
    else:
        if a < b: return 1
        if a < b: return -1

    # The numbers have similar magnitude but different exponents.
    # So we subtract and check the sign of resulting mantissa.
    return cmp(fsub(s, t, 5, ROUND_FLOOR)[0], 0)


#----------------------------------------------------------------------------#
#                                                                            #
#                            Basic arithmetic                                #
#                                                                            #
#----------------------------------------------------------------------------#

def fadd(s, t, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """Floating-point addition. Given two tuples s and t containing the
    components of floating-point numbers, return their sum rounded to 'prec'
    bits using the 'rounding' mode, represented as a tuple of components."""

    #  General algorithm: we set min(s.exp, t.exp) = 0, perform exact integer
    #  addition, and then round the result.
    #                   exp = 0
    #                       |
    #                       v
    #          11111111100000   <-- s.man (padded with zeros from shifting)
    #      +        222222222   <-- t.man (no shifting necessary)
    #          --------------
    #      =   11111333333333

    # We assume that s has the higher exponent. If not, just switch them:
    if t[1] > s[1]:
        s, t = t, s

    sman, sexp, sbc = s
    tman, texp, tbc = t

    # Check if one operand is zero. Float(0) always has exp = 0; if the
    # other operand has a large exponent, its mantissa will unnecessarily
    # be shifted a huge number of bits if we don't check for this case.
    if not tman:
        return normalize(sman, sexp, prec, rounding)
    if not sman:
        return normalize(tman, texp, prec, rounding)

    # More generally, if one number is huge and the other is small,
    # and in particular, if their mantissas don't overlap at all at
    # the current precision level, we can avoid work.

    #       precision
    #    |            |
    #     111111111
    #  +                 222222222
    #     ------------------------
    #  #  1111111110000...

    # However, if the rounding isn't to nearest, correct rounding mandates
    # the result should be adjusted up or down.

    if sexp - texp > 10:
        bitdelta = (sbc + sexp) - (tbc + texp)
        if bitdelta > prec + 5:
            if rounding > 4:     # nearby rounding
                return normalize(sman, sexp, prec, rounding)

            # shift s and add a dummy bit outside the precision range to
            # force rounding up or down
            offset = min(bitdelta + 3, prec+3)
            sman <<= offset
            if tman > 0:
                sman += 1
            else:
                sman -= 1
            return normalize(sman, sexp-offset, prec, rounding)

    # General case
    return normalize(tman+(sman<<(sexp-texp)), texp, prec, rounding)


def fsub(s, t, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """Floating-point subtraction"""
    return fadd(s, (-t[0], t[1], t[2]), prec, rounding)

def fneg(s, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """Floating-point negation. In addition to changing sign, rounds to
    the specified precision."""
    return normalize(-s[0], s[1], prec, rounding)

def fabs(s, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    man, exp, bc = s
    if man < 0:
        return normalize(-man, exp, prec, rounding)
    return normalize(man, exp, prec, rounding)


def fmul(s, t, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """Floating-point multiplication"""

    sman, sexp, sbc = s
    tman, texp, tbc = t

    # This is very simple. A possible optimization would be to throw
    # away some bits when prec is much smaller than sbc+tbc
    return normalize(sman*tman, sexp+texp, prec, rounding)


def fdiv(s, t, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """Floating-point division"""
    sman, sexp, sbc = s
    tman, texp, tbc = t

    # Perform integer division between mantissas. The mantissa of s must
    # be padded appropriately to preserve accuracy.
    
    # Note: this algorithm does produce slightly wrong rounding in corner
    # cases. Padding with a few extra bits makes the chance very small.
    # Changing '12' to something lower will reveal the error in the
    # test_standard_float test case
    extra = prec - sbc + tbc + 12
    if extra < 12:
        extra = 12

    return normalize((sman<<extra)//tman, sexp-texp-extra, prec, rounding)


def fpow(s, n, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """Compute s**n, where n is an integer"""
    n = int(n)
    if n == 0: return fone
    if n == 1: return normalize(s[0], s[1], prec, rounding)
    if n == 2: return fmul(s, s, prec, rounding)
    if n == -1: return fdiv(fone, s, prec, rounding)
    if n < 0:
        return fdiv(fone, fpow(s, -n, prec+3, ROUND_FLOOR), prec, rounding)
    # Now we perform binary exponentiation. Need to estimate precision
    # to avoid rounding from temporary operations. Roughly log_2(n)
    # operations are performed.
    prec2 = prec + int(4*math.log(n, 2) + 4)
    man, exp, bc = normalize(s[0], s[1], prec2, ROUND_FLOOR)
    pm, pe, pbc = fone
    while n:
        if n & 1:
            pm, pe, pbc = normalize(pm*man, pe+exp, prec2, ROUND_FLOOR)
            n -= 1
        man, exp, bc = normalize(man*man, exp+exp, prec2, ROUND_FLOOR)
        n = n // 2
    return normalize(pm, pe, prec, rounding)


"""
Square roots are most efficiently computed with Newton's method.
Two functions are implemented: _sqrt_fixed and _sqrt_fixed2.

  _sqrt_fixed uses the iteration r_{n+1} = (r_n + y/r_n)/2,
  which is just Newton's method applied to the equation r**2 = y.

  _sqrt_fixed2 uses the iteration r_{n+1} = r_n*(3 - y*r_n**2)
  to calculate 1/sqrt(y), and then multiplies by y to obtain
  sqrt(y).

The first iteration is slightly faster at low precision levels, since
it essentially just requires one division at each step, compared to
the three multiplications in the second formula. However, the second
iteration is much better at extremely high precision levels. This is
due to the fact that Python uses the Karatsuba algorithm for integer
multiplication, which is asymptotically faster than its division
algorithm.

Both functions use fixed-point arithmetic and assume that the input y
is a big integer, i.e. given the integer y and precision prec,
they return floor(sqrt(x) * 2**prec) where y = floor(x * 2**prec).

The functions currently assume that x ~= 1. (TODO: make the code
work for x of arbitrary magnitude.) The main fsqrt() function
fiddles with the exponent of the input to reduce it to unit
magnitude before passing it to _sqrt_fixed or _sqrt_fixed2.

"""

def _sqrt_fixed(y, prec):
    # get 50-bit initial guess from regular float math
    if prec < 200:
        r = int(y**0.5 * 2.0**(50-prec*0.5))
    else:
        r = int((y >> (prec-100))**0.5)
    prevp = 50
    for p in giant_steps(50, prec+8):
        # Newton iteration: r_{n+1} = (r_{n} + y/r_{n})/2
        # print "sqrt", p
        r = lshift_quick(r, p-prevp-1) + (rshift_quick(y, prec-p-prevp+1)//r)
        prevp = p
    return r >> 8

def _sqrt_fixed2(y, prec):
    r = float_to_pyfloat(normalize(y, -prec, 64, ROUND_FLOOR)) ** -0.5
    r = int(r * 2**50)
    prevp = 50
    for p in giant_steps(50, prec+8):
        # print "sqrt", p
        r2 = rshift_quick(r*r, 2*prevp - p)
        A = lshift_quick(r, p-prevp)
        T = rshift_quick(y, prec-p)
        S = (T*r2) >> p
        B = (3 << p) - S
        r = (A*B)>>(p+1)
        prevp = p
    r = (r * y) >> prec
    return r >> 8

def fsqrt(s, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """
    If x is a positive Float, sqrt(x) returns the square root of x as a
    Float, rounded to the current working precision.
    """
    man, exp, bc = s
    if not man: return fzero
    if (man, exp) == (1, 0): return fone

    prec2 = prec + 10

    # Convert to a fixed-point number with prec bits. Adjust
    # exponents to be even so that they can be divided in half
    if prec2 & 1:
        prec2 += 1
    if exp & 1:
        exp -= 1
        man <<= 1
    shift = bitcount(man) - prec2
    shift -= shift & 1
    man = rshift_quick(man, shift)

    if prec < 65000:
        man = _sqrt_fixed(man, prec2)
    else:
        man = _sqrt_fixed2(man, prec2)

    return normalize(man, (exp+shift-prec2)//2, prec, rounding)


def fhypot(x, y, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    if y == fzero: return fabs(x, prec, rounding)
    if x == fzero: return fabs(y, prec, rounding)
    RF = ROUND_FLOOR
    hypot2 = fadd(fmul(x,x,prec+4,RF), fmul(y,y,prec+4,RF), prec+4, RF)
    return fsqrt(hypot2, prec, rounding)


#----------------------------------------------------------------------------#
#                                                                            #
#                         Mathematical constants                             #
#                                                                            #
#----------------------------------------------------------------------------#

# Only re-compute a constant if the precision level is raised
def _constmemo(f):
    f.memo_prec = -1
    f.memo_val = None
    def calc(prec):
        if prec == f.memo_prec: return f.memo_val
        if prec < f.memo_prec: return f.memo_val >> (f.memo_prec-prec)
        f.memo_val = f(prec)
        f.memo_prec = prec
        return f.memo_val
    return calc

# Evaluate a Machin-like formula, i.e., a rational combination of
# of acot(n) or acoth(n) for specific integer values of n
def _machin(coefs, prec, hyperbolic=False):
    prec += 10
    def acot(x):
        # Series expansion for atan/acot, optimized for integer arguments
        s = w = (1<<prec)//x; x2 = x*x; n = 3
        while 1:
            w //= x2
            term = w // n
            if not term: break
            if hyperbolic or n & 2 == 0: s += term
            else: s -= term
            n += 2
        return s
    s = 0
    for a, b in coefs:
        s += a * acot(b)
    return (s >> 10)

"""
At low precision, pi can be calculated easily using Machin's formula
pi = 16*acot(5)-4*acot(239). For high precision, we use the Brent-Salamin
algorithm based on the arithmetic-geometric mean. See for example Wikipedia
(http://en.wikipedia.org/wiki/Brent-Salamin_algorithm) or "Pi and the AGM" by
Jonathan and Peter Borwein (Wiley, 1987). The algorithm (as stated in the
Wikipedia article) consists of setting

  a_0 = 1;  b_0 = 1/sqrt(2);  t_0 = 1/4;  p_0 = 1

and computing

  a_{n+1} = (a_n + b_n)/2
  b_{n+1} = sqrt(a_n * b_n)
  t_{n+1} = t_n - p_n*(a_n - a_{n+1})**2
  p_{n+1} = 2*p_n

for n = 0, 1, 2, 3, ..., after which the approximation is given by
pi ~= (a_n + b_n)**2 / (4*t_n). Each step roughly doubles the number of
correct digits.
"""

def _pi_agm(prec, verbose=False, verbose_base=10):
    prec += 50
    a = 1 << prec
    if verbose: print "  computing initial square root..."
    b = _sqrt_fixed2(a>>1, prec)
    t = a >> 2
    p = 1
    step = 1
    while 1:
        an = (a+b)>>1
        adiff = a - an
        if verbose:
            logdiff = math.log(max(1, adiff), verbose_base)
            digits = int(prec/math.log(verbose_base,2) - logdiff)
            print "  iteration", step, ("(accuracy ~= %i base-%i digits)" % \
                (digits, verbose_base))
        if p > 16 and abs(adiff) < 1000:
            break
        prod = (a*b)>>prec
        b = _sqrt_fixed2(prod, prec)
        t = t - p*((adiff**2) >> prec)
        p = 2*p
        a = an
        step += 1
    if verbose: print "  final division"
    return ((((a+b)**2) >> 2) // t) >> 50

@_constmemo
def pi_fixed(prec):
    if prec < 10000:
        return _machin([(16, 5), (-4, 239)], prec)
    else:
        return _pi_agm(prec)

def fpi(prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """Compute a floating-point approximation of pi"""
    return normalize(pi_fixed(prec+5), -prec-5, prec, rounding)


# Logarithms of integers can be computed easily using
# Machin-like formulas

@_constmemo
def log2_fixed(prec):
    return _machin([(18, 26), (-2, 4801), (8, 8749)], prec, True)

def flog2(prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    return normalize(log2_fixed(prec+5), -prec-5, prec, rounding)

@_constmemo
def log10_fixed(prec):
    return _machin([(46, 31), (34, 49), (20, 161)], prec, True)

def flog10(prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    return normalize(log10_fixed(prec+5), -prec-5, prec, rounding)

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

@_constmemo
def gamma_fixed(prec):
    prec += 30
    # choose p such that exp(-4*(2**p)) < 2**-n
    p = int(math.log((prec/4) * math.log(2), 2)) + 1
    n = 1<<p; r=one=1<<prec
    H, A, B, npow, k, d = 0, 0, 0, 1, 1, 1
    while r:
        A += (r * H) >> prec
        B += r
        r = r * (n*n) // (k*k)
        H += one // k
        k += 1
    S = ((A<<prec) // B) - p*log2_fixed(prec)
    return S >> 30

def fgamma(prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    return normalize(gamma_fixed(prec+5), -prec-5, prec, rounding)


#----------------------------------------------------------------------------#
#                                                                            #
#                          Transcendental functions                          #
#                                                                            #
#----------------------------------------------------------------------------#

"""
The exponential function has a rapidly convergent Maclaurin series:

    exp(x) = 1 + x + x**2/2! + x**3/3! + x**4/4! + ...

The series can be summed very easily using fixed-point arithmetic.
The convergence can be improved further, using a trick due to
Richard P. Brent: instead of computing exp(x) directly, we choose a
small integer r (say, r=10) and compute exp(x/2**r)**(2**r).

The optimal value for r depends on the Python platform, the magnitude
of x and the target precision, and has to be estimated from
experimental timings. One test with x ~= 0.3 showed that
r = 2.2*prec**0.42 gave a good fit to the optimal values for r for
prec between 1 and 10000 bits, on one particular machine.

This optimization makes the summation about twice as fast at
low precision levels and much faster at high precision
(roughly five times faster at 1000 decimal digits).

If |x| is very large, we first rewrite it as t + n*log(2) with the
integer n chosen such that |t| <= log(2), and then calculate
exp(x) as exp(t)*(2**n), using the Maclaurin series for exp(t)
(the multiplication by 2**n just amounts to shifting the exponent).
"""

def exp_series(x, prec):
    r = int(2.2 * prec ** 0.42)
    # XXX: more careful calculation of guard bits
    guards = r + 3
    if prec > 60:
        guards += int(math.log(prec))
    prec2 = prec + guards
    x = rshift_quick(x, r - guards)
    s = (1 << prec2) + x
    a = x
    k = 2
    # Sum exp(x/2**r)
    while 1:
        a = ((a*x) >> prec2) // k
        if not a: break
        s += a
        k += 1
    # Calculate s**(2**r) by repeated squaring
    for j in range(r):
        s = (s*s) >> prec2
    return s >> guards

def fexp(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    man, exp, bc = x
    # extra precision needs to be similar in magnitude to log_2(|x|)
    prec2 = prec + 6 + max(0, bc+exp)
    t = make_fixed(x, prec2)
    # abs(x) > 1?
    if exp+bc > 1:  #fcmp(fabs(x), fone) > 0:
        lg2 = log2_fixed(prec2)
        n, t = divmod(t, lg2)
    else:
        n = 0
    return normalize(exp_series(t, prec2), -prec2+n, prec, rounding)


"""
The basic strategy for computing log(x) is to set r = log(x) and use
Newton's method to solve the equation exp(r) = x. We set the initial
value r_0 to math.log(x) and then iterate r_{n+1} = r_n + exp(-r_n) - 1
until convergence. As with square roots, we increase the working
precision dynamically during the process so that only one full-precision
evaluation of exp is required.

log(x) is small for most inputs, so the r values can safely be
computed using fixed-point arithmetic. However, when x has a very
large or small exponent, we can improve performance through the
normalization log(t * 2**n) = log(t) + n*log(2), choosing n such
that 0.5 <= t <= 1 (for example).

There are some caveats: if x is extremely close to 1, the working
precision must be increased to maintain high relative precision in the
output (alternatively, the series approximation for log(1+x) could
be used in that case).
"""

# This function performs the Newton iteration using fixed-point
# arithmetic. x is assumed to have magnitude ~= 1
def _log_newton(x, prec):
    extra = 8
    # 50-bit approximation
    #r = int(_clog(Float((x, -prec), 64)) * 2.0**50)
    fx = math.log(float_to_pyfloat((x, -prec, 1)))
    r = int(fx * 2.0**50)
    prevp = 50
    for p in giant_steps(50, prec+extra):
        rb = lshift_quick(r, p-prevp)
        e = exp_series(-rb, p)
        r = rb + ((rshift_quick(x, prec-p)*e)>>p) - (1 << p)
        prevp = p
    return r >> extra

def flog(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    if x == fzero: raise ValueError, "logarithm of 0"
    if x == fone:  return fzero
    man, exp, bc = x
    if man < 0: raise ValueError, "logarithm of a negative number"
    # Estimated precision needed for log(t) + n*log(2)
    prec2 = prec + int(math.log(1+abs(bc+exp), 2)) + 10
    # Watch out for the case when x is very close to 1
    if -1 < bc + exp < 2:
        near_one = fabs(fsub(x, fone))
        if near_one == 0:
            return fzero
        # estimate how close
        prec2 += -(near_one[1]) - bitcount(near_one[0])
    # Separate mantissa and exponent, calculate, join parts
    t = rshift_quick(man, bc-prec2)
    l = _log_newton(t, prec2)
    a = (exp + bc) * log2_fixed(prec2)
    return normalize(l+a, -prec2, prec, rounding)


## XXX: need to increase working precision here
#def fpow(x, y):
#    return exp(log(x) * y)


"""
We compute sin(x) around 0 from its Taylor series, and cos(x) around 0
from sqrt(1-sin(x)**2). This way we can simultaneously compute sin and
cos, which are often needed together (e.g. for the tangent function or
the complex exponential), with little extra cost compared to computing
just one of them. The main reason for computing sin first (and not cos
from sin) is to obtain high relative accuracy for x extremely close to
0, where the operation sqrt(1-cos(x)**2) can cause huge cancellations.

For any value of x, we can reduce it to the interval A = [-pi/4, pi/4]
(where the Taylor series converges quickly) by translations, changing
signs, and switching the roles of cos and sin:

   A : sin(x) = sin(x)           cos(x) = cos(x)
   B : sin(x) = cos(x-pi/2)      cos(x) = -sin(x-pi/2)
   C : sin(x) = -sin(x-pi)       cos(x) = -cos(x-pi)
   D : sin(x) = -cos(x-3*pi/2)   cos(x) = sin(x-3*pi/2)

|     A      |      B     |      C     |     D     |
v            v            v            v           v

   1 |  ____   ..........                            ____
     |      _..          ..                        __
     |      . __           .                     __
     |    ..    _           ..                  _
     |   .       __           .               __
-----| -.----------_-----------.-------------_-----------
     | .            _           ..          _           .
     |               __           .       __           .
     |                 _           ..    _           ..
     |                  __           . __           .
     |                    __         _..          ..
  -1 |                      _________   ..........
      0                       pi                     2*pi

"""

def _sin_series(x, prec):
    x2 = (x*x) >> prec
    s = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (-k*(k-1))
        s += a
        k += 2
    return s

def _trig_reduce(x, prec):
    pi_ = pi_fixed(prec)
    pi4 = pi_ >> 2
    pi2 = pi_ >> 1
    n, rem = divmod(x + pi4, pi2)
    rem -= pi4
    return n, rem

def cos_sin(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """
    cos_sin(x) calculates both the cosine and the sine of x rounded
    to the nearest Float value, and returns the tuple (cos(x), sin(x)).
    """
    man, exp, bc = x
    bits_from_unit = abs(bc + exp)
    prec2 = prec + bits_from_unit + 15
    xf = make_fixed(x, prec2)
    n, rx = _trig_reduce(xf, prec2)
    case = n % 4
    one = 1 << prec2
    if case == 0:
        s = _sin_series(rx, prec2)
        c = _sqrt_fixed(one - ((s*s)>>prec2), prec2)
    elif case == 1:
        c = -_sin_series(rx, prec2)
        s = _sqrt_fixed(one - ((c*c)>>prec2), prec2)
    elif case == 2:
        s = -_sin_series(rx, prec2)
        c = -_sqrt_fixed(one - ((s*s)>>prec2), prec2)
    elif case == 3:
        c = _sin_series(rx, prec2)
        s = -_sqrt_fixed(one - ((c*c)>>prec2), prec2)
    c = normalize(c, -prec2, prec, rounding)
    s = normalize(s, -prec2, prec, rounding)
    return c, s

def fcos(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    return cos_sin(x, prec, rounding)[0]

def fsin(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    return cos_sin(x, prec, rounding)[1]

def ftan(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    c, s = cos_sin(x, prec+2, rounding)
    return fdiv(s, c, prec, rounding)

#----------------------------------------------------------------------
# Inverse tangent
#

"""
Near x = 0, use atan(x) = x - x**3/3 + x**5/5 - ...
Near x = 1, use atan(x) = y/x * (1 + 2/3*y + 2*4/3/5*y**2 + ...)
where y = x**2/(1+x**2).

TODO: these series are not impressively fast. It is probably better
to calculate atan from tan, using Newton's method or even the
secant method.
"""

def _atan_series_1(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    man, exp, bc = x
    # Increase absolute precision when extremely close to 0
    bc = bitcount(man)
    diff = -(bc + exp)
    prec2 = prec
    if diff > 10:
        if 3*diff - 4 > prec:  # x**3 term vanishes; atan(x) ~x
            return normalize(man, exp, prec, rounding)
        prec2 = prec + diff
    prec2 += 15  # XXX: better estimate for number of guard bits
    x = make_fixed(x, prec2)
    x2 = (x*x)>>prec2; one = 1<<prec2; s=a=x
    for n in xrange(1, 1000000):
        a = (a*x2) >> prec2
        s += a // ((-1)**n * (n+n+1))
        if -100 < a < 100:
            break
    return normalize(s, -prec2, prec, rounding)

def _atan_series_2(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    prec2 = prec + 15
    x = make_fixed(x, prec2)
    one = 1<<prec2; x2 = (x*x)>>prec2; y=(x2<<prec2)//(one+x2)
    s = a = one
    for n in xrange(1, 1000000):
        a = ((a*y)>>prec2) * (2*n) // (2*n+1)
        if a < 100:
            break
        s += a
    return normalize(y*s//x, -prec2, prec, rounding)

_cutoff_1 = (5, -3, 3)   # ~0.6
_cutoff_2 = (3, -1, 2)   # 1.5

def fatan(x, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    if x[0] < 0:
        t = fatan(fneg(x), prec+4, ROUND_FLOOR)
        return normalize(-t[0], t[1], prec, rounding)
    if fcmp(x, _cutoff_1) < 0:
        return _atan_series_1(x, prec, rounding)
    if fcmp(x, _cutoff_2) < 0:
        return _atan_series_2(x, prec, rounding)
    # For large x, use atan(x) = pi/2 - atan(1/x)
    if x[1] > 10*prec:
        pi = fpi(prec, rounding)
        pihalf = pi[0], pi[1]-1, pi[2]
    else:
        pi = fpi(prec+4, ROUND_FLOOR)
        pihalf = pi[0], pi[1]-1, pi[2]
        t = fatan(fdiv(fone, x, prec+4, ROUND_FLOOR), prec+4, ROUND_FLOOR)
        return fsub(pihalf, t, prec, rounding)


#----------------------------------------------------------------------------#
#                                                                            #
#                           Complex functions                                #
#                                                                            #
#----------------------------------------------------------------------------#

def fcabs(re, im, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    return fhypot(re, im, prec, rounding)


# For complex square roots, we have sqrt(a+b*I) = sqrt((r+a)/2) + 
# I*b/sqrt(2*(r+a)) where r = abs(a+b*I), when a+b*I is not a negative
# real number (http://en.wikipedia.org/wiki/Square_root)
def fcsqrt(re, im, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    if re == im == fzero:
        return (re, im)
    if re[0] < 0 and im[0] == 0:
        return (fzero, fsqrt(fneg(re, prec, rounding), prec, rounding))
    RF, prec2 = ROUND_FLOOR, prec+20
    rpx = fadd(fcabs(re, im, prec2, RF), re, prec2, RF)
    are = fsqrt(fdiv(rpx, ftwo, prec2, RF), prec, rounding)
    aim = fdiv(im, fsqrt(fmul(rpx, ftwo, prec2, RF), prec2, RF), prec, rounding)
    return are, aim

def fcexp(re, im, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    mag = fexp(re, prec+4, ROUND_FLOOR)
    are, aim = cos_sin(im, prec+4, ROUND_FLOOR)
    return fmul(mag, are, prec, rounding), fmul(mag, aim, prec, rounding)

def fcsin(re, im, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    prec2 = prec+4
    RF = ROUND_FLOOR
    # sin(x+y*i) = sin(x)*cosh(y)+cos(x)*sinh(y)*i
    c, s = cos_sin(re, prec2, RF)
    expb1 = fexp(im, prec2, RF)
    expb2 = fdiv(fone, expb1, prec2, RF)
    ch = fmul(fadd(expb1, expb2, prec2, RF), fhalf, prec2, RF)
    sh = fmul(fsub(expb1, expb2, prec2, RF), fhalf, prec2, RF)
    return fmul(s, ch, prec, rounding), fmul(c, sh, prec, rounding)

def fccos(re, im, prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    prec2 = prec+4
    RF = ROUND_FLOOR
    # cos(x+y*i) = cos(x)*cosh(y)-sin(x)*sinh(y)*i
    c, s = cos_sin(re, prec2, RF)
    expb1 = fexp(im, prec2, RF)
    expb2 = fdiv(fone, expb1, prec2, RF)
    ch = fmul(fadd(expb1, expb2, prec2, RF), fhalf, prec2, RF)
    sh = fmul(fsub(expb1, expb2, prec2, RF), fhalf, prec2, RF)
    return fmul(c, ch, prec, rounding), \
        fneg(fmul(s, sh, prec, rounding), prec, rounding)

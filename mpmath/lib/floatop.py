"""
This module implements fundamental operations for floating-point
arithmetic (bit fiddling, normalization, arithmetic, ...).
"""

import math
from bisect import bisect
from random import randrange

# Regular (i.e. not inf/nan) numbers are represented in "raw" form as tuples
# (mantissa, exponent, bitcount) [short: man, exp, bc]
fzero = (0, 0, 0)
fone = (1, 0, 1)
ftwo = (1, 1, 1)
fthree = (3, 0, 2)
ften = (5, 1, 3)
fhalf = (1, -1, 1)

# Special numbers can be identified by looking for bc == -1. The choice of
# representation is entirely arbitrary, but fast. A number can be unpacked
# into its components as usual; only a single extra integer comparison needs
# to be expended to determine whether the number is regular.
finf = ('+inf', 0, -1)
fninf = ('-inf', 0, -1)
fnan = ('nan', 0, -1)


#----------------------------------------------------------------------------#
#                                                                            #
#                            Integer bit-level utilities                     #
#                                                                            #
#----------------------------------------------------------------------------#

def trailing(n):
    """Count the number of trailing zero bits in abs(n)."""
    if not n:
        return 0
    t = 0
    while not n & 1:
        n >>= 1
        t += 1
    return t

def bitcount(n, powers=[1<<_ for _ in range(300)]):
    """Calculate the bit size of abs(n)."""
    n = abs(n)
    bc = bisect(powers, n)
    if bc != 300:
        return bc
    bc = int(math.log(n, 2)) - 4
    return bc + bitcount(n>>bc)

# To the extent possible, we will try to avoid the above, slow functions
# and instead look up correction terms from precomputed tables
trailtable = map(trailing, range(256))
bctable = map(bitcount, range(1024))

# Rounding functions are variations of the shift-right operator, i.e., they
# calculate x >> n, but round differently when the discarded part is nonzero.

def round_floor(x, n):
    return x >> n

def round_ceiling(x, n):
    return -((-x) >> n)

def round_down(x, n):
    if x > 0: return x >> n
    else:     return -((-x) >> n)

def round_up(x, n):
    if x > 0: return -((-x) >> n)
    else:     return x >> n

def round_half_up(x, n):
    positive = x > 0
    if positive: t = x >> (n-1)
    else:        t = (-x) >> (n-1)
    if t & 1:
        if positive: return (t>>1)+1
        else:        return -((t>>1)+1)
    if positive: return t>>1
    else:        return -(t>>1)

def round_half_down(x, n):
    positive = x > 0
    if positive: t = x >> (n-1)
    else:        t = (-x) >> (n-1)
    if t & 1 and x & ((1<<(n-1))-1):
        if positive: return (t>>1)+1
        else:        return -((t>>1)+1)
    if positive: return t>>1
    else:        return -(t>>1)

# Half-even rounding is the default rounding mode, and also the slowest due
# to its complexity. Storing small bit masks in a table pays off.

def round_half_even(x, n, masks=[0]+[((1<<(_-1))-1) for _ in range(1, 300)]):
    if n < 300: mask = masks[n]
    else:       mask = ((1<<(n-1))-1)
    if x > 0:
        t = x >> (n-1)
        if t&1 and (t&2 or x&mask):
            return (t>>1)+1
        return t>>1
    else:
        t = (-x) >> (n-1)
        if t&1 and (t&2 or x&mask):
            return -((t>>1)+1)
        return -(t>>1)


#----------------------------------------------------------------------------#
#                                                                            #
#             Basic functions for number conversion / construction           #
#                                                                            #
#----------------------------------------------------------------------------#

# normalize() is the workhorse function in mpmath. All floating-point
# operations are implemented according to the pattern
#
#  1) convert floating-point problem to an equivalent integer problem
#  2) solve integer problem using Python integer arithmetic
#  3) use normalize() to convert the integer solution to a canonical
#     floating-point number
#
# A number of hacks are used in normalize() to reduce the overhead of step
# (3) as far as possible.

def normalize(man, exp, bc, prec, rounding):
    """
    normalize(man, exp, bc, prec, rounding) -> return tuple representing
    a fully rounded and normalized raw mpf with value (man * 2**exp)

    The mantissa is rounded in the specified direction if the number of
    bits exceeds the precision. Trailing zero bits are also stripped
    from the mantissa to ensure that the representation is canonical.

    Note: bc *must* be exact. If the bc is not known, use from_man_exp
    or fpos to normalize a raw mpf.
    """
    if not man:
        return fzero
    delta = bc - prec
    if bc > prec:
        man = rounding(man, delta)
        exp += delta
        bc = prec
    if not man & 1:
        while not man & 0xff:
            man >>= 8; exp += 8; bc -= 8
        t = trailtable[man & 0xff]
        man >>= t; exp += t; bc -= t
    if man == 1 or man == -1:
        bc = 1
    return man, exp, bc

def from_man_exp(man, exp, prec=None, rounding=None):
    """Create raw mpf from (man, exp) pair. If no precision is specified,
    the mantissa is stored exactly."""
    if prec is None:
        t = trailing(man)
        man >>= t
        return (man, exp+t, bitcount(man))
    return normalize(man, exp, bitcount(man), prec, rounding)

def to_man_exp(s):
    """Return (man, exp) of a raw mpf. Raise an error if inf/nan."""
    man, exp, bc = s
    if bc == -1:
        raise ValueError("mantissa and exponent are undefined for %s" % man)
    return man, exp

def from_int(n, prec=None, rounding=None):
    """Create a raw mpf from an integer. If no precision is specified,
    the mantissa is stored exactly."""
    return from_man_exp(n, 0, prec, rounding)

def to_int(s, rounding=round_down):
    """Convert a raw mpf to the nearest int. Rounding is done down by
    default (same as int(float) in Python), but can be changed. If the
    input is inf/nan, an exception is raised."""
    man, exp, bc = s
    if bc == -1:
        raise ValueError("cannot convert %s to int" % man)
    if exp > 0:
        return man << exp
    return rounding(man, -exp)

def fceil(s, prec, rounding):
    """Calculate ceil of a raw mpf, and round the result in the given
    direction (not necessarily ceiling). Note: returns a raw mpf
    representing an integer, not a Python int."""
    sman, sexp, sbc = s
    if sbc == -1:
        return s
    if sexp > 0:
        return fpos(s, prec, rounding)
    return from_int(to_int(s, round_ceiling), prec, rounding)

def ffloor(s, prec, rounding):
    """Calculate floor of a raw mpf, and round the result in the given
    direction (not necessarily floor). Note: returns a raw mpf
    representing an integer, not a Python int."""
    sman, sexp, sbc = s
    if sbc == -1:
        return s
    if sexp > 0:
        return fpos(s, prec, rounding)
    return from_int(to_int(s, round_floor), prec, rounding)

def from_float(x, prec=None, rounding=None):
    """Create a raw mpf from a Python float, rounding if necessary.
    If prec >= 53, the result is guaranteed to represent exactly the
    same number as the input. If prec is not specified, use prec=53."""
    if prec is None:
        prec, rounding = 53, round_floor
    # frexp only raises an exception for nan on some platforms
    if x != x:
        return fnan
    try:
        m, e = math.frexp(x)
    except:
        if x == 1e1000: return finf
        if x == -1e1000: return fninf
        return fnan
    return from_man_exp(int(m*(1<<53)), e-53, prec, rounding)

def to_float(s):
    """Convert a raw mpf to a Python float. The result is exact if the
    bitcount of s is <= 53 and no underflow/overflow occurs. An
    OverflowError is raised if the number is too large to be
    represented as a regular float."""
    man, exp, bc = s
    if bc == -1:
        if s == finf: return 1e1000
        if s == fninf: return -1e1000
        return 1e1000/1e1000
    if bc < 100:
        return math.ldexp(man, exp)
    # Try resizing the mantissa. Overflow may still happen here.
    n = bc - 53
    m = man >> n
    return math.ldexp(m, exp + n)

def from_rational(p, q, prec, rounding):
    """Create a raw mpf from a rational number p/q, rounding if
    necessary."""
    return fdiv(from_int(p), from_int(q), prec, rounding)

def to_rational(s):
    """Convert a raw mpf to a rational number. Return integers (p, q)
    such that s = p/q exactly."""
    man, exp, bc = s
    if bc == -1:
        raise ValueError("cannot convert %s to a rational number" % man)
    if exp >= 0:
        return man * (1<<exp), 1
    else:
        return man, 1<<(-exp)

def frand(prec):
    """Return a raw mpf chosen randomly from [0, 1), with prec bits
    in the mantissa."""
    return from_man_exp(randrange(0, 1<<prec), -prec, prec, round_floor)


#----------------------------------------------------------------------------#
#                                                                            #
#                       Arithmetic operations, etc.                          #
#                                                                            #
#----------------------------------------------------------------------------#

def feq(s, t):
    """Test equality of two raw mpfs. This is simply tuple comparion
    unless either number is nan, in which case the result is False."""
    if fnan in (s, t):
        return False
    return s == t

def fhash(s):
    man, exp, bc = s
    try:
        # Try to be compatible with hash values for floats and ints
        return hash(to_float(s))
    except OverflowError:
        # We must unfortunately sacrifice compatibility with ints here. We
        # could do hash(man << exp) when the exponent is positive, but
        # this would cause unreasonable inefficiency for large numbers.
        return hash(s)

def fcmp(s, t):
    """Compare the raw mpfs s and t. Return -1 if s < t, 0 if s == t,
    and 1 if s > t. (Same convention as Python's cmp() function.)"""

    # In principle, a comparison amounts to determining the sign of s-t.
    # A full subtraction is relatively slow, however, so we first try to
    # look at the components.
    sman, sexp, sbc = s
    tman, texp, tbc = t

    # Handle special numbers
    if sbc == -1 or tbc == -1:
        if s is t: return 0
        # Follow same convention as Python's cmp for float nan
        if t is fnan: return 1
        if s is finf: return 1
        return -1

    # Very easy cases: check for zeros and opposite signs
    if not tman: return cmp(sman, 0)
    if not sman: return cmp(0, tman)
    if sman > 0 and tman < 0: return 1
    if sman < 0 and tman > 0: return -1

    # This reduces to direct integer comparison
    if sexp == texp: return cmp(sman, tman)

    # Check position of the highest set bit in each number. If
    # different, there is certainly an inequality.
    a = sbc + sexp
    b = tbc + texp
    if sman > 0:
        if a < b: return -1
        if a > b: return 1
    else:
        if a < b: return 1
        if a > b: return -1

    # Both numbers have the same highest bit. Subtract to find
    # how the lower bits compare.
    return cmp(fsub(s, t, 5, round_floor)[0], 0)

def fpos(s, prec, rounding):
    """Calculate 0+s for a raw mpf (i.e., just round s to the specified
    precision)."""
    man, exp, bc = s
    if bc == -1:
        return s
    return normalize(man, exp, bc, prec, rounding)

def fneg(s, prec=None, rounding=None):
    """Negate a raw mpf (return -s), rounding the result to the
    specified precision. The prec argument can be omitted to do the
    operation exactly."""
    man, exp, bc = s
    if bc == -1:
        if s is finf: return fninf
        if s is fninf: return finf
        return fnan
    if prec is None:
        return (-man, exp, bc)
    return normalize(-man, exp, bc, prec, rounding)

def fabs(s, prec=None, rounding=None):
    """Return abs(s) of the raw mpf s, rounded to the specified
    precision. The prec argument can be omitted to generate an
    exact result."""
    man, exp, bc = s
    if bc == -1:
        if s is fninf: return finf
        return s
    if prec is None:
        if man < 0:
            return (-man, exp, bc)
        return s
    if man < 0:
        return normalize(-man, exp, bc, prec, rounding)
    return normalize(man, exp, bc, prec, rounding)

def fsign(s):
    """Return -1, 0, or 1 (as a Python int, not a raw mpf) depending on
    whether s is negative, zero, or positive. (Nan is taken to give 0.)"""
    man, exp, bc = s
    if bc == -1:
        if s is finf: return 1
        if s is fninf: return -1
        return 0
    return cmp(man, 0)

def fadd(s, t, prec, rounding):
    # We will assume below that s has the higher exponent. If not, we swap
    # them before unpacking::
    if t[1] > s[1]:
        s, t = t, s
    sman, sexp, sbc = s
    tman, texp, tbc = t

    if sbc == -1 or tbc == -1:
        either = s, t
        if fnan in either: return fnan
        if finf in either and fninf in either: return fnan
        if finf in either: return finf
        return fninf

    # Check if one operand is zero. Zero always has exp = 0; if the
    # other operand has a huge exponent, its mantissa will unnecessarily
    # be shifted into something huge if we don't check for this case.
    if not tman: return normalize(sman, sexp, sbc, prec, rounding)
    if not sman: return normalize(tman, texp, tbc, prec, rounding)

    # More generally, if one number is huge and the other is small,
    # and in particular, if their mantissas don't overlap at all at
    # the current precision level, we can avoid work.
    #         precision
    #      |            |
    #       111111111
    #    +                    222222222
    #       ------------------------
    #       1111111110000... (222)
    delta = sbc + sexp - tbc - texp
    if delta > prec + 4:
        # The result may have to be rounded up or down. So we shift s
        # and add a dummy bit outside the precision range to force
        # rounding.
        offset = min(delta, prec) + 4
        sman <<= offset
        if tman > 0:
            sman += 1
        else:
            sman -= 1
        # TODO: use that bc ~= sbc+offset
        bc = bitcount(sman)
        return normalize(sman, sexp-offset, bc, prec, rounding)

    # General algorithm: we set min(sexp, texp) = 0, perform exact
    # integer addition, and then round the result.
    #                 exp = 0
    #                     v
    #        11111111100000   <-- sman (padded with zeros from shifting)
    #    +        222222222   <-- tman (no shifting necessary)
    #        --------------
    #    =   11111333333333
    offset = sexp - texp
    man = tman + (sman<<offset)

    # At this point we are almost done; what remains is to determine the
    # bitcount of the mantissa.

    # If signs are equal the bitcount is given approximately by the inputs'
    # bitcounts (a quick correction lookup gives the exact value). Todo:
    # this is also doable when signs differ if the high bits have different
    # positions.
    sbc += offset
    if (sman<0) == (tman<0):
        if tbc > sbc: bc = tbc - 4
        else:         bc = sbc - 4
        if bc < 4:    bc = bctable[abs(man)]
        else:         bc += bctable[abs(man)>>bc]
    else:
        # Subtraction might have cancelled many bits; slow case
        bc = bitcount(man)

    return normalize(man, texp, bc, prec, rounding)

def fsub(s, t, prec, rounding):
    """Return the difference of two raw mpfs, s-t. This function is
    simply a wrapper of fadd that changes the sign of t."""
    man, exp, bc = t
    if bc == -1:
        return fadd(s, fneg(t, prec, rounding), prec, rounding)
    return fadd(s, (-man, exp, bc), prec, rounding)

def fmul(s, t, prec, rounding):
    sman, sexp, sbc = s
    tman, texp, tbc = t
    if sbc == -1 or tbc == -1:
        if fnan in (s, t): return fnan
        if tbc == -1: s, t = t, s
        if t == fzero: return fnan
        return {1:finf, -1:fninf}[fsign(s) * fsign(t)]
    # Compared to addition, multiplication is a piece of cake!
    man = sman*tman
    # Happily, the bitcount is extremely simple to estimate for
    # multiplication, allowing us to skip an expensive call to
    # bitcount(). Note: the bitcount actually becomes wrong when sbc
    # or tbc is 0, but that is accidentally ok since normalize checks
    # for a zero mantissa and handles that specially.
    bc = sbc + tbc - 4
    if bc < 4: bc = bctable[abs(man)]
    else:      bc += bctable[abs(man)>>bc]
    return normalize(man, sexp+texp, bc, prec, rounding)

def fshift_exact(s, n):
    """Quickly multiply the raw mpf s by 2**n without rounding."""
    man, exp, bc = s
    if bc == -1:
        return s
    if not man:
        return s
    return man, exp+n, bc

def fdiv(s, t, prec, rounding, bct=bctable):
    """Floating-point division"""
    sman, sexp, sbc = s
    tman, texp, tbc = t
    if sbc == -1 or tbc == -1:
        if fnan in (s, t): return fnan
        if sbc == tbc == -1: return fnan
        if tbc != -1:
            if t == fzero:
                return fnan
            return {1:finf, -1:fninf}[fsign(s) * fsign(t)]
        return fzero
    if not tman:
        return fnan
    # Same strategy as for addition: if there is a remainder, perturb
    # the result a few bits outside the precision range before rounding
    extra = prec-sbc+tbc+5
    if extra < 5:
        extra = 5
    quot, rem = divmod(sman<<extra, tman)
    if rem:
        quot = (quot << 5) + 1
        extra += 5
    bc = sbc+extra-tbc-4
    if bc < 4: bc = bctable[abs(quot)]
    else:      bc += bctable[abs(quot)>>bc]
    return normalize(quot, sexp-texp-extra, bc, prec, rounding)

def fmod(s, t, prec, rounding):
    sman, sexp, sbc = s
    tman, texp, tbc = t
    if sbc == -1 or tbc == -1:
        return fnan
    # Important special case: do nothing if t is larger
    if ((sman < 0) == (tman < 0)) and texp > sexp+sbc:
        return s
    # Another important special case: this allows us to do e.g. x % 1.0
    # to find the fractional part of x, and it'll work when x is huge.
    if tman == 1 and sexp > texp+tbc:
        return fzero
    base = min(sexp, texp)
    man = (sman << (sexp-base)) % (tman << (texp-base))
    return normalize(man, base, bitcount(man), prec, rounding)

reciprocal_rounding = {
  round_down : round_up,
  round_up : round_down,
  round_floor : round_ceiling,
  round_ceiling : round_floor,
  round_half_down : round_half_up,
  round_half_up : round_half_down,
  round_half_even : round_half_even
}

negative_rounding = {
  round_down : round_down,
  round_up : round_up,
  round_floor : round_ceiling,
  round_ceiling : round_floor,
  round_half_down : round_half_down,
  round_half_up : round_half_up,
  round_half_even : round_half_even
}

def fpow(s, n, prec, rounding):
    """Compute s**n, where n is an integer"""
    if s[2] == -1:
        if s == finf:
            if n > 0: return s
            if n == 0: return fnan
            return fzero
        if s == fninf:
            if n > 0: return [finf, fninf][n & 1]
            if n == 0: return fnan
            return fzero
        return fnan

    n = int(n)
    if n == 0: return fone
    if n == 1: return fpos(s, prec, rounding)
    if n == 2: return fmul(s, s, prec, rounding)
    if n == -1: return fdiv(fone, s, prec, rounding)
    if n < 0:
        inverse = fpow(s, -n, prec+5, reciprocal_rounding[rounding])
        return fdiv(fone, inverse, prec, rounding)

    man, exp, bc = s

    # Use exact integer power when the exact mantissa is small
    if bc*n < 5000 or abs(man) == 1:
        return from_man_exp(man**n, exp*n, prec, rounding)

    # Use directed rounding all the way through to maintain rigorous
    # bounds for interval arithmetic
    sign = 1
    rounding2 = rounding
    if man < 0 and n % 2:
        sign = -1
        rounding2 = negative_rounding[rounding]
    man = abs(man)

    # Now we perform binary exponentiation. Need to estimate precision
    # to avoid rounding from temporary operations. Roughly log_2(n)
    # operations are performed.
    prec2 = prec + 4*bitcount(n) + 4
    pm, pe, pbc = fone
    while n:
        if n & 1:
            pm, pe, pbc = from_man_exp(pm*man, pe+exp, prec2, rounding2)
            n -= 1
        man, exp, bc = from_man_exp(man*man, exp+exp, prec2, rounding2)
        n = n // 2

    return from_man_exp(sign*pm, pe, prec, rounding)

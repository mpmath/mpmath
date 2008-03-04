__docformat__ = 'plaintext'

import math
from bisect import bisect
from random import randrange

# Same as standard Python float
STANDARD_PREC = 53

LOG2_10 = math.log(10,2)  # 3.3219...

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

def make_fixed(s, prec):
    """Convert a floating-point number to a fixed-point big integer"""
    man, exp, bc = s
    offset = exp + prec
    if offset >= 0:
        return man << offset
    else:
        return man >> (-offset)

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
    if exp >= 0:
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

def fmuli(s, n, prec, rounding):
    """Multiply by a Python integer."""
    man, exp, bc = s
    if bc == -1:
        return fmul(s, from_int(n), prec, rounding)
    if not n:
        return fzero
    man *= n
    # Generally n will be small
    try:
        bc += bctable[abs(n)] - 4
    except:
        bc += bitcount(n) - 4
    if bc < 4: bc = bctable[abs(man)]
    else:      bc += bctable[abs(man)>>bc]
    return normalize(man, exp, bc, prec, rounding)

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

def fpow(s, t, prec, rounding):
    """Compute s**t. Raise ValueError if s is negative and t is
    fractional."""
    sman, sexp, sbc = s
    tman, texp, tbc = t
    if sman < 0 and texp < 0:
        raise ValueError
    if texp >= 0:
        return fpowi(s, tman << texp, prec, rounding)
    # s**(n/2) = sqrt(s)**n
    if texp == -1:
        if tman == 1:
            return fsqrt(s, prec, rounding)
        return fpowi(fsqrt(s, prec+10, rounding), tman, prec, rounding)
    # General formula: s**t = exp(t*log(s))
    # TODO: handle rounding direction of the logarithm carefully
    c = flog(s, prec+10, rounding)
    return fexp(fmul(t, c, prec+10, rounding), prec, rounding)

def fpowi(s, n, prec, rounding):
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
        inverse = fpowi(s, -n, prec+5, reciprocal_rounding[rounding])
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
    while 1:
        if n & 1:
            pm, pe, pbc = from_man_exp(pm*man, pe+exp, prec2, rounding2)
            n -= 1
            if not n:
                break
        man, exp, bc = from_man_exp(man*man, exp+exp, prec2, rounding2)
        n = n // 2

    return from_man_exp(sign*pm, pe, prec, rounding)



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
        s = fdiv(s, fpowi(ften, b, bitprec, RF), bitprec, RF)
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
        s = fmul(s, fpowi(ften, exp, prec+10, round_floor), prec, rounding)
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


def sqrt_initial(y, prec):
    """Given y = floor(x * 2**prec), compute floor(sqrt(x) * 2**50),
    i.e. calculate a 50-bit approximation of the square root. This is
    done quickly using regular floating-point arithmetic. It is
    assumed that x ~= 1."""

    # Two cases; the second avoids overflow
    if prec < 200: return int(y**0.5 * 2.0**(50 - prec*0.5))
    else:          return int((y >> (prec-100))**0.5)


# XXX: doesn't work
def invsqrt_initial(y, prec):
    """Like sqrt_initial, but computes 1/sqrt(y) instead of sqrt(y)."""
    if prec < 200: return int(y**-0.5 * 2.0**(50 + prec*0.5))
    else:          return int((y >> (prec-100)) ** -0.5)


def sqrt_fixed(y, prec):
    """
    Square root of a fixed-point number.

    Given the big integer y = floor(x * 2**prec), this function returns
    floor(r * 2**prec) where r = sqrt(x).

    We start from a 50-bit estimate for r generated with ordinary
    floating-point arithmetic, and then refines the value to full
    accuracy using the iteration

                 1  /        y  \
        r     = --- | r  +  --- |
         n+1     2  \ n     r_n /

    which is simply Newton's method applied to the equation r**2 = y.

    Newton's method doubles the accuracy with each step. We make use of
    this fact by only using a precision corresponding to the current
    accuracy during intermediate iterations. For example, with a 50-bit
    accurate r_1, r_2 can be computed using 100-bit precision, r_3
    using 200-bit precision, and so on. (In practice, the precision
    levels must be chosen slightly more conservatively to account for
    rounding errors in the last one or two bits.)

    It is assumed that x ~= 1 (the main fsqrt() function fiddles with
    the exponent of the input to reduce it to unit magnitude before
    passing it here.)

    TODO: it would be possible to specify separate precision levels
    for the input and output, which could be useful when calculating
    pure-integer square roots.
    """

    r = sqrt_initial(y, prec)
    extra = 10
    prevp = 50

    for p in giant_steps(50, prec+extra):

        # Explanation: in the first term, we shift by the appropriate number
        # of bits to convert r from the previous precision to the current one.
        # The "-1" divides by two in the same step.

        # In the second term, we do a fixed-point division using the usual
        # formula (y<<precision)//r. The precision term is a bit messy and
        # takes into account the fact that y, r_n and r_{n+1} all have
        # different precision levels. As before, the "-1" divides by two.
        r = lshift_quick(r, p-prevp-1) + (lshift_quick(y, p+prevp-prec-1)//r)

        prevp = p

    return r >> extra


def sqrt_fixed2(y, prec):
    """
    This function is essentially equivalent to sqrt_fixed (see its
    documentation), but uses an asymptotically faster algorithm.

    Instead of using Newton's method to calculate sqrt(y) directly,
    we calculate 1/sqrt(y) with Newton's method and multiply by y to
    obtain sqrt(y). The Newton iteration for 1/sqrt(y) is

                 r
                  n      /            2 \
        r    =  ----  *  | 3  - y * r   |.
         n+1      2      \           n  /

    This is slightly slower at low precision levels since it requires
    three multiplications in each step, as opposed to the single
    division in the Newton iteration for sqrt(y).

    However, since Python uses Karatsuba algorithm for multiplication,
    three multiplications can be performed much more quickly than a
    single division at high precision levels. In practice, the cutoff
    where sqrt_fixed2 becomes faster than sqrt_fixed seems to be around
    60,000 bits.
    """

    # XXX
    r = to_float(from_man_exp(y, -prec, 64, round_floor)) ** -0.5
    r = int(r * 2**50)

    # r = invsqrt_initial(y, prec)

    extra = 10
    prevp = 50

    for p in giant_steps(50, prec+extra):

        # This is even messier than in sqrt_fixed. As many shifts as possible
        # have been combined together for optimal speed, at a slight expense
        # of legibility.

        # Compute r**2 at precision p.
        r2 = rshift_quick(r*r, 2*prevp - p)

        # A = r, converted from precision prevp to p
        A = lshift_quick(r, p-prevp)

        # S = y * r2, computed at precision p. We shift y by '-prec' to
        # account for its initial precision, and by 'p' for the fixed-point
        # multiplication
        S = (lshift_quick(y, p-prec) * r2) >> p

        # B = (3-S) and finally the outer product, both done at precision p
        B = (3<<p) - S
        r = (A*B) >> (p+1)

        prevp = p

    # Finally, sqrt(y) = y * (1/sqrt(y))
    r = (r * y) >> prec

    return r >> extra




def fsqrt(s, prec, rounding):
    """
    Floating-point square root.

    Returns a tuple representing the square root of s, rounded to the
    nearest floating-point number in the specified rounding direction.
    The input must be a tuple representing a nonnegative floating-point
    number.
    """
    if s == fone:
        return fone

    man, exp, bc = s

    if not man:
        return fzero

    # Convert to a fixed-point number with prec2 bits. Adjust
    # exponents to be even so that they can be divided in half
    prec2 = prec + 12 + (prec & 1)

    if exp & 1:
        exp -= 1
        man <<= 1
        bc += 1

    # Mantissa may have more bits than we need. Trim it down.
    shift = bc - prec2
    shift -= shift & 1
    man = rshift_quick(man, shift)

    if prec < 65000:
        man = sqrt_fixed(man, prec2)
    else:
        man = sqrt_fixed2(man, prec2)

    return from_man_exp(man, (exp+shift-prec2)>>1, prec, rounding)



def fhypot(x, y, prec, rounding):
    if y == fzero:
        return fabs(x, prec, rounding)
    if x == fzero:
        return fabs(y, prec, rounding)

    RF = round_floor
    hypot2 = fadd(fmul(x,x,prec+4,RF), fmul(y,y,prec+4,RF), prec+4, RF)

    return fsqrt(hypot2, prec, rounding)


def constant_memo(f):
    """Cache computed values of mathematical constants"""
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


"""Provides real valued functons.

Transcendental functions for real numbers:
* exp
* log
* sin/cos/tan
* sinh/cosh/tanh

"""



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

def fexp(x, prec, rounding):
    man, exp, bc = x
    if bc == -1:
        if x == fninf:
            return fzero
        return x
    # extra precision needs to be similar in magnitude to log_2(|x|)
    prec2 = prec + 6 + max(0, bc+exp)
    t = make_fixed(x, prec2)
    # abs(x) > 1?
    if exp+bc > 1:  #fcmp(fabs(x), fone) > 0:
        lg2 = log2_fixed(prec2)
        n, t = divmod(t, lg2)
    else:
        n = 0
    return from_man_exp(exp_series(t, prec2), -prec2+n, prec, rounding)


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
    fx = math.log(to_float((x, -prec, bitcount(x))))
    r = int(fx * 2.0**50)
    prevp = 50
    for p in giant_steps(50, prec+extra):
        rb = lshift_quick(r, p-prevp)
        e = exp_series(-rb, p)
        r = rb + ((rshift_quick(x, prec-p)*e)>>p) - (1 << p)
        prevp = p
    return r >> extra

def flog(x, prec, rounding):
    if x == fzero: return fnan
    if x == fone:  return fzero
    man, exp, bc = x
    if bc == -1:
        if x == finf: return finf
        return fnan
    if man < 0: return fnan
    # Estimated precision needed for log(t) + n*log(2)
    prec2 = prec + int(math.log(1+abs(bc+exp), 2)) + 10
    # Watch out for the case when x is very close to 1
    if -1 < bc + exp < 2:
        near_one = fabs(fsub(x, fone, STANDARD_PREC, round_floor), STANDARD_PREC, round_floor)
        if near_one == 0:
            return fzero
        # estimate how close
        prec2 += -(near_one[1]) - bitcount(near_one[0])
    # Separate mantissa and exponent, calculate, join parts
    t = rshift_quick(man, bc-prec2)
    l = _log_newton(t, prec2)
    a = (exp + bc) * log2_fixed(prec2)
    return from_man_exp(l+a, -prec2, prec, rounding)



"""
We compute sin(x) around 0 from its Taylor series, and cos(x) around 0
from sqrt(1-sin(x)**2). This way we can simultaneously compute sin and
cos, which are often needed together (e.g. for the tangent function or
the complex exponential), with little extra cost compared to computing
just one of them. The main reason for computing sin first (and not sin
from cos) is to obtain high relative accuracy for x extremely close to
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


TODO: could use cos series too when extremely close to 0
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

def cos_sin(x, prec, rounding):
    """Simultaneously compute (cos(x), sin(x)) for real x."""
    man, exp, bc = x
    if bc == -1:
        return (fnan, fnan)
    bits_from_unit = abs(bc + exp)
    prec2 = prec + bits_from_unit + 15
    xf = make_fixed(x, prec2)
    n, rx = _trig_reduce(xf, prec2)
    case = n % 4
    one = 1 << prec2
    if case == 0:
        s = _sin_series(rx, prec2)
        c = sqrt_fixed(one - ((s*s)>>prec2), prec2)
    elif case == 1:
        c = -_sin_series(rx, prec2)
        s = sqrt_fixed(one - ((c*c)>>prec2), prec2)
    elif case == 2:
        s = -_sin_series(rx, prec2)
        c = -sqrt_fixed(one - ((s*s)>>prec2), prec2)
    elif case == 3:
        c = _sin_series(rx, prec2)
        s = -sqrt_fixed(one - ((c*c)>>prec2), prec2)
    c = from_man_exp(c, -prec2, prec, rounding)
    s = from_man_exp(s, -prec2, prec, rounding)
    return c, s

def fcos(x, prec, rounding):
    return cos_sin(x, prec, rounding)[0]

def fsin(x, prec, rounding):
    return cos_sin(x, prec, rounding)[1]

def ftan(x, prec, rounding):
    c, s = cos_sin(x, prec+6, round_floor)
    return fdiv(s, c, prec, rounding)


#----------------------------------------------------------------------
# Hyperbolic functions
#

def _sinh_series(x, prec):
    x2 = (x*x) >> prec
    s = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (k*(k-1))
        s += a
        k += 2
    return s

def cosh_sinh(x, prec, rounding):
    """Simultaneously compute (cosh(x), sinh(x)) for real x"""
    man, exp, bc = x
    if bc == -1:
        if x == finf: return (finf, finf)
        if x == fninf: return (finf, fninf)
        return fnan

    high_bit = exp + bc
    prec2 = prec + 6

    if high_bit < -3:
        # Extremely close to 0, sinh(x) ~= x and cosh(x) ~= 1
        # TODO: support directed rounding
        if high_bit < -prec-2:
            return (fone, fpos(x, prec, rounding))

        # Avoid cancellation when computing sinh
        # TODO: might be faster to use sinh series directly
        prec2 += (-high_bit) + 4

    # In the general case, we use
    #    cosh(x) = (exp(x) + exp(-x))/2
    #    sinh(x) = (exp(x) - exp(-x))/2
    # and note that the exponential only needs to be computed once.
    ep = fexp(x, prec2, round_floor)
    em = fdiv(fone, ep, prec2, round_floor)
    ch = fshift_exact(fadd(ep, em, prec, rounding), -1)
    sh = fshift_exact(fsub(ep, em, prec, rounding), -1)
    return ch, sh

def fcosh(x, prec, rounding):
    """Compute cosh(x) for a real argument x"""
    return cosh_sinh(x, prec, rounding)[0]

def fsinh(x, prec, rounding):
    """Compute sinh(x) for a real argument x"""
    return cosh_sinh(x, prec, rounding)[1]

def ftanh(x, prec, rounding):
    """Compute tanh(x) for a real argument x"""
    ch, sh = cosh_sinh(x, prec+6, round_floor)
    return fdiv(sh, ch, prec, rounding)


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

def _atan_series_1(x, prec, rounding):
    man, exp, bc = x
    # Increase absolute precision when extremely close to 0
    bc = bitcount(man)
    diff = -(bc + exp)
    prec2 = prec
    if diff > 10:
        if 3*diff - 4 > prec:  # x**3 term vanishes; atan(x) ~x
            return from_man_exp(man, exp, prec, rounding)
        prec2 = prec + diff
    prec2 += 15  # XXX: better estimate for number of guard bits
    x = make_fixed(x, prec2)
    x2 = (x*x)>>prec2; one = 1<<prec2; s=a=x
    for n in xrange(1, 1000000):
        a = (a*x2) >> prec2
        s += a // ((-1)**n * (n+n+1))
        if -100 < a < 100:
            break
    return from_man_exp(s, -prec2, prec, rounding)

def _atan_series_2(x, prec, rounding):
    prec2 = prec + 15
    x = make_fixed(x, prec2)
    one = 1<<prec2; x2 = (x*x)>>prec2; y=(x2<<prec2)//(one+x2)
    s = a = one
    for n in xrange(1, 1000000):
        a = ((a*y)>>prec2) * (2*n) // (2*n+1)
        if a < 100:
            break
        s += a
    return from_man_exp(y*s//x, -prec2, prec, rounding)

_cutoff_1 = (5, -3, 3)   # ~0.6
_cutoff_2 = (3, -1, 2)   # 1.5

def fatan(x, prec, rounding):
    man, exp, bc = x
    if bc == -1:
        if x == finf: return fshift_exact(fpi(prec, round_down), -1)
        if x == fninf: return fneg(fshift_exact(fpi(prec, round_down), -1))
        return fnan
    if man < 0:
        t = fatan(fneg(x), prec+4, round_floor)
        return from_man_exp(-t[0], t[1], prec, rounding)
    if fcmp(x, _cutoff_1) < 0:
        return _atan_series_1(x, prec, rounding)
    if fcmp(x, _cutoff_2) < 0:
        return _atan_series_2(x, prec, rounding)
    # For large x, use atan(x) = pi/2 - atan(1/x)
    if x[1] > 10*prec:
        pi = fpi(prec, rounding)
        pihalf = pi[0], pi[1]-1, pi[2]
    else:
        pi = fpi(prec+4, round_floor)
        pihalf = pi[0], pi[1]-1, pi[2]
        t = fatan(fdiv(fone, x, prec+4, round_floor), prec+4, round_floor)
        return fsub(pihalf, t, prec, rounding)


# Use fastest rounding mode for intermediate calculations
RF = round_floor

def fcabs(a, b, prec, rounding):
    """Absolute value of a complex number, |a+bi|. Returns a single
    real number."""
    return fhypot(a, b, prec, rounding)

def fcmul(a, b, c, d, prec, rounding):
    """Complex multiplication.

    Returns the real and imaginary part of (a+bi)*(c+di), rounded to
    the specified precision. The rounding mode applies to the real and
    imaginary parts separately."""
    # All-real case
    if b == d == fzero:
        return fmul(a, c, prec, rounding), fzero
    ep = prec + 10
    re = fsub(fmul(a,c, ep, RF), fmul(b,d, ep, RF), prec, rounding)
    im = fadd(fmul(a,d, ep, RF), fmul(b,c, ep, RF), prec, rounding)
    return re, im

def fcsqrt(a, b, prec, rounding):
    """Complex square root (principal branch).

    We have sqrt(a+bi) = sqrt((r+a)/2) + b/sqrt(2*(r+a))*i where
    r = abs(a+bi), when a+bi is not a negative real number."""
    if a == b == fzero:
        return (a, b)
    # When a+bi is a negative real number, we get a real sqrt times i
    if a[0] < 0 and b == fzero:
        im = fsqrt(fneg(a), prec, rounding)
        return (fzero, im)
    ep = prec+20
    t  = fadd(fcabs(a, b, ep, RF), a, ep, RF)  # t = abs(a+bi) + a
    u  = fmul(t, fhalf, ep, RF)                # u = t / 2
    re = fsqrt(u, prec, rounding)              # re = sqrt(u)
    v  = fmul(t, ftwo, ep, RF)                 # v = t * 2
    w  = fsqrt(v, ep, RF)                      # w = sqrt(v)
    im = fdiv(b, w, prec, rounding)            # im = b / w
    return re, im

def fcexp(a, b, prec, rounding):
    """
    Complex exponential function.

    We use the direct formula exp(a+bi) = exp(a) * (cos(b) + sin(b)*i)
    for the computation. This formula is very nice because it is
    perfectly stable; since we just do real multiplications, the only
    numerical errors that can creep in are single-ulp rounding errors.

    The formula is efficient since mpmath's real exp is quite fast and
    since we can compute cos and sin simultaneously.

    It is no problem if a and b are large; if the implementations of
    exp/cos/sin are accurate and efficient for all real numbers, then
    so is this function for all complex numbers.
    """
    mag = fexp(a, prec+4, rounding)
    c, s = cos_sin(b, prec+4, rounding)
    re = fmul(mag, c, prec, rounding)
    im = fmul(mag, s, prec, rounding)
    return re, im

def fccos(a, b, prec, rounding):
    """Complex cosine. The formula used is cos(a+bi) = cos(a)*cosh(b) -
    sin(a)*sinh(b)*i.

    The same comments apply as for the complex exp: only real
    multiplications are performed, so no cancellation errors are
    possible. The formula is also efficient since we can compute both
    pairs (cos, sin) and (cosh, sinh) in single steps."""
    ep = prec + 6
    c, s = cos_sin(a, ep, RF)
    ch, sh = cosh_sinh(b, ep, RF)
    re = fmul(c, ch, prec, rounding)
    im = fmul(s, sh, prec, rounding)
    return re, fneg(im)

def fcsin(a, b, prec, rounding):
    """Complex sine. We have sin(a+bi) = sin(a)*cosh(b) +
    cos(a)*sinh(b)*i. See the docstring for fccos for additional
    comments."""
    ep = prec + 6
    c, s = cos_sin(a, ep, RF)
    ch, sh = cosh_sinh(b, ep, RF)
    re = fmul(s, ch, prec, rounding)
    im = fmul(c, sh, prec, rounding)
    return re, im

def fccosh(a, b, prec, rounding):
    """Complex hyperbolic cosine. Computed as cosh(z) = cos(z*i)."""
    return fccos(b, fneg(a), prec, rounding)

def fcsinh(a, b, prec, rounding):
    """Complex hyperbolic sine. Computed as sinh(z) = -i*sin(z*i)."""
    b, a = fcsin(b, a, prec, rounding)
    return a, b

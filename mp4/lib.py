#import os
#os.environ['MPMATH_NOGMPY'] = '1'

from libintmath import isqrt, sqrtrem
from libelefun import (mpf_cos, mpf_sin, mpf_exp, mpf_log, mpf_pow_int,
    mpf_pow, mpf_atan2, cos_sin, mpf_nthroot, mpf_log_hypot)

from libelefun import ln2_fixed, pi_fixed, e_fixed, ln10_fixed, phi_fixed, degree_fixed
from gammazeta import euler_fixed, catalan_fixed, glaisher_fixed, apery_fixed

from libmpf import mpf_sqrt

from libmpf import MP_BASE as MPZ
from libmpf import MODE
from libmpf import bitcount, trailing, trailtable, from_man_exp
from libmpf import from_man_exp
from libintmath import numeral, powers, bctable, isqrt_fast

from math import ldexp as math_ldexp
from math import frexp as math_frexp
from math import log as math_log

from bisect import bisect

if MODE == 'python':
    powers = [1<<_ for _ in range(300)]
    def bitcount(n):
        if n < (1<<300):
            return bisect(powers, n)
        bc = int(math_log(n, 2)) - 4
        return bc + bctable[n>>bc]
    def bitcount2(n):
        if n < 0:
            n = -n
        if n < (1<<300):
            return bisect(powers, n)
        bc = int(math_log(n, 2)) - 4
        return bc + bctable[n>>bc]
else:
    bitcount2 = bitcount

bctable = map(bitcount, range(1024))

MPZ_0 = MPZ(0)
MPZ_1 = MPZ(1)

def i_trim(man, exp=0, prec=0, rounding='f'):
    """
    If *prec > 0* and *man* is nonzero, round *man* to a width of at most
    *prec* bits, in the specified rounding direction, and remove trailing
    zero bits from *man*. The total number of bits by which *man* has
    been right-shifted is added to *exp*, and *(man, exp)* is returned.

    If *prec=0*, trailing zero bits are removed from *man*, but no
    rounding is performed. This effectively factors *man* into an odd
    part and a power of two.

    If *man = 0*, the result is always *(0, 0)*.

    Factoring:

    >>> i_trim(10)
    (5, 1)
    >>> i_trim(10, 100)
    (5, 101)

    Basic directional rounding:

    >>> a = (3<<20) + 1
    >>> assert i_trim(a, 0, 5, 'n') == (3, 20)
    >>> assert i_trim(a, 0, 5, 'u') == (25, 17)
    >>> assert i_trim(a, 0, 5, 'd') == (3, 20)
    >>> assert i_trim(a, 0, 5, 'c') == (25, 17)
    >>> assert i_trim(a, 0, 5, 'f') == (3, 20)
    >>> assert i_trim(-a, 0, 5, 'n') == (-3, 20)
    >>> assert i_trim(-a, 0, 5, 'u') == (-25, 17)
    >>> assert i_trim(-a, 0, 5, 'd') == (-3, 20)
    >>> assert i_trim(-a, 0, 5, 'c') == (-3, 20)
    >>> assert i_trim(-a, 0, 5, 'f') == (-25, 17)

    Hopefully a comprehensive set of tests:

    >>> assert i_trim(0, -4, 4, 'd') == (0, 0)
    >>> assert i_trim(0xf0, -4, 4, 'd') == (15, 0)
    >>> assert i_trim(0xf1, -4, 4, 'd') == (15, 0)
    >>> assert i_trim(0xff, -4, 4, 'd') == (15, 0)
    >>> assert i_trim(-0xf0, -4, 4, 'd') == (-15, 0)
    >>> assert i_trim(-0xf1, -4, 4, 'd') == (-15, 0)
    >>> assert i_trim(-0xff, -4, 4, 'd') == (-15, 0)
    >>> assert i_trim(0, -4, 4, 'u') == (0, 0)
    >>> assert i_trim(0xf0, -4, 4, 'u') == (15, 0)
    >>> assert i_trim(0xf1, -4, 4, 'u') == (1, 4)
    >>> assert i_trim(0xff, -4, 4, 'u') == (1, 4)
    >>> assert i_trim(-0xf0, -4, 4, 'u') == (-15, 0)
    >>> assert i_trim(-0xf1, -4, 4, 'u') == (-1, 4)
    >>> assert i_trim(-0xff, -4, 4, 'u') == (-1, 4)
    >>> assert i_trim(0, -4, 4, 'f') == (0, 0)
    >>> assert i_trim(0xf0, -4, 4, 'f') == (15, 0)
    >>> assert i_trim(0xf1, -4, 4, 'f') == (15, 0)
    >>> assert i_trim(0xff, -4, 4, 'f') == (15, 0)
    >>> assert i_trim(-0xf0, -4, 4, 'f') == (-15, 0)
    >>> assert i_trim(-0xf1, -4, 4, 'f') == (-1, 4)
    >>> assert i_trim(-0xff, -4, 4, 'f') == (-1, 4)
    >>> assert i_trim(0, -4, 4, 'c') == (0, 0)
    >>> assert i_trim(0xf0, -4, 4, 'c') == (15, 0)
    >>> assert i_trim(0xf1, -4, 4, 'c') == (1, 4)
    >>> assert i_trim(0xff, -4, 4, 'c') == (1, 4)
    >>> assert i_trim(-0xf0, -4, 4, 'c') == (-15, 0)
    >>> assert i_trim(-0xf1, -4, 4, 'c') == (-15, 0)
    >>> assert i_trim(-0xff, -4, 4, 'c') == (-15, 0)
    >>> assert i_trim(0, -4, 4, 'n') == (0, 0)
    >>> assert i_trim(0xf0, -4, 4, 'n') == (15, 0)
    >>> assert i_trim(0xf7, -4, 4, 'n') == (15, 0)
    >>> assert i_trim(0xf8, -4, 4, 'n') == (1, 4)    # 1111.1000 -> 10000.0
    >>> assert i_trim(0xf9, -4, 4, 'n') == (1, 4)    # 1111.1001 -> 10000.0
    >>> assert i_trim(0xe8, -4, 4, 'n') == (7, 1)    # 1110.1000 -> 1110.0
    >>> assert i_trim(0xe9, -4, 4, 'n') == (15, 0)     # 1110.1001 -> 1111.0
    >>> assert i_trim(-0xf0, -4, 4, 'n') == (-15, 0)
    >>> assert i_trim(-0xf7, -4, 4, 'n') == (-15, 0)
    >>> assert i_trim(-0xf8, -4, 4, 'n') == (-1, 4)
    >>> assert i_trim(-0xf9, -4, 4, 'n') == (-1, 4)
    >>> assert i_trim(-0xe8, -4, 4, 'n') == (-7, 1)
    >>> assert i_trim(-0xe9, -4, 4, 'n') == (-15, 0)
    >>> assert i_trim(72057594037927935, -56, 53, 'u') == (1, 0)
    >>> assert i_trim(73786976294838205979, -65, 53, 'n') == (1, 1)
    >>> assert i_trim(-73786976294838205979, -65, 53, 'n') == (-1, 1)
    >>> assert i_trim(31, 0, 4, 'u') == (1, 5)
    >>> assert i_trim(-31, 0, 4, 'f') == (-1, 5)
    >>> assert i_trim(255, 0, 7, 'u') == (1, 8)
    >>> assert i_trim(-255, 0, 7, 'f') == (-1, 8)

    """
    # positive
    if man > 0:

        # Round only if necessary
        if prec and man >> prec:

            # Count bits inline
            # n = bitcount(man) - prec
            if man < (1<<300):
                bc = bisect(powers, man)
            else:
                bc = int(math_log(man, 2)) - 4
                bc += bctable[man>>bc]
            n = bc - prec
            # End inline bitcount

            # Round
            if rounding == 'n':
                n1 = n-1
                t = man >> n1
                if t & 1 and ((t & 2) or (man & ((1<<n1)-1))):
                    man = (t >> 1) + 1
                else:
                    man = t >> 1
            elif rounding in 'fd':
                man >>= n
            else:
                man = -((-man) >> n)
            exp += n

        # Remove any trailing zeros
        if man & 1: return man, exp
        if man & 3: return man>>1, exp+1
        t = trailtable[man & 255]
        if t:
            return man>>t, exp+t
        n = (man^(man-1)) >> 1
        if n < 1024:
            t = bctable[n]
            return man>>t, exp+t
        t = int(math_log(n,2)+0.5)
        return man>>t, exp+t

    # negative
    elif man:

        # Work with a positive mantissa and change the sign back on return.
        # This necessary because bitwise operations are much slower on
        # negative numbers
        pman = -man
        if prec and pman >> prec:

            # Inline bitcounting
            #n = bitcount(pman) - prec
            if pman < (1<<300):
                bc = bisect(powers, pman)
            else:
                bc = int(math_log(pman, 2)) - 4
                bc += bctable[pman>>bc]
            n = bc - prec
            # End bitcount

            if rounding == 'n':
                n1 = n-1
                t = pman >> n1
                if t & 1 and ((t & 2) or (pman & ((1<<n1)-1))):
                    pman = (t >> 1) + 1
                else:
                    pman = t >> 1
            elif rounding in 'fu':
                pman = -(man >> n)
            else:
                pman = pman >> n
            exp += n

        # Remove any trailing zeros
        if pman & 1: return -pman, exp
        if pman & 3: return -(pman>>1), exp+1
        t = trailtable[pman & 255]
        if t:
            return -(pman>>t), exp+t
        n = (pman^(pman-1)) >> 1
        if n < 1024:
            t = bctable[n]
            return -(pman>>t), exp+t
        t = int(math_log(n,2)+0.5)
        return -(pman>>t), exp+t

    else:
        return man, 0

if MODE == 'gmpy':
    import gmpy
    def i_trim(man, exp=0, prec=0, rounding='f'):
        sign, man, exp, bc = from_man_exp(man, exp, prec, rounding)
        if sign:
            return -man, exp
        return man, exp
    if hasattr(gmpy, "_mpmath_trim"):
        pass
        #i_trim2 = gmpy._mpmath_trim

        #def i_trim(man, exp=0, prec=0, rounding='f'):
        #    man, exp = i_trim2(man, exp, prec, rounding)
        #    return man, int(exp)

        #def i_trim(*args, **kwargs):
        #    print "gmpy._mpmath_trim(" + ", ".join(repr(x) for x in args) + ")"
        #    v = gmpy._mpmath_trim(*args, **kwargs)
        #    print "ok", repr(v)
        #    return v

        #i_trim = gmpy._mpmath_trim
        #print "modified trim"
    #print "using gmpy"

def i_add(a, b, c, d, prec=0, rounding='f'):
    """
    Form `a2^b + c2^d`, rounding the result to 

    XXX: assumes both nonzero (?)
    """
    offset = b - d
    if offset >= 0:
        if offset > 300:
            if prec:
                if a:
                    abc = bitcount2(a)
                    cbc = bitcount2(c)
                    delta = abc - cbc + offset
                    if delta > prec + 4:
                        offset = prec + 4
                        a <<= offset
                        if c < 0: a -= 1
                        else:     a += 1
                        return i_trim(a, b-offset, prec, rounding)
                else:
                    return i_trim(c, d, prec, rounding)
        return i_trim((a<<offset)+c, d, prec, rounding)
    else:
        if offset < -300:
            if prec:
                if c:
                    abc = bitcount2(a)
                    cbc = bitcount2(c)
                    delta = cbc - abc - offset
                    if delta > prec + 4:
                        offset = prec + 4
                        c <<= offset
                        if a < 0: c -= 1
                        else:     c += 1
                        return i_trim(c, d-offset, prec, rounding)
                else:
                    return i_trim(a, b, prec, rounding)
        return i_trim((c<<(-offset))+a, b, prec, rounding)

def i_div(a, b, c, d, prec=0, rounding='f'):
    """
    Inputs need not be normalized.
    """
    asign = a < 0
    csign = c < 0
    if csign:
        c = -c
    # Fast division by power of two
    if c == 1:
        if csign:
            a = -a
        return i_trim(a, b-d, prec, rounding)
    if asign:
        a = -a

    abc = bitcount(a)
    cbc = bitcount(c)
    extra = max(0, prec-abc+cbc) + 5

    #extra = prec - int((math_log(a)-math_log(c))*1.4426950408889634) + 5

    #if extra < 5:
    #    extra = 5
    quot, rem = divmod(a<<extra, c)
    if rem:
        quot = (quot<<1) + 1
        extra += 1
    if asign ^ csign:
        quot = -quot
    return i_trim(quot, b-d-extra, prec, rounding)

"""
def i_div(a, b, c, d, prec=0, rounding='f'):
    from mpmath.libmpf import mpf_div
    sign, man, exp, bc = mpf_div(from_man_exp(a,b), from_man_exp(c,d), prec, rounding)
    if sign:
        return -man, exp
    return man, exp
"""

def i_sqrt(man, exp, prec, rounding='f'):
    if man <= 0:
        if man:
            raise ValueError
        return MPZ_0, 0
    if exp & 1:
        exp -= 1
        man <<= 1
    elif man == 1:
        return i_trim(man, exp//2, prec, rounding)
    #shift = 2*prec - int(math_log(man, 2.)) + 4
    shift = 2*prec - bitcount2(man) + 4
    if shift < 4:
        shift = 4
    shift += shift & 1
    if rounding in 'fd':
        man = isqrt(man<<shift)
    else:
        man, rem = sqrtrem(man<<shift)
        # Perturb up
        if rem:
            man = (man<<1)+1
            shift += 2
    return i_trim(man, (exp-shift)//2, prec, rounding)


"""
def i_sqrt(a, b, prec=0, rounding='f'):
    sign, man, exp, bc = mpf_sqrt(from_man_exp(a,b), prec, rounding)
    return man, exp
"""

def i_exp(a, b, prec, rounding='f'):
    sign, man, exp, bc = mpf_exp(from_man_exp(a, b), prec, rounding)
    return man, exp

def i_ln(a, b, prec, rounding='f'):
    sign, man, exp, bc = mpf_log(from_man_exp(a, b), prec, rounding)
    if sign:
        return -man, exp
    return man, exp

def i_nthroot(a, b, n, prec, rounding='f'):
    sign, man, exp, bc = mpf_nthroot(from_man_exp(a, b), n, prec, rounding)
    if sign:
        return -man, exp
    return man, exp

# XXX: argument order
def i_atan2(a, b, c, d, prec, rounding='f'):
    x = from_man_exp(a, b)
    y = from_man_exp(c, d)
    sign, man, exp, bc = mpf_atan2(y, x, prec, rounding)
    if sign:
        return -man, exp
    return man, exp

def i_ln_hypot(a, b, c, d, prec, rounding='f'):
    x = from_man_exp(a, b)
    y = from_man_exp(c, d)
    sign, man, exp, bc = mpf_log_hypot(y, x, prec, rounding)
    if sign:
        return -man, exp
    return man, exp

def i_pow(a, b, c, d, prec, rounding='f'):
    x = from_man_exp(a, b)
    y = from_man_exp(c, d)
    sign, man, exp, bc = mpf_pow(x, y, prec, rounding)
    if sign:
        return -man, exp
    return man, exp

def i_pow_n(a, b, n, prec, rounding='f'):
    sign, man, exp, bc = mpf_pow_int(from_man_exp(a, b), n, prec, rounding)
    if sign:
        return -man, exp
    return man, exp

def i_cos_sin(a, b, prec, rounding='f'):
    a, b = cos_sin(from_man_exp(a, b), prec, rounding)
    asign, aman, aexp, abc = a
    bsign, bman, bexp, bbc = b
    if asign: aman = -aman
    if bsign: bman = -bman
    return aman, aexp, bman, bexp


def light_exp_series(x, prec, r):
    x >>= r
    # 1 + x + x^2/2! + x^3/3! + x^4/4! + ... =
    # (1 + x^2/2! + ...) + x * (1 + x^2/3! + ...)
    s0 = s1 = (MPZ_1 << prec)
    k = 2
    a = x2 = (x*x) >> prec
    while a:
        a //= k; s0 += a; k += 1
        a //= k; s1 += a; k += 1
        a = (a*x2) >> prec
    # Calculate s**(2**r) by repeated squaring
    s1 = (s1*x) >> prec
    s = s0 + s1
    while r:
        s = (s*s) >> prec
        r -= 1
    return s

def exponential_series(x, prec, alt=0, J=2, r=0, exp=0):
    """Sums (cosh(x), sinh(x)) (or (cos(x), sin(x)) with alt=1),
    using J concurrent series, r argument halvings."""
    x >>= r
    one = MPZ_1 << prec
    if J <= 2:
        x2 = a = (x*x) >> prec
        x4 = (x2*x2) >> prec
        s0 = s1 = MPZ_0
        k = 2
        while a:
            a //= (k-1)*k; s0 += a; k += 2
            a //= (k-1)*k; s1 += a; k += 2
            a = (a*x4) >> prec
        s1 = (x2*s1) >> prec
        if alt: c = s1 - s0 + one
        else:   c = s1 + s0 + one
    else:
        x2 = a = (x*x) >> prec
        xpowers = [one, x2]
        for i in xrange(1, J):
            xpowers.append((xpowers[-1]*x2)>>prec)
        sums = [MPZ_0] * J
        k = 2
        while a:
            for i in xrange(J):
                a //= (k-1)*k
                if alt and k & 2: sums[i] -= a
                else:             sums[i] += a
                k += 2
            a = (a*xpowers[-1]) >> prec
        for i in xrange(1, J):
            sums[i] = (sums[i]*xpowers[i]) >> prec
        c = sum(sums) + one
    if exp:
        c += isqrt_fast(abs((one<<prec) - c*c))
        for i in xrange(r):
            c = (c*c) >> prec
        return c

    else:
        # Repeatedly apply the double-angle formula
        # cosh(2*x) = 2*cosh(x)^2 - 1
        # cos(2*x) = 2*cos(x)^2 - 1
        pshift = prec-1
        for i in xrange(r):
            c = ((c*c) >> pshift) - one
        # With the abs, this is the same for sinh and sin
        s = isqrt_fast(abs((one<<prec) - c*c))
        return c, s

def logarithmic_series(x, prec, alt=0, J=2):
    """Sums atanh(x) (or atan(x) with alt=1), using J concurrent series."""
    if J <= 2:
        s0 = a = x
        x2 = x*x >> prec
        x4 = x2*x2 >> prec
        s1 = a//3
        a = x*x4 >> prec
        k = 5
        while a:
            s0 += a // k; k += 2
            s1 += a // k; k += 2
            a = a*x4 >> prec
        s1 = x2*s1 >> prec
        if alt:
            return s0 - s1
        else:
            return s0 + s1
    else:
        x2 = x*x >> prec
        a = x
        xpowers = [(MPZ_1 << prec), x2]
        for i in xrange(1, J):
            xpowers.append((xpowers[-1]*x2)>>prec)
        sums = [MPZ_0] * J
        k = 1
        while a:
            for i in xrange(J):
                if alt and k & 2: sums[i] -= a // k
                else:             sums[i] += a // k
                k += 2
            a = (a*xpowers[-1]) >> prec
        for i in xrange(1, J):
            sums[i] = (sums[i]*xpowers[i]) >> prec
        return sum(sums)

ln2_cache = [MPZ_1, 0]

def i_ln2(prec):
    cman, cprec = ln2_cache
    shift = cprec - prec
    if shift >= 0:
        return cman >> shift 
    else:
        ln2_cache[:] = [ln2_fixed(prec), prec]
        return i_ln2(prec)

pi_cache = [MPZ_1, 0]

def i_pi(prec):
    cman, cprec = pi_cache
    shift = cprec - prec
    if shift >= 0:
        return cman >> shift 
    else:
        pi_cache[:] = [pi_fixed(prec), prec]
        return i_pi(prec)


if MODE == 'gmpy':
    EXP_SERIES_CUTOFF = 450
else:
    EXP_SERIES_CUTOFF = 850


def i_exp2(man, exp, prec, rounding='f'):
    if not man:
        return MPZ_1, 0
    bc = bitcount2(man)
    mag = exp + bc
    if mag < -prec-10:
        return MPZ_1, 0
    # extra precision needs to be similar in magnitude to log_2(|x|)
    # for the modulo reduction, plus r for the error from squaring r times
    wp = prec + max(0, mag)
    light_series = wp < EXP_SERIES_CUTOFF
    if light_series:
        r = int(wp ** 0.5)
    else:
        # Todo: tuning / appropriate cutoffs
        J = 0.5 * wp ** 0.3
        r = 2.4 * J
        if wp < 1500:
            J = 2
        J = int(J)
        r = int(r)
    if mag < 0:
        r = max(1, r+mag)

    wp += r + 20
    offset = exp + wp
    if offset >= 0:
        t = man << offset
    else:
        t = man >> (-offset)

    if mag > 1:
        ln2 = i_ln2(wp)
        n, t = divmod(t, ln2)
        n = int(n)
    else:
        n = 0
    if light_series:
        m = light_exp_series(t, wp, r)
    else:
        m = exponential_series(t, wp, 0, J, r, 1)

    if rounding:
        return i_trim(m, n-wp, prec, rounding)

    return m, n-wp

# XXX
def i_cos_sin2(man, exp, prec, rounding=None, mulpi=False):
    if not man:
        return MPZ_1, 0, MPZ_0, 0
    bc = bitcount2(man)
    mag = exp + bc
    wp = prec + mag + 20
    if wp < 1500:
        J = 2
        r = int(0.5 * wp ** 0.4)
    else:
        J = 0.5 * wp ** 0.3
        r = 2.4 * J
        J = int(J)
        r = int(r)
    wp += r
    offset = exp + wp
    if offset >= 0:
        t = man << offset
    else:
        t = man >> (-offset)
    pi4 = i_pi(wp-2)
    n, t = divmod(t, pi4)
    n = int(n)
    cos, sin = exponential_series(t, wp, 1, J, r)
    if rounding:
        cosm, cose = i_trim(cos, -wp, prec, rounding)
        sinm, sine = i_trim(sin, -wp, prec, rounding)
        return cosm, cose, sinm, sine
    return cos, -wp, sin, -wp


from mpmath.lib import bitcount, to_fixed, from_man_exp, round_nearest
from mpmath.libmpc import complex_int_pow
from mpmath.mptypes import mp, mpnumeric, mpf, mpc, make_mpf, make_mpc

from factorials import int_fac

# A Bessel function of the first kind of integer order, J_n(x), is
# given by the power series

#             oo
#             ___         k         2 k + n
#            \        (-1)     / x \
#    J_n(x) = )    ----------- | - |
#            /___  k! (k + n)! \ 2 /
#            k = 0

# This series can be summed effectively using fixed-point arithmetic,
# basing each term on the previous. In fact, it is very similar to
# the Taylor series summation for a sine or cosine; the main difference
# is that the initial term contains a power and a factorial.

# Simplifying the quotient between two successive terms gives the
# ratio x^2 / (-4*k*(k+n)). Hence, we only need one full-precision
# multiplication and one division by a small integer per term.
# The complex version is very similar, the only difference being
# that the multiplication is actually 4 multiplies.

# TODO: for extremely large x, we could use an asymptotic
# trigonometric approximation.

# TODO: recompute at higher precision if the fixed-point mantissa
# is very small

def mpf_jn_series(n, x, prec):
    negate = n < 0 and n & 1
    n = abs(n)
    origprec = prec
    prec += 20 + bitcount(abs(n))
    x = to_fixed(x, prec)
    x2 = (x**2) >> prec
    if not n:
        s = t = 1 << prec
    else:
        s = t = (x**n // int_fac(n)) >> ((n-1)*prec + n)
    k = 1
    while t:
        t = ((t * x2) // (-4*k*(k+n))) >> prec
        s += t
        k += 1
    if negate:
        s = -s
    return from_man_exp(s, -prec, origprec, round_nearest)

def mpc_jn_series(n, z, prec):
    negate = n < 0 and n & 1
    n = abs(n)
    origprec = prec
    prec += 20 + bitcount(abs(n))
    zre, zim = z
    zre = to_fixed(zre, prec)
    zim = to_fixed(zim, prec)
    z2re = (zre**2 - zim**2) >> prec
    z2im = (zre*zim) >> (prec-1)
    if not n:
        sre = tre = 1 << prec
        sim = tim = 0
    else:
        re, im = complex_int_pow(zre, zim, n)
        sre = tre = (re // int_fac(n)) >> ((n-1)*prec + n)
        sim = tim = (im // int_fac(n)) >> ((n-1)*prec + n)
    k = 1
    while abs(tre) + abs(tim) > 3:
        p = -4*k*(k+n)
        tre, tim = tre*z2re - tim*z2im, tim*z2re + tre*z2im
        tre = (tre // p) >> prec
        tim = (tim // p) >> prec
        sre += tre
        sim += tim
        k += 1
    if negate:
        sre = -sre
        sim = -sim
    re = from_man_exp(sre, -prec, origprec, round_nearest)
    im = from_man_exp(sim, -prec, origprec, round_nearest)
    return re, im

def jv(v, x):
    """
    Compute the Bessel function J_v(x).
    """
    assert type(v) is int
    prec = mp.prec
    x = mpnumeric(x)
    if isinstance(x, mpf):
        return make_mpf(mpf_jn_series(v, x._mpf_, prec))
    if isinstance(x, mpc):
        return make_mpc(mpc_jn_series(v, (x.real._mpf_, x.imag._mpf_), prec))

jn = jv

def j0(x):
    """Bessel function J_0(x)."""
    return jv(0, x)

def j1(x):
    """Bessel function J_1(x)."""
    return jv(1, x)

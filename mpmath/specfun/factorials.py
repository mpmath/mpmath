from mpmath.lib import *
from mpmath.libmpc import *
from mpmath.mptypes import mp, mpf, make_mpf, make_mpc, convert_lossless

def int_fac(n, memo={0:1, 1:1}):
    """Return n factorial (for integers n >= 0 only)."""
    f = memo.get(n)
    if f:
        return f
    k = len(memo)
    p = memo[k-1]
    while k <= n:
        p *= k
        if k < 1024:
            memo[k] = p
        k += 1
    return p

spouge_cache = {}

def calc_spouge_coefficients(a, prec):
    rnd = round_floor
    wp = prec + int(a*1.4)
    c = [0] * a
    # b = exp(a-1)
    b = fexp(from_int(a-1), wp, rnd)
    # e = exp(1)
    e = fexp(fone, wp, rnd)
    # sqrt(2*pi)
    sq2pi = fsqrt(fshift(fpi(wp, rnd), 1), wp, rnd)
    c[0] = to_fixed(sq2pi, prec)
    for k in xrange(1, a):
        # c[k] = ((-1)**(k-1) * (a-k)**k) * b / sqrt(a-k)
        term = fmuli(b, ((-1)**(k-1) * (a-k)**k), wp, rnd)
        term = fdiv(term, fsqrt(from_int(a-k), wp, rnd), wp, rnd)
        c[k] = to_fixed(term, prec)
        # b = b / (e * k)
        b = fdiv(b, fmul(e, from_int(k), wp, rnd), wp, rnd)
    return c

# Cached lookup of coefficients
def get_spouge_coefficients(prec):

    # This exact precision has been used before
    if prec in spouge_cache:
        return spouge_cache[prec]

    for p in spouge_cache:
        if 0.8 <= float(p)/prec < 1:
            return spouge_cache[p]

    # Here we estimate the value of a based on Spouge's inequality for
    # the relative error
    a = max(3, int(0.39*prec))  # ~= 1.26*n

    coefs = calc_spouge_coefficients(a, prec)
    spouge_cache[prec] = (prec, a, coefs)
    return spouge_cache[prec]

def spouge_sum_real(x, prec, a, c):
    x = to_fixed(x, prec)
    s = c[0]
    for k in xrange(1, a):
        s += (c[k] << prec) // (x + (k << prec))
    return from_man_exp(s, -prec, prec, round_floor)

# Unused: for fast computation of gamma(p/q)
def spouge_sum_rational(p, q, prec, a, c):
    s = c[0]
    for k in xrange(1, a):
        s += c[k] * q // (p+q*k)
    return from_man_exp(s, -prec, prec, round_floor)

# For a complex number a + b*I, we have
#
#        c_k          (a+k)*c_k     b * c_k
#  -------------  =   ---------  -  ------- * I
#  (a + b*I) + k          M            M
#
#                 2    2      2   2              2
#  where M = (a+k)  + b   = (a + b ) + (2*a*k + k )
#
def spouge_sum_complex(re, im, prec, a, c):
    re = to_fixed(re, prec)
    im = to_fixed(im, prec)
    sre, sim = c[0], 0
    mag = ((re**2)>>prec) + ((im**2)>>prec)
    for k in xrange(1, a):
        M = mag + re*(2*k) + ((k**2) << prec)
        sre += (c[k] * (re + (k << prec))) // M
        sim -= (c[k] * im) // M
    re = from_man_exp(sre, -prec, prec, round_floor)
    im = from_man_exp(sim, -prec, prec, round_floor)
    return re, im

def mpf_gamma(x, prec, rounding, p1=1):
    sign, man, exp, bc = x
    if exp >= 0:
        if sign or (p1 and not man):
            raise ValueError("gamma function pole")
        # A direct factorial is fastest
        if exp + bc <= 10:
            return from_int(int_fac((man<<exp)-p1), prec, rounding)
    wp = prec + 15
    rnd = round_floor
    if p1:
        x = fsub(x, fone, wp, rnd)
    # x < 0.25
    if sign or exp+bc < -1:
        # gamma = pi / (sin(pi*x) * gamma(1-x))
        wp += 15
        pi = fpi(wp, rnd)
        pix = fmul(x, pi, wp, rnd)
        t = fsin(pix, wp, rnd)
        g = mpf_gamma(fsub(fone, x, wp, rnd), wp, rnd)
        return fdiv(pix, fmul(t, g, wp, rnd), prec, rounding)
    sprec, a, c = get_spouge_coefficients(wp)
    s = spouge_sum_real(x, sprec, a, c)
    # gamma = exp(log(x+a)*(x+0.5) - xpa) * s
    xpa = fadd(x, from_int(a), wp, rnd)
    logxpa = flog(xpa, wp, rnd)
    xph = fadd(x, fhalf, wp, rnd)
    t = fsub(fmul(logxpa, xph, wp, rnd), xpa, wp, rnd)
    t = fmul(fexp(t, wp, rnd), s, prec, rounding)
    return t

def mpc_gamma(x, prec, rounding, p1=1):
    re, im = x
    if im == fzero:
        return mpf_gamma(re, prec, rounding, p1), fzero
    wp = prec + 25
    rnd = round_floor
    sign, man, exp, bc = re
    if p1:
        re = fsub(re, fone, wp, rnd)
        x = re, im
    if sign or exp+bc < -1:
        # Reflection formula
        wp += 15
        pi = fpi(wp, rnd), fzero
        pix = mpc_mul(x, pi, wp, rnd)
        t = mpc_sin(pix, wp, rnd)
        u = mpc_sub(mpc_one, x, wp, rnd)
        g = mpc_gamma(u, wp, rnd)
        w = mpc_mul(t, g, wp, rnd)
        return mpc_div(pix, w, wp, rnd)
    sprec, a, c = get_spouge_coefficients(wp)
    s = spouge_sum_complex(re, im, sprec, a, c)
    # gamma = exp(log(x+a)*(x+0.5) - xpa) * s
    repa = fadd(re, from_int(a), wp, rnd)
    logxpa = mpc_log((repa, im), wp, rnd)
    reph = fadd(re, fhalf, wp, rnd)
    t = mpc_sub(mpc_mul(logxpa, (reph, im), wp, rnd), (repa, im), wp, rnd)
    t = mpc_mul(mpc_exp(t, wp, rnd), s, prec, rounding)
    return t

def gamma(x):
    x = convert_lossless(x)
    prec = mp.prec
    if isinstance(x, mpf):
        return make_mpf(mpf_gamma(x._mpf_, prec, round_nearest, 1))
    else:
        return make_mpc(mpc_gamma(x._mpc_, prec, round_nearest, 1))

def factorial(x):
    x = convert_lossless(x)
    prec = mp.prec
    if isinstance(x, mpf):
        return make_mpf(mpf_gamma(x._mpf_, prec, round_nearest, 0))
    else:
        return make_mpc(mpc_gamma(x._mpc_, prec, round_nearest, 0))

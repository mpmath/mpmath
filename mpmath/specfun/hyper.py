"""
Hypergeometric functions.
"""

from ..lib import from_man_exp, to_fixed
from ..mptypes import mp, mpf, mpc, make_mpf, make_mpc, \
    convert_lossless, inf, pi

def mpf_hyp2f1((ap, aq), (bp, bq), (cp, cq), x, prec, rnd):
    wp = prec + 25
    x = to_fixed(x, wp)
    s = p = 1 << wp
    n = 1
    while 1:
        p = (p * (ap*bp*cq*x) // (n*aq*bq*cp)) >> wp
        if -100 < p < 100:
            break
        s += p; n += 1; ap += aq; bp += bq; cp += cq
    return from_man_exp(s, -wp, prec, rnd)

def mpc_hyp2f1((ap, aq), (bp, bq), (cp, cq), z, prec, rnd):
    wp = prec + 25
    zre, zim = z
    zre = to_fixed(zre, wp)
    zim = to_fixed(zim, wp)
    sre = pre = 1 << wp
    sim = pim = 0
    n = 1
    while 1:
        r1 = ap*bp*cq
        r2 = n*aq*bq*cp
        pre, pim = pre*zre - pim*zim, pim*zre + pre*zim
        pre = ((pre * r1) // r2) >> wp
        pim = ((pim * r1) // r2) >> wp
        if -100 < pre < 100 and -100 < pim < 100:
            break
        sre += pre; sim += pim; n += 1; ap += aq; bp += bq; cp += cq
    re = from_man_exp(sre, -wp, prec, rnd)
    im = from_man_exp(sim, -wp, prec, rnd)
    return re, im

def as_frac(x):
    if isinstance(x, (int, long)):
        return x, 1
    p, q = x
    return p, q

def check_convergence(a,b,c,x):
    ap,aq = a
    bp,bq = b
    cp,cq = c
    # c may not be a negative integer
    if cp*cq < 0 and not cp % cq:
        return False
    if abs(x) < 1:
        return True
    return False

def hyp2f1(a,b,c,x):
    """Hypergeometric function 2F1(a,b,c,x). The parameters a, b, c
    must be ints or fractions p/q specified as tuples (p, q). Currently
    x is restricted to |x| < 1, where the function can be computed from
    the hypergeometric series."""
    a = as_frac(a)
    b = as_frac(b)
    c = as_frac(c)
    x = convert_lossless(x)
    if not check_convergence(a,b,c,x):
        raise NotImplementedError("divergent hypergeometric series")
    if isinstance(x, mpf):
        return make_mpf(mpf_hyp2f1(a,b,c,x._mpf_,mp.prec,mp.rounding[0]))
    if isinstance(x, mpc):
        return make_mpc(mpc_hyp2f1(a,b,c,x._mpc_,mp.prec,mp.rounding[0]))

def ellipk(m):
    """Complete elliptic integral of the first kind, K(m). Note that
    the argument is the parameter m = k^2, not the modulus k."""
    m = convert_lossless(m)
    if m == 1:
        return inf
    return pi/2 * hyp2f1((1,2),(1,2),(1,1),m)

def ellipe(m):
    """Complete elliptic integral of the second kind, E(m). Note that
    the argument is the parameter m = k^2, not the modulus k."""
    m = convert_lossless(m)
    if m == 1:
        return m
    return pi/2 * hyp2f1((1,2),(-1,2),(1,1),m)

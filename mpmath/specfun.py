"""
Miscellaneous special functions
"""

from lib import *
from libmpc import *
from mptypes import *
from mptypes import constant, funcwrapper

from gammazeta import *

__docformat__ = 'plaintext'

#---------------------------------------------------------------------------#
#                                                                           #
#                       First some mathematical constants                   #
#                                                                           #
#---------------------------------------------------------------------------#


# The golden ratio is given by phi = (1 + sqrt(5))/2

@constant_memo
def phi_fixed(prec):
    prec += 10
    sqrt = [sqrt_fixed2, sqrt_fixed][prec < 20000]
    a = sqrt(MP_FIVE<<prec, prec) + (MP_ONE << prec)
    return a >> 11

# Catalan's constant is computed using Lupas's rapidly convergent series
# (listed on http://mathworld.wolfram.com/CatalansConstant.html)
#            oo
#            ___       n-1  8n     2                   3    2
#        1  \      (-1)    2   (40n  - 24n + 3) [(2n)!] (n!)
#  K =  ---  )     -----------------------------------------
#       64  /___               3               2
#                             n  (2n-1) [(4n)!]
#           n = 1

@constant_memo
def catalan_fixed(prec):
    prec = prec + 20
    a = one = MP_ONE << prec
    s, t, n = 0, 1, 1
    while t:
        a *= 32 * n**3 * (2*n-1)
        a //= (3-16*n+16*n**2)**2
        t = a * (-1)**(n-1) * (40*n**2-24*n+3) // (n**3 * (2*n-1))
        s += t
        n += 1
    return s >> (20 + 6)

# Khinchin's constant is relatively difficult to compute. Here
# we use the rational zeta series

#                    oo                2*n-1
#                   ___                ___
#                   \   ` zeta(2*n)-1  \   ` (-1)^(k+1)
#  log(K)*log(2) =   )    ------------  )    ----------
#                   /___.      n       /___.      k
#                   n = 1              k = 1

# which adds half a digit per term. The essential trick for achieving
# reasonable efficiency is to recycle both the values of the zeta
# function (essentially Bernoulli numbers) and the partial terms of
# the inner sum.

# An alternative might be to use K = 2*exp[1/log(2) X] where

#      / 1     1       [ pi*x*(1-x^2) ]
#  X = |    ------ log [ ------------ ].
#      / 0  x(1+x)     [  sin(pi*x)   ]

# and integrate numerically. In practice, this seems to be slightly
# slower than the zeta series at high precision.

@constant_memo
def khinchin_fixed(prec):
    orig = mp.prec
    try:
        mp.prec = int(prec + prec**0.5 + 15)
        s = mpf(0)
        t = one = mpf(1)
        fac = mpf(4)
        pipow = twopi2 = (2*pi)**2
        n = 1
        while 1:
            zeta2n = (-1)**(n+1) * bernoulli(2*n) * pipow / fac
            term = ((zeta2n - 1) * t) / n
            # print n, nstr(term)
            if term < eps:
                break
            s += term
            t += (one/(2*n+1) - one/(2*n))
            n += 1
            fac *= (2*n)*(2*n-1)
            pipow *= twopi2
        return to_fixed(exp(s/ln2)._mpf_, prec)
    finally:
        mp.prec = orig

# Glaisher's constant is defined as A = exp(1/2 - zeta'(-1)).
# One way to compute it would be to perform direct numerical
# differentiation, but computing arbitrary Riemann zeta function
# values at high precision is expensive. We instead use the formula

#     A = exp((6 (-zeta'(2))/pi^2 + log 2 pi + gamma)/12)

# and compute zeta'(2) from the series representation

#              oo
#              ___
#             \     log k
#  -zeta'(2) = )    -----
#             /___     2
#                    k
#            k = 2

# This series converges exceptionally slowly, but can be accelerated
# using Euler-Maclaurin formula. The important insight is that the
# E-M integral can be done in closed form and that the high order
# are given by

#    n  /       \
#   d   | log x |   a + b log x
#   --- | ----- | = -----------
#     n |   2   |      2 + n
#   dx  \  x    /     x

# where a and b are integers given by a simple recurrence. Note
# that just one logarithm is needed. However, lots of integer
# logarithms are required for the initial summation.

# This algorithm could possibly be turned into a faster algorithm
# for general evaluation of zeta(s) or zeta'(s); this should be
# looked into.

@constant_memo
def glaisher_fixed(prec):
    orig = mp.prec
    try:
        dps = mp.dps
        mp.prec = prec + 30
        N = int(1.0*dps + 5)
        logs = log_range()
        s = mpf(0)
        # E-M step 1: sum log(k)/k**2 for k = 2..N-1
        for n in range(2, N):
            # print n, N
            logn = logs.next()
            s += logn / n**2
        logN = logs.next()
        # E-M step 2: integral of log(x)/x**2 from N to inf
        s += (1+logN)/N
        # E-M step 3: endpoint correction term f(N)/2
        s += logN/(N**2 * 2)
        # E-M step 4: the series of derivatives
        pN, a, b, j, fac, k = N**3, 1, -2, 3, 2, 1
        while 1:
            # D(2*k-1) * B(2*k) / fac(2*k) [D(n) = nth derivative]
            D = (a+b*logN)/pN
            B = bernoulli(2*k)
            term = B * D / fac
            if abs(term) < eps:
                break
            # print k, nstr(term)
            s -= term
            # Advance derivative twice
            a, b, pN, j = b-a*j, -j*b, pN*N, j+1
            a, b, pN, j = b-a*j, -j*b, pN*N, j+1
            k += 1
            fac *= (2*k) * (2*k-1)
        A = exp((6*s/pi**2 + log(2*pi) + euler)/12)
        return to_fixed(A._mpf_, prec)
    finally:
        mp.prec = orig

# Apery's constant can be computed using the very rapidly convergent
# series
#              oo
#              ___              2                      10
#             \         n  205 n  + 250 n + 77     (n!)
#  zeta(3) =   )    (-1)   -------------------  ----------
#             /___               64                      5
#             n = 0                             ((2n+1)!)

@constant_memo
def apery_fixed(prec):
    prec += 20
    d = MP_ONE << prec
    term = MP_BASE(77) << prec
    n = 1
    s = MP_ZERO
    while term:
        s += term
        d *= (n**10)
        d //= (((2*n+1)**5) * (2*n)**5)
        term = (-1)**n * (205*(n**2) + 250*n + 77) * d
        n += 1
    return s >> (20 + 6)

fme = from_man_exp

def mpf_phi(p, r): return fme(phi_fixed(p+10), -p-10, p, r)
def mpf_khinchin(p, r): return fme(khinchin_fixed(p+10), -p-10, p, r)
def mpf_glaisher(p, r): return fme(glaisher_fixed(p+10), -p-10, p, r)
def mpf_apery(p, r): return fme(apery_fixed(p+10), -p-10, p, r)
def mpf_catalan(p, r): return fme(catalan_fixed(p+10), -p-10, p, r)

phi = constant(mpf_phi, "Golden ratio (phi)")
catalan = constant(mpf_catalan, "Catalan's constant")
khinchin = constant(mpf_khinchin, "Khinchin's constant")
glaisher = constant(mpf_glaisher, "Glaisher's constant")
apery = constant(mpf_apery, "Apery's constant")


def isnpint(x):
    if not x:
        return True
    if isinstance(x, mpf):
        sign, man, exp, bc = x._mpf_
        return sign and exp >= 0
    if isinstance(x, mpc):
        return not x.imag and isnpint(x.real)

def gammaprod(a, b):
    """
    Computes the product / quotient of gamma functions

        G(a_0) G(a_1) ... G(a_p)
        ------------------------
        G(b_0) G(b_1) ... G(a_q)

    with proper cancellation of poles (interpreting the expression as a
    limit). Returns +inf if the limit diverges.
    """
    a = [convert_lossless(x) for x in a]
    b = [convert_lossless(x) for x in b]
    poles_num = []
    poles_den = []
    regular_num = []
    regular_den = []
    for x in a: [regular_num, poles_num][isnpint(x)].append(x)
    for x in b: [regular_den, poles_den][isnpint(x)].append(x)
    # One more pole in numerator or denominator gives 0 or inf
    if len(poles_num) < len(poles_den): return mpf(0)
    if len(poles_num) > len(poles_den): return mpf('+inf')
    # All poles cancel
    # lim G(i)/G(j) = (-1)**(i+j) * gamma(1-j) / gamma(1-i)
    p = mpf(1)
    orig = mp.prec
    try:
        mp.prec = orig + 15
        while poles_num:
            i = poles_num.pop()
            j = poles_den.pop()
            p *= (-1)**(i+j) * gamma(1-j) / gamma(1-i)
        for x in regular_num: p *= gamma(x)
        for x in regular_den: p /= gamma(x)
    finally:
        mp.prec = orig
    return +p

def binomial(n, k):
    """Binomial coefficient, C(n,k) = n!/(k!*(n-k)!)."""
    return gammaprod([n+1], [k+1, n-k+1])

def rf(x, n):
    """Rising factorial (Pochhammer symbol), x^(n)"""
    return gammaprod([x+n], [x])

def ff(x, n):
    """Falling factorial, x_(n)"""
    return gammaprod([x+1], [x-n+1])


#---------------------------------------------------------------------------#
#                                                                           #
#                          Hypergeometric functions                         #
#                                                                           #
#---------------------------------------------------------------------------#

import operator

"""
TODO:
  * By counting the number of multiplications vs divisions,
    the bit size of p can be kept around wp instead of growing
    it to n*wp for some (possibly large) n

  * Due to roundoff error, the series may fail to converge
    when x is negative and the convergence is slow.

"""

def hypsum(ar, af, ac, br, bf, bc, x):
    """
    Generic hypergeometric summation. This function computes:

            1   a_1 a_2 ...     1  (a_1 + 1) (a_2 + 1) ...  2
        1 + --  ----------- x + -- ----------------------- x  + ...
            1!  b_1 b_2 ...     2! (b_1 + 1) (b_2 + 1) ...

    The a_i and b_i sequences are separated by type:

    ar - list of a_i rationals [p,q]
    af - list of a_i mpf value tuples
    ac - list of a_i mpc value tuples
    br - list of b_i rationals [p,q]
    bf - list of b_i mpf value tuples
    bc - list of b_i mpc value tuples

    Note: the rational coefficients will be updated in-place and must
    hence be mutable (lists rather than tuples).

    x must be an mpf or mpc instance.
    """

    have_float = af or bf
    have_complex = ac or bc

    prec = mp.prec
    rnd = round_nearest
    wp = prec + 25

    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        y = MP_ZERO
    else:
        have_complex = 1
        x, y = x._mpc_
        x = to_fixed(x, wp)
        y = to_fixed(y, wp)

    sre = pre = one = MP_ONE << wp
    sim = pim = MP_ZERO

    n = 1

    # Need to shift down by wp once for each fixed-point multiply
    # At minimum, we multiply by once by x each step
    shift = 1

    # Fixed-point real coefficients
    if have_float:
        len_af = len(af)
        len_bf = len(bf)
        range_af = range(len_af)
        range_bf = range(len_bf)
        for i in range_af: af[i] = to_fixed(af[i], wp)
        for i in range_bf: bf[i] = to_fixed(bf[i], wp)
        shift += len_af

    if have_complex:
        len_ac = len(ac)
        len_bc = len(bc)
        range_ac = range(len_ac)
        range_bc = range(len_bc)
        for i in range_ac: ac[i] = [to_fixed(ac[i][0], wp), to_fixed(ac[i][1], wp)]
        for i in range_bc: bc[i] = [to_fixed(bc[i][0], wp), to_fixed(bc[i][1], wp)]
        shift += len_ac

    aqs = [a[1] for a in ar]
    bqs = [b[1] for b in br]
    aqprod = reduce(operator.mul, aqs, 1)
    bqprod = reduce(operator.mul, bqs, 1)

    assert shift >= 0

    while 1:
        # Integer and rational part of product
        mul = bqprod
        div = n * aqprod
        for ap, aq in ar: mul *= ap
        for bp, bq in br: div *= bp

        if have_complex:
            # Multiply by rational factors
            pre *= mul
            pim *= mul
            # Multiply by z
            pre, pim = pre*x - pim*y, pim*x + pre*y
            # Multiply by real factors
            for a in af:
                pre *= a
                pim *= a
            # Multiply by complex factors
            for are, aim in ac:
                pre, pim = pre*are - pim*aim, pim*are + pre*aim
            # Divide by rational factors
            pre //= div
            pim //= div
            # Divide by real factors
            for b in bf:
                pre = (pre << wp) // b
                pim = (pim << wp) // b
            # Divide by complex factors
            for bre, bim in bc:
                mag = bre*bre + bim*bim
                re = pre*bre + pim*bim
                im = pim*bre - pre*bim
                pre = (re << wp) // mag
                pim = (im << wp) // mag
        elif have_float:
            # Multiply and divide by real and rational factors, and x
            for a in af: pre *= a
            for b in bf:
                pre = (pre << wp) // b
            pre = (pre * (mul * x)) // div

        else:
            # Multiply and divide by rational factors and x
            pre = (pre * (mul * x)) // div

        pre >>= (wp*shift)
        sre += pre

        if have_complex:
            pim >>= (wp*shift)
            sim += pim
            if (-100 < pre < 100) and (-100 < pim < 100):
                break
        else:
            if -100 < pre < 100:
                break

        # Add 1 to all as and bs
        n += 1
        for ap_aq in ar: ap_aq[0] += ap_aq[1]
        for bp_bq in br: bp_bq[0] += bp_bq[1]
        if have_float:
            for i in range_af: af[i] += one
            for i in range_bf: bf[i] += one
        if have_complex:
            for i in range_ac: ac[i][0] += one
            for i in range_bc: bc[i][0] += one

    re = from_man_exp(sre, -wp, prec, rnd)
    if have_complex:
        return make_mpc((re, from_man_exp(sim, -wp, prec, rnd)))
    else:
        return make_mpf(re)


#---------------------------------------------------------------------------#
#   Special-case implementation for rational parameters. These are          #
#   about 2x faster at low precision                                        #
#---------------------------------------------------------------------------#

def sum_hyp0f1_rat((bp, bq), x):
    """Sum 0F1 for rational a. x must be mpf or mpc."""
    prec = mp.prec
    rnd = round_nearest
    wp = prec + 25
    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        s = p = MP_ONE << wp
        n = 1
        while 1:
            p = (p * (bq*x) // (n*bp)) >> wp
            if -100 < p < 100:
                break
            s += p; n += 1; bp += bq
        return make_mpf(from_man_exp(s, -wp, prec, rnd))
    else:
        wp = prec + 25
        zre, zim = x._mpc_
        zre = to_fixed(zre, wp)
        zim = to_fixed(zim, wp)
        sre = pre = MP_ONE << wp
        sim = pim = MP_ZERO
        n = 1
        while 1:
            r1 = bq
            r2 = n*bp
            pre, pim = pre*zre - pim*zim, pim*zre + pre*zim
            pre = ((pre * r1) // r2) >> wp
            pim = ((pim * r1) // r2) >> wp
            if -100 < pre < 100 and -100 < pim < 100:
                break
            sre += pre; sim += pim; n += 1; bp += bq
        re = from_man_exp(sre, -wp, prec, rnd)
        im = from_man_exp(sim, -wp, prec, rnd)
        return make_mpc((re, im))


def sum_hyp1f1_rat((ap, aq), (bp, bq), x):
    """Sum 1F1 for rational a, b. x must be mpf or mpc."""
    prec = mp.prec
    rnd = round_nearest
    wp = prec + 25
    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        s = p = MP_ONE << wp
        n = 1
        while 1:
            p = (p * (ap*bq*x) // (n*aq*bp)) >> wp
            if -100 < p < 100:
                break
            s += p; n += 1; ap += aq; bp += bq
        return make_mpf(from_man_exp(s, -wp, prec, rnd))
    else:
        wp = prec + 25
        zre, zim = x._mpc_
        zre = to_fixed(zre, wp)
        zim = to_fixed(zim, wp)
        sre = pre = MP_ONE << wp
        sim = pim = MP_ZERO
        n = 1
        while 1:
            r1 = ap*bq
            r2 = n*aq*bp
            pre, pim = pre*zre - pim*zim, pim*zre + pre*zim
            pre = ((pre * r1) // r2) >> wp
            pim = ((pim * r1) // r2) >> wp
            if -100 < pre < 100 and -100 < pim < 100:
                break
            sre += pre; sim += pim; n += 1; ap += aq; bp += bq
        re = from_man_exp(sre, -wp, prec, rnd)
        im = from_man_exp(sim, -wp, prec, rnd)
        return make_mpc((re, im))

def sum_hyp2f1_rat((ap, aq), (bp, bq), (cp, cq), x):
    """Sum 2F1 for rational a, b, c. x must be mpf or mpc"""
    prec = mp.prec
    rnd = round_nearest
    wp = prec + 25
    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        s = p = MP_ONE << wp
        n = 1
        while 1:
            p = (p * (ap*bp*cq*x) // (n*aq*bq*cp)) >> wp
            if -100 < p < 100:
                break
            s += p; n += 1; ap += aq; bp += bq; cp += cq
        return make_mpf(from_man_exp(s, -wp, prec, rnd))
    else:
        wp = prec + 25
        zre, zim = x._mpc_
        zre = to_fixed(zre, wp)
        zim = to_fixed(zim, wp)
        sre = pre = MP_ONE << wp
        sim = pim = MP_ZERO
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
        return make_mpc((re, im))

def parse_param(x):
    if isinstance(x, tuple):
        p, q = x
        return [[p, q]], [], []
    if isinstance(x, (int, long)):
        return [[x, 1]], [], []
    x = convert_lossless(x)
    if isinstance(x, mpf):
        return [], [x._mpf_], []
    if isinstance(x, mpc):
        return [], [], [x._mpc_]

class _mpq(tuple):
    @property
    def _mpf_(self):
        return (mpf(self[0])/self[1])._mpf_
    def __add__(self, other):
        if isinstance(other, _mpq):
            a, b = self
            c, d = other
            return _mpq((a*d+b*c, b*d))
        return NotImplemented
    def __sub__(self, other):
        if isinstance(other, _mpq):
            a, b = self
            c, d = other
            return _mpq((a*d-b*c, b*d))
        return NotImplemented

_1 = _mpq((1,1))
_0 = _mpq((0,1))

def _as_num(x):
    if isinstance(x, list):
        return _mpq(x)
    return x

def eval_hyp2f1(a,b,c,z):
    ar, af, ac = parse_param(a)
    br, bf, bc = parse_param(b)
    cr, cf, cc = parse_param(c)
    absz = abs(z)
    if absz == 1:
        # TODO: determine whether it actually does, and otherwise
        # return infinity instead
        print "Warning: 2F1 might not converge for |z| = 1"
    if absz <= 1:
        if ar and br and cr:
            return sum_hyp2f1_rat(ar[0], br[0], cr[0], z)
        return hypsum(ar+br, af+bf, ac+bc, cr, cf, cc, z)
    # Use 1/z transformation
    a = (ar and _as_num(ar[0])) or convert_lossless(a)
    b = (br and _as_num(br[0])) or convert_lossless(b)
    c = (cr and _as_num(cr[0])) or convert_lossless(c)
    orig = mp.prec
    try:
        mp.prec = orig + 15
        h1 = eval_hyp2f1(a, _1-c+a, _1-b+a, 1/z)
        h2 = eval_hyp2f1(b, _1-c+b, _1-a+b, 1/z)
        #s1 = G(c)*G(b-a)/G(b)/G(c-a) * (-z)**(-a) * h1
        #s2 = G(c)*G(a-b)/G(a)/G(c-b) * (-z)**(-b) * h2
        f1 = gammaprod([c,b-a],[b,c-a])
        f2 = gammaprod([c,a-b],[a,c-b])
        s1 = f1 * (-z)**(_0-a) * h1
        s2 = f2 * (-z)**(_0-b) * h2
        v = s1 + s2
    finally:
        mp.prec = orig
    return +v

#---------------------------------------------------------------------------#
#                      And now the user-friendly versions                   #
#---------------------------------------------------------------------------#

def hyper(as, bs, z):
    """
    Hypergeometric function pFq,

          [ a_1, a_2, ..., a_p |    ]
      pFq [                    |  z ]
          [ b_1, b_2, ..., b_q |    ]

    The parameter lists as and bs may contain real or complex numbers.
    Exact rational parameters can be given as tuples (p, q).
    """
    p = len(as)
    q = len(bs)
    z = convert_lossless(z)
    degree = p, q
    if degree == (0, 1):
        br, bf, bc = parse_param(bs[0])
        if br:
            return sum_hyp0f1_rat(br[0], z)
        return hypsum([], [], [], br, bf, bc, z)
    if degree == (1, 1):
        ar, af, ac = parse_param(as[0])
        br, bf, bc = parse_param(bs[0])
        if ar and br:
            a, b = ar[0], br[0]
            return sum_hyp1f1_rat(a, b, z)
        return hypsum(ar, af, ac, br, bf, bc, z)
    if degree == (2, 1):
        return eval_hyp2f1(as[0],as[1],bs[0],z)
    ars, afs, acs, brs, bfs, bcs = [], [], [], [], [], []
    for a in as:
        r, f, c = parse_param(a)
        ars += r
        afs += f
        acs += c
    for b in bs:
        r, f, c = parse_param(b)
        brs += r
        bfs += f
        bcs += c
    return hypsum(ars, afs, acs, brs, bfs, bcs, z)

def hyp0f1(a, z):
    """Hypergeometric function 0F1. hyp0f1(a,z) is equivalent
    to hyper([], [a], z); see documentation for hyper() for more
    information."""
    return hyper([], [a], z)

def hyp1f1(a,b,z):
    """Hypergeometric function 1F1. hyp1f1(a,b,z) is equivalent
    to hyper([a], [b], z); see documentation for hyper() for more
    information."""
    return hyper([a], [b], z)

def hyp2f1(a,b,c,z):
    """Hypergeometric function 2F1. hyp2f1(a,b,c,z) is equivalent
    to hyper([a,b], [c], z); see documentation for hyper() for more
    information."""
    return hyper([a,b], [c], z)

@funcwrapper
def lower_gamma(a,z):
    """Lower incomplete gamma function gamma(a, z)"""
    # XXX: may need more precision
    return hyp1f1(1, 1+a, z) * z**a * exp(-z) / a

@funcwrapper
def upper_gamma(a,z):
    """Upper incomplete gamma function Gamma(a, z)"""
    return gamma(a) - lower_gamma(a, z)

def mpf_erf(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if x == fzero: return fzero
        if x == finf: return fone
        if x== fninf: return fnone
        return fnan
    size = exp + bc
    lg = math.log
    # The approximation erf(x) = 1 is accurate to > x^2 * log(e,2) bits
    if size > 3 and 2*(size-1) + 0.528766 > lg(prec,2):
        # TODO: interval rounding
        return [fone, fnone][sign]
    # erf(x) ~ x close to 0
    if size < -prec:
        # TODO: interval rounding
        return fdiv(fshift(x,1), fsqrt(fpi(prec+10), prec+10), prec, rnd)
    wp = prec + abs(size) + 20
    # Taylor series for erf, fixed-point summation
    t = abs(to_fixed(x, wp))
    t2 = (t*t) >> wp
    s, term, k = t, 12345, 1
    while term:
        t = ((t * t2) >> wp) // k
        term = t // (2*k+1)
        if k & 1:
            s -= term
        else:
            s += term
        k += 1
    s = (s << (wp+1)) // sqrt_fixed(pi_fixed(wp), wp)
    if sign:
        s = -s
    # 1-wp multiplies by two
    return from_man_exp(s, -wp, wp, rnd)

# If possible, we use the asymptotic series for erfc.
# This is an alternating divergent asymptotic series, so
# the error is at most equal to the first omitted term.
# Here we check if the smallest term is small enough
# for a given x and precision
def erfc_check_series(x, prec):
    n = to_int(x)
    if n**2 * 1.44 > prec:
        return True
    return False

def mpf_erfc(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if x == fzero: return fone
        if x == finf: return fzero
        if x == fninf: return ftwo
        return fnan
    wp = prec + 20
    regular_erf = sign or man+exp < 2
    if regular_erf or not erfc_check_series(x, wp):
        if regular_erf:
            return fsub(fone, mpf_erf(x, prec+5), prec, rnd)
        # 1-erf(x) ~ exp(-x^2), increase prec to deal with cancellation
        n = to_int(x)
        return fsub(fone, mpf_erf(x, prec + int(n**2*1.44) + 10), prec, rnd)
    s = term = MP_ONE << wp
    term_prev = 0
    t = (2 * to_fixed(x, wp) ** 2) >> wp
    k = 1
    while 1:
        term = ((term * (2*k - 1)) << wp) // t
        if k > 4 and term > term_prev or not term:
            break
        if k & 1:
            s -= term
        else:
            s += term
        term_prev = term
        #print k, to_str(from_man_exp(term, -wp, 50), 10)
        k += 1
    s = (s << wp) // sqrt_fixed(pi_fixed(wp), wp)
    s = from_man_exp(s, -wp, wp)
    z = fexp(fneg(fmul(x,x,wp),wp),wp)
    y = fdiv(fmul(z, s, wp), x, prec, rnd)
    return y

@funcwrapper
def erf(z):
    """Error function, erf(z)"""
    if z.imag:
        # XXX: may have issues with accuracy for large complex z
        return (2/sqrt(pi)*z) * sum_hyp1f1_rat((1,2),(3,2), -z**2)
    v = mpf_erf(z.real._mpf_, *prec_rounding)
    if isinstance(z, mpf):
        return make_mpf(v)
    else:
        return make_mpc((v, fzero))

@funcwrapper
def erfc(z):
    """Complementary error function, erfc(z) = 1-erf(z)"""
    if z.imag:
        # XXX: may have issues with accuracy for large complex z
        return 1-erf(z)
    v = mpf_erfc(z.real._mpf_, *prec_rounding)
    if isinstance(z, mpf):
        return make_mpf(v)
    else:
        return make_mpc((v, fzero))

@funcwrapper
def erfi(z):
    """Imaginary error function, erfi(z)"""
    return (2/sqrt(pi)*z) * sum_hyp1f1_rat((1,2),(3,2), z**2)

@funcwrapper
def npdf(x, mu=0, sigma=1):
    """
    npdf(x, mu=0, sigma=1) -- probability density function of a
    normal distribution with mean value mu and variance sigma^2.
    """
    sigma = convert_lossless(sigma)
    return exp(-(x-mu)**2/(2*sigma**2)) / (sigma*sqrt(2*pi))

@funcwrapper
def ncdf(x, mu=0, sigma=1):
    """
    ncdf(x, mu=0, sigma=1) -- cumulative distribution function of
    a normal distribution with mean value mu and variance sigma^2.
    """
    a = (x-mu)/(sigma*sqrt(2))
    if a < 0:
        return erfc(-a)/2
    else:
        return (1+erf(a))/2

@funcwrapper
def ei(z):
    """Exponential integral, Ei(z)"""
    if z == inf:
        return z
    if z == -inf:
        return -mpf(0)
    v = z*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1]],[],[],z) + \
        (log(z)-log(1/z))/2 + euler
    if isinstance(z, mpf) and z < 0:
        return v.real
    return v

@funcwrapper
def li(z):
    """Logarithmic integral, li(z)"""
    if not z:
        return z
    if z == 1:
        return -inf
    return ei(log(z))

@funcwrapper
def ci(z):
    """Cosine integral, Ci(z)"""
    if z == inf:
        return 1/z
    if not z:
        return -inf
    z2 = -(z/2)**2
    return euler + log(z) + \
        z2*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1],[3,2]],[],[],z2)

@funcwrapper
def si(z):
    """Sine integral, Si(z)"""
    if z == inf:
        return pi/2
    if z == -inf:
        return -pi/2
    z2 = -(z/2)**2
    return z*hypsum([[1,2]],[],[],[[3,2],[3,2]],[],[],z2)

@funcwrapper
def chi(z):
    """Hyperbolic cosine integral, Chi(z)"""
    if not z:
        return -inf
    z2 = (z/2)**2
    return euler + log(z) + \
        z2*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1],[3,2]],[],[],z2)

@funcwrapper
def shi(z):
    """Hyperbolic sine integral, Shi(z)"""
    z2 = (z/2)**2
    return z*hypsum([[1,2]],[],[],[[3,2],[3,2]],[],[],z2)

@funcwrapper
def fresnels(z):
    """Fresnel integral S, S(z)"""
    if z == inf:
        return mpf(0.5)
    if z == -inf:
        return mpf(-0.5)
    return pi*z**3/6*hypsum([[3,4]],[],[],[[3,2],[7,4]],[],[],-pi**2*z**4/16)

@funcwrapper
def fresnelc(z):
    """Fresnel integral C, C(z)"""
    if z == inf:
        return mpf(0.5)
    if z == -inf:
        return mpf(-0.5)
    return z*hypsum([[1,4]],[],[],[[1,2],[5,4]],[],[],-pi**2*z**4/16)

@funcwrapper
def airyai(z):
    """Airy function, Ai(z)"""
    if z == inf:
        return 1/z
    if z == -inf:
        return mpf(0)
    z3 = z**3 / 9
    a = sum_hyp0f1_rat((2,3), z3) / (cbrt(9) * gamma(mpf(2)/3))
    b = z * sum_hyp0f1_rat((4,3), z3) / (cbrt(3) * gamma(mpf(1)/3))
    return a - b

@funcwrapper
def airybi(z):
    """Airy function, Bi(z)"""
    if z == inf:
        return z
    if z == -inf:
        return mpf(0)
    z3 = z**3 / 9
    rt = nthroot(3, 6)
    a = sum_hyp0f1_rat((2,3), z3) / (rt * gamma(mpf(2)/3))
    b = z * rt * sum_hyp0f1_rat((4,3), z3) / gamma(mpf(1)/3)
    return a + b

@funcwrapper
def ellipe(m):
    """Complete elliptic integral of the second kind, E(m). Note that
    the argument is the parameter m = k^2, not the modulus k."""
    if m == 1:
        return m
    return pi/2 * sum_hyp2f1_rat((1,2),(-1,2),(1,1), m)

@funcwrapper
def ellipk(m):
    """Complete elliptic integral of the first kind, K(m). Note that
    the argument is the parameter m = k^2, not the modulus k."""
    # Poor implementation:
    # return pi/2 * sum_hyp2f1_rat((1,2),(1,2),(1,1), m)
    if m == 1:
        return inf
    if isnan(m):
        return m
    if isinf(m):
        return 1/m
    s = sqrt(m)
    a = (1-s)/(1+s)
    v = pi/4*(1+a)/agm(1,a)
    if isinstance(m, mpf) and m < 1:
        return v.real
    return v

# TODO: for complex a, b handle the branch cut correctly
@funcwrapper
def agm(a, b=1):
    """Arithmetic-geometric mean of a and b. Can be called with
    a single argument, computing agm(a,1) = agm(1,a)."""
    if not a or not b:
        return a*b
    weps = eps * 16
    half = mpf(0.5)
    while abs(a-b) > weps:
        a, b = (a+b)*half, (a*b)**half
    return a

@funcwrapper
def jacobi(n, a, b, x):
    """Jacobi polynomial P_n^(a,b)(x)."""
    return binomial(n+a,n) * hyp2f1(-n,1+n+a+b,a+1,(1-x)/2)

@funcwrapper
def legendre(n, x):
    """Legendre polynomial P_n(x)."""
    if isint(n):
        n = int(n)
    if x == -1:
        # TODO: hyp2f1 should handle this
        if x == int(x):
            return (-1)**(n + (n>=0)) * mpf(-1)
        return inf
    return hyp2f1(-n,n+1,1,(1-x)/2)

@funcwrapper
def chebyt(n, x):
    """Chebyshev polynomial of the first kind T_n(x)."""
    return hyp2f1(-n,n,0.5,(1-x)/2)

@funcwrapper
def chebyu(n, x):
    """Chebyshev polynomial of the second kind U_n(x)."""
    return (n+1) * hyp2f1(-n, n+2, 1.5, (1-x)/2)

# A Bessel function of the first kind of integer order, J_n(x), is
# given by the power series

#             oo
#             ___         k         2 k + n
#            \        (-1)     / x \
#    J_n(x) = )    ----------- | - |
#            /___  k! (k + n)! \ 2 /
#            k = 0

# Simplifying the quotient between two successive terms gives the
# ratio x^2 / (-4*k*(k+n)). Hence, we only need one full-precision
# multiplication and one division by a small integer per term.
# The complex version is very similar, the only difference being
# that the multiplication is actually 4 multiplies.

# In the general case, we have
# J_v(x) = (x/2)**v / v! * 0F1(v+1, (-1/4)*z**2)

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
        s = t = MP_ONE << prec
    else:
        s = t = (x**n // int_fac(n)) >> ((n-1)*prec + n)
    k = 1
    while t:
        t = ((t * x2) // (-4*k*(k+n))) >> prec
        s += t
        k += 1
    if negate:
        s = -s
    return make_mpf(from_man_exp(s, -prec, origprec, round_nearest))

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
        sre = tre = MP_ONE << prec
        sim = tim = MP_ZERO
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
    return make_mpc((re, im))

@funcwrapper
def jv(v, x):
    """Bessel function J_v(x)."""
    if isint(v):
        if isinstance(x, mpf):
            return mpf_jn_series(int(v), x._mpf_, mp.prec)
        if isinstance(x, mpc):
            return mpc_jn_series(int(v), (x.real._mpf_, x.imag._mpf_), mp.prec)
    hx = x/2
    return hx**v * hyp0f1(v+1, -hx**2) / factorial(v)

jn = jv

def j0(x):
    """Bessel function J_0(x)."""
    return jv(0, x)

def j1(x):
    """Bessel function J_1(x)."""
    return jv(1, x)

#---------------------------------------------------------------------------#
#                                                                           #
#                               Miscellaneous                               #
#                                                                           #
#---------------------------------------------------------------------------#

# TODO: remove this; use log_int_fixed internally instead
def log_range():
    """Generate log(2), log(3), log(4), ..."""
    prec = mp.prec + 20
    one = 1 << prec
    L = log2_fixed(prec)
    p = 2
    while 1:
        yield mpf((L, -prec))
        s = 0
        u = one
        k = 1
        a = (2*p+1)**2
        while u:
            s += u // k
            u //= a
            k += 2
        L += 2*s//(2*p+1)
        p += 1

@funcwrapper
def lambertw(z, k=0, approx=None):
    """
    lambertw(z,k) gives the kth branch of the Lambert W function W(z),
    defined as the kth solution of z = W(z)*exp(W(z)).

    lambertw(z) == lambertw(z, k=0) gives the principal branch
    value (0th branch solution), which is real for z > -1/e .

    The k = -1 branch is real for -1/e < z < 0. All branches except
    k = 0 have a logarithmic singularity at 0.

    The definition, implementation and choice of branches is based
    on Corless et al, "On the Lambert W function", Adv. Comp. Math. 5
    (1996) 329-359, available online here:
    http://www.apmaths.uwo.ca/~djeffrey/Offprints/W-adv-cm.pdf

    TODO: use a series expansion when extremely close to the branch point
    at -1/e and make sure that the proper branch is chosen there
    """
    if isnan(z):
        return z
    mp.prec += 20
    # We must be extremely careful near the singularities at -1/e and 0
    u = exp(-1)
    if abs(z) <= u:
        if not z:
            # w(0,0) = 0; for all other branches we hit the pole
            if not k:
                return z
            return -inf
        if not k:
            w = z
        # For small real z < 0, the -1 branch behaves roughly like log(-z)
        elif k == -1 and not z.imag and z.real < 0:
            w = log(-z)
        # Use a simple asymptotic approximation.
        else:
            w = log(z)
            # The branches are roughly logarithmic. This approximation
            # gets better for large |k|; need to check that this always
            # works for k ~= -1, 0, 1.
            if k: w += k * 2*pi*j
    else:
        if z == inf: return z
        if z == -inf: return nan
        # Simple asymptotic approximation as above
        w = log(z)
        if k: w += k * 2*pi*j
    # Use Halley iteration to solve w*exp(w) = z
    two = mpf(2)
    weps = ldexp(eps, 15)
    for i in xrange(100):
        ew = exp(w)
        wew = w*ew
        wewz = wew-z
        wn = w - wewz/(wew+ew-(w+two)*wewz/(two*w+two))
        if abs(wn-w) < weps*abs(wn):
            return wn
        else:
            w = wn
    print "Warning: Lambert W iteration failed to converge:", z
    return wn

"""
Miscellaneous special functions
"""

from libmpf import *
from libmpc import *
from mptypes import *
from mptypes import constant, funcwrapper

import libhyper


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

mpq_1 = _mpq((1,1))
mpq_0 = _mpq((0,1))

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

def _as_num(x):
    if isinstance(x, list):
        return _mpq(x)
    return x

def hypsum(ar, af, ac, br, bf, bc, x):
    prec, rnd = prec_rounding
    if hasattr(x, "_mpf_") and not (ac or bc):
        v = libhyper.hypsum_internal(ar, af, ac, br, bf, bc, x._mpf_, None, prec, rnd)
        return make_mpf(v)
    else:
        if hasattr(x, "_mpc_"):
            re, im = x._mpc_
        else:
            re, im = x._mpf_, fzero
        v = libhyper.hypsum_internal(ar, af, ac, br, bf, bc, re, im, prec, rnd)
        return make_mpc(v)

def eval_hyp2f1(a,b,c,z):
    prec, rnd = prec_rounding
    ar, af, ac = parse_param(a)
    br, bf, bc = parse_param(b)
    cr, cf, cc = parse_param(c)
    absz = abs(z)
    if absz == 1:
        # TODO: determine whether it actually does, and otherwise
        # return infinity instead
        print "Warning: 2F1 might not converge for |z| = 1"
    if absz <= 1:
        # All rational
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
        h1 = eval_hyp2f1(a, mpq_1-c+a, mpq_1-b+a, 1/z)
        h2 = eval_hyp2f1(b, mpq_1-c+b, mpq_1-a+b, 1/z)
        #s1 = G(c)*G(b-a)/G(b)/G(c-a) * (-z)**(-a) * h1
        #s2 = G(c)*G(a-b)/G(a)/G(c-b) * (-z)**(-b) * h2
        f1 = gammaprod([c,b-a],[b,c-a])
        f2 = gammaprod([c,a-b],[a,c-b])
        s1 = f1 * (-z)**(mpq_0-a) * h1
        s2 = f2 * (-z)**(mpq_0-b) * h2
        v = s1 + s2
    finally:
        mp.prec = orig
    return +v

def sum_hyp0f1_rat(a, z):
    prec, rnd = prec_rounding
    if hasattr(z, "_mpf_"):
        return make_mpf(libhyper.mpf_hyp0f1_rat(a, z._mpf_, prec, rnd))
    else:
        return make_mpc(libhyper.mpc_hyp0f1_rat(a, z._mpc_, prec, rnd))

def sum_hyp1f1_rat(a, b, z):
    prec, rnd = prec_rounding
    if hasattr(z, "_mpf_"):
        return make_mpf(libhyper.mpf_hyp1f1_rat(a, b, z._mpf_, prec, rnd))
    else:
        return make_mpc(libhyper.mpc_hyp1f1_rat(a, b, z._mpc_, prec, rnd))

def sum_hyp2f1_rat(a, b, c, z):
    prec, rnd = prec_rounding
    if hasattr(z, "_mpf_"):
        return make_mpf(libhyper.mpf_hyp2f1_rat(a, b, c, z._mpf_, prec, rnd))
    else:
        return make_mpc(libhyper.mpc_hyp2f1_rat(a, b, c, z._mpc_, prec, rnd))


#---------------------------------------------------------------------------#
#                      And now the user-friendly versions                   #
#---------------------------------------------------------------------------#

def hyper(a_s, b_s, z):
    """
    Hypergeometric function pFq,

          [ a_1, a_2, ..., a_p |    ]
      pFq [                    |  z ]
          [ b_1, b_2, ..., b_q |    ]

    The parameter lists a_s and b_s may contain real or complex numbers.
    Exact rational parameters can be given as tuples (p, q).
    """
    p = len(a_s)
    q = len(b_s)
    z = convert_lossless(z)
    degree = p, q
    if degree == (0, 1):
        br, bf, bc = parse_param(b_s[0])
        if br:
            return sum_hyp0f1_rat(br[0], z)
        return hypsum([], [], [], br, bf, bc, z)
    if degree == (1, 1):
        ar, af, ac = parse_param(a_s[0])
        br, bf, bc = parse_param(b_s[0])
        if ar and br:
            a, b = ar[0], br[0]
            return sum_hyp1f1_rat(a, b, z)
        return hypsum(ar, af, ac, br, bf, bc, z)
    if degree == (2, 1):
        return eval_hyp2f1(a_s[0], a_s[1], b_s[0], z)
    ars, afs, acs, brs, bfs, bcs = [], [], [], [], [], []
    for a in a_s:
        r, f, c = parse_param(a)
        ars += r
        afs += f
        acs += c
    for b in b_s:
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

erf = mpfunc("erf", libhyper.mpf_erf, libhyper.mpc_erf,
    "Error function, erf(z)")
erfc = mpfunc("erfc", libhyper.mpf_erfc, libhyper.mpc_erfc,
    "Complementary error function, erfc(z) = 1-erf(z)")

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

@funcwrapper
def jv(v, x):
    """Bessel function J_v(x)."""
    if isint(v):
        if isinstance(x, mpf):
            return make_mpf(libhyper.mpf_besseljn(int(v), x._mpf_, mp.prec))
        if isinstance(x, mpc):
            return make_mpc(libhyper.mpc_besseljn(int(v), x._mpc_, mp.prec))
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

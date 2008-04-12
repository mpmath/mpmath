"""
Hypergeometric functions.
"""

from mpmath.lib import from_man_exp, to_fixed
from mpmath.mptypes import mp, mpf, mpc, make_mpf, make_mpc, \
    convert_lossless, inf, pi, extraprec, eps, sqrt

from factorials import gammaquot, binomial

import operator

"""
TODO:
  * By counting the number of multiplications vs divisions,
    the bit size of p can be kept around wp instead of growing
    it to n*wp for some (possible large) n
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
    rnd = mp.rounding[0]
    wp = prec + 25

    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        y = 0
    else:
        have_complex = 1
        x, y = x._mpc_
        x = to_fixed(x, wp)
        y = to_fixed(y, wp)

    sre = pre = one = 1 << wp
    sim = pim = 0

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

def sum_hyp1f1_rat((ap, aq), (bp, bq), x):
    """Sum 1F1 for rational a, b, c. x must be mpf or mpc."""
    prec = mp.prec
    rnd = mp.rounding[0]
    wp = prec + 25
    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        s = p = 1 << wp
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
        sre = pre = 1 << wp
        sim = pim = 0
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
    rnd = mp.rounding[0]
    wp = prec + 25
    if isinstance(x, mpf):
        x = to_fixed(x._mpf_, wp)
        s = p = 1 << wp
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
        f1 = gammaquot([c,b-a],[b,c-a])
        f2 = gammaquot([c,a-b],[a,c-b])
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

def funcwrapper(f):
    def g(z):
        orig = mp.prec
        rnd = mp.rounding[0]
        try:
            z = convert_lossless(z)
            mp.prec = orig + 10
            v = f(z)
        finally:
            mp.prec = orig
        return +v
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g

@funcwrapper
def erf(z):
    """Error function, erf(z)"""
    return (2/sqrt(pi)*z) * sum_hyp1f1_rat((1,2),(3,2), -z**2)

@funcwrapper
def ellipk(m):
    """Complete elliptic integral of the first kind, K(m). Note that
    the argument is the parameter m = k^2, not the modulus k."""
    if m == 1:
        return inf
    return pi/2 * sum_hyp2f1_rat((1,2),(1,2),(1,1), m)

@funcwrapper
def ellipe(m):
    """Complete elliptic integral of the second kind, E(m). Note that
    the argument is the parameter m = k^2, not the modulus k."""
    if m == 1:
        return m
    return pi/2 * sum_hyp2f1_rat((1,2),(-1,2),(1,1), m)

# TODO: for complex a, b handle the branch cut correctly
@extraprec(15, normalize_output=True)
def agm(a, b):
    """Arithmetic-geometric mean of a and b."""
    a = convert_lossless(a)
    b = convert_lossless(b)
    if not a or not b:
        return a*b
    weps = eps * 16
    half = mpf(0.5)
    while abs(a-b) > weps:
        a, b = (a+b)*half, (a*b)**half
    return a

def jacobi(n, a, b, x):
    """Jacobi polynomial P_n^(a,b)(x)."""
    orig = mp.prec
    try:
        mp.prec = orig + 15
        x = convert_lossless(x)
        v = binomial(n+a,n) * hyp2f1(-n,1+n+a+b,a+1,(1-x)/2)
    finally:
        mp.prec = orig
    return +v

def legendre(n, x):
    """Legendre polynomial P_n(x)."""
    orig = mp.prec
    try:
        mp.prec = orig + 15
        x = convert_lossless(x)
        if not isinstance(n, (int, long)):
            n = convert_lossless(n)
        if x == -1:
            # TODO: hyp2f1 should handle this
            if x == int(x):
                return (-1)**(n + (n>=0)) * mpf(-1)
            return inf
        v = hyp2f1(-n,n+1,1,(1-x)/2)
    finally:
        mp.prec = orig
    return +v

def chebyt(n, x):
    """Chebyshev polynomial of the first kind T_n(x)."""
    orig = mp.prec
    try:
        mp.prec = orig + 15
        x = convert_lossless(x)
        v = hyp2f1(-n,n,0.5,(1-x)/2)
    finally:
        mp.prec = orig
    return +v

def chebyu(n, x):
    """Chebyshev polynomial of the second kind U_n(x)."""
    orig = mp.prec
    try:
        mp.prec = orig + 15
        x = convert_lossless(x)
        v = (n+1) * hyp2f1(-n,n+2,1.5,(1-x)/2)
    finally:
        mp.prec = orig
    return +v

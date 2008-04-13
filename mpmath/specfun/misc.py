# XXX: this file needs cleanup

__docformat__ = 'plaintext'

from mpmath.lib import *
from mpmath.libmpc import *
from mpmath.mptypes import *

#----------------------------------------------------------------------
# Factorial related functions
#

# For internal use
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

"""
We compute the gamma function using Spouge's approximation

    x! = (x+a)**(x+1/2) * exp(-x-a) * [c_0 + S(x) + eps]

where S(x) is the sum of c_k/(x+k) from k = 1 to a-1 and the coefficients
are given by

  c_0 = sqrt(2*pi)

         (-1)**(k-1)
  c_k =  ----------- (a-k)**(k-1/2) exp(-k+a),  k = 1,2,...,a-1
          (k - 1)!

Due to an inequality proved by Spouge, if we choose a = int(1.26*n), the
error eps is less than 10**-n for any x in the right complex half-plane
(assuming a > 2). In practice, it seems that a can be chosen quite a bit
lower still (30-50%); this possibility should be investigated.

Reference:
John L. Spouge, "Computation of the gamma, digamma, and trigamma
functions", SIAM Journal on Numerical Analysis 31 (1994), no. 3, 931-944.
"""

spouge_cache = {}

def calc_spouge_coefficients(a, prec):
    wp = prec + int(a*1.4)
    c = [0] * a
    # b = exp(a-1)
    b = fexp(from_int(a-1), wp)
    # e = exp(1)
    e = fexp(fone, wp)
    # sqrt(2*pi)
    sq2pi = fsqrt(fshift(fpi(wp), 1), wp)
    c[0] = to_fixed(sq2pi, prec)
    for k in xrange(1, a):
        # c[k] = ((-1)**(k-1) * (a-k)**k) * b / sqrt(a-k)
        term = fmuli(b, ((-1)**(k-1) * (a-k)**k), wp)
        term = fdiv(term, fsqrt(from_int(a-k), wp), wp)
        c[k] = to_fixed(term, prec)
        # b = b / (e * k)
        b = fdiv(b, fmul(e, from_int(k), wp), wp)
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

def mpf_gamma(x, prec, rounding=round_fast, p1=1):
    sign, man, exp, bc = x
    if exp >= 0:
        if sign or (p1 and not man):
            raise ValueError("gamma function pole")
        # A direct factorial is fastest
        if exp + bc <= 10:
            return from_int(int_fac((man<<exp)-p1), prec, rounding)
    wp = prec + 15
    if p1:
        x = fsub(x, fone, wp)
    # x < 0.25
    if sign or exp+bc < -1:
        # gamma = pi / (sin(pi*x) * gamma(1-x))
        wp += 15
        pi = fpi(wp)
        pix = fmul(x, pi, wp)
        t = fsin(pix, wp)
        g = mpf_gamma(fsub(fone, x, wp), wp)
        return fdiv(pix, fmul(t, g, wp), prec, rounding)
    sprec, a, c = get_spouge_coefficients(wp)
    s = spouge_sum_real(x, sprec, a, c)
    # gamma = exp(log(x+a)*(x+0.5) - xpa) * s
    xpa = fadd(x, from_int(a), wp)
    logxpa = flog(xpa, wp)
    xph = fadd(x, fhalf, wp)
    t = fsub(fmul(logxpa, xph, wp), xpa, wp)
    t = fmul(fexp(t, wp), s, prec, rounding)
    return t

def mpc_gamma(x, prec, rounding=round_fast, p1=1):
    re, im = x
    if im == fzero:
        return mpf_gamma(re, prec, rounding, p1), fzero
    wp = prec + 25
    sign, man, exp, bc = re
    if p1:
        re = fsub(re, fone, wp)
        x = re, im
    if sign or exp+bc < -1:
        # Reflection formula
        wp += 15
        pi = fpi(wp), fzero
        pix = mpc_mul(x, pi, wp)
        t = mpc_sin(pix, wp)
        u = mpc_sub(mpc_one, x, wp)
        g = mpc_gamma(u, wp)
        w = mpc_mul(t, g, wp)
        return mpc_div(pix, w, wp)
    sprec, a, c = get_spouge_coefficients(wp)
    s = spouge_sum_complex(re, im, sprec, a, c)
    # gamma = exp(log(x+a)*(x+0.5) - xpa) * s
    repa = fadd(re, from_int(a), wp)
    logxpa = mpc_log((repa, im), wp)
    reph = fadd(re, fhalf, wp)
    t = mpc_sub(mpc_mul(logxpa, (reph, im), wp), (repa, im), wp)
    t = mpc_mul(mpc_exp(t, wp), s, prec, rounding)
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

def isnpint(x):
    if not x:
        return True
    if isinstance(x, mpf):
        sign, man, exp, bc = x._mpf_
        return sign and exp >= 0
    if isinstance(x, mpc):
        return not x.imag and isnpint(x.real)

def gammaquot(a, b):
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
    return gammaquot([n+1], [k+1, n-k+1])

def rf(x, n):
    """Rising factorial (Pochhammer symbol), x_(n)"""
    return gammaquot([x+n], [x])

def ff(x, n):
    """Falling factorial, x_(n)"""
    return gammaquot([x+1], [x-n+1])



#---------------------------------------------------------------------------#
#                                                                           #
#                         Riemann zeta function                             #
#                                                                           #
#---------------------------------------------------------------------------#

"""
We use zeta(s) = eta(s) * (1 - 2**(1-s)) and Borwein's approximation
                  n-1
                  ___       k
             -1  \      (-1)  (d_k - d_n)
  eta(s) ~= ----  )     ------------------
             d_n /___              s
                 k = 0      (k + 1)
where
             k
             ___                i
            \     (n + i - 1)! 4
  d_k  =  n  )    ---------------.
            /___   (n - i)! (2i)!
            i = 0

If s = a + b*I, the absolute error for eta(s) is bounded by

    3 (1 + 2|b|)
    ------------ * exp(|b| pi/2)
               n
    (3+sqrt(8))

Disregarding the linear term, we have approximately,

  log(err) ~= log(exp(1.58*|b|)) - log(5.8**n)
  log(err) ~= 1.58*|b| - log(5.8)*n
  log(err) ~= 1.58*|b| - 1.76*n
  log2(err) ~= 2.28*|b| - 2.54*n

So for p bits, we should choose n > (p + 2.28*|b|) / 2.54.

Reference:
Peter Borwein, "An Efficient Algorithm for the Riemann Zeta Function"
http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P117.ps

http://en.wikipedia.org/wiki/Dirichlet_eta_function
"""

d_cache = {}

def zeta_coefs(n):
    if n in d_cache:
        return d_cache[n]
    ds = [0] * (n+1)
    d = 1
    s = ds[0] = 1
    for i in range(1, n+1):
        d = d * 4 * (n+i-1) * (n-i+1)
        d //= ((2*i) * ((2*i)-1))
        s += d
        ds[i] = s
    d_cache[n] = ds
    return ds

# Integer logarithms
_log_cache = {}

def _logk(k):
    p = mp.prec
    if k in _log_cache and _log_cache[k][0] >= p:
        return +_log_cache[k][1]
    else:
        x = log(k)
        _log_cache[k] = (p, x)
        return x

@extraprec(10, normalize_output=True)
def zeta(s):
    """Returns the Riemann zeta function of s."""
    s = convert_lossless(s)
    if s.real < 0:
        # Reflection formula (XXX: gets bad around the zeros)
        return 2**s * pi**(s-1) * sin(pi*s/2) * gamma(1-s) * zeta(1-s)
    else:
        p = mp.prec
        n = int((p + 2.28*abs(float(mpc(s).imag)))/2.54) + 3
        d = zeta_coefs(n)
        if isinstance(s, mpf) and s == int(s):
            sint = int(s)
            t = 0
            for k in range(n):
                t += (((-1)**k * (d[k] - d[n])) << p) // (k+1)**sint
            return (mpf((t, -p)) / -d[n]) / (1 - mpf(2)**(1-sint))
        else:
            t = mpf(0)
            for k in range(n):
                t += (-1)**k * mpf(d[k]-d[n]) * exp(-_logk(k+1)*s)
            return (t / -d[n]) / (mpf(1) - exp(log(2)*(1-s)))
    

@extraprec(5, normalize_output=True)
def bernoulli(n):
    """nth Bernoulli number, B_n"""
    if n == 1:
        return mpf(-0.5)
    if n & 1:
        return mpf(0)
    m = n // 2
    return (-1)**(m-1) * 2 * factorial(n) / (2*pi)**n * zeta(n)

from mpmath.lib import fmuli, fmul, fpos, fdiv, fadd, fdivi, fsub, from_int, round_down, fone, fzero
from mpmath.mptypes import make_mpf

# For sequential computation of Bernoulli numbers, we use Ramanujan's formula

#                            / n + 3 \
#   B   =  (A(n) - S(n))  /  |       |
#    n                       \   n   /

# where A(n) = (n+3)/3 when n = 0 or 2 (mod 6), A(n) = -(n+3)/6
# when n = 4 (mod 6), and

#          [n/6]
#           ___
#          \      /  n + 3  \
#   S(n) =  )     |         | * B
#          /___   \ n - 6*k /    n-6*k
#          k = 1

def bernoulli2n():
    """Generates B(2), B(4), B(6), ..."""
    oprec = mp.prec
    rounding = mp.rounding[0]
    prec = oprec + 30
    computed = {0:fone}
    m, bin1, bin = 2, 1, 10
    f3 = from_int(3)
    f6 = from_int(6)
    while 1:
        case = m % 6
        s = fzero
        if m < 6: a = 0
        else:     a = bin1
        for j in xrange(1, m//6+1):
            s = fadd(s, fmuli(computed[m-6*j], a, prec), prec)
            # Inner binomial coefficient
            j6 = 6*j
            a *= ((m-5-j6)*(m-4-j6)*(m-3-j6)*(m-2-j6)*(m-1-j6)*(m-j6))
            a //= ((4+j6)*(5+j6)*(6+j6)*(7+j6)*(8+j6)*(9+j6))
        if case == 0: b = fdivi(m+3, f3, prec)
        if case == 2: b = fdivi(m+3, f3, prec)
        if case == 4: b = fdivi(-m-3, f6, prec)
        b = fdiv(fsub(b, s, prec), from_int(bin), prec)
        computed[m] = b
        yield make_mpf(fpos(b, oprec, rounding))
        m += 2
        bin = bin * ((m+2)*(m+3)) // (m*(m-1))
        if m > 6: bin1 = bin1 * ((2+m)*(3+m)) // ((m-7)*(m-6))

from mpmath.lib import log2_fixed

def logk():
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

@extraprec(30, normalize_output=True)
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
    z = convert_lossless(z)
    if isnan(z):
        return z
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

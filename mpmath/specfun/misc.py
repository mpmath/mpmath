# XXX: this file needs cleanup

__docformat__ = 'plaintext'

from mpmath.mptypes import mpnumeric, mpf, mpc, pi, exp, log, sqrt, sin,\
    power, extraprec, mp
from mpmath.lib import bitcount, to_fixed, from_man_exp, round_nearest, fmuli

from factorials import gamma, factorial

def _fix(x, prec):
    return to_fixed(x._mpf_, prec)

def _re(s):
    if isinstance(s, mpf):
        return s
    return s.real

def _im(s):
    if isinstance(s, mpf):
        return s
    return s.imag

#---------------------------------------------------------------------------#
#                                                                           #
#                       Incomplete gamma functions                          #
#                                                                           #
#---------------------------------------------------------------------------#

"""
We compute the lower incomplete gamma function g(a,z) using the formula
g(a,z) = z**a * exp(-z) * S(a,z) / a, where
                 oo
                ___            k
               \              z
  S(a,z) = 1 +  )     ------------------.
               /___   (a+1)(a+2)...(a+k)
               k = 1

Then, in turn, various functions such as erf and exponential integrals
can be computed from the incomplete gamma function.
"""

def _lower_gamma_series(are, aim, zre, zim, prec):
    are = _fix(are, prec)
    aim = _fix(aim, prec)
    zre = _fix(zre, prec)
    zim = _fix(zim, prec)
    one = 1 << prec
    cre = sre = one
    cim = sim = 0
    while abs(cre) > 3 or abs(cim) > 3:
        # c = (c * z) << prec
        cre, cim = (cre*zre-cim*zim)>>prec, (cim*zre+cre*zim)>>prec
        # c = c / (a+k)
        are += one
        mag = ((are**2 + aim**2) >> prec)
        cre, cim = (cre*are + cim*aim)//mag, (cim*are - cre*aim)//mag
        sre += cre
        sim += cim
    sre = mpf((sre, -prec))
    sim = mpf((sim, -prec))
    return mpc(sre, sim)

@extraprec(15, normalize_output=True)
def lower_gamma(a, z):
    """Returns the lower incomplete gamma function gamma(a, z)"""
    # XXX: may need more precision
    a = mpc(a)
    z = mpc(z)
    s = _lower_gamma_series(a.real, a.imag, z.real, z.imag, mp.prec)
    return exp(log(z)*a) * exp(-z) * s / a

@extraprec(10, normalize_output=True)
def upper_gamma(a, z):
    """Returns the upper incomplete gamma function Gamma(a, z)"""
    return gamma(a) - lower_gamma(a, z)

'''
# Using hyper implementation instead

def erf(x):
    """Returns the error function of x."""
    x = mpnumeric(x)
    if x == 0:
        return x
    w = mpc(x)
    if w.real < 0:
        if isinstance(x, mpf):
            return -erf(-x)
        return -erf(-w)
    oldprec = mp.prec
    mp.prec += 10

    y = lower_gamma(0.5, w**2) / sqrt(pi)
    if _re(x) == 0 and _im(x) < 0:
        y = -y

    if isinstance(x, mpf):
        y = y.real

    mp.prec = oldprec
    return +y
'''


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

_d_cache = {}

def _zeta_coefs(n):
    if n in _d_cache:
        return _d_cache[n]
    ds = [0] * (n+1)
    d = 1
    s = ds[0] = 1
    for i in range(1, n+1):
        d = d * 4 * (n+i-1) * (n-i+1)
        d //= ((2*i) * ((2*i)-1))
        s += d
        ds[i] = s
    _d_cache[n] = ds
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

def zeta(s):
    """Returns the Riemann zeta function of s."""
    oldprec = mp.prec
    mp.prec += 10
    s = mpnumeric(s)
    if _re(s) < 0:
        # Reflection formula (XXX: gets bad around the zeros)
        y = power(2, s) * power(pi, s-1) * sin(pi*s/2) * gamma(1-s) * zeta(1-s)
    else:
        p = mp.prec
        n = int((p + 2.28*abs(float(mpc(s).imag)))/2.54) + 3
        d = _zeta_coefs(n)
        if isinstance(s, mpf) and s == int(s):
            sint = int(s)
            t = 0
            for k in range(n):
                t += (((-1)**k * (d[k] - d[n])) << p) // (k+1)**sint
            y = (mpf((t, -p)) / -d[n]) / (mpf(1) - mpf(2)**(1-sint))
        else:
            t = mpf(0)
            for k in range(n):
                t += (-1)**k * mpf(d[k]-d[n]) * exp(-_logk(k+1)*s)
            y = (t / -d[n]) / (mpf(1) - exp(log(2)*(1-s)))
    mp.prec = oldprec
    return +y

@extraprec(5, normalize_output=True)
def bernoulli(n):
    """Returns the nth Bernoulli number"""
    if n == 1:
        return mpf(-0.5)
    if n & 1:
        return mpf(0)
    m = n // 2
    return (-1)**(m-1) * 2 * factorial(n) / (2*pi)**(n) * zeta(n)

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
    """Generate B(2), B(4), B(6), ..."""
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

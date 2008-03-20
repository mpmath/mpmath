__docformat__ = 'plaintext'

from mpmath.mptypes import mpnumeric, mpf, mpc, pi, euler, exp, log, sqrt, sin,\
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

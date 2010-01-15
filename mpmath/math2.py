"""
This module complements the math and cmath builtin modules by providing
fast machine precision versions of some additional functions (gamma, ...).
"""

import operator
import math
import cmath

# Irrational constants
pi = 3.1415926535897932385
e = 2.7182818284590452354
sqrt2 = 1.4142135623730950488
sqrt5 = 2.2360679774997896964
phi = 1.6180339887498948482
ln2 = 0.69314718055994530942
ln10 = 2.302585092994045684
euler = 0.57721566490153286061
catalan = 0.91596559417721901505
khinchin = 2.6854520010653064453
apery = 1.2020569031595942854

def _mathfun_real(f_real, f_complex):
    def f(x, **kwargs):
        if type(x) is float:
            return f_real(x)
        if type(x) is complex:
            return f_complex(x)
        try:
            x = float(x)
            return f_real(x)
        except (TypeError, ValueError):
            x = complex(x)
            return f_complex(x)
    f.__name__ = f_real.__name__
    return f

def _mathfun(f_real, f_complex):
    def f(x, **kwargs):
        if type(x) is complex:
            return f_complex(x)
        try:
            return f_real(float(x))
        except (TypeError, ValueError):
            return f_complex(complex(x))
    f.__name__ = f_real.__name__
    return f

def _mathfun_n(f_real, f_complex):
    def f(*args, **kwargs):
        try:
            return f_real(*(float(x) for x in args))
        except (TypeError, ValueError):
            return f_complex(*(complex(x) for x in args))
    f.__name__ = f_real.__name__
    return f

pow = _mathfun_n(operator.pow, lambda x, y: complex(x)**y)
log = _mathfun_n(math.log, cmath.log)
sqrt = _mathfun(math.sqrt, cmath.sqrt)
exp = _mathfun_real(math.exp, cmath.exp)

cos = _mathfun_real(math.cos, cmath.cos)
sin = _mathfun_real(math.sin, cmath.sin)
tan = _mathfun_real(math.tan, cmath.tan)

acos = _mathfun(math.acos, cmath.acos)
asin = _mathfun(math.asin, cmath.asin)
atan = _mathfun_real(math.atan, cmath.atan)

cosh = _mathfun_real(math.cosh, cmath.cosh)
sinh = _mathfun_real(math.sinh, cmath.sinh)
tanh = _mathfun_real(math.tanh, cmath.tanh)

floor = _mathfun_real(math.floor,
    lambda z: complex(math.floor(z.real), math.floor(z.imag)))
ceil = _mathfun_real(math.ceil,
    lambda z: complex(math.ceil(z.real), math.ceil(z.imag)))

def nthroot(x, n):
    r = 1./n
    try:
        return x ** r
    except ValueError:
        return complex(x) ** r

def cos_sin(x, **kwargs):
    if type(x) is complex:
        return cmath.cos(x), cmath.sin(x)
    else:
        return math.cos(x), math.sin(x)

def _sinpi_real(x):
    n, r = divmod(x, 0.5)
    r *= pi
    n %= 4
    if n == 0: return math.sin(r)
    if n == 1: return math.cos(r)
    if n == 2: return -math.sin(r)
    if n == 3: return -math.cos(r)

def _cospi_real(x):
    n, r = divmod(x, 0.5)
    r *= pi
    n %= 4
    if n == 0: return math.cos(r)
    if n == 1: return -math.sin(r)
    if n == 2: return -math.cos(r)
    if n == 3: return math.sin(r)

def _sinpi_complex(z):
    n, r = divmod(z.real, 0.5)
    z = pi*complex(r, z.imag)
    n %= 4
    if n == 0: return cmath.sin(z)
    if n == 1: return cmath.cos(z)
    if n == 2: return -cmath.sin(z)
    if n == 3: return -cmath.cos(z)

def _cospi_complex(z):
    n, r = divmod(z.real, 0.5)
    z = pi*complex(r, z.imag)
    n %= 4
    if n == 0: return cmath.cos(z)
    if n == 1: return -cmath.sin(z)
    if n == 2: return -cmath.cos(z)
    if n == 3: return cmath.sin(z)

cospi = _mathfun_real(_cospi_real, _cospi_complex)
sinpi = _mathfun_real(_sinpi_real, _sinpi_complex)

def tanpi(x):
    return sinpi(x) / cospi(x)

def cotpi(x):
    return cospi(x) / sinpi(x)

INF = 1e300*1e300
NINF = -INF
NAN = INF-INF
EPS = 2.2204460492503131e-16

_exact_gamma = (INF, 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0,
  362880.0, 3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0,
  1307674368000.0, 20922789888000.0, 355687428096000.0, 6402373705728000.0,
  121645100408832000.0, 2432902008176640000.0)

_max_exact_gamma = len(_exact_gamma)-1

# Lanczos coefficients used by the GNU Scientific Library
_lanczos_g = 7
_lanczos_p = (0.99999999999980993, 676.5203681218851, -1259.1392167224028,
     771.32342877765313, -176.61502916214059, 12.507343278686905,
     -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7)
 
def _gamma_real(x):
    _intx = int(x)
    if _intx == x:
        if _intx <= 0:
            #return (-1)**_intx * INF
            raise ZeroDivisionError("gamma function pole")
        if _intx <= _max_exact_gamma:
            return _exact_gamma[_intx]
    if x < 0.5:
        # TODO: sinpi
        return pi / (_sinpi_real(x)*_gamma_real(1-x))
    else:
        x -= 1.0
        r = _lanczos_p[0]
        for i in range(1, _lanczos_g+2):
            r += _lanczos_p[i]/(x+i)
        t = x + _lanczos_g + 0.5
        return 2.506628274631000502417 * t**(x+0.5) * math.exp(-t) * r

def _gamma_complex(x):
    if not x.imag:
        return complex(_gamma_real(x.real))
    if x.real < 0.5:
        # TODO: sinpi
        return pi / (_sinpi_complex(x)*_gamma_complex(1-x))
    else:
        x -= 1.0
        r = _lanczos_p[0]
        for i in range(1, _lanczos_g+2):
            r += _lanczos_p[i]/(x+i)
        t = x + _lanczos_g + 0.5
        return 2.506628274631000502417 * t**(x+0.5) * cmath.exp(-t) * r

gamma = _mathfun_real(_gamma_real, _gamma_complex)

def factorial(x):
    return gamma(x+1.0)

def arg(x):
    if type(x) is float:
        return math.atan2(0.0,x)
    return math.atan2(x.imag,x.real)

# XXX: broken for negatives
def loggamma(x):
    if type(x) not in (float, complex):
        try:
            x = float(x)
        except (ValueError, TypeError):
            x = complex(x)
    p = 0.
    while abs(x) < 11:
        p -= log(x)
        x += 1.0
    s = 0.918938533204672742 + (x-0.5)*log(x) - x
    r = 1./x
    r2 = r*r
    s += 0.083333333333333333333*r; r *= r2
    s += -0.0027777777777777777778*r; r *= r2
    s += 0.00079365079365079365079*r; r *= r2
    s += -0.0005952380952380952381*r; r *= r2
    s += 0.00084175084175084175084*r; r *= r2
    s += -0.0019175269175269175269*r; r *= r2
    s += 0.0064102564102564102564*r; r *= r2
    s += -0.02955065359477124183*r
    return s + p

_psi_coeff = [
0.083333333333333333333,
-0.0083333333333333333333,
0.003968253968253968254,
-0.0041666666666666666667,
0.0075757575757575757576,
-0.021092796092796092796,
0.083333333333333333333,
-0.44325980392156862745,
3.0539543302701197438,
-26.456212121212121212]

def _digamma_real(x):
    _intx = int(x)
    if _intx == x:
        if _intx <= 0:
            raise ZeroDivisionError("polygamma pole")
    if x < 0.5:
        x = 1.0-x
        s = pi*cotpi(x)
    else:
        s = 0.0
    while x < 10.0:
        s -= 1.0/x
        x += 1.0
    x2 = x**-2
    t = x2
    for c in _psi_coeff:
        s -= c*t
        if t < 1e-20:
            break
        t *= x2
    return s + math.log(x) - 0.5/x

def _digamma_complex(x):
    if not x.imag:
        return complex(_digamma_real(x.real))
    if x.real < 0.5:
        x = 1.0-x
        s = pi*cotpi(x)
    else:
        s = 0.0
    while abs(x) < 10.0:
        s -= 1.0/x
        x += 1.0
    x2 = x**-2
    t = x2
    for c in _psi_coeff:
        s -= c*t
        if abs(t) < 1e-20:
            break
        t *= x2
    return s + cmath.log(x) - 0.5/x

digamma = _mathfun_real(_digamma_real, _digamma_complex)

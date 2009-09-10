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

def nthroot(x, n):
    r = 1./n
    try:
        return x ** r
    except ValueError:
        return complex(x) ** r

def cos_sin(x):
    if type(x) is complex:
        return cmath.cos(x), cmath.sin(x)
    else:
        return math.cos(x), math.sin(x)

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
    try:
        _intx = int(x)
        if _intx == x:
            if _intx <= 0:
                return (-1)**_intx * INF
            if _intx <= _max_exact_gamma:
                return _exact_gamma[_intx]
    except:
        pass
    if x < 0.5:
        return pi / (math.sin(pi*x)*_gamma_real(1-x))
    else:
        x -= 1.0
        r = _lanczos_p[0]
        for i in range(1, _lanczos_g+2):
            r += _lanczos_p[i]/(x+i)
        t = x + _lanczos_g + 0.5
        return 2.506628274631000502417 * t**(x+0.5) * math.exp(-t) * r

def _gamma_complex(x):
    if not x.imag:
        return complex(gamma(x.real))
    if x.real < 0.5:
        return pi / (cmath.sin(pi*x)*_gamma_complex(1-x))
    else:
        x -= 1.0
        r = _lanczos_p[0]
        for i in range(1, _lanczos_g+2):
            r += _lanczos_p[i]/(x+i)
        t = x + _lanczos_g + 0.5
        return 2.506628274631000502417 * t**(x+0.5) * cmath.exp(-t) * r

gamma = _mathfun_real(_gamma_real, _gamma_complex)

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


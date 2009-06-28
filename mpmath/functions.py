"""
This module defines most special functions and mathematical constants
provided by mpmath. [Exception: elliptic functions are currently
in elliptic.py]

Most of the actual computational code is located in the lib* modules
(libelefun, libhyper, ...); this module simply wraps this code to
handle precision management in a user friendly way, provide type
conversions, etc.

In addition, this module defines a number of functions that would
be inconvenient to define in the lib* modules, due to requiring
high level operations (e.g. numerical quadrature) for the computation,
or the need to support multiple arguments of mixed types.

"""

import libmpf
import libelefun
import libmpc
import libmpi
import gammazeta
import libhyper
import libintmath

from mptypes import (\
    MultiPrecisionArithmetic,
    def_mp_builtin,
    defun_wrapped,
    defun,
    defun_static,
    mp,
    constant,
    ComplexResult,
)

def_mp_builtin = def_mp_builtin
libelefun = libelefun
libmpc = libmpc
libhyper = libhyper
gammazeta = gammazeta

# The following multiprecision functions are implemented entirely with
# low-level code
sqrt = def_mp_builtin('sqrt', libelefun.mpf_sqrt, libmpc.mpc_sqrt, libmpi.mpi_sqrt, "principal square root")
cbrt = def_mp_builtin('cbrt', libelefun.mpf_cbrt, libmpc.mpc_cbrt, None, "principal cubic root")
exp = def_mp_builtin('exp', libelefun.mpf_exp, libmpc.mpc_exp, libmpi.mpi_exp, "exponential function")
ln = def_mp_builtin('ln', libelefun.mpf_log, libmpc.mpc_log, libmpi.mpi_log, "natural logarithm")
cos = def_mp_builtin('cos', libelefun.mpf_cos, libmpc.mpc_cos, libmpi.mpi_cos, "cosine")
sin = def_mp_builtin('sin', libelefun.mpf_sin, libmpc.mpc_sin, libmpi.mpi_sin, "sine")
tan = def_mp_builtin('tan', libelefun.mpf_tan, libmpc.mpc_tan, libmpi.mpi_tan, "tangent")
cosh = def_mp_builtin('cosh', libelefun.mpf_cosh, libmpc.mpc_cosh, None, "hyperbolic cosine")
sinh = def_mp_builtin('sinh', libelefun.mpf_sinh, libmpc.mpc_sinh, None, "hyperbolic sine")
tanh = def_mp_builtin('tanh', libelefun.mpf_tanh, libmpc.mpc_tanh, None, "hyperbolic tangent")
acos = def_mp_builtin('acos', libelefun.mpf_acos, libmpc.mpc_acos, None, "inverse cosine")
asin = def_mp_builtin('asin', libelefun.mpf_asin, libmpc.mpc_asin, None, "inverse sine")
atan = def_mp_builtin('atan', libelefun.mpf_atan, libmpc.mpc_atan, None, "inverse tangent")
asinh = def_mp_builtin('asinh', libelefun.mpf_asinh, libmpc.mpc_asinh, None, "inverse hyperbolic sine")
acosh = def_mp_builtin('acosh', libelefun.mpf_acosh, libmpc.mpc_acosh, None, "inverse hyperbolic cosine")
atanh = def_mp_builtin('atanh', libelefun.mpf_atanh, libmpc.mpc_atanh, None, "inverse hyperbolic tangent")
cospi = def_mp_builtin('cospi', libelefun.mpf_cos_pi, libmpc.mpc_cos_pi, None, "")
sinpi = def_mp_builtin('sinpi', libelefun.mpf_sin_pi, libmpc.mpc_sin_pi, None, "")
floor = def_mp_builtin('floor', libmpf.mpf_floor, libmpc.mpc_floor, None, "")
ceil = def_mp_builtin('ceil', libmpf.mpf_ceil, libmpc.mpc_ceil, None, "")
fibonacci = def_mp_builtin('fibonacci', libelefun.mpf_fibonacci, libmpc.mpc_fibonacci, None, "")
zeta = def_mp_builtin('zeta', gammazeta.mpf_zeta, gammazeta.mpc_zeta, None, "Riemann zeta function")
altzeta = def_mp_builtin('altzeta', gammazeta.mpf_altzeta, gammazeta.mpc_altzeta, None, "Dirichlet eta function")
gamma = def_mp_builtin('gamma', gammazeta.mpf_gamma, gammazeta.mpc_gamma, None, "gamma function")
factorial = def_mp_builtin('factorial', gammazeta.mpf_factorial, gammazeta.mpc_factorial, None, "factorial")
harmonic = def_mp_builtin('harmonic', gammazeta.mpf_harmonic, gammazeta.mpc_harmonic, None, "nth harmonic number")
ci = def_mp_builtin('ci', libhyper.mpf_ci, libhyper.mpc_ci, None, "")
si = def_mp_builtin('si', libhyper.mpf_si, libhyper.mpc_si, None, "")
ellipk = def_mp_builtin('ellipk', libhyper.mpf_ellipk, libhyper.mpc_ellipk, None, "")
ellipe = def_mp_builtin('ellipe', libhyper.mpf_ellipe, libhyper.mpc_ellipe, None, "")
agm1 = def_mp_builtin('agm1', libhyper.mpf_agm1, libhyper.mpc_agm1, None, "Fast alias for agm(1,a) = agm(a,1)")

_erf = def_mp_builtin('_erf', libhyper.mpf_erf, None, None, "Error function, erf(z)")
_erfc = def_mp_builtin('_erfc', libhyper.mpf_erfc, None, None, "Complementary error function, erfc(z) = 1-erf(z)")


fac = MultiPrecisionArithmetic.fac = factorial
fib = MultiPrecisionArithmetic.fib = fibonacci


# The main reason why each constant is a class and not just an instance
# is that Sphinx won't show docstrings for single instances

def defconst(name, func, descr):
    MultiPrecisionArithmetic._constants.append((name, func, descr))

defconst("pi", libelefun.mpf_pi, "pi")
defconst("degree", libelefun.mpf_degree, "degree")
defconst("e", libelefun.mpf_e, "e")
defconst("ln2", libelefun.mpf_ln2, "ln(2)")
defconst("ln10", libelefun.mpf_ln10, "ln(10)")
defconst("phi", libelefun.mpf_phi, "Golden ratio (phi)")
defconst("euler", gammazeta.mpf_euler, "Euler's constant (gamma)")
defconst("catalan", gammazeta.mpf_catalan, "Catalan's constant")
defconst("khinchin", gammazeta.mpf_khinchin, "Khinchin's constant")
defconst("glaisher", gammazeta.mpf_glaisher, "Glaisher's constant")
defconst("apery", gammazeta.mpf_apery, "Apery's constant")
defconst("mertens", gammazeta.mpf_mertens, "Mertens' constant")
defconst("twinprime", gammazeta.mpf_twinprime, "Twin prime constant")

mp._create_constants(globals())

def funcwrapper(f):
    def g(*args, **kwargs):
        orig = mp.prec
        try:
            args = [mp.convert(z) for z in args]
            mp.prec = orig + 10
            v = f(*args, **kwargs)
        finally:
            mp.prec = orig
        return +v
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g

def altfunc(f, name, desc):
    def g(self, x):
        orig = self.prec
        try:
            self.prec = orig + 10
            return self.one/f(x)
        finally:
            self.prec = orig
    g.__name__ = name
    g.__doc__ = "Computes the %s of x, 1/%s(x)" % (desc, f.__name__)
    return defun(g)

def altinvfunc(f, name, desc):
    def g(self, x):
        orig = self.prec
        try:
            self.prec = orig + 10
            return f(self.one/x)
        finally:
            self.prec = orig
    setattr(MultiPrecisionArithmetic, name, g)
    g.__name__ = name
    g.__doc__ = "Computes the inverse %s of x, %s(1/x)" % (desc, f.__name__)
    return defun(g)

sec = altfunc(cos, 'sec', 'secant')
csc = altfunc(sin, 'csc', 'cosecant')
cot = altfunc(tan, 'cot', 'cotangent')
sech = altfunc(cosh, 'sech', 'hyperbolic secant')
csch = altfunc(sinh, 'csch', 'hyperbolic cosecant')
coth = altfunc(tanh, 'coth', 'hyperbolic cotangent')
asec = altinvfunc(acos, 'asec', 'secant')
acsc = altinvfunc(asin, 'acsc', 'cosecant')
acot = altinvfunc(atan, 'acot', 'cotangent')
asech = altinvfunc(acosh, 'asech', 'hyperbolic secant')
acsch = altinvfunc(asinh, 'acsch', 'hyperbolic cosecant')
acoth = altinvfunc(atanh, 'acoth', 'hyperbolic cotangent')


@defun_wrapped
def sinc(ctx, x):
    if ctx.isinf(x):
        return 1/x
    if not x:
        return x+1
    return ctx.sin(x)/x

@defun_wrapped
def sincpi(ctx, x):
    if ctx.isinf(x):
        return 1/x
    if not x:
        return x+1
    return ctx.sinpi(x)/(ctx.pi*x)

# TODO: tests; improve implementation
@defun_wrapped
def expm1(ctx, x):
    if not x:
        return type(x)(0)
    # exp(x) - 1 ~ x
    if ctx.mag(x) < -ctx.prec:
        return x + 0.5*x**2
    # TODO: accurately eval the smaller of the real/imag parts
    return sum_accurately(ctx, lambda: iter([ctx.exp(x),-1]),1)

@defun_wrapped
def powm1(ctx, x, y):
    mag = ctx.mag
    one = ctx.one
    w = x**y - one
    M = mag(w)
    # Only moderate cancellation
    if M > -8:
        return w
    # Check for the only possible exact cases
    if not w:
        if (not y) or (x in (1, -1, 1j, -1j) and ctx.isint(y)):
            return w
    x1 = x - one
    magy = mag(y)
    lnx = ctx.ln(x)
    # Small y: x^y - 1 ~ log(x)*y + O(log(x)^2 * y^2)
    if magy + mag(lnx) < -ctx.prec:
        return lnx*y + (lnx*y)**2/2
    # TODO: accurately eval the smaller of the real/imag part
    return sum_accurately(ctx, lambda: iter([x**y, -1]), 1)

@defun
def _rootof1(ctx, k, n):
    k = int(k)
    n = int(n)
    k %= n
    if not k:
        return ctx.one
    elif 2*k == n:
        return -ctx.one
    elif 4*k == n:
        return ctx.j
    elif 4*k == 3*n:
        return -ctx.j
    return ctx.exp(2*ctx.pi*k/n*ctx.j)

@defun
def root(ctx, x, n, k=0):
    n = int(n)
    x = ctx.convert(x)
    if k:
        # Special case: there is an exact real root
        if (n & 1 and 2*k == n-1) and (not ctx.im(x)) and (ctx.re(x) < 0):
            return -ctx.root(-x, n)
        # Multiply by root of unity
        prec = ctx.prec
        try:
            ctx.prec += 10
            v = ctx.root(x, n, 0) * ctx._rootof1(k, n)
        finally:
            ctx.prec = prec
        return +v
    if hasattr(x, '_mpf_'):
        try:
            return ctx.make_mpf(libelefun.mpf_nthroot(x._mpf_, n, *ctx._prec_rounding))
        except ComplexResult:
            if ctx.trap_complex:
                raise
            x = (x._mpf_, libmpf.fzero)
    else:
        x = x._mpc_
    return ctx.make_mpc(libmpc.mpc_nthroot(x, n, *ctx._prec_rounding))

nthroot = MultiPrecisionArithmetic.nthroot = root

@defun
def unitroots(ctx, n, primitive=False):
    gcd = libintmath.gcd
    prec = ctx.prec
    try:
        ctx.prec += 10
        if primitive:
            v = [ctx._rootof1(k,n) for k in range(n) if gcd(k,n) == 1]
        else:
            # TODO: this can be done *much* faster
            v = [ctx._rootof1(k,n) for k in range(n)]
    finally:
        ctx.prec = prec
    return [+x for x in v]

@defun
def hypot(ctx, x, y):
    r"""
    Computes the Euclidean norm of the vector `(x, y)`, equal
    to `\sqrt{x^2 + y^2}`. Both `x` and `y` must be real."""
    x = ctx.convert(x)
    y = ctx.convert(y)
    return ctx.make_mpf(libmpf.mpf_hypot(x._mpf_, y._mpf_, *ctx._prec_rounding))

@defun
def ldexp(ctx, x, n):
    r"""
    Computes `x 2^n` efficiently. No rounding is performed.
    The argument `x` must be a real floating-point number (or
    possible to convert into one) and `n` must be a Python ``int``.

        >>> from mpmath import *
        >>> ldexp(1, 10)
        mpf('1024.0')
        >>> ldexp(1, -3)
        mpf('0.125')

    """
    x = ctx.convert(x)
    return ctx.make_mpf(libmpf.mpf_shift(x._mpf_, n))

@defun
def frexp(ctx, x):
    r"""
    Given a real number `x`, returns `(y, n)` with `y \in [0.5, 1)`,
    `n` a Python integer, and such that `x = y 2^n`. No rounding is
    performed.

        >>> from mpmath import *
        >>> frexp(7.5)
        (mpf('0.9375'), 3)

    """
    x = ctx.convert(x)
    y, n = libmpf.mpf_frexp(x._mpf_)
    return ctx.make_mpf(y), n

@defun
def sign(ctx, x):
    r"""
    Returns the sign of `x`, defined as `\mathrm{sign}(x) = x / |x|`
    (with the special case `\sign(0) = 0`)::

        >>> from mpmath import *
        >>> sign(10)
        mpf('1.0')
        >>> sign(-10)
        mpf('-1.0')
        >>> sign(0)
        mpf('0.0')

    Note that the sign function is also defined for complex numbers,
    for which it gives the projection onto the unit circle::

        >>> mp.dps = 15
        >>> print sign(1+j)
        (0.707106781186547 + 0.707106781186547j)

    """
    x = ctx.convert(x)
    if not x or ctx.isnan(x):
        return x
    if ctx.is_real_type(x):
        return ctx.mpf(cmp(x, 0))
    return x / abs(x)

@defun
def arg(ctx, x):
    r"""
    Computes the complex argument (phase) of `x`, defined as the
    signed angle between the positive real axis and `x` in the
    complex plane::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print arg(3)
        0.0
        >>> print arg(3+3j)
        0.785398163397448
        >>> print arg(3j)
        1.5707963267949
        >>> print arg(-3)
        3.14159265358979
        >>> print arg(-3j)
        -1.5707963267949

    The angle is defined to satisfy `-\pi < \arg(x) \le \pi` and
    with the sign convention that a nonnegative imaginary part
    results in a nonnegative argument.

    The value returned by :func:`arg` is an ``mpf`` instance.
    """
    x = ctx.convert(x)
    return ctx.atan2(x.imag, x.real)

@defun
def fabs(ctx, x):
    r"""
    Returns the absolute value of `x`, `|x|`. Unlike :func:`abs`,
    :func:`fabs` converts non-mpmath numbers (such as ``int``)
    into mpmath numbers::

        >>> from mpmath import *
        >>> fabs(3)
        mpf('3.0')
        >>> fabs(-3)
        mpf('3.0')
        >>> fabs(3+4j)
        mpf('5.0')

    """
    return abs(ctx.convert(x))

@defun
def re(ctx, x):
    r"""
    Returns the real part of `x`, `\Re(x)`. Unlike ``x.real``,
    :func:`re` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> re(3)
        mpf('3.0')
        >>> re(-1+4j)
        mpf('-1.0')
    """
    return ctx.convert(x).real

@defun
def im(ctx, x):
    r"""
    Returns the imaginary part of `x`, `\Im(x)`. Unlike ``x.imag``,
    :func:`im` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> im(3)
        mpf('0.0')
        >>> im(-1+4j)
        mpf('4.0')

    """
    return ctx.convert(x).imag

@defun
def conj(ctx, x):
    r"""
    Returns the complex conjugate of `x`, `\overline{x}`. Unlike
    ``x.conjugate()``, :func:`im` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> conj(3)
        mpf('3.0')
        >>> conj(-1+4j)
        mpc(real='-1.0', imag='-4.0')

    """
    return ctx.convert(x).conjugate()


@defun
def log(ctx, x, b=None):
    if b is None:
        return ln(x)
    wp = ctx.prec + 20
    return ctx.ln(x, prec=wp) / ctx.ln(b, prec=wp)

@defun
def log10(ctx, x):
    r"""
    Computes the base-10 logarithm of `x`, `\log_{10}(x)`. ``log10(x)``
    is equivalent to ``log(x, 10)``.
    """
    return ctx.log(x, 10)

@defun
def power(ctx, x, y):
    return ctx.convert(x) ** ctx.convert(y)

@defun
def modf(ctx,x,y):
    return ctx.convert(x) % ctx.convert(y)

@defun
def degrees(ctx,x):
    return x / ctx.degree

@defun
def radians(ctx,x):
    return x * ctx.degree

@defun
def atan2(ctx, y, x):
    x = ctx.convert(x)
    y = ctx.convert(y)
    return ctx.make_mpf(libelefun.mpf_atan2(y._mpf_, x._mpf_, *ctx._prec_rounding))

@defun
def psi(ctx, m, z):
    z = ctx.convert(z)
    m = int(m)
    if ctx.is_real_type(z):
        return ctx.make_mpf(gammazeta.mpf_psi(m, z._mpf_, *ctx._prec_rounding))
    else:
        return ctx.make_mpc(gammazeta.mpc_psi(m, z._mpc_, *ctx._prec_rounding))

@defun
def psi0(ctx, z):
    """Shortcut for psi(0,z) (the digamma function)"""
    return ctx.psi(0, z)

@defun
def psi1(ctx, z):
    """Shortcut for psi(1,z) (the trigamma function)"""
    return ctx.psi(1, z)

@defun
def psi2(ctx, z):
    """Shortcut for psi(2,z) (the tetragamma function)"""
    return ctx.psi(2, z)

@defun
def psi3(ctx, z):
    """Shortcut for psi(3,z) (the pentagamma function)"""
    return ctx.psi(3, z)

polygamma = MultiPrecisionArithmetic.polygamma = psi
digamma = MultiPrecisionArithmetic.digamma = psi0
trigamma = MultiPrecisionArithmetic.trigamma = psi1
tetragamma = MultiPrecisionArithmetic.tetragamma = psi2
pentagamma = MultiPrecisionArithmetic.pentagamma = psi3

@defun
def bernoulli(ctx, n):
    return ctx.make_mpf(gammazeta.mpf_bernoulli(int(n), *ctx._prec_rounding))

bernfrac = defun_static(gammazeta.bernfrac)

@defun
def stieltjes(ctx, n, a=1):
    n = ctx.convert(n)
    a = ctx.convert(a)
    if n < 0:
        return ctx.bad_domain("Stieltjes constants defined for n >= 0")
    if hasattr(ctx, "stieltjes_cache"):
        stieltjes_cache = ctx.stieltjes_cache
    else:
        stieltjes_cache = ctx.stieltjes_cache = {}
    if a == 1:
        if n == 0:
            return +ctx.euler
        if n in stieltjes_cache:
            prec, s = stieltjes_cache[n]
            if prec >= ctx.prec:
                return +s
    mag = 1
    def f(x):
        xa = x/a
        v = (xa-ctx.j)*ctx.ln(a-ctx.j*x)**n/(1+xa**2)/(ctx.exp(2*ctx.pi*x)-1)
        return v.real / mag
    orig = ctx.prec
    try:
        # Normalize integrand by approx. magnitude to
        # speed up quadrature (which uses absolute error)
        if n > 50:
            ctx.prec = 20
            mag = ctx.quad(f, [0,ctx.inf], maxdegree=3)
        ctx.prec = orig + 10 + int(n**0.5)
        s = ctx.quad(f, [0,ctx.inf], maxdegree=20)
        v = ctx.ln(a)**n/(2*a) - ctx.ln(a)**(n+1)/(n+1) + 2*s/a*mag
    finally:
        ctx.prec = orig
    if a == 1 and ctx.isint(n):
        stieltjes_cache[n] = (ctx.prec, v)
    return +v

@defun
def gammaprod(ctx, a, b, _infsign=False):
    a = [ctx.convert(x) for x in a]
    b = [ctx.convert(x) for x in b]
    poles_num = []
    poles_den = []
    regular_num = []
    regular_den = []
    for x in a: [regular_num, poles_num][ctx.isnpint(x)].append(x)
    for x in b: [regular_den, poles_den][ctx.isnpint(x)].append(x)
    # One more pole in numerator or denominator gives 0 or inf
    if len(poles_num) < len(poles_den): return ctx.zero
    if len(poles_num) > len(poles_den):
        # Get correct sign of infinity for x+h, h -> 0 from above
        # XXX: hack, this should be done properly
        if _infsign:
            a = [x and x*(1+ctx.eps) or x+ctx.eps for x in poles_num]
            b = [x and x*(1+ctx.eps) or x+ctx.eps for x in poles_den]
            return sign(gammaprod(a+regular_num,b+regular_den)) * ctx.inf
        else:
            return ctx.inf
    # All poles cancel
    # lim G(i)/G(j) = (-1)**(i+j) * gamma(1-j) / gamma(1-i)
    p = ctx.one
    orig = ctx.prec
    try:
        ctx.prec = orig + 15
        while poles_num:
            i = poles_num.pop()
            j = poles_den.pop()
            p *= (-1)**(i+j) * ctx.gamma(1-j) / ctx.gamma(1-i)
        for x in regular_num: p *= ctx.gamma(x)
        for x in regular_den: p /= ctx.gamma(x)
    finally:
        ctx.prec = orig
    return +p

@defun
def beta(ctx, x, y):
    x = ctx.convert(x)
    y = ctx.convert(y)
    if ctx.isinf(y):
        x, y = y, x
    if ctx.isinf(x):
        if x == ctx.inf and not y.imag:
            if y == ctx.ninf:
                return ctx.nan
            if y > 0:
                return ctx.zero
            if ctx.isint(y):
                return ctx.nan
            if y < 0:
                return ctx.sign(ctx.gamma(y)) * ctx.inf
        return ctx.nan
    return ctx.gammaprod([x, y], [x+y])

@defun_wrapped
def betainc(ctx, a, b, x1=0, x2=1, regularized=False):
    if x1 == x2:
        v = 0
    elif not x1:
        if x1 == 0 and x2 == 1:
            v = ctx.beta(a, b)
        else:
            v = x2**a * ctx.hyp2f1(a, 1-b, a+1, x2) / a
    else:
        m, d = ctx.nint_distance(a)
        if m <= 0:
            if d < -ctx.prec:
                h = +ctx.eps
                ctx.prec *= 2
                a += h
            elif d < -4:
                ctx.prec -= d
        s1 = x2**a * ctx.hyp2f1(a,1-b,a+1,x2)
        s2 = x1**a * ctx.hyp2f1(a,1-b,a+1,x1)
        v = (s1 - s2) / a
    if regularized:
        v /= ctx.beta(a,b)
    return v

@defun
def binomial(ctx, n, k):
    return ctx.gammaprod([n+1], [k+1, n-k+1])

@defun
def rf(ctx, x, n):
    return ctx.gammaprod([x+n], [x])

@defun
def ff(ctx, x, n):
    return ctx.gammaprod([x+1], [x-n+1])

@defun_wrapped
def fac2(ctx, x):
    if ctx.isinf(x):
        if x == ctx.inf:
            return x
        return ctx.nan
    return 2**(x/2)*(ctx.pi/2)**((ctx.cospi(x)-1)/4)*ctx.gamma(x/2+1)


#---------------------------------------------------------------------------#
#                                                                           #
#                          Hypergeometric functions                         #
#                                                                           #
#---------------------------------------------------------------------------#

from libmpf import from_rational

class _mpq(tuple):

    def _mpmath_(self, prec, rounding):
        # XXX
        return mp.make_mpf(from_rational(self[0], self[1], prec, rounding))
        #(mpf(self[0])/self[1])._mpf_

    @property
    def _mpq_(self):
        return self

    def __int__(self):
        a, b = self
        return a // b

    def __abs__(self):
        a, b = self
        return _mpq((abs(a), b))

    def __neg__(self):
        a, b = self
        return _mpq((-a, b))

    def __nonzero__(self):
        return bool(self[0])

    def __cmp__(self, other):
        if type(other) is int and self[1] == 1:
            return cmp(self[0], other)
        return cmp(mp.mpf(self), other)

    def __add__(self, other):
        if isinstance(other, _mpq):
            a, b = self
            c, d = other
            return _mpq((a*d+b*c, b*d))
        if isinstance(other, (int, long)):
            a, b = self
            return _mpq((a+b*other, b))
        return NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, _mpq):
            a, b = self
            c, d = other
            return _mpq((a*d-b*c, b*d))
        if isinstance(other, (int, long)):
            a, b = self
            return _mpq((a-b*other, b))
        return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, _mpq):
            a, b = self
            c, d = other
            return _mpq((b*c-a*d, b*d))
        if isinstance(other, (int, long)):
            a, b = self
            return _mpq((b*other-a, b))
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, _mpq):
            a, b = self
            c, d = other
            return _mpq((a*c, b*d))
        if isinstance(other, (int, long)):
            a, b = self
            return _mpq((a*other, b))
        return NotImplemented

    def __pow__(self, other):
        if type(other) is int:
            a, b = self
            return _mpq((a**other, b**other))
        return NotImplemented

    __rmul__ = __mul__


mpq_1 = _mpq((1,1))
mpq_0 = _mpq((0,1))
mpq_1_2 = _mpq((1,2))
mpq_3_2 = _mpq((3,2))
mpq_1_4 = _mpq((1,4))
mpq_1_16 = _mpq((1,16))
mpq_3_16 = _mpq((3,16))
mpq_5_2 = _mpq((5,2))

@defun
def _hyp_parse_param(ctx, x):
    if isinstance(x, tuple):
        p, q = x
        return [[p, q]], [], [], _mpq(x)
    if isinstance(x, (int, long)):
        return [[x, 1]], [], [], x
    x = ctx.convert(x)
    if hasattr(x, '_mpf_'):
        sign, man, exp, bc = _mpf_ = x._mpf_
        # Recognize simple rationals
        if exp >= -4:
            if sign:
                man = -man
            if exp >= 0:
                return [[int(man)<<exp, 1]], [], [], x
            return [[int(man), 2**(-exp)]], [], [], x
        else:
            return [], [_mpf_], [], x
    if hasattr(x, '_mpc_'):
        return [], [], [x._mpc_], x

def _as_num(x):
    if isinstance(x, list):
        return _mpq(x)
    return x

@defun
def hypsum(ctx, ar, af, ac, br, bf, bc, x):
    prec, rnd = ctx._prec_rounding
    if hasattr(x, '_mpf_') and not (ac or bc):
        v = libhyper.hypsum_internal(ar, af, ac, br, bf, bc, x._mpf_, None, prec, rnd)
        return ctx.make_mpf(v)
    else:
        if hasattr(x, '_mpc_'):
            re, im = x._mpc_
        else:
            re, im = x._mpf_, libmpf.fzero
        v = libhyper.hypsum_internal(ar, af, ac, br, bf, bc, re, im, prec, rnd)
        return ctx.make_mpc(v)

@defun
def hyper(ctx, a_s, b_s, z, **kwargs):
    """
    Hypergeometric function, general case.
    """
    p = len(a_s)
    q = len(b_s)
    z = ctx.convert(z)
    degree = p, q
    # Handle common cases
    if p == 0:
        if q == 1:
            return ctx.hyp0f1(b_s[0], z, **kwargs)
        elif q == 0:
            return ctx.exp(z)
    elif p == 1:
        if q == 1:
            return ctx.hyp1f1(a_s[0], b_s[0], z, **kwargs)
        elif q == 2:
            return ctx.hyp1f2(a_s[0], b_s[0], b_s[1], z, **kwargs)
        elif q == 0:
            return (1-z)**(-a_s[0])
    elif p == 2:
        if q == 1:
            return ctx.hyp2f1(a_s[0], a_s[1], b_s[0], z, **kwargs)
        elif q == 2:
            return ctx.hyp2f2(a_s[0], a_s[1], b_s[0], b_s[1], z, **kwargs)
        elif q == 3:
            return ctx.hyp2f3(a_s[0], a_s[1], b_s[0], b_s[1], b_s[2], z, **kwargs)
        elif q == 0:
            return ctx.hyp2f0(a_s[0], a_s[1], z, **kwargs)
    # TODO: do convergence tests here
    z = ctx.convert(z)
    ars, afs, acs, brs, bfs, bcs = [], [], [], [], [], []
    for a in a_s:
        r, f, c, a = ctx._hyp_parse_param(a)
        ars += r
        afs += f
        acs += c
    for b in b_s:
        r, f, c, b = ctx._hyp_parse_param(b)
        brs += r
        bfs += f
        bcs += c
    return ctx.hypsum(ars, afs, acs, brs, bfs, bcs, z)

@defun
def sum_hyp0f1_rat(ctx, a, z):
    prec, rnd = ctx._prec_rounding
    if hasattr(z, "_mpf_"):
        return ctx.make_mpf(libhyper.mpf_hyp0f1_rat(a, z._mpf_, prec, rnd))
    else:
        return ctx.make_mpc(libhyper.mpc_hyp0f1_rat(a, z._mpc_, prec, rnd))

@defun
def sum_hyp1f1_rat(ctx, a, b, z):
    prec, rnd = ctx._prec_rounding
    if hasattr(z, "_mpf_"):
        return ctx.make_mpf(libhyper.mpf_hyp1f1_rat(a, b, z._mpf_, prec, rnd))
    else:
        return ctx.make_mpc(libhyper.mpc_hyp1f1_rat(a, b, z._mpc_, prec, rnd))

@defun
def sum_hyp2f1_rat(ctx, a, b, c, z):
    prec, rnd = ctx._prec_rounding
    if hasattr(z, "_mpf_"):
        return ctx.make_mpf(libhyper.mpf_hyp2f1_rat(a, b, c, z._mpf_, prec, rnd))
    else:
        return ctx.make_mpc(libhyper.mpc_hyp2f1_rat(a, b, c, z._mpc_, prec, rnd))

@defun
def _hyp_check_convergence(ctx, a_s, b_s, z, prec, n=None):
    p = len(a_s)
    q = len(b_s)
    a = max([1] + map(abs, a_s))
    b = min([1] + map(abs, b_s))
    z = abs(z)
    amag = ctx.mag(a)
    bmag = ctx.mag(b)
    zmag = ctx.mag(z)
    tol = -prec
    # z extremely tiny
    if p*amag - q*bmag + zmag < tol:
        return True
    a = int(a)
    b = int(b)
    # XXX: better choice of n?

    if n:
        n = int(n)
    else:
        n = prec
    nmag = ctx.mag(n)

    #n = int(1/z)
    #nmag = -zmag

    # Accurately estimate size of nth term using Stirling's formula
    t = n*zmag
    t -= (n*nmag-n)
    t += p*((a+n)*max(amag,nmag)-(a*amag-a))
    t -= q*((b+n)*max(bmag,nmag)-(b*bmag-b))

    u = z**n / ctx.fac(n) * rf(a,n)**p / rf(b,n)**q

    #print "debug", prec, "kek", n, t, tol, t < tol, "real", int(log(abs(u),2))
    return t < tol

@defun
def hyp0f1(ctx, b, z, **kwargs):
    """
    Hypergeometric 0F1.
    """
    br, bf, bc, b = ctx._hyp_parse_param(b)
    z = ctx.convert(z)
    magz = z and ctx.mag(z)
    if magz >= 8 and not kwargs.get('force_series'):
        #if ctx._hyp2f0_check_convergence(b, b, 1/abs(z)**0.5, ctx.prec+40+magz//2):
        if ctx._hyp_check_convergence([b, b], [], 1/abs(z)**0.5, ctx.prec+40+magz//2):
            # http://functions.wolfram.com/HypergeometricFunctions/
            # Hypergeometric0F1/06/02/03/0004/
            # We don't need hypercomb because the only possible singularity
            # occurs when the value is undefined. However, we should perhaps
            # still check for cancellation...
            # TODO: handle the all-real case more efficiently!
            # TODO: figure out how much precision is needed (exponential growth)
            orig = ctx.prec
            try:
                ctx.prec += 12 + magz//2
                w = ctx.sqrt(-z)
                jw = ctx.j*w
                u = 1/(4*jw)
                c = mpq_1_2 - b
                E = exp(2*jw)
                H1 = (-jw)**c/E*ctx.hyp2f0(b-mpq_1_2, mpq_3_2-b, -u,
                    force_series=True)
                H2 = (jw)**c*E*ctx.hyp2f0(b-mpq_1_2, mpq_3_2-b, u,
                    force_series=True)
                v = gamma(b)/(2*sqrt(pi))*(H1 + H2)
            finally:
                ctx.prec = orig
            if ctx.is_real_type(b) and ctx.is_real_type(z):
                v = v.real
            return +v
    if br:
        return ctx.sum_hyp0f1_rat(br[0], z)
    return ctx.hypsum([], [], [], br, bf, bc, z)

@defun
def hyp1f1(ctx, a, b, z, **kwargs):
    """
    Hypergeometric 1F1.
    """
    ar, af, ac, a = ctx._hyp_parse_param(a)
    br, bf, bc, b = ctx._hyp_parse_param(b)
    z = ctx.convert(z)
    if not z:
        return ctx.one+z
    if ctx.mag(z) >= 7:
        if ctx.isinf(z):
            if ctx.sign(a) == ctx.sign(b) == ctx.sign(z) == 1:
                return ctx.inf
            return ctx.nan * z
        # TODO: extra precision?
        # Check with twice the precision because a limit could be invoked
        #if ctx._hyp2f0_check_convergence(a, a-b, 1/z, 2*ctx.prec+40):
        checked = 2*ctx.prec+40
        if ctx._hyp_check_convergence([a, a-b], [], 1/z, 2*ctx.prec+40):
            sector = z.imag < 0 and z.real <= 0
            def h(a,b):
                if sector:
                    E = ctx.exp(-ctx.j * ctx.pi*a)
                else:
                    E = ctx.exp(ctx.j * ctx.pi*a)
                rz = 1/z
                T1 = ([E,z], [1,-a], [b], [b-a], [a, 1+a-b], [], -rz)
                T2 = ([ctx.exp(z),z], [1,a-b], [b], [a], [b-a, 1-a], [], rz)
                return T1, T2
            v = ctx.hypercomb(h, [a,b], force_series=True)
            if ctx.is_real_type(a) and ctx.is_real_type(b) and ctx.is_real_type(z):
                v = v.real
            return +v
    if ar and br:
        a, b = ar[0], br[0]
        return ctx.sum_hyp1f1_rat(a, b, z)
    return ctx.hypsum(ar, af, ac, br, bf, bc, z)

def _hyp2f1_gosper(ctx,a,b,c,z):
    # Use Gosper's recurrence
    # See http://www.math.utexas.edu/pipermail/maxima/2006/000126.html
    a = ctx.convert(a)
    b = ctx.convert(b)
    c = ctx.convert(c)
    d = ctx.mpf(0)
    e = ctx.mpf(1)
    f = ctx.mpf(0)
    k = 0
    # Common subexpression elimination, unfortunately making
    # things a bit unreadable. The formula is quite messy to begin
    # with, though...
    abz = a*b*z
    ch = c/2
    c1h = (c+1)/2
    nz = 1-z
    g = z/nz
    abg = a*b*g
    cba = c-b-a
    z2 = z-2
    tol = -ctx.prec - 10
    while 1:
        kch = k+ch
        kakbz = (k+a)*(k+b)*z / (4*(k+1)*kch*(k+c1h))
        d1 = kakbz*(e-(k+cba)*d*g)
        e1 = kakbz*(d*abg+(k+c)*e)
        f1 = f + e - d*(k*(cba*z+k*z2-c)-abz)/(2*kch*nz)
        if ctx.mag(f1-f) < tol:
            break
        d, e, f = d1, e1, f1
        k += 1
    return f1

@defun
def hyp2f1(ctx,a,b,c,z,**kwargs):
    """
    Hypergeometric 2F1.
    """
    prec, rnd = ctx._prec_rounding
    ar, af, ac, a = ctx._hyp_parse_param(a)
    br, bf, bc, b = ctx._hyp_parse_param(b)
    cr, cf, cc, c = ctx._hyp_parse_param(c)
    z = ctx.convert(z)

    if z == 1:
        # TODO: the following logic can be simplified
        convergent = ctx.re(c-a-b) > 0
        finite = (ctx.isint(a) and a <= 0) or (ctx.isint(b) and b <= 0)
        zerodiv = ctx.isint(c) and c <= 0 and not \
            ((ctx.isint(a) and c <= a <= 0) or (ctx.isint(b) and c <= b <= 0))
        #print "bz", a, b, c, z, convergent, finite, zerodiv
        # Gauss's theorem gives the value if convergent
        if (convergent or finite) and not zerodiv:
            return ctx.gammaprod([c, c-a-b], [c-a, c-b], _infsign=True)
        # Otherwise, there is a pole and we take the
        # sign to be that when approaching from below
        # XXX: this evaluation is not necessarily correct in all cases
        return ctx.hyp2f1(a,b,c,1-ctx.eps*2) * ctx.inf

    # Equal to 1 (first term), unless there is a subsequent
    # division by zero
    if not z:
        # Division by zero but power of z is higher than
        # first order so cancels
        if c or a == 0 or b == 0:
            return 1+z
        # Indeterminate
        return ctx.nan

    # Hit zero denominator unless numerator goes to 0 first
    if ctx.isint(c) and c <= 0:
        if (ctx.isint(a) and c <= a <= 0) or \
           (ctx.isint(b) and c <= b <= 0):
            pass
        else:
            # Pole in series
            return ctx.inf
    elif a == c:
        return (1-z)**(-b)
    elif b == c:
        return (1-z)**(-a)

    absz = abs(z)

    # Fast case: standard series converges rapidly,
    # possibly in finitely many terms
    if absz <= 0.8 or (ctx.isint(a) and a <= 0 and a >= -1000) or \
                      (ctx.isint(b) and b <= 0 and b >= -1000):
        # All rational
        if ar and br and cr:
            return ctx.sum_hyp2f1_rat(ar[0], br[0], cr[0], z)
        return ctx.hypsum(ar+br, af+bf, ac+bc, cr, cf, cc, z)

    a = (ar and _as_num(ar[0])) or ctx.convert(a)
    b = (br and _as_num(br[0])) or ctx.convert(b)
    c = (cr and _as_num(cr[0])) or ctx.convert(c)

    orig = ctx.prec
    try:
        ctx.prec += 10

        # Use 1/z transformation
        if absz >= 1.3:
            def h(a,b):
                t = mpq_1-c; ab = a-b; rz = 1/z
                T1 = ([-z],[-a], [c,-ab],[b,c-a], [a,t+a],[mpq_1+ab],  rz)
                T2 = ([-z],[-b], [c,ab],[a,c-b], [b,t+b],[mpq_1-ab],  rz)
                return T1, T2
            v = ctx.hypercomb(h, [a,b])

        # Use 1-z transformation
        elif abs(1-z) <= 0.75:
            def h(a,b):
                t = c-a-b; ca = c-a; cb = c-b; rz = 1-z
                T1 = [], [], [c,t], [ca,cb], [a,b], [1-t], rz
                T2 = [rz], [t], [c,a+b-c], [a,b], [ca,cb], [1+t], rz
                return T1, T2
            v = ctx.hypercomb(h, [a,b])

        # Remaining part of unit circle
        else:
            v = _hyp2f1_gosper(ctx,a,b,c,z)
    finally:
        ctx.prec = orig
    return +v

@defun
def hypercomb(ctx, function, params=[], **kwargs):
    orig = ctx.prec
    sumvalue = 0
    dist = ctx.nint_distance
    ninf = ctx.ninf
    convert = ctx.convert
    orig_params = params[:]
    try:
        while 1:
            ctx.prec += 10
            orig2 = ctx.prec
            params = orig_params[:]
            terms = function(*params)
            perturb = False
            recompute = False
            for term in terms:
                w_s, c_s, alpha_s, beta_s, a_s, b_s, z = term
                extraprec = 0
                try:
                    # Avoid division by zero in leading factors (TODO:
                    # also check for near division by zero?)
                    for k, w in enumerate(w_s):
                        if not w:
                            if ctx.re(c_s[k]) <= 0:
                                perturb = recompute = True
                                raise StopIteration
                    # Check for gamma and series poles and near-poles
                    for data in alpha_s, beta_s, b_s:
                        for i, x in enumerate(data):
                            n, d = dist(x)
                            # Poles
                            if n > 0:
                                continue
                            if d < -orig:
                                perturb = recompute = True
                                raise StopIteration
                            elif d < -4:
                                ctx.prec += (-d)
                                recompute = True
                except StopIteration:
                    pass
            if perturb:
                """
                # Should check for poles far from 0 and ensure that
                # the summation includes those terms
                minterms = 0
                for term in terms:
                    for b in term[-2]:
                        n, d = dist(b)
                        if d < -4 and n <= 0:
                            minterms = max(minterms, -n)
                """
                h = ctx.ldexp(1,-orig-10)
                ctx.prec = (orig2+10)*2
                for k in range(len(params)):
                    params[k] += h
                    # Heuristically ensure that the perturbations
                    # are "independent" so that two perturbations
                    # don't accidentally cancel each other out
                    # in a subtraction.
                    h += h/(k+1)
            if recompute:
                terms = function(*params)
            evaluated_terms = []
            for w_s, c_s, alpha_s, beta_s, a_s, b_s, z in terms:
                v = ctx.hyper(a_s, b_s, z, **kwargs)
                for a in alpha_s: v *= ctx.gamma(a)
                for b in beta_s: v /= ctx.gamma(b)
                for w, c in zip(w_s, c_s):
                    v *= convert(w) ** c
                evaluated_terms.append(v)
            if len(terms) == 1:
                sumvalue = evaluated_terms[0]
                break
            # Problem: reconcile this with intentional cancellation
            elif kwargs.get('check_cancellation'):
                sumvalue = sum(evaluated_terms)
                c = max(ctx.mag(x) for x in evaluated_terms) - ctx.mag(sumvalue)
                if c < ctx.prec - orig:
                    break
                else:
                    ctx.prec += min(c, orig)
                    if ctx.prec > 50*orig:
                        raise ValueError("hypercomb() was called with "
                            "check_cancellation=True but failed to converge "
                            "within a reasonable number of steps. The function "
                            "value is probably either zero or extremely small.")
                    continue
            else:
                sumvalue = sum(evaluated_terms)
                break
    finally:
        ctx.prec = orig
    return +sumvalue


class NoConvergence(Exception):
    pass


@defun
def hyp2f2(ctx,a1,a2,b1,b2,z,**kwargs):
    a1r, a1f, a1c, a1 = ctx._hyp_parse_param(a1)
    a2r, a2f, a2c, a2 = ctx._hyp_parse_param(a2)
    b1r, b1f, b1c, b1 = ctx._hyp_parse_param(b1)
    b2r, b2f, b2c, b2 = ctx._hyp_parse_param(b2)
    z = ctx.convert(z)

    absz = abs(z)
    magz = ctx.mag(z)
    orig = ctx.prec

    # Asymptotic expansion is ~ exp(z)
    asymp_extraprec = magz

    # Asymptotic series is in terms of 3F1
    can_use_asymptotic = (not kwargs.get('force_series')) and \
        (ctx.mag(absz) > 7) and \
        (ctx.sqrt(absz) > 1.5*orig) and \
        ctx._hyp_check_convergence([a1, a1-b1+1, a1-b2+1], [a1-a2+1],
            1/absz, orig+40+asymp_extraprec)

    # TODO: much of the following could be shared with 2F3 instead of
    # copypasted
    if can_use_asymptotic:
        #print "using asymp"
        try:
            ctx.prec += asymp_extraprec
            # http://functions.wolfram.com/HypergeometricFunctions/
            # Hypergeometric2F2/06/02/02/0002/
            def h(a1,a2,b1,b2):
                X = a1+a2-b1-b2
                A2 = a1+a2
                B2 = b1+b2
                c = {}
                c[0] = ctx.one
                c[1] = (A2-1)*X+b1*b2-a1*a2
                s1 = 0
                k = 0
                tprev = 0
                while 1:
                    if k not in c:
                        uu1 = 1-B2+2*a1+a1**2+2*a2+a2**2-A2*B2+a1*a2+b1*b2+(2*B2-3*(A2+1))*k+2*k**2
                        uu2 = (k-A2+b1-1)*(k-A2+b2-1)*(k-X-2)
                        c[k] = ctx.one/k * (uu1*c[k-1]-uu2*c[k-2])
                    t1 = c[k] * z**(-k)
                    if abs(t1) < 0.1*ctx.eps:
                        #print "Convergence :)"
                        break
                    # Quit if the series doesn't converge quickly enough
                    if k > 5 and abs(tprev) / abs(t1) < 1.5:
                        #print "No convergence :("
                        raise NoConvergence
                    s1 += t1
                    tprev = t1
                    k += 1
                S = ctx.exp(z)*s1
                T1 = [z,S], [X,1], [b1,b2],[a1,a2],[],[],0
                T2 = [-z],[-a1],[b1,b2,a2-a1],[a2,b1-a1,b2-a1],[a1,a1-b1+1,a1-b2+1],[a1-a2+1],-1/z
                T3 = [-z],[-a2],[b1,b2,a1-a2],[a1,b1-a2,b2-a2],[a2,a2-b1+1,a2-b2+1],[-a1+a2+1],-1/z
                return T1, T2, T3
            v = hypercomb(h, [a1,a2,b1,b2])
            if sum(ctx.is_real_type(u) for u in [a1,a2,b1,b2,z]) == 5:
                v = ctx.re(v)
            return v
        except NoConvergence:
            pass
        finally:
            ctx.prec = orig

    #print "not using asymp"
    return ctx.hypsum(a1r+a2r, a1f+a2f, a1c+a2c, b1r+b2r, b1f+b2f, b1c+b2c, z)




@defun
def hyp1f2(ctx,a1,b1,b2,z,**kwargs):
    a1r, a1f, a1c, a1 = ctx._hyp_parse_param(a1)
    b1r, b1f, b1c, b1 = ctx._hyp_parse_param(b1)
    b2r, b2f, b2c, b2 = ctx._hyp_parse_param(b2)
    z = ctx.convert(z)

    absz = abs(z)
    magz = ctx.mag(z)
    orig = ctx.prec

    # Asymptotic expansion is ~ exp(sqrt(z))
    asymp_extraprec = z and magz//2

    # Asymptotic series is in terms of 3F0
    can_use_asymptotic = (not kwargs.get('force_series')) and \
        (ctx.mag(absz) > 19) and \
        (ctx.sqrt(absz) > 1.5*orig) and \
        ctx._hyp_check_convergence([a1, a1-b1+1, a1-b2+1], [],
            1/absz, orig+40+asymp_extraprec)

    # TODO: much of the following could be shared with 2F3 instead of
    # copypasted
    if can_use_asymptotic:
        #print "using asymp"
        try:
            ctx.prec += asymp_extraprec
            # http://functions.wolfram.com/HypergeometricFunctions/
            # Hypergeometric1F2/06/02/03/
            def h(a1,b1,b2):
                X = mpq_1_2*(a1-b1-b2+mpq_1_2)
                c = {}
                c[0] = ctx.one
                c[1] = 2*(mpq_1_4*(3*a1+b1+b2-2)*(a1-b1-b2)+b1*b2-mpq_3_16)
                c[2] = 2*(b1*b2+mpq_1_4*(a1-b1-b2)*(3*a1+b1+b2-2)-mpq_3_16)**2+\
                    mpq_1_16*(-16*(2*a1-3)*b1*b2 + \
                    4*(a1-b1-b2)*(-8*a1**2+11*a1+b1+b2-2)-3)
                s1 = 0
                s2 = 0
                k = 0
                tprev = 0
                while 1:
                    if k not in c:
                        uu1 = (3*k**2+(-6*a1+2*b1+2*b2-4)*k + 3*a1**2 - \
                            (b1-b2)**2 - 2*a1*(b1+b2-2) + mpq_1_4)
                        uu2 = (k-a1+b1-b2-mpq_1_2)*(k-a1-b1+b2-mpq_1_2)*\
                            (k-a1+b1+b2-mpq_5_2)
                        c[k] = ctx.one/(2*k)*(uu1*c[k-1]-uu2*c[k-2])
                    w = c[k] * (-z)**(-0.5*k)
                    t1 = (-ctx.j)**k * ctx.mpf(2)**(-k) * w
                    t2 = ctx.j**k * ctx.mpf(2)**(-k) * w
                    if abs(t1) < 0.1*ctx.eps:
                        #print "Convergence :)"
                        break
                    # Quit if the series doesn't converge quickly enough
                    if k > 5 and abs(tprev) / abs(t1) < 1.5:
                        #print "No convergence :("
                        raise NoConvergence
                    s1 += t1
                    s2 += t2
                    tprev = t1
                    k += 1
                S = ctx.exp(ctx.j*(ctx.pi*X+2*ctx.sqrt(-z)))*s1 + \
                    ctx.exp(-ctx.j*(ctx.pi*X+2*ctx.sqrt(-z)))*s2
                T1 = [0.5*S, ctx.pi, -z], [1, -0.5, X], [b1, b2], [a1],\
                    [], [], 0
                T2 = [-z], [-a1], [b1,b2],[b1-a1,b2-a1], \
                    [a1,a1-b1+1,a1-b2+1], [], 1/z
                return T1, T2
            v = hypercomb(h, [a1,b1,b2])
            if sum(ctx.is_real_type(u) for u in [a1,b1,b2,z]) == 4:
                v = ctx.re(v)
            return v
        except NoConvergence:
            pass
        finally:
            ctx.prec = orig

    #print "not using asymp"
    return ctx.hypsum(a1r, a1f, a1c, b1r+b2r, b1f+b2f, b1c+b2c, z)



@defun
def hyp2f3(ctx,a1,a2,b1,b2,b3,z,**kwargs):
    prec, rnd = ctx._prec_rounding
    a1r, a1f, a1c, a1 = ctx._hyp_parse_param(a1)
    a2r, a2f, a2c, a2 = ctx._hyp_parse_param(a2)
    b1r, b1f, b1c, b1 = ctx._hyp_parse_param(b1)
    b2r, b2f, b2c, b2 = ctx._hyp_parse_param(b2)
    b3r, b3f, b3c, b3 = ctx._hyp_parse_param(b3)
    z = ctx.convert(z)

    #return ctx.hypsum(a1r+a2r, a1f+a2f, a1c+a2c, b1r+b2r+b3r,
    #    b1f+b2f+b3f, b1c+b2c+b3c, z)

    absz = abs(z)
    magz = ctx.mag(z)

    # Asymptotic expansion is ~ exp(sqrt(z))
    asymp_extraprec = z and magz//2
    orig = ctx.prec

    # Asymptotic series is in terms of 4F1
    # The square root below empirically provides a plausible criterion
    # for the leading series to converge
    can_use_asymptotic = (not kwargs.get('force_series')) and \
        (ctx.mag(absz) > 19) and \
        (ctx.sqrt(absz) > 1.5*orig) and \
        ctx._hyp_check_convergence([a2, a2-b1+1, a2-b2+1, a2-b3+1], [-a1+a2+1],
            1/absz, orig+40+asymp_extraprec)

    if can_use_asymptotic:
        #print "using asymp"
        try:
            ctx.prec += asymp_extraprec
            # http://functions.wolfram.com/HypergeometricFunctions/
            # Hypergeometric2F3/06/02/03/01/0002/
            def h(a1,a2,b1,b2,b3):
                X = mpq_1_2*(a1+a2-b1-b2-b3+mpq_1_2)
                A2 = a1+a2
                B3 = b1+b2+b3
                A = a1*a2
                B = b1*b2+b3*b2+b1*b3
                R = b1*b2*b3
                c = {}
                c[0] = ctx.one
                c[1] = 2*(B - A + mpq_1_4*(3*A2+B3-2)*(A2-B3) - mpq_3_16)
                c[2] = mpq_1_2*c[1]**2 + mpq_1_16*(-16*(2*A2-3)*(B-A) + 32*R +\
                    4*(-8*A2**2 + 11*A2 + 8*A + B3 - 2)*(A2-B3)-3)
                s1 = 0
                s2 = 0
                k = 0
                tprev = 0
                while 1:
                    if k not in c:
                        uu1 = (k-2*X-3)*(k-2*X-2*b1-1)*(k-2*X-2*b2-1)*\
                            (k-2*X-2*b3-1)
                        uu2 = (4*(k-1)**3 - 6*(4*X+B3)*(k-1)**2 + \
                            2*(24*X**2+12*B3*X+4*B+B3-1)*(k-1) - 32*X**3 - \
                            24*B3*X**2 - 4*B - 8*R - 4*(4*B+B3-1)*X + 2*B3-1)
                        uu3 = (5*(k-1)**2+2*(-10*X+A2-3*B3+3)*(k-1)+2*c[1])
                        c[k] = ctx.one/(2*k)*(uu1*c[k-3]-uu2*c[k-2]+uu3*c[k-1])
                    w = c[k] * (-z)**(-0.5*k)
                    t1 = (-ctx.j)**k * ctx.mpf(2)**(-k) * w
                    t2 = ctx.j**k * ctx.mpf(2)**(-k) * w
                    #print k, abs(t1)
                    if abs(t1) < 0.1*ctx.eps:
                        #print "Convergence :)"
                        break
                    # Quit if the series doesn't converge quickly enough
                    if k > 5 and abs(tprev) / abs(t1) < 1.5:
                        #print "No convergence :("
                        raise NoConvergence
                    s1 += t1
                    s2 += t2
                    tprev = t1
                    k += 1
                S = ctx.exp(ctx.j*(ctx.pi*X+2*ctx.sqrt(-z)))*s1 + \
                    ctx.exp(-ctx.j*(ctx.pi*X+2*ctx.sqrt(-z)))*s2
                T1 = [0.5*S, ctx.pi, -z], [1, -0.5, X], [b1, b2, b3], [a1, a2],\
                    [], [], 0
                T2 = [-z], [-a1], [b1,b2,b3,a2-a1],[a2,b1-a1,b2-a1,b3-a1], \
                    [a1,a1-b1+1,a1-b2+1,a1-b3+1], [a1-a2+1], 1/z
                T3 = [-z], [-a2], [b1,b2,b3,a1-a2],[a1,b1-a2,b2-a2,b3-a2], \
                    [a2,a2-b1+1,a2-b2+1,a2-b3+1],[-a1+a2+1], 1/z
                return T1, T2, T3
            v = hypercomb(h, [a1,a2,b1,b2,b3])
            if sum(ctx.is_real_type(u) for u in [a1,a2,b1,b2,b3,z]) == 6:
                v = ctx.re(v)
            return v
        except NoConvergence:
            pass
        finally:
            ctx.prec = orig

    #print "not using asymp"
    return ctx.hypsum(a1r+a2r, a1f+a2f, a1c+a2c, b1r+b2r+b3r, b1f+b2f+b3f,
        b1c+b2c+b3c, z)


@defun
def hyp2f0(ctx, a, b, z, **kwargs):
    """
    Hypergeometric 2F0.
    """
    ar, af, ac, a = ctx._hyp_parse_param(a)
    br, bf, bc, b = ctx._hyp_parse_param(b)
    z = ctx.convert(z)
    # TODO: also use series when it terminates
    if kwargs.get('force_series') or \
        ctx._hyp_check_convergence([a, b], [], z, ctx.prec+40) or \
        (ctx.isint(a) and -100 < a <= 0) or \
        (ctx.isint(b) and -100 < b <= 0):
        return ctx.hypsum(ar+br, af+bf, ac+bc, [], [], [], z)
    def h(a, b):
        w = ctx.sinpi(b)
        rz = -1/z
        T1 = ([ctx.pi,w,rz],[1,-1,a],[],[a-b+1,b],[a],[b],rz)
        T2 = ([-ctx.pi,w,rz],[1,-1,1+a-b],[],[a,2-b],[a-b+1],[2-b],rz)
        return T1, T2
    return ctx.hypercomb(h, [a, 1+a-b], check_cancellation=True)

@defun
def hyperu(ctx, a,b,z):
    ar, af, ac, a = ctx._hyp_parse_param(a)
    br, bf, bc, b = ctx._hyp_parse_param(b)
    z = ctx.convert(z)
    if not z:
        if ctx.re(b) <= 1:
            return ctx.gammaprod([1-b],[a-b+1])
        else:
            return ctx.inf + z
    bb = 1+a-b
    if ctx._hyp_check_convergence([a, bb], [], 1/z, ctx.prec+40):
        def h(a,b):
            rz = -1/z
            return [([z],[-a],[],[],[a,b],[],rz)]
        return hypercomb(h, [a,bb], force_series=True)
    else:
        def h(a,b):
            w = sinpi(b)
            T1 = ([pi,w],[1,-1],[],[a-b+1,b],[a],[b],z)
            T2 = ([-pi,w,z],[1,-1,1-b],[],[a,2-b],[a-b+1],[2-b],z)
            return T1, T2
        return hypercomb(h, [a,b], check_cancellation=True)

@defun
def _lower_gamma(ctx, z, b):
    return ctx.hyp1f1(1, 1+z, b) * b**z * ctx.exp(-b) / z

def _check_pos(x):
    try:
        return x > 0
    except TypeError:
        return False

@defun_wrapped
def gammainc(ctx, z, a=0, b=None, regularized=False):
    if b is None:
        b = ctx.inf
    ln = ctx.ln
    ei = ctx.ei
    if b == ctx.inf:
        if not a:
            v = ctx.gamma(z)
        else:
            if not z:
                # Reduces to exponential integral. Mind branch cuts.
                if _check_pos(a):
                    return -ei(-a)
                else:
                    return -ei(-a) + (ln(-a)-ln(-1/a))/2-ln(a)
            # XXX: avoid poles
            v = ctx.gamma(z) - ctx._lower_gamma(z, a)
    elif not a:
        v = ctx._lower_gamma(z, b)
    else:
        if not z:
            # Reduces to exponential integral
            if _check_pos(a) and _check_pos(b):
                return ei(-b) - ei(-a)
            else:
                return ei(-b)-ei(-a) + \
                    (ln(-a)-ln(-1/a))/2-ln(a) + \
                    (ln(-1/b)-ln(-b))/2+ln(b)
        # XXX: avoid poles
        v = ctx._lower_gamma(z, b) - ctx._lower_gamma(z, a)
    if regularized:
        return v / ctx.gamma(z)
    else:
        return v

@defun_wrapped
def _erf_complex(ctx, z):
    v = (2/ctx.sqrt(ctx.pi))*z * ctx.hyp1f1((1,2),(3,2), -z**2)
    if not z.real:
        v = v.imag*ctx.j
    return v

@defun_wrapped
def _erfc_complex(ctx, z):
    if ctx.re(z) > 2:
        v = ctx.exp(-z*z)/ctx.sqrt(ctx.pi)*ctx.hyperu((1,2),(1,2),z**2)
    else:
        v = 1 - ctx._erf_complex(z)
    if not z.real:
        v = 1+v.imag*ctx.j
    return v

@defun
def erf(ctx, z):
    z = ctx.convert(z)
    if hasattr(z, "_mpf_"):
        return ctx._erf(z)
    elif hasattr(z, "_mpc_"):
        if z.imag:
            return ctx._erf_complex(z)
        else:
            return ctx.mpc(ctx._erf(z.real))

@defun
def erfc(ctx, z):
    z = ctx.convert(z)
    if hasattr(z, "_mpf_"):
        return ctx._erfc(z)
    elif hasattr(z, "_mpc_"):
        if z.imag:
            return ctx._erfc_complex(z)
        else:
            return ctx.mpc(ctx._erfc(z.real))

@defun_wrapped
def erfi(ctx, z):
    if not z:
        return z
    v = (2/ctx.sqrt(ctx.pi)*z) * ctx.hyp1f1((1,2), (3,2), z**2)
    if not z.real:
        v = v.imag*ctx.j
    return v

@defun_wrapped
def erfinv(ctx, x):
    if x.imag or (x < -1) or (x > 1):
        return ctx.bad_domain("erfinv(x) is defined only for -1 <= x <= 1")
    if ctx.isnan(x): return x
    if not x: return x
    if x == 1: return ctx.inf
    if x == -1: return ctx.ninf
    if abs(x) < 0.9:
        a = 0.53728*x**3 + 0.813198*x
    else:
        # An asymptotic formula
        u = ctx.ln(2/ctx.pi/(abs(x)-1)**2)
        a = ctx.sign(x) * ctx.sqrt(u - ctx.ln(u))/ctx.sqrt(2)
    ctx.prec += 10
    return ctx.findroot(lambda t: ctx.erf(t)-x, a)

@defun_wrapped
def npdf(ctx, x, mu=0, sigma=1):
    sigma = ctx.convert(sigma)
    return ctx.exp(-(x-mu)**2/(2*sigma**2)) / (sigma*ctx.sqrt(2*ctx.pi))

@defun_wrapped
def ncdf(ctx, x, mu=0, sigma=1):
    a = (x-mu)/(sigma*ctx.sqrt(2))
    if a < 0:
        return ctx.erfc(-a)/2
    else:
        return (1+ctx.erf(a))/2

def ei_as(ctx, a):
    extra = 10
    ctx.dps += extra
    s = k = p = 1
    while abs(p) > ctx.eps:
        p = (p*k)/a
        s += p
        k += 1
    s = (s * ctx.exp(a))/a
    ctx.dps -= extra
    return s

@defun_wrapped
def ei(ctx, z):
    if z == ctx.inf:
        return z
    if z == ctx.ninf:
        return -ctx.zero
    if not z:
        return ctx.ninf
    if abs(z) > ctx.prec * 0.7 + 50:
        r = ei_as(ctx, z)
        if z.imag > 0:
            r += ctx.j*ctx.pi
        elif z.imag < 0:
            r -= ctx.j*ctx.pi
        return r
    v = z*hyp2f2(1,1,2,2,z) + ctx.euler
    if z.imag:
        v += (ctx.ln(z)-ctx.ln(1/z))/2
    else:
        v += ctx.ln(abs(z))
    return v

@defun_wrapped
def expint(ctx, *args):
    if len(args) == 1:
        n = ctx.one
        z = args[0]
    else:
        n, z = args
    if ctx.isnan(n) or ctx.isnan(z):
        return z*n
    if z == ctx.inf:
        return type(z)(0)
    if z == 0:
        # integral from 1 to infinity of t^n
        if ctx.re(n) <= 1:
            # TODO: reasonable sign of infinity
            return type(z)(ctx.inf)
        else:
            return ctx.one/(n-1)
    if n == 0:
        return ctx.exp(-z)/z
    if n == -1:
        return ctx.exp(-z)*(z+1)/z**2
    # Perturb if at pole
    m, d = ctx.nint_distance(n)
    if ctx.re(n) > 0:
        if d < -ctx.prec:
            h = +ctx.eps
            ctx.prec *= 2
            n += h
            z += h
        elif d < -4:
            ctx.prec -= d
    # XXX: this is entirely arbitrary
    ctx.prec += 4*(int(abs(n)) + int(abs(z)))
    # Main formula
    # TODO: use asymptotic expansions, either here or in gammainc
    return z**(n-1) * ctx.gammainc(1-n, z)

@defun_wrapped
def li(ctx, z):
    if not z:
        return z
    if z == 1:
        return ctx.ninf
    return ctx.ei(ctx.ln(z))

@defun_wrapped
def chi(ctx, z):
    if not z:
        return ctx.ninf
    if z == ctx.inf or z == ctx.ninf:
        return ctx.inf
    z2 = (z/2)**2
    return ctx.euler + ctx.ln(z) + z2*ctx.hyp2f3(1,1,2,2,(3,2),z2)

@defun_wrapped
def shi(ctx, z):
    if z == ctx.inf:
        return z
    if z == ctx.ninf:
        return z
    z2 = (z/2)**2
    return z*ctx.hyp1f2((1,2),(3,2),(3,2),z2)

@defun_wrapped
def fresnels(ctx, z):
    if z == ctx.inf:
        return ctx.mpf(0.5)
    if z == ctx.ninf:
        return ctx.mpf(-0.5)
    return ctx.pi*z**3/6*ctx.hypsum([[3,4]],[],[],[[3,2],[7,4]],[],[],-ctx.pi**2*z**4/16)

@defun_wrapped
def fresnelc(ctx, z):
    if z == ctx.inf:
        return ctx.mpf(0.5)
    if z == ctx.ninf:
        return ctx.mpf(-0.5)
    return z*ctx.hypsum([[1,4]],[],[],[[1,2],[5,4]],[],[],-ctx.pi**2*z**4/16)

@defun_wrapped
def airyai(ctx, z):
    if z == ctx.inf or z == ctx.ninf:
        return 1/z
    if z:
        # Account for exponential scaling
        ctx.prec += max(0, int(1.5*ctx.mag(z)))
    if z.real > 4:
        # We could still use 1F1, but it results in huge cancellation;
        # the following expansion is better
        w = z**1.5
        r = -ctx.mpf(3)/(4*w)
        v = ctx.exp(-2*w/3)/(2*ctx.sqrt(ctx.pi)*ctx.nthroot(z,4))
        v *= ctx.hyp2f0((1,6),(5,6),r)
        return v
    elif z.real > 1:
        # If not using asymptotic series:
        # cancellation: both terms are ~ 2^(z^1.5),
        # result is ~ 2^(-z^1.5), so need ~2*z^1.5 extra bits
        ctx.prec += 2*int(z.real**1.5)
    z3 = z**3 / 9
    a = ctx.hyp0f1((2,3), z3) / (ctx.cbrt(9) * ctx.gamma(ctx.mpf(2)/3))
    b = z * ctx.hyp0f1((4,3), z3) / (ctx.cbrt(3) * ctx.gamma(ctx.mpf(1)/3))
    return a - b

@defun_wrapped
def airybi(ctx, z):
    if z == ctx.inf:
        return z
    if z == ctx.ninf:
        return 1/z
    if z:
        # Account for exponential scaling
        ctx.prec += max(0, int(1.5*ctx.mag(z)))
    z3 = z**3 / 9
    rt = ctx.nthroot(3, 6)
    a = ctx.hyp0f1((2,3), z3) / (rt * ctx.gamma(ctx.mpf(2)/3))
    b = z * rt * ctx.hyp0f1((4,3), z3) / ctx.gamma(ctx.mpf(1)/3)
    return a + b

@defun
def agm(ctx, a, b=1):
    if b == 1:
        return ctx.agm1(a)
    a = ctx.convert(a)
    b = ctx.convert(b)
    prec, rounding = ctx._prec_rounding
    if hasattr(a, '_mpf_') and hasattr(b, '_mpf_'):
        try:
            v = libhyper.mpf_agm(a._mpf_, b._mpf_, prec, rounding)
            return ctx.make_mpf(v)
        except ComplexResult:
            pass
    if hasattr(a, '_mpf_'): a = (a._mpf_, libmpf.fzero)
    else: a = a._mpc_
    if hasattr(b, '_mpf_'): b = (b._mpf_, libmpf.fzero)
    else: b = b._mpc_
    return ctx.make_mpc(libhyper.mpc_agm(a, b, prec, rounding))

@defun_wrapped
def jacobi(ctx, n, a, b, x):
    return ctx.binomial(n+a,n) * ctx.hyp2f1(-n,1+n+a+b,a+1,(1-x)/2)

@defun_wrapped
def legendre(ctx, n, x):
    if ctx.isint(n):
        n = int(n)
    return ctx.hyp2f1(-n,n+1,1,(1-x)/2)

@defun_wrapped
def chebyt(ctx, n, x):
    return ctx.hyp2f1(-n,n,(1,2),(1-x)/2)

@defun_wrapped
def chebyu(ctx, n, x):
    return (n+1) * ctx.hyp2f1(-n, n+2, (3,2), (1-x)/2)

@defun
def j0(ctx, x):
    """Computes the Bessel function `J_0(x)`. See :func:`besselj`."""
    return ctx.besselj(0, x)

@defun
def j1(ctx, x):
    """Computes the Bessel function `J_1(x)`.  See :func:`besselj`."""
    return ctx.besselj(1, x)

@defun
def besselj(ctx, n, z, derivative=0):
    if type(n) is int:
        n_isint = True
    else:
        n = ctx.convert(n)
        n_isint = ctx.isint(n)
        if n_isint:
            n = int(n)
    if n_isint and n < 0:
        return (-1)**n * ctx.besselj(-n, z, derivative)
    z = ctx.convert(z)
    M = ctx.mag(z)
    if derivative:
        d = ctx.convert(derivative)
        # TODO: the integer special-casing shouldn't be necessary.
        # However, the hypergeometric series gets inaccurate for large d
        # because of inaccurate pole cancellation at a pole far from
        # zero (needs to be fixed in hypercomb or hypsum)
        if ctx.isint(d) and d >= 0:
            d = int(d)
            orig = ctx.prec
            try:
                ctx.prec += 15
                v = ctx.fsum((-1)**k * ctx.binomial(d,k) * ctx.besselj(2*k+n-d,z)
                    for k in range(d+1))
            finally:
                ctx.prec = orig
            v *= ctx.mpf(2)**(-d)
        else:
            def h(n,d):
                r = ctx.fmul(ctx.fmul(z, z, prec=ctx.prec+M), -0.25, exact=True)
                B = [0.5*(n-d+1), 0.5*(n-d+2), n+1]
                T = [([2,ctx.pi,z],[d-2*n,0.5,n-d],[n+1],B,[(n+1)*0.5,(n+2)*0.5],B,r)]
                return T
            v = ctx.hypercomb(h, [n,d])
    # Fast case: J_n(x), n int, appropriate magnitude for fixed-point calculation
    elif (not derivative) and n_isint and abs(M) < 10 and abs(n) < 20:
        prec, rounding = ctx._prec_rounding
        if hasattr(z, '_mpf_'):
            v = ctx.make_mpf(libhyper.mpf_besseljn(n, z._mpf_, prec, rounding))
        elif hasattr(z, '_mpc_'):
            v = ctx.make_mpc(libhyper.mpc_besseljn(n, z._mpc_, prec, rounding))
        else:
            raise TypeError
    elif not z:
        if not n:
            v = ctx.one + n+z
        elif ctx.re(n) > 0:
            v = n*z
        else:
            v = ctx.inf + z + n
    else:
        v = 0
        orig = ctx.prec
        try:
            # XXX: workaround for accuracy in low level hypergeometric series
            # when alternating, large arguments
            ctx.prec += min(3*abs(M), ctx.prec)
            w = ctx.fmul(z, 0.5, exact=True)
            def h(n):
                r = ctx.fneg(ctx.fmul(w, w, prec=ctx.prec+M), exact=True)
                return [([w], [n], [], [n+1], [], [n+1], r)]
            v = ctx.hypercomb(h, [n])
        finally:
            ctx.prec = orig
        v = +v
    return v

@defun
def besseli(ctx, n, z, derivative=0):
    n = ctx.convert(n)
    z = ctx.convert(z)
    if not z:
        if derivative:
            raise ValueError
        if not n:
            # I(0,0) = 1
            return 1+n+z
        if ctx.isint(n):
            return 0*(n+z)
        r = ctx.re(n)
        if r == 0:
            return ctx.nan*(n+z)
        elif r > 0:
            return 0*(n+z)
        else:
            return ctx.inf+(n+z)
    M = ctx.mag(z)
    if derivative:
        d = ctx.convert(derivative)
        def h(n,d):
            r = ctx.fmul(ctx.fmul(z, z, prec=ctx.prec+M), 0.25, exact=True)
            B = [0.5*(n-d+1), 0.5*(n-d+2), n+1]
            T = [([2,ctx.pi,z],[d-2*n,0.5,n-d],[n+1],B,[(n+1)*0.5,(n+2)*0.5],B,r)]
            return T
        v = ctx.hypercomb(h, [n,d])
    else:
        def h(n):
            w = ctx.fmul(z, 0.5, exact=True)
            r = ctx.fmul(w, w, prec=ctx.prec+M)
            return [([w], [n], [], [n+1], [], [n+1], r)]
        v = ctx.hypercomb(h, [n])
    return v

@defun_wrapped
def bessely(ctx, n, z, derivative=0):
    if not z:
        if derivative:
            # Not implemented
            raise ValueError
        if not n:
            # ~ log(z/2)
            return -ctx.inf + (n+z)
        if ctx.im(n):
            return nan * (n+z)
        r = ctx.re(n)
        q = n+0.5
        if ctx.isint(q):
            if n > 0:
                return -ctx.inf + (n+z)
            else:
                return 0 * (n+z)
        if r < 0 and int(ctx.floor(q)) % 2:
            return ctx.inf + (n+z)
        else:
            return ctx.ninf + (n+z)
    ctx.prec += 10
    m, d = ctx.nint_distance(n)
    if d < -ctx.prec:
        h = +ctx.eps
        ctx.prec *= 2
        n += h
    elif d < 0:
        ctx.prec -= d
    # TODO: avoid cancellation for imaginary arguments
    return (ctx.besselj(n,z,derivative)*ctx.cospi(n) - \
        ctx.besselj(-n,z,derivative))/ctx.sinpi(n)

@defun_wrapped
def besselk(ctx, n, z):
    if not z:
        return ctx.inf
    M = ctx.mag(z)
    if M < 1:
        # Represent as limit definition
        def h(n):
            r = (z/2)**2
            T1 = [z, 2], [-n, n-1], [n], [], [], [1-n], r
            T2 = [z, 2], [n, -n-1], [-n], [], [], [1+n], r
            return T1, T2
    # We could use the limit definition always, but it leads
    # to very bad cancellation (of exponentially large terms)
    # for large real z
    # Instead represent in terms of 2F0
    else:
        ctx.prec += M
        def h(n):
            return [([pi/2, z, exp(-z)], [0.5,-0.5,1], [], [], \
                [n+0.5, 0.5-n], [], -1/(2*z))]
    return ctx.hypercomb(h, [n], check_cancellation=True)

@defun_wrapped
def hankel1(ctx,n,x):
    return ctx.besselj(n,x) + ctx.j*ctx.bessely(n,x)

@defun_wrapped
def hankel2(ctx,n,x):
    return ctx.besselj(n,x) - ctx.j*ctx.bessely(n,x)

@defun
def whitm(ctx,k,m,z):
    k = ctx.convert(k)
    m = ctx.convert(m)
    z = ctx.convert(z)
    def h(k,m):
        return [([ctx.exp(-z/2), z], [1, m+0.5], [], [], [0.5+m-k], [1+2*m], z)]
    return ctx.hypercomb(h, [k,m])

@defun
def whitw(ctx,k,m,z):
    k = ctx.convert(k)
    m = ctx.convert(m)
    z = ctx.convert(z)
    def h(k,m):
        w = ctx.exp(-z/2)
        T1 = [w, z], [1, 0.5-m], [2*m], [0.5+m-k], [0.5-m-k], [1-2*m], z
        T2 = [w, z], [1, 0.5+m], [-2*m], [0.5-m-k], [0.5+m-k], [1+2*m], z
        return T1, T2
    return ctx.hypercomb(h, [k,m])

@defun
def struveh(ctx,n,z):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/StruveH/26/01/02/
    def h(n):
        return [([z/2, 0.5*ctx.sqrt(ctx.pi)], [n+1, -1], [], [n+1.5], [1], [1.5, n+1.5], -(z/2)**2)]
    return ctx.hypercomb(h, [n])

@defun
def struvel(ctx,n,z):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/StruveL/26/01/02/
    def h(n):
        return [([z/2, 0.5*ctx.sqrt(ctx.pi)], [n+1, -1], [], [n+1.5], [1], [1.5, n+1.5], (z/2)**2)]
    return ctx.hypercomb(h, [n])

@defun
def ber(ctx, n, z):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinBer2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [ctx.cospi(0.75*n), z/2], [1, n], [], [n+1], [], [0.5, 0.5*(n+1), 0.5*n+1], r
        T2 = [-ctx.sinpi(0.75*n), z/2], [1, n+2], [], [n+2], [], [1.5, 0.5*(n+3), 0.5*n+1], r
        return T1, T2
    return ctx.hypercomb(h, [n])

@defun
def bei(ctx, n, z):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinBei2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [ctx.cospi(0.75*n), z/2], [1, n+2], [], [n+2], [], [1.5, 0.5*(n+3), 0.5*n+1], r
        T2 = [ctx.sinpi(0.75*n), z/2], [1, n], [], [n+1], [], [0.5, 0.5*(n+1), 0.5*n+1], r
        return T1, T2
    return ctx.hypercomb(h, [n])

@defun
def ker(ctx, n, z):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinKer2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [2, z, 4*cospi(0.25*n)], [-n-3, n, 1], [-n], [], [], [0.5, 0.5*(1+n), 0.5*(n+2)], r
        T2 = [2, z, -sinpi(0.25*n)], [-n-3, 2+n, 1], [-n-1], [], [], [1.5, 0.5*(3+n), 0.5*(n+2)], r
        T3 = [2, z, 4*cospi(0.75*n)], [n-3, -n, 1], [n], [], [], [0.5, 0.5*(1-n), 1-0.5*n], r
        T4 = [2, z, -sinpi(0.75*n)], [n-3, 2-n, 1], [n-1], [], [], [1.5, 0.5*(3-n), 1-0.5*n], r
        return T1, T2, T3, T4
    return ctx.hypercomb(h, [n])

@defun
def kei(ctx, n, z):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinKei2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [-cospi(0.75*n), 2, z], [1, n-3, 2-n], [n-1], [], [], [1.5, 0.5*(3-n), 1-0.5*n], r
        T2 = [-sinpi(0.75*n), 2, z], [1, n-1, -n], [n], [], [], [0.5, 0.5*(1-n), 1-0.5*n], r
        T3 = [-sinpi(0.25*n), 2, z], [1, -n-1, n], [-n], [], [], [0.5, 0.5*(n+1), 0.5*(n+2)], r
        T4 = [-cospi(0.25*n), 2, z], [1, -n-3, n+2], [-n-1], [], [], [1.5, 0.5*(n+3), 0.5*(n+2)], r
        return T1, T2, T3, T4
    return ctx.hypercomb(h, [n])

@defun_wrapped
def lambertw(ctx, z, k=0, approx=None):
    if ctx.isnan(z):
        return z
    ctx.prec += 20
    # We must be extremely careful near the singularities at -1/e and 0
    u = ctx.exp(-1)
    if abs(z) <= u:
        if not z:
            # w(0,0) = 0; for all other branches we hit the pole
            if not k:
                return z
            return ctx.ninf
        if not k:
            w = z
        # For small real z < 0, the -1 branch behaves roughly like log(-z)
        elif k == -1 and not ctx.im(z) and ctx.re(z) < 0:
            w = ctx.ln(-z)
        # Use a simple asymptotic approximation.
        else:
            w = ctx.ln(z)
            # The branches are roughly logarithmic. This approximation
            # gets better for large |k|; need to check that this always
            # works for k ~= -1, 0, 1.
            if k: w += k * 2*ctx.pi*ctx.j
    elif k == 0 and ctx.im(z) and abs(z) <= 0.6:
        w = z
    else:
        if z == ctx.inf:
            if k == 0:
                return z
            else:
                return z + 2*k*ctx.pi*ctx.j
        if z == ctx.ninf:
            return (-z) + (2*k+1)*ctx.pi*ctx.j
        # Simple asymptotic approximation as above
        w = ctx.ln(z)
        if k: w += k * 2*ctx.pi*ctx.j
    # Use Halley iteration to solve w*exp(w) = z
    two = ctx.mpf(2)
    weps = ctx.ldexp(ctx.eps, 15)
    for i in xrange(100):
        ew = ctx.exp(w)
        wew = w*ew
        wewz = wew-z
        wn = w - wewz/(wew+ew-(w+two)*wewz/(two*w+two))
        if abs(wn-w) < weps*abs(wn):
            return wn
        else:
            w = wn
    print "Warning: Lambert W iteration failed to converge:", z
    return wn

@defun_wrapped
def barnesg(ctx, z):
    if ctx.isinf(z):
        if z == ctx.inf:
            return z
        return ctx.nan
    if ctx.isnan(z):
        return z
    if (not z.imag) and z.real <= 0 and ctx.isint(z.real):
        return z*0
    # Account for size (would not be needed if computing log(G))
    if abs(z) > 5:
        ctx.dps += 2*ctx.log(abs(z),2)
    # Estimate terms for asymptotic expansion
    N = ctx.dps // 2 + 5
    G = 1
    while ctx.re(z) < N:
        G /= ctx.gamma(z)
        z += 1
    z -= 1
    s = ctx.mpf(1)/12
    s -= ctx.log(ctx.glaisher)
    s += z*ctx.log(2*ctx.pi)/2
    s += (z**2/2-ctx.mpf(1)/12)*ctx.log(z)
    s -= 3*z**2/4
    z2k = z2 = z**2
    for k in xrange(1, N+1):
        t = ctx.bernoulli(2*k+2) / (4*k*(k+1)*z2k)
        if abs(t) < ctx.eps:
            #print k, N      # check how many terms were needed
            break
        z2k *= z2
        s += t
    #if k == N:
    #    print "warning: series for barnesg failed to converge"
    return G*ctx.exp(s)

@defun
def superfac(ctx, z):
    return ctx.barnesg(z+2)

@defun_wrapped
def hyperfac(ctx, z):
    # XXX: estimate needed extra bits accurately
    if z == ctx.inf:
        return z
    if abs(z) > 5:
        extra = 4*int(ctx.log(abs(z),2))
    else:
        extra = 0
    ctx.prec += extra
    if not z.imag and z.real < 0 and ctx.isint(z.real):
        n = int(ctx.re(z))
        h = ctx.hyperfac(-n-1)
        if ((n+1)//2) & 1:
            h = -h
        if ctx.is_complex_type(z):
            return h + 0j
        return h
    zp1 = z+1
    # Wrong branch cut
    #v = ctx.gamma(zp1)**z
    #ctx.prec -= extra
    #return v / ctx.barnesg(zp1)
    v = ctx.exp(z*ctx.loggamma(zp1))
    ctx.prec -= extra
    return v / ctx.barnesg(zp1)

@defun_wrapped
def loggamma(ctx, z):
    a = z.real
    b = z.imag
    if not b and a > 0:
        return ctx.ln(ctx.gamma(z))
    u = ctx.arg(z)
    w = ctx.ln(ctx.gamma(z))
    if b:
        gi = -b - u/2 + a*u + b*ctx.ln(abs(z))
        n = ctx.floor((gi-w.imag)/(2*ctx.pi)+0.5) * (2*ctx.pi)
        return w + n*ctx.j
    elif a < 0:
        n = int(ctx.floor(a))
        w += (n-(n%2))*ctx.pi*ctx.j
    return w

@defun_wrapped
def siegeltheta(ctx, t):
    if t.imag:
        # XXX: cancellation occurs
        a = ctx.loggamma(0.25+0.5j*t)
        b = ctx.loggamma(0.25-0.5j*t)
        return -ctx.ln(ctx.pi)/2*t - 0.5j*(a-b)
    else:
        if ctx.isinf(t):
            return t
        return ctx.loggamma(0.25+0.5j*t).imag - ctx.ln(ctx.pi)/2*t

@defun_wrapped
def grampoint(ctx, n):
    # ctxsymptotic expansion, from
    # http://mathworld.wolfram.com/GramPoint.html
    g = 2*ctx.pi*ctx.exp(1+ctx.lambertw((8*n+1)/(8*ctx.e)))
    return ctx.findroot(lambda t: ctx.siegeltheta(t)-ctx.pi*n, g)

@defun_wrapped
def siegelz(ctx, t):
    v = ctx.exp(ctx.j*ctx.siegeltheta(t))*ctx.zeta(0.5+ctx.j*t)
    if ctx.is_real_type(t):
        return v.real
    return v

_zeta_zeros = [
14.134725142,21.022039639,25.010857580,30.424876126,32.935061588,
37.586178159,40.918719012,43.327073281,48.005150881,49.773832478,
52.970321478,56.446247697,59.347044003,60.831778525,65.112544048,
67.079810529,69.546401711,72.067157674,75.704690699,77.144840069,
79.337375020,82.910380854,84.735492981,87.425274613,88.809111208,
92.491899271,94.651344041,95.870634228,98.831194218,101.317851006,
103.725538040,105.446623052,107.168611184,111.029535543,111.874659177,
114.320220915,116.226680321,118.790782866,121.370125002,122.946829294,
124.256818554,127.516683880,129.578704200,131.087688531,133.497737203,
134.756509753,138.116042055,139.736208952,141.123707404,143.111845808,
146.000982487,147.422765343,150.053520421,150.925257612,153.024693811,
156.112909294,157.597591818,158.849988171,161.188964138,163.030709687,
165.537069188,167.184439978,169.094515416,169.911976479,173.411536520,
174.754191523,176.441434298,178.377407776,179.916484020,182.207078484,
184.874467848,185.598783678,187.228922584,189.416158656,192.026656361,
193.079726604,195.265396680,196.876481841,198.015309676,201.264751944,
202.493594514,204.189671803,205.394697202,207.906258888,209.576509717,
211.690862595,213.347919360,214.547044783,216.169538508,219.067596349,
220.714918839,221.430705555,224.007000255,224.983324670,227.421444280,
229.337413306,231.250188700,231.987235253,233.693404179,236.524229666,
]

def _load_zeta_zeros(url):
    import urllib
    d = urllib.urlopen(url)
    L = [float(x) for x in d.readlines()]
    # Sanity check
    assert round(L[0]) == 14
    _zeta_zeros[:] = L

@defun
def zetazero(ctx, n, url='http://www.dtc.umn.edu/~odlyzko/zeta_tables/zeros1'):
    n = int(n)
    if n < 0:
        return zetazero(-n).conjugate()
    if n == 0:
        raise ValueError("n must be nonzero")
    if n > len(_zeta_zeros) and n <= 100000:
        _load_zeta_zeros(url)
    if n > len(_zeta_zeros):
        raise NotImplementedError("n too large for zetazeros")
    return ctx.mpc(0.5, ctx.findroot(ctx.siegelz, _zeta_zeros[n-1]))

@defun_wrapped
def riemannr(ctx, x):
    if x == 0:
        return ctx.zero
    # Check if a simple asymptotic estimate is accurate enough
    if abs(x) > 1000:
        a = ctx.li(x)
        b = 0.5*ctx.li(ctx.sqrt(x))
        if abs(b) < abs(a)*ctx.eps:
            return a
    if abs(x) < 0.01:
        # XXX
        ctx.prec += int(-ctx.log(abs(x),2))
    # Sum Gram's series
    s = t = ctx.one
    u = ctx.ln(x)
    k = 1
    while abs(t) > abs(s)*ctx.eps:
        t = t * u / k
        s += t / (k * ctx.zeta(k+1))
        k += 1
    return s

@defun_static
def primepi(x):
    x = int(x)
    if x < 2:
        return 0
    from gammazeta import list_primes
    return len(list_primes(x))

@defun_wrapped
def primepi2(ctx, x):
    x = int(x)
    if x < 2:
        return ctx.mpi(0,0)
    if x < 2657:
        return ctx.mpi(ctx.primepi(x))
    mid = ctx.li(x)
    # Schoenfeld's estimate for x >= 2657, assuming RH
    err = ctx.sqrt(x,rounding='u')*ctx.ln(x,rounding='u')/8/ctx.pi(rounding='d')
    a = ctx.floor((ctx.mpi(mid)-err).a, rounding='d')
    b = ctx.ceil((ctx.mpi(mid)+err).b, rounding='u')
    return ctx.mpi(a, b)

@defun_wrapped
def primezeta(ctx, s):
    if ctx.isnan(s):
        return s
    if ctx.re(s) <= 0:
        raise ValueError("prime zeta function defined only for re(s) > 0")
    if s == 1:
        return ctx.inf
    if s == 0.5:
        return ctx.mpc(ctx.ninf, ctx.pi)
    r = ctx.re(s)
    if r > ctx.prec:
        return 0.5**s
    else:
        wp = ctx.prec + int(r)
        def terms():
            orig = ctx.prec
            # zeta ~ 1+eps; need to set precision
            # to get logarithm accurately
            k = 0
            while 1:
                k += 1
                u = libintmath.moebius(k)
                if not u:
                    continue
                ctx.prec = wp
                t = u*ctx.ln(ctx.zeta(k*s))/k
                if not t:
                    return
                #print ctx.prec, ctx.nstr(t)
                ctx.prec = orig
                yield t
    return sum_accurately(ctx, terms)

@defun_wrapped
def bernpoly(ctx, n, z):
    n = int(n)
    assert n >= 0
    # XXX: optimize further
    if abs(z) > 2:
        s = t = ctx.one
        r = ctx.one/z
        for k in xrange(1,n+1):
            t = t*(n+1-k)/k*r
            if not (k > 2 and k & 1):
                u = t*ctx.bernoulli(k)
                s += u
                if abs(u) < ctx.eps:
                    break
        return z**n * s
    return sum(ctx.binomial(n,k)*ctx.bernoulli(k)*z**(n-k) for k in xrange(0,n+1))

# TODO: this should be implemented low-level
def polylog_series(ctx, s, z):
    tol = +ctx.eps
    l = ctx.zero
    k = 1
    zk = z
    while 1:
        term = zk / k**s
        l += term
        if abs(term) < tol:
            break
        zk *= z
        k += 1
    return l

def polylog_continuation(ctx, n, z):
    if n < 0:
        return z*0
    twopij = 2j * ctx.pi
    a = -twopij**n/ctx.fac(n) * ctx.bernpoly(n, ctx.ln(z)/twopij)
    if ctx.is_real_type(z) and z < 0:
        a = a.real
    if z.imag < 0 or (z.imag == 0 and z.real >= 1):
        a -= twopij*ctx.ln(z)**(n-1)/ctx.fac(n-1)
    return a

def polylog_unitcircle(ctx, n, z):
    tol = +ctx.eps
    if n > 1:
        l = ctx.zero
        logz = ctx.ln(z)
        logmz = ctx.one
        m = 0
        while 1:
            if (n-m) != 1:
                term = ctx.zeta(n-m) * logmz / ctx.fac(m)
                if term and abs(term) < tol:
                    break
                l += term
            logmz *= logz
            m += 1
        l += ctx.ln(z)**(n-1)/ctx.fac(n-1)*(ctx.harmonic(n-1)-ctx.ln(-ctx.ln(z)))
    elif n < 1:  # else
        l = ctx.fac(-n)*(-ctx.ln(z))**(n-1)
        logz = ctx.ln(z)
        logkz = ctx.one
        k = 0
        while 1:
            b = ctx.bernoulli(k-n+1)
            if b:
                term = b*logkz/(ctx.fac(k)*(k-n+1))
                if abs(term) < tol:
                    break
                l -= term
            logkz *= logz
            k += 1
    else:
        raise ValueError
    if ctx.is_real_type(z) and z < 0:
        l = l.real
    return l

@defun_wrapped
def polylog(ctx, s, z):
    if z == 1:
        return ctx.zeta(s)
    if z == -1:
        return -ctx.altzeta(s)
    if s == 0:
        return z/(1-z)
    if s == 1:
        return -ctx.ln(1-z)
    if s == -1:
        return z/(1-z)**2
    if abs(z) <= 0.75 or (not ctx.isint(s) and abs(z) < 0.99):
        return polylog_series(ctx, s, z)
    if abs(z) >= 1.4 and ctx.isint(s):
        return (-1)**(s+1)*polylog_series(ctx, s, 1/z) + polylog_continuation(ctx, s, z)
    if ctx.isint(s):
        return polylog_unitcircle(ctx, int(s), z)
    raise NotImplementedError("polylog for arbitrary s and z")
    # This could perhaps be used in some cases
    #from quadrature import quad
    #return quad(lambda t: t**(s-1)/(exp(t)/z-1),[0,inf])/gamma(s)

# Experimental code; could be used elsewhere
def sum_accurately(ctx, terms, check_step=1):
    orig = ctx.prec
    extra = 10
    while 1:
        ctx.prec = orig + extra
        max_term = ctx.ninf
        s = 0
        k = 0
        for term in terms():
            s += term
            if not k % check_step and term:
                abs_term = abs(term)
                abs_sum = abs(s)
                max_term = max(max_term, abs_term)
                if abs_term <= ctx.eps*abs_sum:
                    break
            k += 1
        if abs_sum:
            cancellation = int(max(0,ctx.log(max_term/abs_sum,2)))
        else:
            cancellation = ctx.prec
        if cancellation < extra:
            break
        else:
            extra += cancellation
    return s

@defun_wrapped
def bell(ctx, n, x=1):
    x = ctx.convert(x)
    if not n:
        if ctx.isnan(x):
            return x
        return type(x)(1)
    if ctx.isinf(x) or ctx.isinf(n) or ctx.isnan(x) or ctx.isnan(n):
        return x**n
    if n == 1: return x
    if n == 2: return x*(x+1)
    if x == 0: return ctx.sincpi(n)
    return _polyexp(ctx, n, x, True) / ctx.exp(x)

def _polyexp(ctx, n, x, extra=False):
    def _terms():
        if extra:
            yield ctx.sincpi(n)
        t = x
        k = 1
        while 1:
            yield k**n * t
            k += 1
            t = t*x/k
    return sum_accurately(ctx, _terms, check_step=4)

@defun_wrapped
def polyexp(ctx, s, z):
    if ctx.isinf(z) or ctx.isinf(s) or ctx.isnan(z) or ctx.isnan(s):
        return z**s
    if z == 0: return z*s
    if s == 0: return ctx.expm1(z)
    if s == 1: return ctx.exp(z)*z
    if s == 2: return ctx.exp(z)*z*(z+1)
    return _polyexp(ctx, s, z)

@defun_wrapped
def cyclotomic(ctx, n, z):
    n = int(n)
    assert n >= 0
    p = ctx.one
    if n == 0:
        return p
    if n == 1:
        return z - p
    if n == 2:
        return z + p
    # Use divisor product representation. Unfortunately, this sometimes
    # includes singularities for roots of unity, which we have to cancel out.
    # Matching zeros/poles pairwise, we have (1-z^a)/(1-z^b) ~ a/b + O(z-1).
    moebius = libintmath.moebius
    a_prod = 1
    b_prod = 1
    num_zeros = 0
    num_poles = 0
    for d in range(1,n+1):
        if not n % d:
            w = moebius(n//d)
            # Use powm1 because it is important that we get 0 only
            # if it really is exactly 0
            b = -ctx.powm1(z, d)
            if b:
                p *= b**w
            else:
                if w == 1:
                    a_prod *= d
                    num_zeros += 1
                elif w == -1:
                    b_prod *= d
                    num_poles += 1
    #print n, num_zeros, num_poles
    if num_zeros:
        if num_zeros > num_poles:
            p *= 0
        else:
            p *= a_prod
            p /= b_prod
    return p


if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    try:
        import psyco; psyco.full()
    except ImportError:
        pass
    import sys
    filter = []
    for i, arg in enumerate(sys.argv):
        if 'functions.py' in arg:
            filter = sys.argv[i+1:]
            break
    import doctest
    globs = globals().copy()
    for obj in globs: #sorted(globs.keys()):
        if filter:
            if not sum([pat in obj for pat in filter]):
                continue
        print obj
        doctest.run_docstring_examples(globs[obj], {})

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
erf = def_mp_builtin("erf", libhyper.mpf_erf, libhyper.mpc_erf, None, "Error function, erf(z)")
erfc = def_mp_builtin("erfc", libhyper.mpf_erfc, libhyper.mpc_erfc, None, "Complementary error function, erfc(z) = 1-erf(z)")
ci = def_mp_builtin('ci', libhyper.mpf_ci, libhyper.mpc_ci, None, "")
si = def_mp_builtin('si', libhyper.mpf_si, libhyper.mpc_si, None, "")
ellipk = def_mp_builtin('ellipk', libhyper.mpf_ellipk, libhyper.mpc_ellipk, None, "")
ellipe = def_mp_builtin('ellipe', libhyper.mpf_ellipe, libhyper.mpc_ellipe, None, "")
agm1 = def_mp_builtin('agm1', libhyper.mpf_agm1, libhyper.mpc_agm1, None, "Fast alias for agm(1,a) = agm(a,1)")

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
            args = [mpmathify(z) for z in args]
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
def sinc(A, x):
    if A.isinf(x):
        return 1/x
    if not x:
        return x+1
    return A.sin(x)/x

@defun_wrapped
def sincpi(A, x):
    if A.isinf(x):
        return 1/x
    if not x:
        return x+1
    return A.sinpi(x)/(A.pi*x)

@defun
def nthroot(A, x, n):
    x = A.convert(x)
    n = int(n)
    if hasattr(x, '_mpf_'):
        try:
            return A.make_mpf(libelefun.mpf_nthroot(x._mpf_, n, *A._prec_rounding))
        except ComplexResult:
            if A.trap_complex:
                raise
            x = (x._mpf_, libmpf.fzero)
    else:
        x = x._mpc_
    return A.make_mpc(libmpc.mpc_nthroot(x, n, *A._prec_rounding))

@defun
def hypot(A, x, y):
    r"""
    Computes the Euclidean norm of the vector `(x, y)`, equal
    to `\sqrt{x^2 + y^2}`. Both `x` and `y` must be real."""
    x = A.convert(x)
    y = A.convert(y)
    return A.make_mpf(libmpf.mpf_hypot(x._mpf_, y._mpf_, *A._prec_rounding))

@defun
def ldexp(A, x, n):
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
    x = A.convert(x)
    return A.make_mpf(libmpf.mpf_shift(x._mpf_, n))

@defun
def frexp(A, x):
    r"""
    Given a real number `x`, returns `(y, n)` with `y \in [0.5, 1)`,
    `n` a Python integer, and such that `x = y 2^n`. No rounding is
    performed.

        >>> from mpmath import *
        >>> frexp(7.5)
        (mpf('0.9375'), 3)

    """
    x = A.convert(x)
    y, n = libmpf.mpf_frexp(x._mpf_)
    return A.make_mpf(y), n

@defun
def sign(A, x):
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
    x = A.convert(x)
    if not x or A.isnan(x):
        return x
    if A.is_real_type(x):
        return A.mpf(cmp(x, 0))
    return x / abs(x)

@defun
def arg(A, x):
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
    x = A.convert(x)
    return A.atan2(x.imag, x.real)

@defun
def fabs(A, x):
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
    return abs(A.convert(x))

@defun
def re(A, x):
    r"""
    Returns the real part of `x`, `\Re(x)`. Unlike ``x.real``,
    :func:`re` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> re(3)
        mpf('3.0')
        >>> re(-1+4j)
        mpf('-1.0')
    """
    return A.convert(x).real

@defun
def im(A, x):
    r"""
    Returns the imaginary part of `x`, `\Im(x)`. Unlike ``x.imag``,
    :func:`im` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> im(3)
        mpf('0.0')
        >>> im(-1+4j)
        mpf('4.0')

    """
    return A.convert(x).imag

@defun
def conj(A, x):
    r"""
    Returns the complex conjugate of `x`, `\overline{x}`. Unlike
    ``x.conjugate()``, :func:`im` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> conj(3)
        mpf('3.0')
        >>> conj(-1+4j)
        mpc(real='-1.0', imag='-4.0')

    """
    return A.convert(x).conjugate()


@defun
def log(A, x, b=None):
    if b is None:
        return ln(x)
    wp = A.prec + 20
    return A.ln(x, prec=wp) / A.ln(b, prec=wp)

@defun
def log10(A, x):
    r"""
    Computes the base-10 logarithm of `x`, `\log_{10}(x)`. ``log10(x)``
    is equivalent to ``log(x, 10)``.
    """
    return A.log(x, 10)

@defun
def power(A, x, y):
    return A.convert(x) ** A.convert(y)

@defun
def modf(A,x,y):
    return A.convert(x) % A.convert(y)

@defun
def degrees(A,x):
    return x / A.degree

@defun
def radians(A,x):
    return x * A.degree

@defun
def atan2(A, y, x):
    x = A.convert(x)
    y = A.convert(y)
    return A.make_mpf(libelefun.mpf_atan2(y._mpf_, x._mpf_, *A._prec_rounding))

@defun
def psi(A, m, z):
    z = A.convert(z)
    m = int(m)
    if A.is_real_type(z):
        return A.make_mpf(gammazeta.mpf_psi(m, z._mpf_, *A._prec_rounding))
    else:
        return A.make_mpc(gammazeta.mpc_psi(m, z._mpc_, *A._prec_rounding))

@defun
def psi0(A, z):
    """Shortcut for psi(0,z) (the digamma function)"""
    return A.psi(0, z)

@defun
def psi1(A, z):
    """Shortcut for psi(1,z) (the trigamma function)"""
    return A.psi(1, z)

@defun
def psi2(A, z):
    """Shortcut for psi(2,z) (the tetragamma function)"""
    return A.psi(2, z)

@defun
def psi3(A, z):
    """Shortcut for psi(3,z) (the pentagamma function)"""
    return A.psi(3, z)

polygamma = MultiPrecisionArithmetic.polygamma = psi
digamma = MultiPrecisionArithmetic.digamma = psi0
trigamma = MultiPrecisionArithmetic.trigamma = psi1
tetragamma = MultiPrecisionArithmetic.tetragamma = psi2
pentagamma = MultiPrecisionArithmetic.pentagamma = psi3

@defun
def bernoulli(A, n):
    return A.make_mpf(gammazeta.mpf_bernoulli(int(n), *A._prec_rounding))

bernfrac = defun_static(gammazeta.bernfrac)

@defun
def stieltjes(A, n, a=1):
    n = A.convert(n)
    a = A.convert(a)
    if n < 0:
        return A.bad_domain("Stieltjes constants defined for n >= 0")
    if hasattr(A, "stieltjes_cache"):
        stieltjes_cache = A.stieltjes_cache
    else:
        stieltjes_cache = A.stieltjes_cache = {}
    if a == 1:
        if n == 0:
            return +A.euler
        if n in stieltjes_cache:
            prec, s = stieltjes_cache[n]
            if prec >= A.prec:
                return +s
    mag = 1
    def f(x):
        xa = x/a
        v = (xa-A.j)*A.ln(a-A.j*x)**n/(1+xa**2)/(A.exp(2*A.pi*x)-1)
        return v.real / mag
    orig = A.prec
    try:
        # Normalize integrand by approx. magnitude to
        # speed up quadrature (which uses absolute error)
        if n > 50:
            A.prec = 20
            mag = A.quad(f, [0,A.inf], maxdegree=3)
        A.prec = orig + 10 + int(n**0.5)
        s = A.quad(f, [0,A.inf], maxdegree=20)
        v = A.ln(a)**n/(2*a) - A.ln(a)**(n+1)/(n+1) + 2*s/a*mag
    finally:
        A.prec = orig
    if a == 1 and A.isint(n):
        stieltjes_cache[n] = (A.prec, v)
    return +v

@defun
def gammaprod(A, a, b):
    a = [A.convert(x) for x in a]
    b = [A.convert(x) for x in b]
    poles_num = []
    poles_den = []
    regular_num = []
    regular_den = []
    for x in a: [regular_num, poles_num][A.isnpint(x)].append(x)
    for x in b: [regular_den, poles_den][A.isnpint(x)].append(x)
    # One more pole in numerator or denominator gives 0 or inf
    if len(poles_num) < len(poles_den): return A.zero
    if len(poles_num) > len(poles_den): return A.inf
    # All poles cancel
    # lim G(i)/G(j) = (-1)**(i+j) * gamma(1-j) / gamma(1-i)
    p = A.one
    orig = A.prec
    try:
        A.prec = orig + 15
        while poles_num:
            i = poles_num.pop()
            j = poles_den.pop()
            p *= (-1)**(i+j) * A.gamma(1-j) / A.gamma(1-i)
        for x in regular_num: p *= A.gamma(x)
        for x in regular_den: p /= A.gamma(x)
    finally:
        A.prec = orig
    return +p

@defun
def beta(A, x, y):
    r"""
    Computes the beta function,
    `B(x,y) = \Gamma(x) \Gamma(y) / \Gamma(x+y)`.
    The beta function is also commonly defined by the integral
    representation

    .. math ::

        B(x,y) = \int_0^1 t^{x-1} (1-t)^{y-1} \, dt

    **Examples**

    For integer and half-integer arguments where all three gamma
    functions are finite, the beta function becomes either rational
    number or a rational multiple of `\pi`::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print beta(5, 2)
        0.0333333333333333
        >>> print beta(1.5, 2)
        0.266666666666667
        >>> print 16*beta(2.5, 1.5)
        3.14159265358979

    Where appropriate, :func:`beta` evaluates limits. A pole
    of the beta function is taken to result in ``+inf``::

        >>> print beta(-0.5, 0.5)
        0.0
        >>> print beta(-3, 3)
        -0.333333333333333
        >>> print beta(-2, 3)
        +inf
        >>> print beta(inf, 1)
        0.0
        >>> print beta(inf, 0)
        nan

    :func:`beta` supports complex numbers and arbitrary precision
    evaluation::

        >>> print beta(1, 2+j)
        (0.4 - 0.2j)
        >>> mp.dps = 25
        >>> print beta(j,0.5)
        (1.079424249270925780135675 - 1.410032405664160838288752j)
        >>> mp.dps = 50
        >>> print beta(pi, e)
        0.037890298781212201348153837138927165984170287886464

    Various integrals can be computed by means of the
    beta function::

        >>> mp.dps = 15
        >>> print quad(lambda t: t**2.5*(1-t)**2, [0, 1])
        0.0230880230880231
        >>> print beta(3.5, 3)
        0.0230880230880231
        >>> print quad(lambda t: sin(t)**4 * sqrt(cos(t)), [0, pi/2])
        0.319504062596158
        >>> print beta(2.5, 0.75)/2
        0.319504062596158

    """
    x = A.convert(x)
    y = A.convert(y)
    if A.isinf(y):
        x, y = y, x
    if A.isinf(x):
        if x == A.inf and not y.imag:
            if y == A.ninf:
                return A.nan
            if y > 0:
                return A.zero
            if A.isint(y):
                return A.nan
            if y < 0:
                return A.sign(A.gamma(y)) * A.inf
        return A.nan
    return A.gammaprod([x, y], [x+y])

@defun
def binomial(A, n, k):
    return A.gammaprod([n+1], [k+1, n-k+1])

@defun
def rf(A, x, n):
    return A.gammaprod([x+n], [x])

@defun
def ff(A, x, n):
    return A.gammaprod([x+1], [x-n+1])

@defun_wrapped
def fac2(A, x):
    if A.isinf(x):
        if x == A.inf:
            return x
        return A.nan
    return 2**(x/2)*(A.pi/2)**((A.cospi(x)-1)/4)*A.gamma(x/2+1)


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

@defun
def _hyp_parse_param(A, x):
    if isinstance(x, tuple):
        p, q = x
        return [[p, q]], [], []
    if isinstance(x, (int, long)):
        return [[x, 1]], [], []
    x = A.convert(x)
    if hasattr(x, '_mpf_'):
        return [], [x._mpf_], []
    if hasattr(x, '_mpc_'):
        return [], [], [x._mpc_]

def _as_num(x):
    if isinstance(x, list):
        return _mpq(x)
    return x

@defun
def hypsum(A, ar, af, ac, br, bf, bc, x):
    prec, rnd = A._prec_rounding
    if hasattr(x, '_mpf_') and not (ac or bc):
        v = libhyper.hypsum_internal(ar, af, ac, br, bf, bc, x._mpf_, None, prec, rnd)
        return A.make_mpf(v)
    else:
        if hasattr(x, '_mpc_'):
            re, im = x._mpc_
        else:
            re, im = x._mpf_, libmpf.fzero
        v = libhyper.hypsum_internal(ar, af, ac, br, bf, bc, re, im, prec, rnd)
        return A.make_mpc(v)

@defun
def eval_hyp2f1(A,a,b,c,z):
    prec, rnd = A._prec_rounding
    ar, af, ac = A._hyp_parse_param(a)
    br, bf, bc = A._hyp_parse_param(b)
    cr, cf, cc = A._hyp_parse_param(c)
    absz = abs(z)
    if absz == 1:
        # TODO: determine whether it actually does, and otherwise
        # return infinity instead
        print "Warning: 2F1 might not converge for |z| = 1"
    if absz <= 1:
        # All rational
        if ar and br and cr:
            return A.sum_hyp2f1_rat(ar[0], br[0], cr[0], z)
        return A.hypsum(ar+br, af+bf, ac+bc, cr, cf, cc, z)
    # Use 1/z transformation
    a = (ar and _as_num(ar[0])) or A.convert(a)
    b = (br and _as_num(br[0])) or A.convert(b)
    c = (cr and _as_num(cr[0])) or A.convert(c)
    orig = A.prec
    try:
        A.prec = orig + 15
        h1 = A.eval_hyp2f1(a, mpq_1-c+a, mpq_1-b+a, 1/z)
        h2 = A.eval_hyp2f1(b, mpq_1-c+b, mpq_1-a+b, 1/z)
        #s1 = G(c)*G(b-a)/G(b)/G(c-a) * (-z)**(-a) * h1
        #s2 = G(c)*G(a-b)/G(a)/G(c-b) * (-z)**(-b) * h2
        f1 = A.gammaprod([c,b-a],[b,c-a])
        f2 = A.gammaprod([c,a-b],[a,c-b])
        s1 = f1 * (-z)**(mpq_0-a) * h1
        s2 = f2 * (-z)**(mpq_0-b) * h2
        v = s1 + s2
    finally:
        A.prec = orig
    return +v

@defun
def sum_hyp0f1_rat(A, a, z):
    prec, rnd = A._prec_rounding
    if hasattr(z, "_mpf_"):
        return A.make_mpf(libhyper.mpf_hyp0f1_rat(a, z._mpf_, prec, rnd))
    else:
        return A.make_mpc(libhyper.mpc_hyp0f1_rat(a, z._mpc_, prec, rnd))

@defun
def sum_hyp1f1_rat(A, a, b, z):
    prec, rnd = A._prec_rounding
    if hasattr(z, "_mpf_"):
        return A.make_mpf(libhyper.mpf_hyp1f1_rat(a, b, z._mpf_, prec, rnd))
    else:
        return A.make_mpc(libhyper.mpc_hyp1f1_rat(a, b, z._mpc_, prec, rnd))

@defun
def sum_hyp2f1_rat(A, a, b, c, z):
    prec, rnd = A._prec_rounding
    if hasattr(z, "_mpf_"):
        return A.make_mpf(libhyper.mpf_hyp2f1_rat(a, b, c, z._mpf_, prec, rnd))
    else:
        return A.make_mpc(libhyper.mpc_hyp2f1_rat(a, b, c, z._mpc_, prec, rnd))


#---------------------------------------------------------------------------#
#                      And now the user-friendly versions                   #
#---------------------------------------------------------------------------#

@defun
def hyper(A, a_s, b_s, z):
    p = len(a_s)
    q = len(b_s)
    z = A.convert(z)
    degree = p, q
    if degree == (0, 1):
        br, bf, bc = A._hyp_parse_param(b_s[0])
        if br:
            return A.sum_hyp0f1_rat(br[0], z)
        return A.hypsum([], [], [], br, bf, bc, z)
    if degree == (1, 1):
        ar, af, ac = A._hyp_parse_param(a_s[0])
        br, bf, bc = A._hyp_parse_param(b_s[0])
        if ar and br:
            a, b = ar[0], br[0]
            return A.sum_hyp1f1_rat(a, b, z)
        return A.hypsum(ar, af, ac, br, bf, bc, z)
    if degree == (2, 1):
        return A.eval_hyp2f1(a_s[0], a_s[1], b_s[0], z)
    ars, afs, acs, brs, bfs, bcs = [], [], [], [], [], []
    for a in a_s:
        r, f, c = A._hyp_parse_param(a)
        ars += r
        afs += f
        acs += c
    for b in b_s:
        r, f, c = A._hyp_parse_param(b)
        brs += r
        bfs += f
        bcs += c
    return A.hypsum(ars, afs, acs, brs, bfs, bcs, z)

@defun
def hyp0f1(A, a, z):
    r"""Hypergeometric function `\,_0F_1`. ``hyp0f1(a,z)`` is equivalent
    to ``hyper([],[a],z)``; see documentation for :func:`hyper` for more
    information."""
    return A.hyper([], [a], z)

@defun
def hyp1f1(A,a,b,z):
    r"""Hypergeometric function `\,_1F_1`. ``hyp1f1(a,b,z)`` is equivalent
    to ``hyper([a],[b],z)``; see documentation for :func:`hyper` for more
    information."""
    return A.hyper([a], [b], z)

@defun
def hyp2f1(A,a,b,c,z):
    r"""Hypergeometric function `\,_2F_1`. ``hyp2f1(a,b,c,z)`` is equivalent
    to ``hyper([a,b],[c],z)``; see documentation for :func:`hyper` for more
    information."""
    return A.hyper([a,b], [c], z)

@defun
def _lower_gamma(A, z, b):
    return A.hyp1f1(1, 1+z, b) * b**z * A.exp(-b) / z

def _check_pos(x):
    try:
        return x > 0
    except TypeError:
        return False

@defun_wrapped
def gammainc(A, z, a=0, b=None, regularized=False):
    if b is None:
        b = A.inf
    ln = A.ln
    ei = A.ei
    if b == A.inf:
        if not a:
            v = A.gamma(z)
        else:
            if not z:
                # Reduces to exponential integral. Mind branch cuts.
                if _check_pos(a):
                    return -ei(-a)
                else:
                    return -ei(-a) + (ln(-a)-ln(-1/a))/2-ln(a)
            # XXX: avoid poles
            v = A.gamma(z) - A._lower_gamma(z, a)
    elif not a:
        v = A._lower_gamma(z, b)
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
        v = A._lower_gamma(z, b) - A._lower_gamma(z, a)
    if regularized:
        return v / A.gamma(z)
    else:
        return v

@defun_wrapped
def erfi(A, z):
    return (2/A.sqrt(A.pi)*z) * A.sum_hyp1f1_rat((1,2),(3,2), z**2)

@defun_wrapped
def erfinv(A, x):
    if x.imag or (x < -1) or (x > 1):
        return A.bad_domain("erfinv(x) is defined only for -1 <= x <= 1")
    if A.isnan(x): return x
    if not x: return x
    if x == 1: return A.inf
    if x == -1: return A.ninf
    if abs(x) < 0.9:
        a = 0.53728*x**3 + 0.813198*x
    else:
        # An asymptotic formula
        u = A.ln(2/A.pi/(abs(x)-1)**2)
        a = A.sign(x) * A.sqrt(u - A.ln(u))/A.sqrt(2)
    return A.findroot(lambda t: A.erf(t)-x, a)

@defun_wrapped
def npdf(A, x, mu=0, sigma=1):
    sigma = A.convert(sigma)
    return A.exp(-(x-mu)**2/(2*sigma**2)) / (sigma*A.sqrt(2*A.pi))

@defun_wrapped
def ncdf(A, x, mu=0, sigma=1):
    a = (x-mu)/(sigma*A.sqrt(2))
    if a < 0:
        return A.erfc(-a)/2
    else:
        return (1+A.erf(a))/2

def ei_as(A, a):
    extra = 10
    A.dps += extra
    s = k = p = 1
    while abs(p) > A.eps:
        p = (p*k)/a
        s += p
        k += 1
    s = (s * A.exp(a))/a
    A.dps -= extra
    return s

@defun_wrapped
def ei(A, z):
    if z == A.inf:
        return z
    if z == A.ninf:
        return -A.zero
    if not z:
        return A.ninf
    if abs(z) > A.prec * 0.7 + 50:
        r = ei_as(A, z)
        if z.imag > 0:
            r += A.j*A.pi
        elif z.imag < 0:
            r -= A.j*A.pi
        return r
    v = z*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1]],[],[],z) + \
        (A.ln(z)-A.ln(1/z))/2 + A.euler
    if A.is_real_type(z) and z < 0:
        return v.real
    return v

@defun_wrapped
def li(A, z):
    if not z:
        return z
    if z == 1:
        return A.ninf
    return A.ei(A.ln(z))

@defun_wrapped
def chi(A, z):
    if not z:
        return A.ninf
    z2 = (z/2)**2
    return A.euler + A.ln(z) + \
        z2*A.hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1],[3,2]],[],[],z2)

@defun_wrapped
def shi(A, z):
    z2 = (z/2)**2
    return z*A.hypsum([[1,2]],[],[],[[3,2],[3,2]],[],[],z2)

@defun_wrapped
def fresnels(A, z):
    if z == A.inf:
        return A.mpf(0.5)
    if z == A.ninf:
        return A.mpf(-0.5)
    return A.pi*z**3/6*A.hypsum([[3,4]],[],[],[[3,2],[7,4]],[],[],-A.pi**2*z**4/16)

@defun_wrapped
def fresnelc(A, z):
    if z == A.inf:
        return A.mpf(0.5)
    if z == A.ninf:
        return A.mpf(-0.5)
    return z*A.hypsum([[1,4]],[],[],[[1,2],[5,4]],[],[],-A.pi**2*z**4/16)

@defun_wrapped
def airyai(A, z):
    if z == A.inf or z == A.ninf:
        return 1/z
    if z.real > 2:
        # cancellation: both terms are ~ 2^(z^1.5),
        # result is ~ 2^(-z^1.5), so need ~2*z^1.5 extra bits
        A.prec += 2*int(z.real**1.5)
    z3 = z**3 / 9
    a = A.sum_hyp0f1_rat((2,3), z3) / (A.cbrt(9) * A.gamma(A.mpf(2)/3))
    b = z * A.sum_hyp0f1_rat((4,3), z3) / (A.cbrt(3) * A.gamma(A.mpf(1)/3))
    return a - b

@defun_wrapped
def airybi(A, z):
    if z == A.inf:
        return z
    if z == A.ninf:
        return 1/z
    z3 = z**3 / 9
    rt = A.nthroot(3, 6)
    a = A.sum_hyp0f1_rat((2,3), z3) / (rt * A.gamma(A.mpf(2)/3))
    b = z * rt * A.sum_hyp0f1_rat((4,3), z3) / A.gamma(A.mpf(1)/3)
    return a + b

@defun
def agm(A, a, b=1):
    if b == 1:
        return A.agm1(a)
    a = A.convert(a)
    b = A.convert(b)
    prec, rounding = A._prec_rounding
    if hasattr(a, '_mpf_') and hasattr(b, '_mpf_'):
        try:
            v = libhyper.mpf_agm(a._mpf_, b._mpf_, prec, rounding)
            return A.make_mpf(v)
        except ComplexResult:
            pass
    if hasattr(a, '_mpf_'): a = (a._mpf_, libmpf.fzero)
    else: a = a._mpc_
    if hasattr(b, '_mpf_'): b = (b._mpf_, libmpf.fzero)
    else: b = b._mpc_
    return A.make_mpc(libhyper.mpc_agm(a, b, prec, rounding))

@defun_wrapped
def jacobi(A, n, a, b, x):
    return A.binomial(n+a,n) * A.hyp2f1(-n,1+n+a+b,a+1,(1-x)/2)

@defun_wrapped
def legendre(A, n, x):
    if A.isint(n):
        n = int(n)
    if x == -1:
        # TODO: hyp2f1 should handle this
        if A.isint(n):
            return (-1)**(n + (n>=0)) * A.mpf(-1)
        if not int(A.floor(A.re(n))) % 2:
            return A.ninf
        return A.inf
    return A.hyp2f1(-n,n+1,1,(1-x)/2)

@defun_wrapped
def chebyt(A, n, x):
    return A.hyp2f1(-n,n,0.5,(1-x)/2)

@defun_wrapped
def chebyu(A, n, x):
    return (n+1) * A.hyp2f1(-n, n+2, 1.5, (1-x)/2)

@defun_wrapped
def _besselj(A, v, x):
    hx = x/2
    return hx**v * A.hyp0f1(v+1, -hx**2) / A.factorial(v)

@defun
def besselj(A, v, x):
    if A.isint(v):
        x = A.convert(x)
        v = int(v)
        prec, rounding = A._prec_rounding
        if hasattr(x, '_mpf_'):
            return A.make_mpf(libhyper.mpf_besseljn(v, x._mpf_, prec, rounding))
        if hasattr(x, '_mpc_'):
            return A.make_mpc(libhyper.mpc_besseljn(v, x._mpc_, prec, rounding))
    return A._besselj(v, x)

@defun
def j0(A, x):
    """Computes the Bessel function `J_0(x)`. See :func:`besselj`."""
    return A.besselj(0, x)

@defun
def j1(A, x):
    """Computes the Bessel function `J_1(x)`.  See :func:`besselj`."""
    return A.besselj(1, x)

@defun_wrapped
def bessely(A,n,x):
    intdist = abs(n.imag) + abs(n.real-A.floor(n.real+0.5))
    if not intdist:
        h = +A.eps
        A.prec *= 2
        n += h
    else:
        A.prec += -int(A.log(intdist, 2)+1)
    return (A.besselj(n,x)*A.cospi(n) - A.besselj(-n,x))/A.sinpi(n)

@defun_wrapped
def besseli(A,n,x):
    if A.isint(n):
        n = abs(int(n))
    hx = x/2
    return hx**n * A.hyp0f1(n+1, hx**2) / A.factorial(n)

@defun_wrapped
def besselk(A,n,x):
    intdist = abs(n.imag) + abs(n.real-A.floor(n.real+0.5))
    if not intdist:
        h = +A.eps
        A.prec *= 2
        n += h
    else:
        A.prec += -int(A.log(intdist, 2)+1)
    return A.pi*(A.besseli(-n,x)-A.besseli(n,x))/(2*A.sinpi(n))

@defun_wrapped
def hankel1(A,n,x):
    return A.besselj(n,x) + A.j*A.bessely(n,x)

@defun_wrapped
def hankel2(A,n,x):
    return A.besselj(n,x) - A.j*A.bessely(n,x)

@defun_wrapped
def lambertw(A, z, k=0, approx=None):
    if A.isnan(z):
        return z
    A.prec += 20
    # We must be extremely careful near the singularities at -1/e and 0
    u = A.exp(-1)
    if abs(z) <= u:
        if not z:
            # w(0,0) = 0; for all other branches we hit the pole
            if not k:
                return z
            return A.ninf
        if not k:
            w = z
        # For small real z < 0, the -1 branch behaves roughly like log(-z)
        elif k == -1 and not A.im(z) and A.re(z) < 0:
            w = A.ln(-z)
        # Use a simple asymptotic approximation.
        else:
            w = A.ln(z)
            # The branches are roughly logarithmic. This approximation
            # gets better for large |k|; need to check that this always
            # works for k ~= -1, 0, 1.
            if k: w += k * 2*A.pi*A.j
    elif k == 0 and A.im(z) and abs(z) <= 0.6:
        w = z
    else:
        if z == A.inf:
            if k == 0:
                return z
            else:
                return z + 2*k*A.pi*A.j
        if z == A.ninf:
            return (-z) + (2*k+1)*A.pi*A.j
        # Simple asymptotic approximation as above
        w = A.ln(z)
        if k: w += k * 2*A.pi*A.j
    # Use Halley iteration to solve w*exp(w) = z
    two = A.mpf(2)
    weps = A.ldexp(A.eps, 15)
    for i in xrange(100):
        ew = A.exp(w)
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
def barnesg(A, z):
    if A.isinf(z):
        if z == A.inf:
            return z
        return A.nan
    if A.isnan(z):
        return z
    if (not z.imag) and z.real <= 0 and A.isint(z.real):
        return z*0
    # Account for size (would not be needed if computing log(G))
    if abs(z) > 5:
        A.dps += 2*A.log(abs(z),2)
    # Estimate terms for asymptotic expansion
    N = A.dps // 2 + 5
    G = 1
    while A.re(z) < N:
        G /= A.gamma(z)
        z += 1
    z -= 1
    s = A.mpf(1)/12
    s -= A.log(A.glaisher)
    s += z*A.log(2*A.pi)/2
    s += (z**2/2-A.mpf(1)/12)*A.log(z)
    s -= 3*z**2/4
    z2k = z2 = z**2
    for k in xrange(1, N+1):
        t = A.bernoulli(2*k+2) / (4*k*(k+1)*z2k)
        if abs(t) < A.eps:
            #print k, N      # check how many terms were needed
            break
        z2k *= z2
        s += t
    #if k == N:
    #    print "warning: series for barnesg failed to converge"
    return G*A.exp(s)

@defun
def superfac(A, z):
    return A.barnesg(z+2)

@defun_wrapped
def hyperfac(A, z):
    # XXX: estimate needed extra bits accurately
    if z == A.inf:
        return z
    if abs(z) > 5:
        extra = 4*int(A.log(abs(z),2))
    else:
        extra = 0
    A.prec += extra
    if not z.imag and z.real < 0 and A.isint(z.real):
        n = int(A.re(z))
        h = A.hyperfac(-n-1)
        if ((n+1)//2) & 1:
            h = -h
        if A.is_complex_type(z):
            return h + 0j
        return h
    zp1 = z+1
    # Wrong branch cut
    #v = A.gamma(zp1)**z
    #A.prec -= extra
    #return v / A.barnesg(zp1)
    v = A.exp(z*A.loggamma(zp1))
    A.prec -= extra
    return v / A.barnesg(zp1)

@defun_wrapped
def loggamma(A, z):
    a = z.real
    b = z.imag
    if not b and a > 0:
        return A.ln(A.gamma(z))
    u = A.arg(z)
    w = A.ln(A.gamma(z))
    if b:
        gi = -b - u/2 + a*u + b*A.ln(abs(z))
        n = A.floor((gi-w.imag)/(2*A.pi)+0.5) * (2*A.pi)
        return w + n*A.j
    elif a < 0:
        n = int(A.floor(a))
        w += (n-(n%2))*A.pi*A.j
    return w

@defun_wrapped
def siegeltheta(A, t):
    if t.imag:
        # XXX: cancellation occurs
        a = A.loggamma(0.25+0.5j*t)
        b = A.loggamma(0.25-0.5j*t)
        return -A.ln(A.pi)/2*t - 0.5j*(a-b)
    else:
        if A.isinf(t):
            return t
        return A.loggamma(0.25+0.5j*t).imag - A.ln(A.pi)/2*t

@defun_wrapped
def grampoint(A, n):
    # Asymptotic expansion, from
    # http://mathworld.wolfram.com/GramPoint.html
    g = 2*A.pi*A.exp(1+A.lambertw((8*n+1)/(8*e)))
    return A.findroot(lambda t: A.siegeltheta(t)-A.pi*n, g)

@defun_wrapped
def siegelz(A, t):
    v = A.exp(A.j*A.siegeltheta(t))*A.zeta(0.5+A.j*t)
    if A.is_real_type(t):
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
def zetazero(A, n, url='http://www.dtc.umn.edu/~odlyzko/zeta_tables/zeros1'):
    n = int(n)
    if n < 0:
        return zetazero(-n).conjugate()
    if n == 0:
        raise ValueError("n must be nonzero")
    if n > len(_zeta_zeros) and n <= 100000:
        _load_zeta_zeros(url)
    if n > len(_zeta_zeros):
        raise NotImplementedError("n too large for zetazeros")
    return A.mpc(0.5, A.findroot(A.siegelz, _zeta_zeros[n-1]))

@defun_wrapped
def riemannr(A, x):
    if x == 0:
        return A.zero
    # Check if a simple asymptotic estimate is accurate enough
    if abs(x) > 1000:
        a = A.li(x)
        b = 0.5*A.li(A.sqrt(x))
        if abs(b) < abs(a)*A.eps:
            return a
    if abs(x) < 0.01:
        # XXX
        A.prec += int(-A.log(abs(x),2))
    # Sum Gram's series
    s = t = A.one
    u = A.ln(x)
    k = 1
    while abs(t) > abs(s)*A.eps:
        t = t * u / k
        s += t / (k * A.zeta(k+1))
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
def primepi2(A, x):
    x = int(x)
    if x < 2:
        return A.mpi(0,0)
    if x < 2657:
        return A.mpi(A.primepi(x))
    mid = A.li(x)
    # Schoenfeld's estimate for x >= 2657, assuming RH
    err = A.sqrt(x,rounding='u')*A.ln(x,rounding='u')/8/A.pi(rounding='d')
    a = A.floor((A.mpi(mid)-err).a, rounding='d')
    b = A.ceil((A.mpi(mid)+err).b, rounding='u')
    return A.mpi(a, b)

@defun_wrapped
def primezeta(A, s):
    if A.isnan(s):
        return s
    if A.re(s) <= 0:
        raise ValueError("prime zeta function defined only for re(s) > 0")
    if s == 1:
        return A.inf
    if s == 0.5:
        return A.mpc(A.ninf, A.pi)
    r = A.re(s)
    if r > A.prec:
        return 0.5**s
    else:
        wp = A.prec + int(r)
        def terms():
            orig = A.prec
            # zeta ~ 1+eps; need to set precision
            # to get logarithm accurately
            k = 0
            while 1:
                k += 1
                u = libintmath.moebius(k)
                if not u:
                    continue
                A.prec = wp
                t = u*A.ln(A.zeta(k*s))/k
                if not t:
                    return
                #print A.prec, A.nstr(t)
                A.prec = orig
                yield t
    return sum_accurately(A, terms)

@defun_wrapped
def bernpoly(A, n, z):
    n = int(n)
    assert n >= 0
    # XXX: optimize
    return sum(A.binomial(n,k)*A.bernoulli(k)*z**(n-k) for k in xrange(0,n+1))

# TODO: this should be implemented low-level
def polylog_series(A, s, z):
    tol = +A.eps
    l = A.zero
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

def polylog_continuation(A, n, z):
    if n < 0:
        return z*0
    twopij = 2j * A.pi
    a = -twopij**n/A.fac(n) * A.bernpoly(n, A.ln(z)/twopij)
    if A.is_real_type(z) and z < 0:
        a = a.real
    if z.imag < 0 or (z.imag == 0 and z.real >= 1):
        a -= twopij*A.ln(z)**(n-1)/A.fac(n-1)
    return a

def polylog_unitcircle(A, n, z):
    tol = +A.eps
    if n > 1:
        l = A.zero
        logz = A.ln(z)
        logmz = A.one
        m = 0
        while 1:
            if (n-m) != 1:
                term = A.zeta(n-m) * logmz / A.fac(m)
                if term and abs(term) < tol:
                    break
                l += term
            logmz *= logz
            m += 1
        l += A.ln(z)**(n-1)/A.fac(n-1)*(A.harmonic(n-1)-A.ln(-A.ln(z)))
    elif n < 1:  # else
        l = A.fac(-n)*(-A.ln(z))**(n-1)
        logz = A.ln(z)
        logkz = A.one
        k = 0
        while 1:
            b = A.bernoulli(k-n+1)
            if b:
                term = b*logkz/(A.fac(k)*(k-n+1))
                if abs(term) < tol:
                    break
                l -= term
            logkz *= logz
            k += 1
    else:
        raise ValueError
    if A.is_real_type(z) and z < 0:
        l = l.real
    return l

@defun_wrapped
def polylog(A, s, z):
    if z == 1:
        return A.zeta(s)
    if z == -1:
        return -A.altzeta(s)
    if s == 0:
        return z/(1-z)
    if s == 1:
        return -A.ln(1-z)
    if s == -1:
        return z/(1-z)**2
    if abs(z) <= 0.75 or (not A.isint(s) and abs(z) < 0.99):
        return polylog_series(A, s, z)
    if abs(z) >= 1.4 and A.isint(s):
        return (-1)**(s+1)*polylog_series(A, s, 1/z) + polylog_continuation(A, s, z)
    if A.isint(s):
        return polylog_unitcircle(A, int(s), z)
    raise NotImplementedError("polylog for arbitrary s and z")
    # This could perhaps be used in some cases
    #from quadrature import quad
    #return quad(lambda t: t**(s-1)/(exp(t)/z-1),[0,inf])/gamma(s)

# Experimental code; could be used elsewhere
def sum_accurately(A, terms, check_step=1):
    orig = A.prec
    extra = 10
    while 1:
        A.prec = orig + extra
        max_term = A.ninf
        s = 0
        k = 0
        for term in terms():
            s += term
            if not k % check_step and term:
                abs_term = abs(term)
                abs_sum = abs(s)
                max_term = max(max_term, abs_term)
                if abs_term <= A.eps*abs_sum:
                    break
            k += 1
        if abs_sum:
            cancellation = int(max(0,A.log(max_term/abs_sum,2)))
        else:
            cancellation = A.prec
        if cancellation < extra:
            break
        else:
            extra += cancellation
    return s

@defun_wrapped
def bell(A, n, x=1):
    x = A.convert(x)
    if not n:
        if A.isnan(x):
            return x
        return type(x)(1)
    if A.isinf(x) or A.isinf(n) or A.isnan(x) or A.isnan(n):
        return x**n
    if n == 1: return x
    if n == 2: return x*(x+1)
    if x == 0: return A.sincpi(n)
    return _polyexp(A, n, x, True) / A.exp(x)

def _polyexp(A, n, x, extra=False):
    def _terms():
        if extra:
            yield A.sincpi(n)
        t = x
        k = 1
        while 1:
            yield k**n * t
            k += 1
            t = t*x/k
    return sum_accurately(A, _terms, check_step=4)

@defun_wrapped
def polyexp(A, s, z):
    if A.isinf(z) or A.isinf(s) or A.isnan(z) or A.isnan(s):
        return z**s
    if z == 0: return z*s
    if s == 0: return A.expm1(z)
    if s == 1: return A.exp(z)*z
    if s == 2: return A.exp(z)*z*(z+1)
    return _polyexp(A, s, z)

# TODO: tests; improve implementation
@defun_wrapped
def expm1(A, x):
    """
    Accurately computes exp(x)-1.
    """
    if not x:
        return type(x)(1)
    return sum_accurately(A, lambda: iter([A.exp(x),-1]),1)

# hack
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    import doctest
    globs = globals().copy()
    for obj in globs: #sorted(globs.keys()):
        print obj
        doctest.run_docstring_examples(globs[obj], {})

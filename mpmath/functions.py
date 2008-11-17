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

from mptypes import (\
    mpnumeric, mpmathify,
    mpf, make_mpf,
    mpc, make_mpc,
    mpi, make_mpi,
    constant,
    prec_rounding, mp,
    extraprec,
    zero, one, inf, ninf, nan, j, isnan, isinf, isint, eps,
    ComplexResult,
)

class _pi(constant):
    r"""
    `\pi`, roughly equal to 3.141592654, represents the area of the unit
    circle, the half-period of trigonometric functions, and many other
    things in mathematics.

    Mpmath can evaluate `\pi` to arbitrary precision:

        >>> from mpmath import *
        >>> mp.dps = 50
        >>> print pi
        3.1415926535897932384626433832795028841971693993751

    This shows digits 99991-100000 of `\pi`:

        >>> mp.dps = 100000
        >>> str(pi)[-10:]
        '5549362464'

    **Possible issues**

    :data:`pi` always rounds to the nearest floating-point
    number when used. This means that exact mathematical identities
    involving `\pi` will generally not be preserved in floating-point
    arithmetic. In particular, multiples of :data:`pi` (except for 
    the trivial case ``0*pi``) are *not* the exact roots of
    :func:`sin`, but differ roughly by the current epsilon:

        >>> mp.dps = 15
        >>> sin(pi)
        mpf('1.2246467991473532e-16')

    One solution is to use the :func:`sinpi` function instead:

        >>> sinpi(1)
        mpf('0.0')

    See the documentation of trigonometric functions for additional
    details.
    """


class _degree(constant): pass

class _e(constant):
    """
    The transcendental number `e` = 2.718281828... is the base of the
    natural logarithm (:func:`ln`) and of the exponential function
    (:func:`exp`).

    Mpmath can be evaluate `e` to arbitrary precision:

        >>> mp.dps = 50
        >>> print e
        2.7182818284590452353602874713526624977572470937

    This shows digits 99991-100000 of `e`:

        >>> mp.dps = 100000
        >>> str(e)[-10:]
        '2100427165'

    **Possible issues**

    :data:`e` always rounds to the nearest floating-point number
    when used, and mathematical identities involving `e` may not
    hold in floating-point arithmetic. For example, ``ln(e)``
    might not evaluate exactly to 1.

    In particular, don't use ``e**x`` to compute the exponential
    function. Use ``exp(x)`` instead; this is both faster and more
    accurate.
    """
    pass


class _ln2(constant): pass

class _ln10(constant): pass

class _phi(constant): pass

class _euler(constant): pass

class _catalan(constant): pass

class _khinchin(constant): pass

class _glaisher(constant):
    """
    Glaisher's constant `A`, also known as the Glaisher-Kinkelin
    constant, is a number approximately equal to 1.282427129 that
    sometimes appears in formulas related to gamma and zeta functions.

    It is defined as `A = \exp(1/12-\zeta'(-1))`. Here `\zeta'(s)`
    denotes the derivative of the Riemann zeta function (see
    :func:`zeta`).

    Mpmath can evaluate Glaisher's constant to arbitrary precision:

        >>> from mpmath import *
        >>> mp.dps = 50
        >>> print glaisher
        1.282427129100622636875342568869791727767688927325

    We can verify that the value computed by :data:`glaisher` is
    correct using mpmath's facilities for numerical
    differentiation and arbitrary evaluation of the zeta function:

        >>> print exp(mpf(1)/12 - diff(zeta, -1))
        1.282427129100622636875342568869791727767688927325

    Here is an example of an integral that can be evaluated in
    terms of Glaisher's constant:

        >>> mp.dps = 15
        >>> print quad(lambda x: log(gamma(x)), [1, 1.5])
        -0.0428537406502909
        >>> print -0.5 - 7*log(2)/24 + log(pi)/4 + 3*log(glaisher)/2
        -0.042853740650291

    Mpmath computes Glaisher's constant by applying Euler-Maclaurin
    summation to a slowly convergent series. The implementation is
    reasonably efficient up to about 10,000 digits. See the source
    code for additional details.

    References:
    http://mathworld.wolfram.com/Glaisher-KinkelinConstant.html

    """

class _apery(constant): pass

# Mathematical constants
pi = _pi(libelefun.mpf_pi, "pi")
degree = _degree(libelefun.mpf_degree, "degree")
e = _e(libelefun.mpf_e, "e")
ln2 = _ln2(libelefun.mpf_ln2, "ln(2)")
ln10 = _ln10(libelefun.mpf_ln10, "ln(10)")
phi = _phi(libelefun.mpf_phi, "Golden ratio (phi)")
euler = _euler(gammazeta.mpf_euler, "Euler's constant (gamma)")
catalan = _catalan(gammazeta.mpf_catalan, "Catalan's constant")
khinchin = _khinchin(gammazeta.mpf_khinchin, "Khinchin's constant")
glaisher = _glaisher(gammazeta.mpf_glaisher, "Glaisher's constant")
apery = _apery(gammazeta.mpf_apery, "Apery's constant")


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

def mpfunc(name, real_f, complex_f, doc, interval_f=None):
    def f(x, **kwargs):
        if not isinstance(x, mpnumeric):
            x = mpmathify(x)
        prec, rounding = prec_rounding
        if kwargs:
            prec = kwargs.get('prec', prec)
            if 'dps' in kwargs:
                prec = dps_to_prec(kwargs['dps'])
            rounding = kwargs.get('rounding', rounding)
        if isinstance(x, mpf):
            try:
                return make_mpf(real_f(x._mpf_, prec, rounding))
            except ComplexResult:
                # Handle propagation to complex
                if mp.trap_complex:
                    raise
                return make_mpc(complex_f((x._mpf_, libmpf.fzero), prec, rounding))
        elif isinstance(x, mpc):
            return make_mpc(complex_f(x._mpc_, prec, rounding))
        elif isinstance(x, mpi):
            if interval_f:
                return make_mpi(interval_f(x._val, prec))
        raise NotImplementedError("%s of a %s" % (name, type(x)))

    f.__name__ = name
    f.__doc__ = "Returns the %s of x" % doc
    return f

def altfunc(f, name, desc):
    def g(x):
        orig = mp.prec
        try:
            mp.prec = orig + 10
            return one/f(x)
        finally:
            mp.prec = orig
    g.__name__ = name
    g.__doc__ = "Returns the %s of x, 1/%s(x)" % (desc, f.__name__)
    return g

def altinvfunc(f, name, desc):
    def g(x):
        orig = mp.prec
        try:
            mp.prec = orig + 10
            return f(one/x)
        finally:
            mp.prec = orig
    g.__name__ = name
    g.__doc__ = "Returns the inverse %s of x, %s(1/x)" % (desc, f.__name__)
    return g

sqrt = mpfunc('sqrt', libelefun.mpf_sqrt, libmpc.mpc_sqrt, "principal square root", libmpi.mpi_sqrt)
sqrt.__doc__ = r"""
``sqrt(x)`` computes the principal square root of `x`, `\sqrt x`.
For positive real numbers, the principal root is simply the
positive square root. For arbitrary complex numbers, the principal
square root is defined to satisfy `\sqrt x = \exp(\log(x)/2)`.
The function thus has a branch cut along the negative half real axis.

For all mpmath numbers ``x``, calling ``sqrt(x)`` is equivalent to
performing ``x**0.5``.

**Examples**

Basic examples and limits::

    >>> from mpmath import *
    >>> mp.dps = 15
    >>> print sqrt(10)
    3.16227766016838
    >>> print sqrt(100)
    10.0
    >>> print sqrt(-4)
    (0.0 + 2.0j)
    >>> print sqrt(1+1j)
    (1.09868411346781 + 0.455089860562227j)
    >>> print sqrt(inf)
    +inf

Square root evaluation is fast at huge precision::

    >>> mp.dps = 50000
    >>> a = sqrt(3)
    >>> str(a)[-10:]
    '9329332814'

:func:`sqrt` supports interval arguments::

    >>> mp.dps = 15
    >>> print sqrt(mpi(16, 100))
    [4.0, 10.0]
    >>> print sqrt(mpi(2))
    [1.4142135623730949234, 1.4142135623730951455]
    >>> print sqrt(mpi(2)) ** 2
    [1.9999999999999995559, 2.0000000000000004441]

"""

cbrt = mpfunc('cbrt', libelefun.mpf_cbrt, libmpc.mpc_cbrt, "principal cubic root")
cbrt.__doc__ = """
``cbrt(x)`` computes the cube root of `x`, `x^{1/3}`. This
function is faster and more accurate than raising to a floating-point
fraction::

    >>> from mpmath import *
    >>> mp.dps = 15
    >>> 125**(mpf(1)/3)
    mpf('4.9999999999999991')
    >>> cbrt(125)
    mpf('5.0')

Every nonzero complex number has three cube roots. This function
returns the cube root defined by `\exp(\log(x)/3)` where the
principal branch of the natural logarithm is used. Note that this
does not give a real cube root for negative real numbers::

    >>> print cbrt(-1)
    (0.5 + 0.866025403784439j)

"""

exp = mpfunc('exp', libelefun.mpf_exp, libmpc.mpc_exp, "exponential function", libmpi.mpi_exp)
ln = mpfunc('ln', libelefun.mpf_log, libmpc.mpc_log, "natural logarithm", libmpi.mpi_log)

ln.__doc__ = r"""Computes the natural logarithm of `x`, `\ln x`.
See :func:`log` for additional documentation."""

cos = mpfunc('cos', libelefun.mpf_cos, libmpc.mpc_cos, "cosine", libmpi.mpi_cos)
sin = mpfunc('sin', libelefun.mpf_sin, libmpc.mpc_sin, "sine", libmpi.mpi_sin)
tan = mpfunc('tan', libelefun.mpf_tan, libmpc.mpc_tan, "tangent", libmpi.mpi_tan)
cosh = mpfunc('cosh', libelefun.mpf_cosh, libmpc.mpc_cosh, "hyperbolic cosine")
sinh = mpfunc('sinh', libelefun.mpf_sinh, libmpc.mpc_sinh, "hyperbolic sine")
tanh = mpfunc('tanh', libelefun.mpf_tanh, libmpc.mpc_tanh, "hyperbolic tangent")

acos = mpfunc('acos', libelefun.mpf_acos, libmpc.mpc_acos, "inverse cosine")
asin = mpfunc('asin', libelefun.mpf_asin, libmpc.mpc_asin, "inverse sine")
atan = mpfunc('atan', libelefun.mpf_atan, libmpc.mpc_atan, "inverse tangent")
asinh = mpfunc('asinh', libelefun.mpf_asinh, libmpc.mpc_asinh, "inverse hyperbolic sine")
acosh = mpfunc('acosh', libelefun.mpf_acosh, libmpc.mpc_acosh, "inverse hyperbolic cosine")
atanh = mpfunc('atanh', libelefun.mpf_atanh, libmpc.mpc_atanh, "inverse hyperbolic tangent")

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

cospi = mpfunc('cospi', libelefun.mpf_cos_pi, libmpc.mpc_cos_pi, "")
sinpi = mpfunc('sinpi', libelefun.mpf_sin_pi, libmpc.mpc_sin_pi, "")

sinpi.__doc__ = """
Computes `\sin(\pi x)`, more accurately than the expression
``sin(pi*x)``::

    >>> from mpmath import *
    >>> mp.dps = 15
    >>> print sinpi(10**10), sin(pi*(10**10))
    0.0 -2.23936276195592e-6
    >>> print sinpi(10**10+0.5), sin(pi*(10**10+0.5))
    1.0 0.999999999998721
"""

cospi.__doc__ = """
Computes `\cos(\pi x)`, more accurately than the expression
``cos(pi*x)``::

    >>> from mpmath import *
    >>> mp.dps = 15
    >>> print cospi(10**10), cos(pi*(10**10))
    1.0 0.999999999997493
    >>> print cospi(10**10+0.5), cos(pi*(10**10+0.5))
    0.0 1.59960492420134e-6
"""

@funcwrapper
def sinc(x):
    r"""
    ``sinc(x)`` computes the unnormalized sinc function, defined as

    .. math ::

        \mathrm{sinc}(x) = \begin{cases}
            \sin(x)/x, & \mbox{if } x \ne 0 \\
            1,         & \mbox{if } x = 0.
        \end{cases}

    See :func:`sincpi` for the normalized sinc function.

    Simple values and limits include::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print sinc(0)
        1.0
        >>> print sinc(1)
        0.841470984807897
        >>> print sinc(inf)
        0.0

    The integral of the sinc function is the sine integral Si::

        >>> print quad(sinc, [0, 1])
        0.946083070367183
        >>> print si(1)
        0.946083070367183
    """
    if isinf(x):
        return 1/x
    if not x:
        return x+1
    return sin(x)/x

@funcwrapper
def sincpi(x):
    r"""
    ``sincpi(x)`` computes the normalized sinc function, defined as

    .. math ::

        \mathrm{sinc}_{\pi}(x) = \begin{cases}
            \sin(\pi x)/(\pi x), & \mbox{if } x \ne 0 \\
            1,                   & \mbox{if } x = 0.
        \end{cases}

    Equivalently, we have
    `\mathrm{sinc}_{\pi}(x) = \mathrm{sinc}(\pi x)`.

    The normalization entails that the function integrates
    to unity over the entire real line::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print quadosc(sincpi, [-inf, inf], period=2.0)
        1.0

    Like, :func:`sinpi`, :func:`sincpi` is evaluated accurately
    at its roots::

        >>> print sincpi(10)
        0.0
    """
    if isinf(x):
        return 1/x
    if not x:
        return x+1
    return sinpi(x)/(pi*x)

floor = mpfunc('floor', libmpf.mpf_floor, libmpc.mpc_floor, "")
floor.__doc__ = r"""
Computes the floor of `x`, `\lfloor x \rfloor`, defined as
the largest integer less than or equal to `x`::

    >>> from mpmath import *
    >>> print floor(3.5)
    3.0

Note: :func:`floor` returns a floating-point number, not a
Python ``int``. If `\lfloor x \rfloor` is too large to be
represented exactly at the present working precision, the
result will be rounded, not necessarily in the floor
direction."""

ceil = mpfunc('ceil', libmpf.mpf_ceil, libmpc.mpc_ceil, "")
ceil.__doc__ = r"""
Computes the ceiling of `x`, `\lceil x \rceil`, defined as
the smallest integer greater than or equal to `x`::

    >>> from mpmath import *
    >>> print ceil(3.5)
    4.0

Note: :func:`ceil` returns a floating-point number, not a
Python ``int``. If `\lceil x \rceil` is too large to be
represented exactly at the present working precision, the
result will be rounded, not necessarily in the ceiling
direction."""

@funcwrapper
def nthroot(x, n):
    r"""
    ``nthroot(x, n)`` computes the principal `n`-th root of `x`,
    `x^{1/n}`. Here `n` must be an integer, and can be negative
    (`x^{-1/n}` is `1/x^{1/n}`).

    For `n = 2` or `n = 3`, using this function is equivalent to
    calling :func:`sqrt` or :func:`cbrt`. In general,
    ``nthroot(x, n)`` is defined to compute `\exp(\log(x)/n)`.

    :func:`nthroot` is implemented to use Newton's method for small
    `n`. At high precision, this makes `x^{1/n}` not much more
    expensive than the regular exponentiation, `x^n`. For very large
    `n`, :func:`nthroot` falls back to use the exponential function.

    :func:`nthroot` is faster and more accurate than raising to a
    floating-point fraction::
    
        >>> from mpmath import *
        >>> mp.dps = 15
        >>> 16807 ** (mpf(1)/5)
        mpf('7.0000000000000009')
        >>> nthroot(16807, 5)
        mpf('7.0')

    """
    n = int(n)
    if isinstance(x, mpf):
        try:
            return make_mpf(libelefun.mpf_nthroot(x._mpf_, n, *prec_rounding))
        except ComplexResult:
            if mp.trap_complex:
                raise
            x = (x._mpf_, libmpf.fzero)
    else:
        x = x._mpc_
    return make_mpc(libmpc.mpc_nthroot(x, n, *prec_rounding))

def hypot(x, y):
    r"""
    Computes the Euclidean norm of the vector `(x, y)`, equal
    to `\sqrt{x^2 + y^2}`. Both `x` and `y` must be real."""
    x = mpmathify(x)
    y = mpmathify(y)
    return make_mpf(libmpf.mpf_hypot(x._mpf_, y._mpf_, *prec_rounding))

def ldexp(x, n):
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
    x = mpmathify(x)
    return make_mpf(libmpf.mpf_shift(x._mpf_, n))

def frexp(x):
    r"""
    Given a real number `x`, returns `(y, n)` with `y \in [0.5, 1)`,
    `n` a Python integer, and such that `x = y 2^n`. No rounding is
    performed.

        >>> from mpmath import *
        >>> frexp(7.5)
        (mpf('0.9375'), 3)

    """
    x = mpmathify(x)
    y, n = libmpf.mpf_frexp(x._mpf_)
    return make_mpf(y), n

def sign(x):
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
    x = mpmathify(x)
    if not x or isnan(x):
        return x
    if isinstance(x, mpf):
        return mpf(cmp(x, 0))
    return x / abs(x)

@extraprec(5)
def arg(x):
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
    x = mpc(x)
    return atan2(x.imag, x.real)

def fabs(x):
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
    return abs(mpmathify(x))

def re(x):
    r"""
    Returns the real part of `x`, `\Re(x)`. Unlike ``x.real``,
    :func:`re` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> re(3)
        mpf('3.0')
        >>> re(-1+4j)
        mpf('-1.0')

    """
    return mpmathify(x).real

def im(x):
    r"""
    Returns the imaginary part of `x`, `\Im(x)`. Unlike ``x.imag``,
    :func:`im` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> im(3)
        mpf('0.0')
        >>> im(-1+4j)
        mpf('4.0')

    """
    return mpmathify(x).imag

def conj(x):
    r"""
    Returns the complex conjugate of `x`, `\overline{x}`. Unlike
    ``x.conjugate()``, :func:`im` converts `x` to a mpmath number::

        >>> from mpmath import *
        >>> conj(3)
        mpf('3.0')
        >>> conj(-1+4j)
        mpc(real='-1.0', imag='-4.0')

    """
    return mpmathify(x).conjugate()

def log(x, b=None):
    r"""
    Computes the base-`b` logarithm of `x`, `\log_b(x)`. If `b` is
    unspecified, :func:`log` computes the natural (base `e`) logarithm
    and is equivalent to :func:`ln`. In general, the base `b` logarithm
    is defined in terms of the natural logarithm as
    `\log_b(x) = \ln(x)/\ln(b)`.

    By convention, we take `\log(0) = -\infty`.

    The natural logarithm is real if `x > 0` and complex if `x < 0` or if
    `x` is complex. The principal branch of the complex logarithm is
    used, meaning that `\Im(\ln(x)) = -\pi < \arg(x) \le \pi`.

    **Examples**

    Some basic values and limits::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print log(1)
        0.0
        >>> print log(2)
        0.693147180559945
        >>> print log(1000,10)
        3.0
        >>> print log(4, 16)
        0.5
        >>> print log(j)
        (0.0 + 1.5707963267949j)
        >>> print log(-1)
        (0.0 + 3.14159265358979j)
        >>> print log(0)
        -inf
        >>> print log(inf)
        +inf

    The natural logarithm is the antiderivative of `1/x`::

        >>> print quad(lambda x: 1/x, [1, 5])
        1.6094379124341
        >>> print log(5)
        1.6094379124341
        >>> print diff(log, 10)
        0.1

    The Taylor series expansion of the natural logarithm around
    `x = 1` has coefficients `(-1)^{n+1}/n`::

        >>> nprint(taylor(log, 1, 7))
        [0.0, 1.0, -0.5, 0.333333, -0.25, 0.2, -0.166667, 0.142857]

    :func:`log` supports arbitrary precision evaluation::

        >>> mp.dps = 50
        >>> print log(pi)
        1.1447298858494001741434273513530587116472948129153
        >>> print log(pi, pi**3)
        0.33333333333333333333333333333333333333333333333333
        >>> mp.dps = 25
        >>> print log(3+4j)
        (1.609437912434100374600759 + 0.9272952180016122324285125j)

    """
    if b is None:
        return ln(x)
    wp = mp.prec + 20
    return ln(x, prec=wp) / ln(b, prec=wp)

def log10(x):
    r"""
    Computes the base-10 logarithm of `x`, `\log_{10}(x)`. ``log10(x)``
    is equivalent to ``log(x, 10)``.
    """
    return log(x, 10)

def power(x, y):
    r"""
    Converts `x` and `y` to mpmath numbers and evaluates
    `x^y = \exp(y \log(x))`::

        >>> from mpmath import *
        >>> mp.dps = 30
        >>> print power(2, 0.5)
        1.41421356237309504880168872421

    This shows the leading few digits of a large Mersenne prime
    (performing the exact calculation ``2**43112609-1`` and
    displaying the result in Python would be very slow)::

        >>> print power(2, 43112609)-1
        3.16470269330255923143453723949e+12978188

    """
    return mpmathify(x) ** mpmathify(y)

def modf(x,y):
    r"""
    Converts `x` and `y` to mpmath numbers and returns `x \mod y`.
    For mpmath numbers, this is equivalent to ``x % y``.

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print modf(100, pi)
        2.61062773871641

    You can use :func:`modf` to compute fractional parts of numbers::

        >>> print modf(10.25, 1)
        0.25

    """
    x = mpmathify(x)
    y = mpmathify(y)
    return x % y

def degrees(x):
    r"""
    Converts the radian angle `x` to a degree angle::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print degrees(pi/3)
        60.0
    """
    return x / degree

def radians(x):
    r"""
    Converts the degree angle `x` to radians::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print radians(60)
        1.0471975511966
    """
    return x * degree

def atan2(y, x):
    r"""
    Computes the two-argument arctangent, `\mathrm{atan2}(y, x)`,
    giving the signed angle between the positive `x`-axis and the
    point `(x, y)` in the 2D plane. This function is defined for
    real `x` and `y` only.

    The two-argument arctangent essentially computes
    `\mathrm{atan}(y/x)`, but accounts for the signs of both
    `x` and `y` to give the angle for the correct quadrant. The
    following examples illustrate the difference::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print atan2(1,1), atan(1/1.)
        0.785398163397448 0.785398163397448
        >>> print atan2(1,-1), atan(1/-1.)
        2.35619449019234 -0.785398163397448
        >>> print atan2(-1,1), atan(-1/1.)
        -0.785398163397448 -0.785398163397448
        >>> print atan2(-1,-1), atan(-1/-1.)
        -2.35619449019234 0.785398163397448

    The angle convention is the same as that used for the complex
    argument; see :func:`arg`.
    """
    x = mpmathify(x)
    y = mpmathify(y)
    return make_mpf(libelefun.mpf_atan2(y._mpf_, x._mpf_, *prec_rounding))


fib = fibonacci = mpfunc('fibonacci', libelefun.mpf_fibonacci, libmpc.mpc_fibonacci, "")

fibonacci.__doc__ = r"""
``fibonacci(n)`` computes the `n`-th Fibonacci number, `F(n)`. The
Fibonacci numbers are defined by the recurrence `F(n) = F(n-1) + F(n-2)`
with the initial values `F(0) = 0`, `F(1) = 1`. :func:`fibonacci`
extends this definition to arbitrary real and complex arguments
using the formula

.. math ::

  F(z) = \frac{\phi^z - \cos(\pi z) \phi^{-z}}{\sqrt 5}

where `\phi` is the golden ratio. :func:`fibonacci` also uses this
continuous formula to compute `F(n)` for extremely large `n`, where
calculating the exact integer would be wasteful.

For convenience, :func:`fib` is available as an alias for
:func:`fibonacci`.

**Basic examples**

Some small Fibonacci numbers are::

    >>> from mpmath import *
    >>> mp.dps = 15
    >>> for i in range(10):
    ...     print fibonacci(i),
    ...
    0.0 1.0 1.0 2.0 3.0 5.0 8.0 13.0 21.0 34.0

    >>> print fibonacci(50)
    12586269025.0

The recurrence for `F(n)` extends backwards to negative `n`::

    >>> for i in range(10):
    ...     print fibonacci(-i),
    ...
    0.0 1.0 -1.0 2.0 -3.0 5.0 -8.0 13.0 -21.0 34.0

Large Fibonacci numbers will be computed approximately unless
the precision is set high enough::

    >>> print fib(200)
    2.8057117299251e+41
    >>> mp.dps = 45
    >>> print fib(200)
    280571172992510140037611932413038677189525.0

:func:`fibonacci` can compute approximate Fibonacci numbers
of stupendous size::

    >>> mp.dps = 15
    >>> print fibonacci(10**25)
    3.49052338550226e+2089876402499787337692720

**Real and complex arguments**

The extended Fibonacci function is an analytic function. The
property `F(z) = F(z-1) + F(z-2)` holds for arbitrary `z`::

    >>> mp.dps = 15
    >>> print fib(pi)
    2.1170270579161
    >>> print fib(pi-1) + fib(pi-2)
    2.1170270579161
    >>> print fib(3+4j)
    (-5248.51130728372 - 14195.962288353j)
    >>> print fib(2+4j) + fib(1+4j)
    (-5248.51130728372 - 14195.962288353j)

The Fibonacci function has infinitely many roots on the
negative half-real axis. The first root is at 0, the second is
close to -0.18, and then there are infinitely many roots that
asymptotically approach `-n+1/2`::

    >>> print findroot(fib, -0.2)
    -0.183802359692956
    >>> print findroot(fib, -2)
    -1.57077646820395
    >>> print findroot(fib, -17)
    -16.4999999596115
    >>> print findroot(fib, -24)
    -23.5000000000479

**Mathematical relationships**

For large `n`, `F(n+1)/F(n)` approaches the golden ratio::

    >>> mp.dps = 50
    >>> print fibonacci(101)/fibonacci(100)
    1.6180339887498948482045868343656381177203127439638
    >>> print phi
    1.6180339887498948482045868343656381177203091798058

The sum of reciprocal Fibonacci numbers converges to an irrational
number for which no closed form expression is known::

    >>> mp.dps = 15
    >>> print nsum(lambda n: 1/fib(n), [1, inf])
    3.35988566624318

Amazingly, however, the sum of odd-index reciprocal Fibonacci
numbers can be expressed in terms of a Jacobi theta function::

    >>> print nsum(lambda n: 1/fib(2*n+1), [0, inf])
    1.82451515740692
    >>> print sqrt(5)*jtheta(2,0,(3-sqrt(5))/2)**2/4
    1.82451515740692

Some related sums can be done in closed form::

    >>> print nsum(lambda k: 1/(1+fib(2*k+1)), [0, inf])
    1.11803398874989
    >>> print phi - 0.5
    1.11803398874989
    >>> f = lambda k:(-1)**(k+1) / sum(fib(n)**2 for n in range(1,k+1))
    >>> print nsum(f, [1, inf])
    0.618033988749895
    >>> print phi-1
    0.618033988749895

**References**

1. http://mathworld.wolfram.com/FibonacciNumber.html
"""


zeta = mpfunc('zeta', gammazeta.mpf_zeta, gammazeta.mpc_zeta, 'Riemann zeta function')
altzeta = mpfunc('zeta', gammazeta.mpf_altzeta, gammazeta.mpc_altzeta, 'Dirichlet eta function')

zeta.__doc__ = r"""
    ``zeta(s)`` computes the Riemann zeta function, `\zeta(s)`.
    The Riemann zeta function is defined for `\Re(s) > 1` by

    .. math ::

      \zeta(s) = 1+\frac{1}{2^s}+\frac{1}{3^s}+\frac{1}{4^s}+\ldots

    and for `\Re(s) \le 1` by analytic continuation. It has a pole
    at `s = 1`.

    **Examples**

    Some exact values of the zeta function are::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print zeta(2)
        1.64493406684823
        >>> print pi**2 / 6
        1.64493406684823
        >>> print zeta(0)
        -0.5
        >>> print zeta(-1)
        -0.0833333333333333
        >>> print zeta(-2)
        0.0

    :func:`zeta` supports arbitrary precision evaluation and
    complex arguments::

        >>> mp.dps = 50
        >>> print zeta(pi)
        1.1762417383825827588721504519380520911697389900217
        >>> print zeta(1+2j)  # doctest: +NORMALIZE_WHITESPACE
        (0.5981655697623817367034568491742186771747764868876 -
        0.35185474521784529049653859679690026505229177886045j)

    The Riemann zeta function has so-called nontrivial zeros on
    the critical line `s = 1/2 + it`::

        >>> mp.dps = 15
        >>> print findroot(zeta, 0.5+14j)
        (0.5 + 14.1347251417347j)
        >>> print findroot(zeta, 0.5+21j)
        (0.5 + 21.0220396387716j)
        >>> print findroot(zeta, 0.5+25j)
        (0.5 + 25.0108575801457j)

    For large positive `s`, `\zeta(s)` rapidly approaches 1::

        >>> print zeta(30)
        1.00000000093133
        >>> print zeta(100)
        1.0
        >>> print zeta(inf)
        1.0

    The following series converges and in fact has a simple
    closed form value::

        >>> print nsum(lambda k: zeta(k)-1, [2, inf])
        1.0

    **Algorithm**

    The primary algorithm is Borwein's algorithm for the Dirichlet
    eta function. Three separate implementations are used: for general
    real arguments, general complex arguments, and for integers. The
    reflection formula is applied to arguments in the negative
    half-plane. For very large real arguments, either direct
    summation or the Euler prime product is used.

    It should be noted that computation of `\zeta(s)` gets very slow
    when `s` is far away from the real axis.

    **References**

    1. http://mathworld.wolfram.com/RiemannZetaFunction.html

    2. http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P155.pdf
"""

altzeta.__doc__ = r"""
    Computes the Dirichlet eta function, `\eta(s)`, also known as the
    alternating zeta function. This function is defined in analogy
    with the Riemann zeta function as providing the sum of the
    alternating series

    .. math ::

        \eta(s) = 1-\frac{1}{2^s}+\frac{1}{3^s}-\frac{1}{4^s}+\ldots

    Note that `\eta(1) = \log(2)` is the alternating harmonic series.
    The eta function unlike the Riemann zeta function is an entire
    function, having a finite value for all complex `s`.

    The alternating and non-alternating zeta functions are related
    via the simple formula

    .. math ::

        \eta(s) = (1 - 2^{1-s}) \zeta(s).

    This formula can be used to define `\eta(s)` for `\Re(s) \le 0`,
    where the series diverges.

    **Examples**

    Some special values are::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print altzeta(1)
        0.693147180559945
        >>> print altzeta(0)
        0.5
        >>> print altzeta(-1)
        0.25
        >>> print altzeta(-2)
        0.0

    An example of a sum that can be computed more accurately and
    efficiently via :func:`altzeta` than via numerical summation::

        >>> sum(-(-1)**n / n**2.5 for n in range(1, 100))
        0.86720495150398402
        >>> print altzeta(2.5)
        0.867199889012184

    At positive even integers, the Dirichlet eta function
    evaluates to a rational multiple of a power of `\pi`::

        >>> print altzeta(2)
        0.822467033424113
        >>> print pi**2/12
        0.822467033424113

    Like the Riemann zeta function, `\eta(s)`, approaches 1
    as `s` approaches positive infinity, although it does
    so from below rather than from above::

        >>> print altzeta(30)
        0.999999999068682
        >>> print altzeta(inf)
        1.0
        >>> altzeta(1000, rounding='d')
        mpf('0.99999999999999989')
        >>> altzeta(1000, rounding='u')
        mpf('1.0')

    **References**

    1. http://mathworld.wolfram.com/DirichletEtaFunction.html

    2. http://en.wikipedia.org/wiki/Dirichlet_eta_function

"""

gamma = mpfunc('gamma', gammazeta.mpf_gamma, gammazeta.mpc_gamma, "gamma function")
factorial = mpfunc('factorial', gammazeta.mpf_factorial, gammazeta.mpc_factorial, "factorial")
fac = factorial

def psi(m, z):
    """
    Gives the polygamma function of order m of z, psi^(m)(z). Special
    cases are the digamma function (psi0), trigamma function (psi1),
    tetragamma (psi2) and pentagamma (psi4) functions.

    The parameter m should be a nonnegative integer.
    """
    z = mpmathify(z)
    m = int(m)
    if isinstance(z, mpf):
        return make_mpf(gammazeta.mpf_psi(m, z._mpf_, *prec_rounding))
    else:
        return make_mpc(gammazeta.mpc_psi(m, z._mpc_, *prec_rounding))

def psi0(z):
    """Shortcut for psi(0,z) (the digamma function)"""
    return psi(0, z)

def psi1(z):
    """Shortcut for psi(1,z) (the trigamma function)"""
    return psi(1, z)

def psi2(z):
    """Shortcut for psi(2,z) (the tetragamma function)"""
    return psi(2, z)

def psi3(z):
    """Shortcut for psi(3,z) (the pentagamma function)"""
    return psi(3, z)

polygamma = psi
digamma = psi0
trigamma = psi1
tetragamma = psi2
pentagamma = psi3

harmonic = mpfunc('harmonic', gammazeta.mpf_harmonic, gammazeta.mpc_harmonic,
    "nth harmonic number")

harmonic.__doc__ = r"""
    If `n` is an integer, ``harmonic(n)`` gives a floating-point
    approximation of the `n`-th harmonic number `H(n)`, defined as

    .. math ::

        H(n) = 1 + \frac{1}{2} + \frac{1}{3} + \ldots + \frac{1}{n}

    The firrst few harmonic numbers are::

        >>> from mpmath import *
        >>> mp.dps = 15    
        >>> for n in range(8):
        ...     print n, harmonic(n)
        ...
        0 0.0
        1 1.0
        2 1.5
        3 1.83333333333333
        4 2.08333333333333
        5 2.28333333333333
        6 2.45
        7 2.59285714285714

    The infinite harmonic series `1 + 1/2 + 1/3 + \ldots` diverges::

        >>> print harmonic(inf)
        +inf

    :func:`harmonic` is evaluated using the digamma function rather
    than by summing the harmonic series term by term. It can therefore
    be computed quickly for arbitrarily large `n`, and even for
    nonintegral arguments::

        >>> print harmonic(10**100)
        230.835724964306
        >>> print harmonic(0.5)
        0.613705638880109
        >>> print harmonic(3+4j)
        (2.24757548223494 + 0.850502209186044j)

    :func:`harmonic` supports arbitrary precision evaluation::

        >>> mp.dps = 50
        >>> print harmonic(11)
        3.0198773448773448773448773448773448773448773448773
        >>> print harmonic(pi)
        1.8727388590273302654363491032336134987519132374152

    The harmonic series diverges, but at a glacial pace. It is possible
    to calculate the exact number of terms required before the sum
    exceeds a given amount, say 100::

        >>> mp.dps = 50
        >>> v = 10**findroot(lambda x: harmonic(10**x) - 100, 10)
        >>> print v
        15092688622113788323693563264538101449859496.864101
        >>> v = int(ceil(v))
        >>> print v
        15092688622113788323693563264538101449859497
        >>> print harmonic(v-1)
        99.999999999999999999999999999999999999999999942747
        >>> print harmonic(v)
        100.000000000000000000000000000000000000000000009

    """

def bernoulli(n):
    r"""
    Computes the nth Bernoulli number, `B_n`, for any integer `n \ge 0`.

    The Bernoulli numbers are rational numbers, but this function
    returns a floating-point approximation. To obtain an exact
    fraction, use :func:`bernfrac` instead.

    **Examples**

    Numerical values of the first few Bernoulli numbers::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> for n in range(15):
        ...     print n, bernoulli(n)
        ...
        0 1.0
        1 -0.5
        2 0.166666666666667
        3 0.0
        4 -0.0333333333333333
        5 0.0
        6 0.0238095238095238
        7 0.0
        8 -0.0333333333333333
        9 0.0
        10 0.0757575757575758
        11 0.0
        12 -0.253113553113553
        13 0.0
        14 1.16666666666667

    Bernoulli numbers can be approximated with arbitrary precision::

        >>> mp.dps = 50
        >>> print bernoulli(100)
        -2.8382249570693706959264156336481764738284680928013e+78

    Arbitrarily large `n` are supported::

        >>> mp.dps = 15
        >>> print bernoulli(10**20 + 2)
        3.09136296657021e+1876752564973863312327

    The Bernoulli numbers are related to the Riemann zeta function
    at integer arguments::

        >>> print -bernoulli(8) * (2*pi)**8 / (2*fac(8))
        1.00407735619794
        >>> print zeta(8)
        1.00407735619794

    **Algorithm**

    For small `n` (`n < 3000`) :func:`bernoulli` uses a recurrence
    formula due to Ramanujan. All results in this range are cached,
    so sequential computation of small Bernoulli numbers is
    guaranteed to be fast.

    For larger `n`, `B_n` is evaluated in terms of the Riemann zeta
    function.
    """
    return make_mpf(gammazeta.mpf_bernoulli(int(n), *prec_rounding))

bernfrac = gammazeta.bernfrac

stieltjes_cache = {}

def stieltjes(n):
    r"""
    For a nonnegative integer `n`, ``stieltjes(n)`` computes the
    `n`-th Stieltjes constant `\gamma_n`, defined as the
    `n`-th coefficient in the Laurent series expansion of the
    Riemann zeta function around the pole at `s = 1`. That is,
    we have:

    .. math ::

      \zeta(s) = \frac{1}{s-1} \sum_{n=0}^{\infty}
          \frac{(-1)^n}{n!} \gamma_n (s-1)^n

    **Examples**

    The zeroth Stieltjes constant is just Euler's constant `\gamma`::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print stieltjes(0)
        0.577215664901533

    Some more values are::

        >>> print stieltjes(1)
        -0.0728158454836767
        >>> print stieltjes(10)
        0.000205332814909065
        >>> print stieltjes(30)
        0.00355772885557316

    An alternative way to compute `\gamma_1`::

        >>> print diff(extradps(25)(lambda x: 1/(x-1) - zeta(x)), 1)
        -0.0728158454836767

    :func:`stieltjes` supports arbitrary precision evaluation,
    and caches computed results::

        >>> mp.dps = 50
        >>> print stieltjes(2)
        -0.0096903631928723184845303860352125293590658061013408

    **Algorithm**

    The calculation is done using numerical differentiation
    for very small `n` (currently `n = 1,2,3`).

    For larger `n`, integration of the Riemann zeta function is
    used. The method should work for any `n` and precision, but
    soon becomes quite slow in practice. The code has been tested
    with `n = 50` and 100 digit precision; that computation took
    about 2 minutes.

    **References**

    1. http://mathworld.wolfram.com/StieltjesConstants.html

    2. http://pi.lacim.uqam.ca/piDATA/stieltjesgamma.txt

    """
    n = int(n)
    if n == 0:
        return +euler
    if n < 0:
        raise ValueError("Stieltjes constants defined for n >= 0")
    if n in stieltjes_cache:
        prec, s = stieltjes_cache[n]
        if prec >= mp.prec:
            return +s
    orig = mp.prec
    try:
        if n <= 3:
            from calculus import diff
            # The Stieltjes constants appear as derivatives of the Dirichlet
            # eta function at s = 1 (although a lot of cruft needs to be
            # included in the formulas). The following special cases were
            # derived using Mathematica.
            #
            # It is possible that a simple recursion formula could be
            # determined. However, this would not necessarily be worthwhile.
            #
            # Note that the nth derivative requires (n+1)*prec, so this
            # is only faster than the integration method for really small n
            # anyway.
            if n == 1:
                v = -(diff(altzeta, 1, 2, direction=1)*3 + 3*euler*log(2)**2-\
                    log(2)**3) / (3*log(4))
            if n == 2:
                v = (diff(altzeta, 1, 3, direction=1) - euler*log(2)**3 + \
                    log(2)**4/4 - 3*log(2)**2*stieltjes(1))/log(8)
            if n == 3:
                v = -(diff(altzeta, 1, 4, direction=1) + euler*log(2)**4 - \
                    log(2)**5/5 + 4*log(2)**3*stieltjes(1) + 6*log(2)**2* \
                    stieltjes(2)) / (4*log(2))
        else:
            from quadrature import quadgl
            def f(x):
                r = exp(pi*j*x)
                return (zeta(r+1) / r**n).real
            p = int(log(factorial(n), 2) + 35)
            mp.prec += p
            u = quadgl(f, [-1, 1])
            v = mpf(-1)**n * factorial(n) * u / 2
    finally:
        mp.prec = orig
    stieltjes_cache[n] = (mp.prec, v)
    return +v

def isnpint(x):
    if not x:
        return True
    if isinstance(x, mpf):
        sign, man, exp, bc = x._mpf_
        return sign and exp >= 0
    if isinstance(x, mpc):
        return not x.imag and isnpint(x.real)

def gammaprod(a, b):
    r"""
    Given iterables `a` and `b`, ``gammaprod(a, b)`` computes the
    product / quotient of gamma functions:

    .. math ::

        \frac{\Gamma(a_0) \Gamma(a_1) \cdots \Gamma(a_p)}
             {\Gamma(b_0) \Gamma(b_1) \cdots \Gamma(b_q)}

    **Handling of poles**

    Unlike direct calls to :func:`gamma`, :func:`gammaprod` considers
    the entire product as a limit and evaluates this limit properly if
    any of the numerator or denominator arguments are nonpositive
    integers such that poles of the gamma function are encountered.

    In particular:

    * If there are equally many poles in the numerator and the
      denominator, the limit is a rational number times the remaining,
      regular part of the product.

    * If there are more poles in the numerator, :func:`gammaprod`
      returns ``+inf``.

    * If there are more poles in the denominator, :func:`gammaprod`
      returns 0.

    **Examples**

    The reciprocal gamma function `1/\Gamma(x)` evaluated at `x = 0`::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> gammaprod([], [0])
        mpf('0.0')

    A limit::

        >>> gammaprod([-4], [-3])
        mpf('-0.25')
        >>> limit(lambda x: gamma(x-1)/gamma(x), -3, direction=1)
        mpf('-0.25')
        >>> limit(lambda x: gamma(x-1)/gamma(x), -3, direction=-1)
        mpf('-0.25')

    """
    a = [mpmathify(x) for x in a]
    b = [mpmathify(x) for x in b]
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

def beta(x, y):
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
    x = mpmathify(x)
    y = mpmathify(y)
    if isinf(y):
        x, y = y, x
    if isinf(x):
        if x == inf and not y.imag:
            if y == -inf:
                return nan
            if y > 0:
                return zero
            if isint(y):
                return nan
            if y < 0:
                return sign(gamma(y)) * inf
        return nan
    return gammaprod([x, y], [x+y])

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
    x = mpmathify(x)
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
            re, im = x._mpf_, libmpf.fzero
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
    a = (ar and _as_num(ar[0])) or mpmathify(a)
    b = (br and _as_num(br[0])) or mpmathify(b)
    c = (cr and _as_num(cr[0])) or mpmathify(c)
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
    r"""
    Evaluates the generalized hypergeometric function

    .. math ::

        \,_pF_q(a_1,\ldots,a_p; b_1,\ldots,b_q; z) =
        \sum_{n=0}^\infty \frac{(a_1)_n (a_2)_n \ldots (a_p)_n}
           {(b_1)_n(b_2)_n\ldots(b_q)_n} \frac{z^n}{n!}

    where `(x)_n` denotes the rising factorial (see :func:`rf`).

    The parameters lists ``a_s`` and ``b_s`` may contain integers,
    real numbers, complex numbers, as well as exact fractions given in
    the form of tuples `(p, q)`. :func:`hyper` is optimized to handle
    integers and fractions more efficiently than arbitrary
    floating-point parameters (since rational parameters are by
    far the most common).

    **Examples**

    We can compare the output of :func:`hyper` with :func:`nsum`::

        >>> from mpmath import *
        >>> mp.dps = 25
        >>> a,b,c,d = 2,3,4,5
        >>> x = 0.25
        >>> print hyper([a,b],[c,d],x)
        1.078903941164934876086237
        >>> fn = lambda n: rf(a,n)*rf(b,n)/rf(c,n)/rf(d,n)*x**n/fac(n)
        >>> print nsum(fn, [0, inf])
        1.078903941164934876086237

    The parameters can be any combination of integers, fractions,
    floats and complex numbers::

        >>> a, b, c, d, e = 1, (-1,2), pi, 3+4j, (2,3)
        >>> x = 0.2j
        >>> print hyper([a,b],[c,d,e],x)
        (0.9923571616434024810831887 - 0.005753848733883879742993122j)
        >>> b, e = -0.5, mpf(2)/3
        >>> fn = lambda n: rf(a,n)*rf(b,n)/rf(c,n)/rf(d,n)/rf(e,n)*x**n/fac(n)
        >>> print nsum(fn, [0, inf])
        (0.9923571616434024810831887 - 0.005753848733883879742993122j)

    """
    p = len(a_s)
    q = len(b_s)
    z = mpmathify(z)
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
    r"""Hypergeometric function `\,_0F_1`. ``hyp0f1(a,z)`` is equivalent
    to ``hyper([],[a],z)``; see documentation for :func:`hyper` for more
    information."""
    return hyper([], [a], z)

def hyp1f1(a,b,z):
    r"""Hypergeometric function `\,_1F_1`. ``hyp1f1(a,b,z)`` is equivalent
    to ``hyper([a],[b],z)``; see documentation for :func:`hyper` for more
    information."""
    return hyper([a], [b], z)

def hyp2f1(a,b,c,z):
    r"""Hypergeometric function `\,_2F_1`. ``hyp2f1(a,b,c,z)`` is equivalent
    to ``hyper([a,b],[c],z)``; see documentation for :func:`hyper` for more
    information."""
    return hyper([a,b], [c], z)

def _lower_gamma(z, b):
    return hyp1f1(1, 1+z, b) * b**z * exp(-b) / z

def _check_pos(x):
    return isinstance(x, mpf) and x > 0

@funcwrapper
def gammainc(z, a=0, b=inf, regularized=False):
    r"""
    ``gammainc(z, a=0, b=inf)`` computes the (generalized) incomplete
    gamma function with integration limits `[a, b]`:

    .. math ::

      \Gamma(z,a,b) = \int_a^b t^{z-1} e^{-t} \, dt

    The generalized incomplete gamma function reduces to the
    following special cases when one or both endpoints are fixed:

    * `\Gamma(z,0,\infty)` is the standard ("complete")
      gamma function, `\Gamma(z)` (available directly
      as the mpmath function :func:`gamma`)
    * `\Gamma(z,a,\infty)` is the "upper" incomplete gamma
      function, `\Gamma(z,a)`
    * `\Gamma(z,0,b)` is the "lower" incomplete gamma
      function, `\gamma(z,b)`.

    Of course, we have
    `\Gamma(z,0,x) + \Gamma(z,x,\infty) = \Gamma(z)`
    for all `z` and `x`.

    Note however that some authors reverse the order of the
    arguments when defining the lower and upper incomplete
    gamma function, so one should be careful to get the correct
    definition.

    If also given the keyword argument ``regularized=True``,
    :func:`gammainc` computes the "regularized" incomplete gamma
    function

    .. math ::

      P(z,a,b) = \frac{\Gamma(z,a,b)}{\Gamma(z)}.

    **Examples**

    We can compare with numerical quadrature to verify that
    :func:`gammainc` computes the integral in the definition::

        >>> from mpmath import *
        >>> mp.dps = 20
        >>> print gammainc(2+3j, 4, 10)
        (0.009772126686277051606 - 0.077063730631298989245j)
        >>> print quad(lambda t: t**(2+3j-1) * exp(-t), [4, 10])
        (0.009772126686277051606 - 0.077063730631298989245j)

    The incomplete gamma functions satisfy simple recurrence
    relations::

        >>> mp.dps = 15
        >>> z = 3.5
        >>> a = 2
        >>> print gammainc(z+1, a), z*gammainc(z,a) + a**z*exp(-a)
        10.6013029693353 10.6013029693353
        >>> print gammainc(z+1,0,a), z*gammainc(z,0,a) - a**z*exp(-a)
        1.03042542723211 1.03042542723211

    If `z` is an integer, the recurrence reduces the incomplete gamma
    function to `P(a) \exp(-a) + Q(b) \exp(-b)` where `P` and
    `Q` are polynomials::

        >>> mp.dps = 15
        >>> print gammainc(1, 2), exp(-2)
        0.135335283236613 0.135335283236613
        >>> mp.dps = 50
        >>> identify(gammainc(6, 1, 2), ['exp(-1)', 'exp(-2)'])
        '(326*exp(-1) + (-872)*exp(-2))'

    The incomplete gamma functions reduce to functions such as
    the exponential integral Ei and the error function for special
    arguments::

        >>> mp.dps = 15
        >>> print gammainc(0, 4), -ei(-4)
        0.00377935240984891 0.00377935240984891
        >>> print gammainc(0.5, 0, 2), sqrt(pi)*erf(sqrt(2))
        1.6918067329452 1.6918067329452

    """
    if b == inf:
        if a == 0:
            v = gamma(z)
        else:
            if z == 0:
                # Reduces to exponential integral. Mind branch cuts.
                if _check_pos(a):
                    return -ei(-a)
                else:
                    return -ei(-a) + (log(-a)-log(-1/a))/2-log(a)
            # XXX: avoid poles
            v = gamma(z) - _lower_gamma(z, a)
    elif a == 0:
        v = _lower_gamma(z, b)
    else:
        if z == 0:
            # Reduces to exponential integral
            if _check_pos(a) and _check_pos(b):
                return ei(-b) - ei(-a)
            else:
                return ei(-b)-ei(-a) + \
                    (log(-a)-log(-1/a))/2-log(a) + \
                    (log(-1/b)-log(-b))/2+log(b)
        # XXX: avoid poles
        v = _lower_gamma(z, b) - _lower_gamma(z, a)
    if regularized:
        return v / gamma(z)
    else:
        return v


erf = mpfunc("erf", libhyper.mpf_erf, libhyper.mpc_erf,
    "Error function, erf(z)")

erf.__doc__ = r"""
Computes the error function, `\mathrm{erf}(x)`. The error
function is the normalized antiderivative of the Gaussian function
`\exp(-t^2)`. More precisely,

.. math::

  \mathrm{erf}(x) = \frac{2}{\sqrt \pi} \int_0^x \exp(-t^2) \,dt

**Basic examples**

Simple values and limits include::

    >>> from mpmath import *
    >>> mp.dps = 15
    >>> print erf(0)
    0.0
    >>> print erf(1)
    0.842700792949715
    >>> print erf(-1)
    -0.842700792949715
    >>> print erf(inf)
    1.0
    >>> print erf(-inf)
    -1.0

For large real `x`, `\mathrm{erf}(x)` approaches 1 very
rapidly::

    >>> print erf(3)
    0.999977909503001
    >>> print erf(5)
    0.999999999998463

The error function is an odd function::

    >>> nprint(taylor(erf, 0, 5))
    [0.0, 1.12838, 0.0, -0.376126, 0.0, 0.112838]

:func:`erf` implements arbitrary-precision evaluation and
supports complex numbers::

    >>> mp.dps = 50
    >>> print erf(0.5)
    0.52049987781304653768274665389196452873645157575796
    >>> mp.dps = 25
    >>> print erf(1+j)
    (1.316151281697947644880271 + 0.1904534692378346862841089j)

**Related functions**

See also :func:`erfc`, which is more accurate for large `x`,
and :func:`erfi` which gives the antiderivative of
`\exp(t^2)`.

The Fresnel integrals :func:`fresnels` and :func:`fresnelc`
are also related to the error function.
"""


erfc = mpfunc("erfc", libhyper.mpf_erfc, libhyper.mpc_erfc,
    "Complementary error function, erfc(z) = 1-erf(z)")

erfc.__doc__ = r"""
Computes the complementary error function,
`\mathrm{erfc}(x) = 1-\mathrm{erf}(x)`.
This function avoids cancellation that occurs when naively
computing the complementary error function as ``1-erf(x)``::

    >>> from mpmath import *
    >>> mp.dps = 15
    >>> print 1 - erf(10)
    0.0
    >>> print erfc(10)
    2.08848758376254e-45

:func:`erfc` works accurately even for ludicrously large
arguments::

    >>> print erfc(10**10)
    4.3504398860243e-43429448190325182776

"""

@funcwrapper
def erfi(z):
    r"""
    Computes the imaginary error function, `\mathrm{erfi}(x)`.
    The imaginary error function is defined in analogy with the
    error function, but with a positive sign in the integrand:

    .. math ::

      \mathrm{erfi}(x) = \frac{2}{\sqrt \pi} \int_0^x \exp(t^2) \,dt

    Whereas the error function rapidly converges to 1 as `x` grows,
    the imaginary error function rapidly diverges to infinity.
    The functions are related as
    `\mathrm{erfi}(x) = -i\,\mathrm{erf}(ix)` for all complex
    numbers `x`.

    **Examples**

    Basic values and limits::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print erfi(0)
        0.0
        >>> print erfi(1)
        1.65042575879754
        >>> print erfi(-1)
        -1.65042575879754
        >>> print erfi(inf)
        +inf
        >>> print erfi(-inf)
        -inf

    Note the symmetry between erf and erfi::

        >>> print erfi(3j)
        (0.0 + 0.999977909503001j)
        >>> print erf(3)
        0.999977909503001
        >>> print erf(1+2j)
        (-0.536643565778565 - 5.04914370344703j)
        >>> print erfi(2+1j)
        (-5.04914370344703 - 0.536643565778565j)

    **Possible issues**

    The current implementation of :func:`erfi` is much less efficient
    and accurate than the one for erf.

    """
    return (2/sqrt(pi)*z) * sum_hyp1f1_rat((1,2),(3,2), z**2)

@funcwrapper
def erfinv(x):
    r"""
    Computes the inverse error function, satisfying

    .. math ::

        \mathrm{erf}(\mathrm{erfinv}(x)) =
        \mathrm{erfinv}(\mathrm{erf}(x)) = x.

    This function is defined only for `-1 \le x \le 1`.

    **Examples**

    Special values include::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print erfinv(0)
        0.0
        >>> print erfinv(1)
        +inf
        >>> print erfinv(-1)
        -inf

    The domain is limited to the standard interval::

        >>> erfinv(2)
        Traceback (most recent call last):
          ...
        ValueError: erfinv(x) is defined only for -1 <= x <= 1

    It is simple to check that :func:`erfinv` computes inverse values of
    :func:`erf` as promised::

        >>> print erf(erfinv(0.75))
        0.75
        >>> print erf(erfinv(-0.995))
        -0.995

    :func:`erfinv` supports arbitrary-precision evaluation::

        >>> mp.dps = 50
        >>> erf(3)
        mpf('0.99997790950300141455862722387041767962015229291260075')
        >>> erfinv(_)
        mpf('3.0')

    A definite integral involving the inverse error function::

        >>> mp.dps = 15
        >>> print quad(erfinv, [0, 1])
        0.564189583547756
        >>> print 1/sqrt(pi)
        0.564189583547756

    The inverse error function can be used to generate random numbers
    with a Gaussian distribution (although this is a relatively
    inefficient algorithm)::

        >>> nprint([erfinv(2*rand()-1) for n in range(6)]) # doctest: +SKIP
        [-0.586747, 1.10233, -0.376796, 0.926037, -0.708142, -0.732012]

    """
    if x.imag or (x < -1) or (x > 1):
        raise ValueError("erfinv(x) is defined only for -1 <= x <= 1")
    if isnan(x): return x
    if not x: return x
    if x == 1: return inf
    if x == -1: return -inf
    if abs(x) < 0.9:
        a = 0.53728*x**3 + 0.813198*x
    else:
        # An asymptotic formula
        u = log(2/pi/(abs(x)-1)**2)
        a = sign(x) * sqrt(u - log(u))/sqrt(2)
    from optimization import findroot
    return findroot(lambda t: erf(t)-x, a)

@funcwrapper
def npdf(x, mu=0, sigma=1):
    r"""
    ``npdf(x, mu=0, sigma=1)`` evaluates the probability density
    function of a normal distribution with mean value `\mu`
    and variance `\sigma^2`.

    Elementary properties of the probability distribution can
    be verified using numerical integration::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print quad(npdf, [-inf, inf])
        1.0
        >>> print quad(lambda x: npdf(x, 3), [3, inf])
        0.5
        >>> print quad(lambda x: npdf(x, 3, 2), [3, inf])
        0.5

    See also :func:`ncdf`, which gives the cumulative
    distribution.
    """
    sigma = mpmathify(sigma)
    return exp(-(x-mu)**2/(2*sigma**2)) / (sigma*sqrt(2*pi))

@funcwrapper
def ncdf(x, mu=0, sigma=1):
    r"""
    ``ncdf(x, mu=0, sigma=1)`` evaluates the cumulative distribution
    function of a normal distribution with mean value `mu`
    and variance `\sigma^2`.

    See also :func:`npdf`, which gives the probability density.

    Elementary properties include::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print ncdf(pi, mu=pi)
        0.5
        >>> print ncdf(-inf)
        0.0
        >>> print ncdf(+inf)
        1.0

    The cumulative distribution is the integral of the density
    function having identical mu and sigma::

        >>> mp.dps = 15
        >>> print diff(ncdf, 2)
        0.053990966513188
        >>> print npdf(2)
        0.053990966513188
        >>> print diff(lambda x: ncdf(x, 1, 0.5), 0)
        0.107981933026376
        >>> print npdf(0, 1, 0.5)
        0.107981933026376

    """
    a = (x-mu)/(sigma*sqrt(2))
    if a < 0:
        return erfc(-a)/2
    else:
        return (1+erf(a))/2

@funcwrapper
def ei(z):
    r"""
    Computes the exponential integral or Ei-function, `\mathrm{Ei}(x)`.
    The exponential integral is defined as

    .. math ::

      \mathrm{Ei}(x) = \int_{-\infty\,}^x \frac{e^t}{t} \, dt.

    When the integration range includes `t = 0`, the exponential
    integral is interpreted as providing the Cauchy principal value.

    For real `x`, the Ei-function behaves roughly like
    `\mathrm{Ei}(x) \approx \exp(x) + \log(|x|)`.

    This function should not be confused with the family of related
    functions denoted by `E_n` which are also called "exponential
    integrals".

    **Basic examples**

    Some basic values and limits are::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print ei(0)
        -inf
        >>> print ei(1)
        1.89511781635594
        >>> print ei(inf)
        +inf
        >>> print ei(-inf)
        0.0

    For `x < 0`, the defining integral can be evaluated
    numerically as a reference::

        >>> print ei(-4)
        -0.00377935240984891
        >>> print quad(lambda t: exp(t)/t, [-inf, -4])
        -0.00377935240984891

    :func:`ei` supports complex arguments and arbitrary
    precision evaluation::

        >>> mp.dps = 50
        >>> mp.dps = 50
        >>> print ei(pi)
        10.928374389331410348638445906907535171566338835056
        >>> mp.dps = 25
        >>> print ei(3+4j)
        (-4.154091651642689822535359 + 4.294418620024357476985535j)

    **Related functions**

    The exponential integral is closely related to the logarithmic
    integral. See :func:`li` for additional information.

    The exponential integral is related to the hyperbolic
    and trigonometric integrals (see :func:`chi`, :func:`shi`,
    :func:`ci`, :func:`si`) similarly to how the ordinary
    exponential function is related to the hyperbolic and
    trigonometric functions::

        >>> mp.dps = 15
        >>> print ei(3)
        9.93383257062542
        >>> print chi(3) + shi(3)
        9.93383257062542
        >>> print ci(3j) - j*si(3j) - pi*j/2
        (9.93383257062542 + 0.0j)

    Beware that logarithmic corrections, as in the last example
    above, are required to obtain the correct branch in general.
    For details, see [1].

    The exponential integral is also a special case of the
    hypergeometric function `\,_2F_2`::

        >>> z = 0.6
        >>> print z*hyper([1,1],[2,2],z) + (ln(z)-ln(1/z))/2 + euler
        0.769881289937359
        >>> print ei(z)
        0.769881289937359

    **References**

    1. Relations between Ei and other functions:
       http://functions.wolfram.com/GammaBetaErf/ExpIntegralEi/27/01/

    2. Abramowitz & Stegun, section 5:
       http://www.math.sfu.ca/~cbm/aands/page_228.htm

    """
    if z == inf:
        return z
    if z == -inf:
        return -mpf(0)
    if not z:
        return -inf
    v = z*hypsum([[1,1],[1,1]],[],[],[[2,1],[2,1]],[],[],z) + \
        (log(z)-log(1/z))/2 + euler
    if isinstance(z, mpf) and z < 0:
        return v.real
    return v

@funcwrapper
def li(z):
    r"""
    Computes the logarithmic integral or li-function
    `\mathrm{li}(x)`, defined by

    .. math ::

        \mathrm{li}(x) = \int_0^x \frac{1}{\log t} \, dt

    The logarithmic integral has a singularity at `x = 1`.

    Note that there is a second logarithmic integral, the Li
    function, defined by

    .. math ::

        \mathrm{Li}(x) = \int_2^x \frac{1}{\log t} \, dt

    This "offset logarithmic integral" can be computed via
    :func:`li` using the simple identity
    `\mathrm{Li}(x) = \mathrm{li}(x) - \mathrm{li}(2)`.

    **Examples**

    Some basic values and limits::

        >>> from mpmath import *
        >>> mp.dps = 30
        >>> print li(0)
        0.0
        >>> print li(1)
        -inf
        >>> print li(1)
        -inf
        >>> print li(2)
        1.04516378011749278484458888919
        >>> print findroot(li, 2)
        1.45136923488338105028396848589
        >>> print li(inf)
        +inf

    The logarithmic integral can be evaluated for arbitrary
    complex arguments::

        >>> mp.dps = 20
        >>> print li(3+4j)
        (3.1343755504645775265 + 2.6769247817778742392j)

    The logarithmic integral is related to the exponential integral::

        >>> print ei(log(3))
        2.1635885946671919729
        >>> print li(3)
        2.1635885946671919729

    The logarithmic integral grows like `O(x/\log(x))`::

        >>> mp.dps = 15
        >>> x = 10**100
        >>> print x/log(x)
        4.34294481903252e+97
        >>> print li(x)
        4.3619719871407e+97

    The prime number theorem states that the number of primes less
    than `x` is asymptotic to `\mathrm{li}(x)`. For example,
    it is known that there are exactly 1,925,320,391,606,803,968,923
    prime numbers less than 10^23 [1]. The logarithmic integral
    provides a very accurate estimate::

        >>> print li(2) + li(10**23)
        1.92532039161405e+21

    A definite integral is::

        >>> print quad(li, [0, 1])
        -0.693147180559945
        >>> print -ln(2)
        -0.693147180559945

    **References**

    1. http://mathworld.wolfram.com/PrimeCountingFunction.html

    2. http://mathworld.wolfram.com/LogarithmicIntegral.html

    """
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
    r"""
    Computes the Airy function `\mathrm{Ai}(x)`, which is
    a solution of the Airy differential equation `y''-xy=0`.
    The Ai-function behaves roughly like a slowly decaying
    sine wave for `x < 0` and like a decreasing exponential for
    `x > 0`.

    Limits and values include::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print airyai(0), 1/(3**(2/3.)*gamma(2/3.))
        0.355028053887817 0.355028053887817
        >>> print airyai(1)
        0.135292416312881
        >>> print airyai(-1)
        0.535560883292352
        >>> print airyai(inf)
        0.0
        >>> print airyai(-inf)
        0.0

    :func:`airyai` uses a series expansion around `x = 0`,
    so it is slow for extremely large arguments. Here are
    some evaluations for moderately large arguments::

        >>> print airyai(-100)
        0.176753393239553
        >>> print airyai(100)
        2.63448215208818e-291
        >>> print airyai(50+50j)
        (-5.31790195707456e-68 - 1.16358800377071e-67j)
        >>> print airyai(-50+50j)
        (1.04124253736317e+158 + 3.3475255449236e+157j)

    The first negative root is::

        >>> print findroot(airyai, -2)
        -2.33810741045977

    We can verify the differential equation::

        >>> for x in [-3.4, 0, 2.5, 1+2j]:
        ...     print abs(diff(airyai, x, 2) - x*airyai(x)) < eps
        ...
        True
        True
        True
        True

    The Taylor series expansion around `x = 0` starts with
    the following coefficients (note that every third term
    is zero)::

        >>> nprint([(abs(c)>eps) and c or 0 for c in taylor(airyai, 0, 5)])
        [0.355028, -0.258819, 0, 5.91713e-2, -2.15683e-2, 0]

    The Airy functions are a special case of Bessel functions.
    For `x < 0`, we have::

        >>> x = 3
        >>> print airyai(-x)
        -0.378814293677658
        >>> p = 2*(x**1.5)/3
        >>> print sqrt(x)*(jv(1/3.,p) + jv(-1/3.,p))/3
        -0.378814293677658

    """
    if z == inf or z == -inf:
        return 1/z
    if z.real > 2:
        # cancellation: both terms are ~ 2^(z^1.5),
        # result is ~ 2^(-z^1.5), so need ~2*z^1.5 extra bits
        mp.prec += 2*int(z.real**1.5)
    z3 = z**3 / 9
    a = sum_hyp0f1_rat((2,3), z3) / (cbrt(9) * gamma(mpf(2)/3))
    b = z * sum_hyp0f1_rat((4,3), z3) / (cbrt(3) * gamma(mpf(1)/3))
    return a - b

@funcwrapper
def airybi(z):
    r"""
    Computes the Airy function `\mathrm{Bi}(x)`, which is
    a solution of the Airy differential equation `y''-xy=0`.
    The Bi-function behaves roughly like a slowly decaying
    sine wave for `x < 0` and like an increasing exponential
    for `x > 0`.

    Limits and values include::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print airybi(0), 1/(3**(1/6.)*gamma(2/3.))
        0.614926627446001 0.614926627446001
        >>> print airybi(1)
        1.20742359495287
        >>> print airybi(-1)
        0.103997389496945
        >>> print airybi(inf)
        +inf
        >>> print airybi(-inf)
        0.0

    :func:`airyai` uses a series expansion around `x = 0`,
    so it is slow for extremely large arguments. Here are
    some evaluations for moderately large arguments::

        >>> print airybi(-100)
        0.0242738876801601
        >>> print airybi(100)
        6.0412239966702e+288
        >>> print airybi(50+50j)
        (-5.32207626732144e+63 + 1.47845029116524e+65j)
        >>> print airybi(-50+50j)
        (-3.3475255449236e+157 + 1.04124253736317e+158j)

    The first negative root is::

        >>> print findroot(airybi, -1)
        -1.17371322270913

    We can verify the differential equation::

        >>> for x in [-3.4, 0, 2.5, 1+2j]:
        ...     print abs(diff(airybi, x, 2) - x*airybi(x)) < eps
        ...
        True
        True
        True
        True

    The Taylor series expansion around `x = 0` starts with
    the following coefficients (note that every third term
    is zero)::

        >>> nprint([(abs(c)>eps) and c or 0 for c in taylor(airybi, 0, 5)])
        [0.614927, 0.448288, 0, 0.102488, 3.73574e-2, 0]

    The Airy functions are a special case of Bessel functions.
    For `x < 0`, we have::

        >>> x = 3
        >>> print airybi(-x)
        -0.198289626374927
        >>> p = 2*(x**1.5)/3
        >>> print sqrt(x/3)*(jv(-1/3.,p) - jv(1/3.,p))
        -0.198289626374926
    """
    if z == inf:
        return z
    if z == -inf:
        return 1/z
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

@funcwrapper
def agm(a, b=1):
    r"""
    ``agm(a, b)`` computes the arithmetic-geometric mean of `a` and
    `b`, defined as the limit of the following iteration:

    .. math ::

        a_0 = a

        b_0 = b

        a_{n+1} = \frac{a_n+b_n}{2}

        b_{n+1} = \sqrt{a_n b_n}

    This function can be called with a single argument, computing
    `\mathrm{agm}(a,1) = \mathrm{agm}(1,a)`.

    **Examples**

    It is a well-known theorem that the geometric mean of
    two distinct positive numbers is less than the arithmetic
    mean. It follows that the arithmetic-geometric mean lies
    between the two means::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> a = mpf(3)
        >>> b = mpf(4)
        >>> print sqrt(a*b)
        3.46410161513775
        >>> print agm(a,b)
        3.48202767635957
        >>> print (a+b)/2
        3.5

    The arithmetic-geometric mean is scale-invariant::

        >>> print agm(10*e, 10*pi)
        29.261085515723
        >>> print 10*agm(e, pi)
        29.261085515723

    As an order-of-magnitude estimate, `\mathrm{agm}(1,x) \approx x`
    for large `x`::

        >>> print agm(10**10)
        643448704.760133
        >>> print agm(10**50)
        1.34814309345871e+48

    The arithmetic-geometric mean can also be computed for complex
    numbers::

        >>> print agm(3, 2+j)
        (2.51055133276184 + 0.547394054060638j)

    The AGM iteration converges very quickly (each step doubles
    the number of correct digits), so :func:`agm` supports efficient
    high-precision evaluation::

        >>> mp.dps = 10000
        >>> a = agm(1,2)
        >>> str(a)[-10:]
        '1679581912'

    **Mathematical relations**

    The arithmetic-geometric mean may be used to evaluate the
    following two parametric definite integrals:

    .. math ::

      I_1 = \int_0^{\infty}
        \frac{1}{\sqrt{(x^2+a^2)(x^2+b^2)}} \,dx

      I_2 = \int_0^{\pi/2}
        \frac{1}{\sqrt{a^2 \cos^2(x) + b^2 \sin^2(x)}} \,dx

    We have::

        >>> mp.dps = 15
        >>> a = 3
        >>> b = 4
        >>> f1 = lambda x: ((x**2+a**2)*(x**2+b**2))**-0.5
        >>> f2 = lambda x: ((a*cos(x))**2 + (b*sin(x))**2)**-0.5
        >>> print quad(f1, [0, inf])
        0.451115405388492
        >>> print quad(f2, [0, pi/2])
        0.451115405388492
        >>> print pi/(2*agm(a,b))
        0.451115405388492

    A formula for `\Gamma(1/4)`::

        >>> print gamma(0.25)
        3.62560990822191
        >>> print sqrt(2*sqrt(2*pi**3)/agm(1,sqrt(2)))
        3.62560990822191

    **Possible issues**

    The branch cut chosen for complex `a` and `b` is somewhat
    arbitrary.

    """
    if not a or not b:
        return a*b
    weps = eps * 16 * max(abs(a), abs(b))
    half = mpf(0.5)
    while abs(a-b) > weps:
        a, b = (a+b)*half, (a*b)**half
    return a

@funcwrapper
def jacobi(n, a, b, x):
    r"""
    ``jacobi(n, a, b, x)`` evaluates the Jacobi polynomial
    `P_n^{(a,b)}(x)`. The Jacobi polynomials are a special
    case of the hypergeometric function `\,_2F_1` given by:

    .. math ::

        P_n^{(a,b)}(x) = {n+a \choose n}
          \,_2F_1\left(-n,1+a+b+n,a+1,\frac{1-x}{2}\right).

    Note that this definition generalizes to nonintegral values
    of `n`. When `n` is an integer, the hypergeometric series
    terminates after a finite number of terms, giving
    a polynomial in `x`.

    **Evaluation of Jacobi polynomials**

    A special evaluation is `P_n^{(a,b)}(1) = {n+a \choose n}`::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print jacobi(4, 0.5, 0.25, 1)
        2.4609375
        >>> print binomial(4+0.5, 4)
        2.4609375

    A Jacobi polynomial of degree `n` is equal to its
    Taylor polynomial of degree `n`. The explicit
    coefficients of Jacobi polynomials can therefore
    be recovered easily using :func:`taylor`::

        >>> for n in range(5):
        ...     nprint(taylor(lambda x: jacobi(n,1,2,x), 0, n))
        ...
        [1.0]
        [-0.5, 2.5]
        [-0.75, -1.5, 5.25]
        [0.5, -3.5, -3.5, 10.5]
        [0.625, 2.5, -11.25, -7.5, 20.625]

    For nonintegral `n`, the Jacobi "polynomial" is no longer
    a polynomial::

        >>> nprint(taylor(lambda x: jacobi(0.5,1,2,x), 0, 4))
        [0.309983, 1.84119, -1.26933, 1.26699, -1.34808]

    **Orthogonality**

    The Jacobi polynomials are orthogonal on the interval
    `[-1, 1]` with respect to the weight function
    `w(x) = (1-x)^a (1+x)^b`. That is,
    `w(x) P_n^{(a,b)}(x) P_m^{(a,b)}(x)` integrates to
    zero if `m \ne n` and to a nonzero number if `m = n`.

    The orthogonality is easy to verify using numerical
    quadrature::

        >>> P = jacobi
        >>> f = lambda x: (1-x)**a * (1+x)**b * P(m,a,b,x) * P(n,a,b,x)
        >>> a = 2
        >>> b = 3
        >>> m, n = 3, 4
        >>> nprint(quad(f, [-1, 1]), 1)
        4.0e-23
        >>> m, n = 4, 4
        >>> print quad(f, [-1, 1])
        1.9047619047619

    **Differential equation**

    The Jacobi polynomials are solutions of the differential
    equation

    .. math ::

      (1-x^2) y'' + (b-a-(a+b+2)x) y' + n (n+a+b+1) y = 0.

    We can verify that :func:`jacobi` approximately satisfies
    this equation::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> a = 2.5
        >>> b = 4
        >>> n = 3
        >>> y = lambda x: jacobi(n,a,b,x)
        >>> x = pi
        >>> A0 = n*(n+a+b+1)*y(x)
        >>> A1 = (b-a-(a+b+2)*x)*diff(y,x)
        >>> A2 = (1-x**2)*diff(y,x,2)
        >>> nprint(A2 + A1 + A0, 1)
        4.0e-12

    The difference of order `10^{-12}` is as close to zero as
    it could be at 15-digit working precision, since the terms
    are large::

        >>> print A0, A1, A2
        26560.2328981879 -21503.7641037294 -5056.46879445852

    """
    return binomial(n+a,n) * hyp2f1(-n,1+n+a+b,a+1,(1-x)/2)

@funcwrapper
def legendre(n, x):
    r"""
    ``legendre(n, x)`` evaluates the Legendre polynomial `P_n(x)`.
    The Legendre polynomials are given by the formula

    .. math ::

        P_n(x) = \frac{1}{2^n n!} \frac{d^n}{dx^n} (x^2 -1)^n.

    Alternatively, they can be computed recursively using

    .. math ::

        P_0(x) = 1

        P_1(x) = x

        (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x).

    A third definition is in terms of the hypergeometric function
    `\,_2F_1`, whereby they can be generalized to arbitrary `n`:

    .. math ::

        P_n(x) = \,_2F_1\left(-n, n+1, 1, \frac{1-x}{2}\right)

    **Basic evaluation**

    The Legendre polynomials assume fixed values at the points
    `x = -1` and `x = 1`::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> nprint([legendre(n, 1) for n in range(6)])
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        >>> nprint([legendre(n, -1) for n in range(6)])
        [1.0, -1.0, 1.0, -1.0, 1.0, -1.0]

    The coefficients of Legendre polynomials can be recovered
    using degree-`n` Taylor expansion::

        >>> for n in range(5):
        ...     nprint(taylor(lambda x: legendre(n, x), 0, n))
        ...
        [1.0]
        [0.0, 1.0]
        [-0.5, 0.0, 1.5]
        [0.0, -1.5, 0.0, 2.5]
        [0.375, 0.0, -3.75, 0.0, 4.375]

    The roots of Legendre polynomials are located symmetrically
    on the interval `[-1, 1]`::

        >>> for n in range(5):
        ...     nprint(polyroots(taylor(lambda x: legendre(n, x), 0, n)[::-1]))
        ...
        []
        [0.0]
        [-0.57735, 0.57735]
        [-0.774597, 0.0, 0.774597]
        [-0.861136, -0.339981, 0.339981, 0.861136]

    An example of an evaluation for arbitrary `n`::

        >>> print legendre(0.75, 2+4j)
        (1.94952805264875 + 2.1071073099422j)

    **Orthogonality**

    The Legendre polynomials are orthogonal on `[-1, 1]` with respect
    to the trivial weight `w(x) = 1`. That is, `P_m(x) P_n(x)`
    integrates to zero if `m \ne n` and to `2/(2n+1)` if `m = n`::

        >>> m, n = 3, 4
        >>> print quad(lambda x: legendre(m,x)*legendre(n,x), [-1, 1])
        0.0
        >>> m, n = 4, 4
        >>> print quad(lambda x: legendre(m,x)*legendre(n,x), [-1, 1])
        0.222222222222222

    **Differential equation**

    The Legendre polynomials satisfy the differential equation

    .. math ::

        ((1-x^2) y')' + n(n+1) y' = 0.

    We can verify this numerically::

        >>> n = 3.6
        >>> x = 0.73
        >>> P = legendre
        >>> A = diff(lambda t: (1-t**2)*diff(lambda u: P(n,u), t), x)
        >>> B = n*(n+1)*P(n,x)
        >>> nprint(A+B,1)
        9.0e-16

    """
    if isint(n):
        n = int(n)
    if x == -1:
        # TODO: hyp2f1 should handle this
        if n == int(n):
            return (-1)**(n + (n>=0)) * mpf(-1)
        return inf
    return hyp2f1(-n,n+1,1,(1-x)/2)

@funcwrapper
def chebyt(n, x):
    r"""
    ``chebyt(n, x)`` evaluates the Chebyshev polynomial of the first
    kind `T_n(x)`, defined by the identity

    .. math ::

        T_n(\cos x) = \cos(n x).

    The Chebyshev polynomials of the first kind are a special
    case of the Jacobi polynomials, and by extension of the
    hypergeometric function `\,_2F_1`. They can thus also be
    evaluated for nonintegral `n`.

    **Basic evaluation**

    The coefficients of the `n`-th polynoimal can be recovered
    using using degree-`n` Taylor expansion::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> for n in range(5):
        ...     nprint(taylor(lambda x: chebyt(n, x), 0, n))
        ...
        [1.0]
        [0.0, 1.0]
        [-1.0, 0.0, 2.0]
        [0.0, -3.0, 0.0, 4.0]
        [1.0, 0.0, -8.0, 0.0, 8.0]

    **Orthogonality**

    The Chebyshev polynomials of the first kind are orthogonal
    on the interval `[-1, 1]` with respect to the weight
    function `w(x) = 1/\sqrt{1-x^2}`::

        >>> f = lambda x: chebyt(m,x)*chebyt(n,x)/sqrt(1-x**2)
        >>> m, n = 3, 4
        >>> nprint(quad(f, [-1, 1]),1)
        -7.0e-28
        >>> m, n = 4, 4
        >>> print quad(f, [-1, 1])
        1.57079632596448

    """
    return hyp2f1(-n,n,0.5,(1-x)/2)

@funcwrapper
def chebyu(n, x):
    r"""
    ``chebyu(n, x)`` evaluates the Chebyshev polynomial of the second
    kind `U_n(x)`, defined by the identity

    .. math ::

        U_n(\cos x) = \frac{\sin((n+1)x)}{\sin(x)}.

    The Chebyshev polynomials of the second kind are a special
    case of the Jacobi polynomials, and by extension of the
    hypergeometric function `\,_2F_1`. They can thus also be
    evaluated for nonintegral `n`.

    **Basic evaluation**

    The coefficients of the `n`-th polynoimal can be recovered
    using using degree-`n` Taylor expansion::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> for n in range(5):
        ...     nprint(taylor(lambda x: chebyu(n, x), 0, n))
        ...
        [1.0]
        [0.0, 2.0]
        [-1.0, 0.0, 4.0]
        [0.0, -4.0, 0.0, 8.0]
        [1.0, 0.0, -12.0, 0.0, 16.0]

    **Orthogonality**

    The Chebyshev polynomials of the second kind are orthogonal
    on the interval `[-1, 1]` with respect to the weight
    function `w(x) = \sqrt{1-x^2}`::

        >>> f = lambda x: chebyu(m,x)*chebyu(n,x)*sqrt(1-x**2)
        >>> m, n = 3, 4
        >>> print quad(f, [-1, 1])
        0.0
        >>> m, n = 4, 4
        >>> print quad(f, [-1, 1])
        1.5707963267949

    """
    return (n+1) * hyp2f1(-n, n+2, 1.5, (1-x)/2)

@funcwrapper
def jv(v, x):
    r"""
    ``jv(n,x)`` computes the Bessel function of the first kind
    `J_n(x)`. Bessel functions of the first kind are defined as
    solutions of the differential equation

    .. math ::

        x^2 y'' + x y' + (x^2 - n^2) y = 0

    which is Laplace's equation in cylindrical coordinates. This
    equation has two solutions for given `n`, where the
    `J_n`-function is the solution that is nonsingular at `x = 0`.
    For positive integer `n`, `J_n(x)` behaves roughly like
    a sine (odd `n`) or cosine (even `n`) multiplied by
    a magnitude factor that decays slowly as `x \to \pm\infty`.

    Generally, `J_n` is a special case of the hypergeometric
    function `\,_0F_1`:

    .. math ::

        J_n(x) = \frac{x^n}{2^n \Gamma(n+1)}
                 \,_0F_1\left(n+1,-\frac{x^2}{4}\right)

    **Examples**

    Evaluation is supported for arbitrary arguments, and at
    arbitrary precision::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print jv(2, 1000)
        -0.024777229528606
        >>> print jv(4, 0.75)
        0.000801070086542314
        >>> print jv(2, 1000j)
        (-2.48071721019185e+432 + 0.0j)
        >>> mp.dps = 25
        >>> print jv(0.75j, 3+4j)
        (-2.778118364828153309919653 - 1.5863603889018621585533j)
        >>> mp.dps = 50
        >>> print jv(1, pi)
        0.28461534317975275734531059968613140570981118184947

    The Bessel functions of the first kind satisfy simple
    symmetries around `x = 0`::

        >>> mp.dps = 15
        >>> nprint([jv(n,0) for n in range(5)])
        [1.0, 0.0, 0.0, 0.0, 0.0]
        >>> nprint([jv(n,pi) for n in range(5)])
        [-0.304242, 0.284615, 0.485434, 0.333458, 0.151425]
        >>> nprint([jv(n,-pi) for n in range(5)])
        [-0.304242, -0.284615, 0.485434, -0.333458, 0.151425]

    Roots of Bessel functions are often used::

        >>> nprint([findroot(j0, k) for k in [2, 5, 8, 11, 14]])
        [2.40483, 5.52008, 8.65373, 11.7915, 14.9309]
        >>> nprint([findroot(j1, k) for k in [3, 7, 10, 13, 16]])
        [3.83171, 7.01559, 10.1735, 13.3237, 16.4706]

    The roots are not periodic, but the distance between successive
    roots asymptotically approaches `2 \pi`. Bessel functions of
    the first kind have the following normalization::

        >>> print quadosc(j0, [0, inf], period=2*pi)
        1.0
        >>> print quadosc(j1, [0, inf], period=2*pi)
        1.0

    For `n = 1/2` or `n = -1/2`, the Bessel function reduces to a
    trigonometric function::

        >>> x = 10
        >>> print jv(0.5, x), sqrt(2/(pi*x))*sin(x)
        -0.13726373575505 -0.13726373575505
        >>> print jv(-0.5, x), sqrt(2/(pi*x))*cos(x)
        -0.211708866331398 -0.211708866331398

    """
    if isint(v):
        v = int(v)
        if isinstance(x, mpf):
            return make_mpf(libhyper.mpf_besseljn(v, x._mpf_, mp.prec))
        if isinstance(x, mpc):
            return make_mpc(libhyper.mpc_besseljn(v, x._mpc_, mp.prec))
    hx = x/2
    return hx**v * hyp0f1(v+1, -hx**2) / factorial(v)

jn = jv

def j0(x):
    """Computes the Bessel function `J_0(x)`. See :func:`jv`."""
    return jv(0, x)

def j1(x):
    """Computes the Bessel function `J_1(x)`.  See :func:`jv`."""
    return jv(1, x)

@funcwrapper
def lambertw(z, k=0, approx=None):
    r"""
    The Lambert W function `W(z)` is defined as the inverse function
    of `w \exp(w)`. In other words, the value of `W(z)` is such that
    `z = W(z) \exp(W(z))` for any complex number `z`.

    The Lambert W function is a multivalued function with infinitely
    many branches. Each branch gives a separate solution of the
    equation `w \exp(w)`. All branches are supported by
    :func:`lambertw`:

    * ``lambertw(z)`` gives the principal solution (branch 0)

    * ``lambertw(z, k)`` gives the solution on branch `k`

    The Lambert W function has two partially real branches: the
    principal branch (`k = 0`) is real for real `z > -1/e`, and the
    `k = -1` branch is real for `-1/e < z < 0`. All branches except
    `k = 0` have a logarithmic singularity at `z = 0`.

    **Basic examples**

    The Lambert W equation is the inverse of `w \exp(w)`::

        >>> from mpmath import *
        >>> mp.dps = 35
        >>> w = lambertw(1)
        >>> print w
        0.56714329040978387299996866221035555
        >>> print w*exp(w)
        1.0

    Any branch gives a valid inverse::

        >>> w = lambertw(1, k=3)
        >>> print w    # doctest: +NORMALIZE_WHITESPACE
        (-2.8535817554090378072068187234910812 + 
          17.113535539412145912607826671159289j)
        >>> print w*exp(w)
        (1.0 + 3.5075477124212226194278700785075126e-36j)

    **Applications to equation-solving**

    The Lambert W function can give the value of the infinite power
    tower `z^{z^{z^{\ldots}}}`::

        >>> def tower(z, n):
        ...     if n == 0:
        ...         return z
        ...     return z ** tower(z, n-1)
        ...
        >>> tower(0.5, 100)
        0.641185744504986
        >>> mp.dps = 50
        >>> print -lambertw(-log(0.5))/log(0.5)
        0.6411857445049859844862004821148236665628209571911

    **Properties**

    The Lambert W function grows roughly like the natural logarithm
    for large arguments::

        >>> mp.dps = 15
        >>> print lambertw(1000)
        5.2496028524016
        >>> print log(1000)
        6.90775527898214
        >>> print lambertw(10**100)
        224.843106445119
        >>> print log(10**100)
        230.258509299405

    The principal branch of the Lambert W function has a rational
    Taylor series expansion around `z = 0`::

        >>> nprint(taylor(lambertw, 0, 6), 10)
        [0.0, 1.0, -1.0, 1.5, -2.666666667, 5.208333333, -10.8]

    Some special values and limits are::

        >>> mp.dps = 15
        >>> print lambertw(0)
        0.0
        >>> print lambertw(1)
        0.567143290409784
        >>> print lambertw(e)
        1.0
        >>> print lambertw(inf)
        +inf
        >>> print lambertw(0, k=-1)
        -inf
        >>> print lambertw(0, k=3)
        -inf
        >>> print lambertw(inf, k=3)
        +inf

    The `k = 0` and `k = -1` branches join at `z = -1/e` where
    `W(z) = -1` for both branches. Since `-1/e` can only be represented
    approximately with mpmath numbers, evaluating the Lambert W function
    at this point only gives `-1` approximately::

        >>> mp.dps = 25
        >>> print lambertw(-1/e, 0)
        -0.999999999999837133022867
        >>> print lambertw(-1/e, -1)
        -1.00000000000016286697718

    If `-1/e` happens to round in the negative direction, there might be
    a small imaginary part::

        >>> mp.dps = 15
        >>> print lambertw(-1/e)
        (-1.0 + 8.22007971511612e-9j)

    **Possible issues**

    The evaluation can become inaccurate very close to the branch point
    at `-1/e`. In some corner cases, :func:`lambertw` might currently
    fail to converge, or can end up on the wrong branch.

    **Algorithm**

    Halley's iteration is used to invert `w \exp(w)`, using a first-order
    asymptotic approximation (`O(\log(w))` or `O(w)`) as the initial
    estimate.

    The definition, implementation and choice of branches is based
    on Corless et al, "On the Lambert W function", Adv. Comp. Math. 5
    (1996) 329-359, available online here:
    http://www.apmaths.uwo.ca/~djeffrey/Offprints/W-adv-cm.pdf

    TODO: use a series expansion when extremely close to the branch point
    at `-1/e` and make sure that the proper branch is chosen there
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
    elif k == 0 and z.imag and abs(z) <= 0.6:
        w = z
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

if __name__ == '__main__':
    import doctest
    doctest.testmod()

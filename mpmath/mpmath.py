from lib import *

class mpnumeric(object):
    """Base class for mpf and mpc. Calling mpnumeric(x) returns an mpf
    if x can be converted to an mpf (if it is a float, int, mpf, ...),
    and an mpc if x is complex."""

    def __new__(cls, val):
        # TODO: should maybe normalize here
        if isinstance(val, cls):
            return val
        if isinstance(val, complex):
            return mpc(val)
        return mpf(val)

def convert_lossless(x):
    """Attempt to convert x to an mpf or mpc losslessly. If x is an
    mpf or mpc, return it unchanged. If x is an int, create an mpf with
    sufficient precision to represent it exactly.

    If x is a str, just convert it to an mpf with the current working
    precision (perhaps this should be done differently...)"""
    if isinstance(x, mpnumeric):
        return x
    if isinstance(x, float):
        return make_mpf(from_float(x, 53, round_floor))
    if isinstance(x, int_types):
        return make_mpf(from_int(x, bitcount(x), round_floor))
    if isinstance(x, complex):
        return mpc(x)
    if isinstance(x, basestring):
        return make_mpf(from_str(x, mpf._prec, mpf._rounding))
    raise TypeError("cannot create mpf from " + repr(x))


def mpf_convert_operand(x):
    if isinstance(x, int_types):
        return make_mpf(from_int(x, bitcount(x), round_floor))
    if isinstance(x, float):
        return make_mpf(from_float(x, 53, round_floor))
    return NotImplemented

def mpf_convert_lhs(x):
    if isinstance(x, complex_types):
        return mpc(x)
    if isinstance(x, int_types):
        return make_mpf(from_int(x, bitcount(x), round_floor))
    if isinstance(x, float):
        return make_mpf(from_float(x, 53, round_floor))
    return NotImplemented



# Global settings
g_prec = 53
g_rounding = round_half_even
g_dps = 15

class context(type):

    def set_rounding(cls, val):
        global g_rounding; g_rounding = val

    _rounding = property(lambda cls: g_rounding, set_rounding)

    def set_prec(self, n):
        global g_prec, g_dps
        g_prec = max(1, int(n))
        g_dps = max(1, int(round(int(n)/LOG2_10)-1))

    prec = property(lambda cls: g_prec, set_prec)
    _prec = prec

    def set_dps(self, n):
        global g_prec, g_dps
        g_prec = max(1, int(round((int(n)+1)*LOG2_10)))
        g_dps = max(1, int(n))

    dps = property(lambda cls: g_dps, set_dps)

    def round_up(cls): cls._rounding = round_up
    def round_down(cls): cls._rounding = round_down
    def round_floor(cls): cls._rounding = round_floor
    def round_ceiling(cls): cls._rounding = round_ceiling
    def round_half_down(cls): cls._rounding = round_half_down
    def round_half_up(cls): cls._rounding = round_half_up
    def round_half_even(cls): cls._rounding = round_half_even
    def round_default(cls): cls._rounding = round_half_even


int_types = (int, long)


def _convert(x):
    """Convet x to mpf data"""
    if isinstance(x, float):
        return from_float(x, mpf._prec, mpf._rounding)
    if isinstance(x, int_types):
        return from_int(x, mpf._prec, mpf._rounding)
    if isinstance(x, basestr):
        return from_str(x, mpf._prec, mpf._rounding)
    raise TypeError("cannot create mpf from " + repr(x))


class mpf(mpnumeric):
    """An mpf instance holds a real-valued floating-point number. mpf:s
    work analogously to Python floats, but support arbitrary-precision
    arithmetic. The mpf class has two properties 'dps' and 'prec' which
    respectively hold the working precision as measured in decimal
    digits and in bits. (The default is 15 digits / 53 bits, the same
    as Python floats.) For example, this increases precision by 10
    bits:

        mpf.prec += 10

    The global working precision controls the precision at which all
    arithmetic operations on mpf:s is carried out. Directed rounding is
    also (partially) implemented; all calculations will be rounded up
    after calling

        mpf.round_up()

    mpf.round_half_even() is the default rounding.

    An mpf is represented internally as a tuple of integers (man, exp,
    bc) where man is the mantissa, exp is the exponent and bc is the
    number of bits in the mantissa (bc <= mpf.prec if the number is
    normalized to the current working precision). Mathematically, that
    means the mpf x has the value x = man * 2**exp. The components can
    be accessed using the .val property of an mpf.

    A useful difference between mpf:s and Python floats is that
    operations on mpf:s that mathematically produce complex numbers
    (like mpf(-1)**0.5) return mpc:s instead of raising exceptions.
    """

    __metaclass__ = context

    special = ''

    def __new__(cls, val=fzero):
        """A new mpf can be created from a Python float, an int, a
        or a decimal string representing a number in floating-point
        format. Examples:

            mpf(25)
            mpf(2.5)
            mpf('2.5')
            mpf('1.6e1000')

        An mpf can also be created from a tuple (man, exp) or
        (man, exp, bc):

            >>> mpf((3, -1))
            1.5
        """
        if isinstance(val, mpf):
            man, exp, bc = val.val
            if bc == -1:
                return val
            return make_mpf(normalize(man, exp, bc, g_prec, g_rounding))
        elif isinstance(val, tuple):
            return make_mpf(from_man_exp(val[0], val[1], g_prec, g_rounding))
        return +convert_lossless(val)

    man_exp = property(lambda self: self.val[0:2])
    man = property(lambda self: self.val[0])
    exp = property(lambda self: self.val[1])
    bc = property(lambda self: self.val[2])

    def __repr__(s): return "mpf('%s')" % to_str(s.val, g_dps+2)
    def __str__(s): return to_str(s.val, g_dps)
    def __hash__(s): return fhash(s.val)
    def __int__(s): return to_int(s.val)
    def __float__(s): return to_float(s.val)
    def __complex__(s): return float(s) + 0j
    def __nonzero__(s): return bool(s.man)

    def __abs__(s): return make_mpf(fabs(s.val, g_prec, g_rounding))
    def __pos__(s): return make_mpf(fpos(s.val, g_prec, g_rounding))
    def __neg__(s): return make_mpf(fneg(s.val, g_prec, g_rounding))

    def __eq__(s, t):
        if isinstance(t, mpf):
            return feq(s.val, t.val)
        if isinstance(t, complex_types):
            return mpc(s) == t
        t = mpf_convert_operand(t)
        if t is NotImplemented:
            return False
        return feq(s.val, t.val)

    def __ne__(s, t):
        return not s.__eq__(t)

    def __cmp__(s, t):
        if not isinstance(t, mpf):
            t = mpf_convert_operand(t)
            if t is NotImplemented:
                return t
        return fcmp(s.val, t.val)

    def binop(s, t, f):
        if isinstance(t, complex_types):
            s = mpc(s)
            if f is feq: return s == t
            if f is fmul: return s * t
            if f is fadd: return s + t
            if f is fsub: return s - t
            if f is fdiv: return s / t
            raise ValueError("bad operation")
        t = mpf_convert_operand(t)
        if t is NotImplemented:
            return t
        return make_mpf(f(s.val, t.val, g_prec, g_rounding))

    def __add__(s, t):
        if isinstance(t, mpf):
            return make_mpf(fadd(s.val, t.val, g_prec, g_rounding))
        return s.binop(t, fadd)

    def __sub__(s, t):
        if isinstance(t, mpf):
            return make_mpf(fsub(s.val, t.val, g_prec, g_rounding))
        return s.binop(t, fsub)

    def __mul__(s, t):
        if isinstance(t, mpf):
            return make_mpf(fmul(s.val, t.val, g_prec, g_rounding))
        return s.binop(t, fmul)

    def __div__(s, t):
        if isinstance(t, mpf):
            return make_mpf(fdiv(s.val, t.val, g_prec, g_rounding))
        return s.binop(t, fdiv)

    def __mod__(s, t):
        if isinstance(t, mpf):
            return make_mpf(fmod(s.val, t.val, g_prec, g_rounding))
        return s.binop(t, fmod)

    def __pow__(s, t):
        if isinstance(t, int_types):
            return make_mpf(fpow(s.val, t, g_prec, g_rounding))
        if not isinstance(t, mpf):
            if isinstance(t, complex_types):
                return power(s, t)
            t = mpf_convert_operand(t)
            if t is NotImplemented:
                return t
        if t.val == fhalf:
            return sqrt(s)
        man, exp, bc = t.val
        if exp >= 0:
            return make_mpf(fpow(s.val, man<<exp, g_prec, g_rounding))
        return power(s, t)

    __radd__ = __add__
    __rmul__ = __mul__
    def __rsub__(s, t): return mpf_convert_lhs(t) - s
    def __rdiv__(s, t): return mpf_convert_lhs(t) / s
    def __rpow__(s, t): return mpf_convert_lhs(t) ** s
    def __rmod__(s, t): return mpf_convert_lhs(t) % s

    def sqrt(s):
        return sqrt(s)

    def ae(s, t, rel_eps=None, abs_eps=None):
        """
        Determine whether the difference between s and t is smaller
        than a given epsilon ("ae" is short for "almost equal").

        Both a maximum relative difference and a maximum difference
        ('epsilons') may be specified. The absolute difference is
        defined as |s-t| and the relative difference is defined
        as |s-t|/max(|s|, |t|).

        If only one epsilon is given, both are set to the same value.
        If none is given, both epsilons are set to 2**(-prec+m) where
        prec is the current working precision and m is a small integer.
        """
        if not isinstance(t, mpf):
            t = mpf(t)
        if abs_eps is None and rel_eps is None:
            rel_eps = abs_eps = make_mpf((1, -mpf._prec+4, 1))
        if abs_eps is None:
            abs_eps = rel_eps
        elif rel_eps is None:
            rel_eps = abs_eps
        diff = abs(s-t)
        if diff <= abs_eps:
            return True
        abss = abs(s)
        abst = abs(t)
        if abss < abst:
            err = diff/abst
        else:
            err = diff/abss
        return err <= rel_eps

    def almost_zero(s, prec):
        """Quick check if |s| < 2**-prec. May return a false negative
        if s is very close to the threshold."""
        return s.bc + s.exp < prec


def make_mpf(tpl, construct=object.__new__, cls=mpf):
    """Create mpf verbatim from a given tuple of data."""
    a = construct(cls)
    a.val = tpl
    return a


inf = make_mpf(finf)
ninf = make_mpf(fninf)
nan = make_mpf(fnan)

def isnan(x):
    if not isinstance(x, mpf):
        return False
    return x.val == fnan

class mpc(mpnumeric):
    """An mpc represents a complex number using a pair of mpf:s (one
    for the real part and another for the imaginary part.) The mpc
    class behaves fairly similarly to Python's complex type."""

    def __new__(cls, real=0, imag=0):
        s = object.__new__(cls)
        if isinstance(real, complex_types):
            real, imag = real.real, real.imag
        s.real = mpf(real)
        s.imag = mpf(imag)
        return s

    def __repr__(s):
        r = repr(s.real)[4:-1]
        i = repr(s.imag)[4:-1]
        return "mpc(real=%s, imag=%s)" % (r, i)

    def __str__(s):
        return "(%s + %sj)" % (s.real, s.imag)

    def __complex__(s):
        return complex(float(s.real), float(s.imag))

    def __pos__(s):
        return mpc(s.real, s.imag)

    def __abs__(s):
        return hypot(s.real, s.imag)

    def __eq__(s, t):
        if not isinstance(t, mpc):
            if isinstance(t, str):
                return False
            t = mpc(t)
        return s.real == t.real and s.imag == t.imag

    def _compare(*args):
        raise TypeError("no ordering relation is defined for complex numbers")

    __gt__ = _compare
    __le__ = _compare
    __gt__ = _compare
    __ge__ = _compare

    def __nonzero__(s):
        return bool(s.real) or bool(s.imag)

    def conjugate(s):
        return mpc(s.real, -s.imag)

    def __add__(s, t):
        if not isinstance(t, mpc):
            t = mpc(t)
        return mpc(s.real+t.real, s.imag+t.imag)

    __radd__ = __add__

    def __neg__(s):
        return mpc(-s.real, -s.imag)

    def __sub__(s, t):
        if not isinstance(t, mpc):
            t = mpc(t)
        return mpc(s.real-t.real, s.imag-t.imag)

    def __rsub__(s, t):
        return (-s) + t

    def __mul__(s, t):
        if not isinstance(t, mpc):
            t = mpc(t)
        return mpc(*fcmul(s.real.val, s.imag.val, t.real.val, t.imag.val,
            g_prec, g_rounding))

    __rmul__ = __mul__

    def __div__(s, t):
        if not isinstance(t, mpc):
            t = mpc(t)
        a = s.real; b = s.imag; c = t.real; d = t.imag
        mag = c*c + d*d
        return mpc((a*c+b*d)/mag, (b*c-a*d)/mag)

    def __rdiv__(s, t):
        return mpc(t) / s

    def __pow__(s, n):
        if n == 0: return mpc(1)
        if n == 1: return +s
        if n == -1: return 1/s
        if n == 2: return s*s
        if isinstance(n, (int, long)) and n > 0:
            # TODO: should increase working precision here
            w = mpc(1)
            while n:
                if n & 1:
                    w = w*s
                    n -= 1
                s = s*s
                n //= 2
            return w
        if n == 0.5:
            return sqrt(s)
        return power(s, n)

    def __rpow__(s, t):
        return convert_lossless(t) ** s

    # TODO: refactor and merge with mpf.ae
    def ae(s, t, rel_eps=None, abs_eps=None):
        if not isinstance(t, mpc):
            t = mpc(t)
        if abs_eps is None and rel_eps is None:
            abs_eps = rel_eps = make_mpf((1, -mpf._prec+4, 1))
        if abs_eps is None:
            abs_eps = rel_eps
        elif rel_eps is None:
            rel_eps = abs_eps
        diff = abs(s-t)
        if diff <= abs_eps:
            return True
        abss = abs(s)
        abst = abs(t)
        if abss < abst:
            err = diff/abst
        else:
            err = diff/abss
        return err <= rel_eps


complex_types = (complex, mpc)

def make_mpc(tpl, construct=object.__new__, cls=mpc):
    a = construct(cls)
    a.real, a.imag = map(make_mpf, tpl)
    return a

j = mpc(0,1)


class constant(mpf):
    """Represents a mathematical constant with dynamic precision.
    When printed or used in an arithmetic operation, a constant
    is converted to a regular mpf at the working precision. A
    regular mpf can also be obtained using the operation +x."""

    def __new__(cls, func, name):
        a = object.__new__(cls)
        a.name = name
        a.func = func
        return a

    @property
    def val(self):
        return self.func(g_prec, g_rounding)

    #def __repr__(self):
    #    return "<%s: %s~>" % (self.name, mpf.__str__(self))


_180 = from_int(180, 10, round_floor)

pi = constant(fpi, "pi")
degree = constant(lambda p, r: fdiv(fpi(p+4, round_floor), _180, p, r), "degree")
e = constant(lambda p, r: fexp(fone, p, r), "e")
euler = constant(fgamma, "Euler's constant gamma")
clog2 = constant(flog2, "log(2)")
clog10 = constant(flog10, "log(10)")


def sqrt(x):
    """For real x >= 0, returns the square root of x. For negative or
    complex x, returns the principal branch of the complex square root
    of x."""
    x = convert_lossless(x)
    if isinstance(x, mpf) and x.val[0] >= 0:
        return make_mpf(fsqrt(x.val, g_prec, g_rounding))
    x = mpc(x)
    return make_mpc(fcsqrt(x.real.val, x.imag.val, g_prec, g_rounding))

def hypot(x, y):
    """Returns the Euclidean distance sqrt(x*x + y*y). Both x and y
    must be real."""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return make_mpf(fhypot(x.val, y.val, g_prec, g_rounding))

def floor(x):
    """Computes the floor function of x. Note: returns an mpf, not a
    Python int. If x is larger than the precision, it will be rounded,
    not necessarily in the floor direction."""
    x = convert_lossless(x)
    return make_mpf(ffloor(x.val, g_prec, g_rounding))

def ceil(x):
    """Computes the ceiling function of x. Note: returns an mpf, not a
    Python int. If x is larger than the precision, it will be rounded,
    not necessarily in the ceiling direction."""
    x = convert_lossless(x)
    return make_mpf(fceil(x.val, g_prec, g_rounding))

# Since E-functions simply map reals to reals and complexes to complexes, we
# can construct all of them the same way (unlike log, sqrt, etc)
def ef(name, real_f, complex_f, doc):
    def f(x):
        x = convert_lossless(x)
        if isinstance(x, mpf):
            return make_mpf(real_f(x.val, g_prec, g_rounding))
        else:
            return make_mpc(complex_f(x.real.val, x.imag.val, g_prec, g_rounding))
    f.__name__ = name
    f.__doc__ = doc
    return f

exp = ef('exp', fexp, fcexp, "Returns the exponential function of x.")
cos = ef('cos', fcos, fccos, "Returns the cosine of x.")
sin = ef('sin', fsin, fcsin, "Returns the sine of x.")
cosh = ef('cosh', fcosh, fccosh, "Returns the hyperbolic cosine of x.")
sinh = ef('sinh', fsinh, fcsinh, "Returns the hyperbolic sine of x.")

# TODO: implement tanh and complex tan in lib instead
def tan(x):
    """Returns the tangent of x."""
    x = convert_lossless(x)
    if isinstance(x, mpf):
        return make_mpf(ftan(x.val, g_prec, g_rounding))
    # the complex division can cause enormous cancellation.
    # TODO: handle more robustly
    mpf._prec += 20
    t = sin(x) / cos(x)
    mpf._prec -= 20
    return +t

def tanh(x):
    """Returns the hyperbolic tangent of x."""
    x = convert_lossless(x)
    oldprec = mpf._prec
    a = abs(x)
    mpf._prec += 10
    high = a.exp + a.bc
    if high < -10:
        if high < (-(mpf._prec-10) * 0.3):
            return x - (x**3)/3 + 2*(x**5)/15
        mpf._prec += (-high)
    a = exp(2*x)
    t = (a-1)/(a+1)
    mpf._prec = oldprec
    return +t

def arg(x):
    """Returns the complex argument (phase) of x. The returned value is
    an mpf instance. The argument is here defined to satisfy
    -pi < arg(x) <= pi. On the negative real half-axis, it is taken to
    be +pi."""
    x = mpc(x)
    mpf._prec += 5
    t = atan2(x.imag, x.real)
    mpf._prec -= 5
    return +t

def log(x, b=None):
    """Returns the base-b logarithm of x. If b is unspecified, return
    the natural (base-e) logarithm. log(x, b) is defined as
    log(x)/log(b). log(0) raises ValueError.

    The natural logarithm is real if x > 0 and complex if x < 0 or if x
    is complex. The principal branch of the complex logarithm is chosen,
    for which Im(log(x)) = -pi < arg(x) <= pi. """
    if b is not None:
        mpf.prec += 3
        a = log(x) / log(b)
        mpf.prec -= 3
        return +a
    x = convert_lossless(x)
    if not x:
        raise ValueError, "logarithm of 0"
    if isinstance(x, mpf) and x.val[0] > 0:
        return make_mpf(flog(x.val, g_prec, g_rounding))
    else:
        return mpc(log(abs(x)), arg(x))

def power(x, y):
    """Returns x**y = exp(y*log(x)) for real or complex x and y."""
    # TODO: better estimate for extra precision needed
    mpf._prec += 10
    t = exp(y * log(x))
    mpf._prec -= 10
    return +t

def atan(x):
    """Returns the inverse tangent of x."""
    x = convert_lossless(x)
    if isinstance(x, mpf):
        return make_mpf(fatan(x.val, g_prec, g_rounding))
    # TODO: maybe get this to agree with Python's cmath atan about the
    # branch to choose on the imaginary axis
    # TODO: handle cancellation robustly
    mpf._prec += 10
    t = (0.5j)*(log(1-1j*x) - log(1+1j*x))
    mpf._prec -= 10
    return +t

def atan2(y,x):
    """atan2(y, x) has the same magnitude as atan(y/x) but accounts for
    the signs of y and x. (Defined for real x and y only.)"""
    x = mpf(x)
    y = mpf(y)
    if y < 0:
        return -atan2(-y, x)
    if not x and not y:
        return mpf(0)
    if y > 0 and x == 0:
        mpf._prec += 2
        t = pi/2
        mpf._prec -= 2
        return t
    mpf._prec += 2
    if x > 0:
        a = atan(y/x)
    else:
        a = pi - atan(-y/x)
    mpf._prec -= 2
    return +a

# TODO: robustly deal with cancellation in all of the following functions

def _asin_complex(z):
    mpf._prec += 10
    t = -1j * log(1j * z + sqrt(1 - z*z))
    mpf._prec -= 10
    return +t

def asin(x):
    """Returns the inverse sine of x. Outside the range [-1, 1], the
    result is complex and defined as the principal branch value of
    -i * log(i * x + sqrt(1 - x**2))."""
    x = convert_lossless(x)
    if isinstance(x, mpf) and abs(x) <= 1:
        return _asin_complex(x).real
    return _asin_complex(x)

def _acos_complex(z):
    mpf._prec += 10
    t = pi/2 + 1j * log(1j * z + sqrt(1 - z*z))
    mpf._prec -= 10
    return +t

def acos(x):
    """Returns the inverse cosine of x. Outside the range [-1, 1], the
    result is complex and defined as the principal branch value of
    pi/2 + i * log(i * x + sqrt(1 - x**2))."""
    x = convert_lossless(x)
    if isinstance(x, mpf) and abs(x) <= 1:
        return _acos_complex(x).real
    return _acos_complex(x)

def asinh(x):
    """Returns the inverse hyperbolic sine of x. For complex x, the
    result is the principal branch value of log(x + sqrt(1 + x**2))."""
    x = convert_lossless(x)
    oldprec = mpf._prec
    a = abs(x)
    mpf._prec += 10
    high = a.exp + a.bc
    if high < -10:
        if high < (-(mpf._prec-10) * 0.3):
            return x - (x**3)/6 + 3*(x**5)/40
        mpf._prec += (-high)
    t = log(x + sqrt(x**2 + 1))
    mpf._prec = oldprec
    return +t

def acosh(x):
    """Returns the inverse hyperbolic cosine of x. The value is
    given by log(x + sqrt(1 + x**2)), where the principal branch is
    used when the result is complex."""
    x = convert_lossless(x)
    mpf._prec += 10
    t = log(x + sqrt(x-1)*sqrt(x+1))
    mpf._prec -= 10
    return +t

def atanh(x):
    """Returns the inverse hyperbolic tangent of x. Outside the range
    [-1, 1], the result is complex and defined as the principal branch
    value of (log(1+x) - log(1-x))/2."""
    x = convert_lossless(x)
    oldprec = mpf._prec
    a = abs(x)
    mpf._prec += 10
    high = a.exp + a.bc
    if high < -10:
        #print mpf._prec, x, x-(x**3)/3+(x**5)/5
        if high < (-(mpf._prec-10) * 0.3):
            return x - (x**3)/3 + (x**5)/5
        mpf._prec += (-high)
    t = 0.5*(log(1+x)-log(1-x))
    mpf._prec = oldprec
    return +t

def rand():
    """Return an mpf chosen randomly from [0, 1)."""
    return make_mpf(frand(mpf._prec))


__all__ = ["mpnumeric", "mpf", "mpc", "pi", "e", "euler", "clog2", "clog10",
  "j", "sqrt", "hypot", "exp", "log", "cos", "sin", "tan", "atan", "atan2",
  "power", "asin", "acos", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
  "arg", "degree", "rand", "inf", "nan", "floor", "ceil", "isnan"]

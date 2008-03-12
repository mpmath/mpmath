"""
This module defines the mpf, mpc classes, and standard functions for
operating with them.
"""
__docformat__ = 'plaintext'

__all__ = ["mpnumeric", "mpf", "mpc", "pi", "e", "euler", "ln2", "ln10",
  "j", "sqrt", "hypot", "exp", "log", "cos", "sin", "tan", "atan", "atan2",
  "power", "asin", "acos", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
  "arg", "degree", "rand", "inf", "nan", "floor", "ceil", "isnan", "almosteq",
  "ldexp", "catalan", "fraction", "nstr", "nprint", "mp", "extraprec",
  "extradps", "workprec", "workdps", "eps"]

from lib import *

int_types = (int, long)

rounding_table = {
  'floor' : round_floor,
  'ceiling' : round_ceiling,
  'down' : round_down,
  'up' : round_up,
  'nearest' : round_nearest,
  'default' : round_nearest
}

reverse_rounding_table = {
  round_floor : 'floor',
  round_ceiling : 'ceiling',
  round_down : 'down',
  round_up : 'up',
  round_nearest : 'nearest'
}

class Context(object):

    __slots__ = []

    def __repr__(self):
        lines = ["Mpmath settings:",
            ("  mp.prec = %s" % self.prec).ljust(30) + "[default: 53]",
            ("  mp.dps = %s" % self.dps).ljust(30) + "[default: 15]",
            ("  mp.rounding = '%s'" % self.rounding).ljust(30) + "[default: 'nearest']"]
        return "\n".join(lines)

    def default(self):
        global g_prec, g_dps, g_rounding
        g_prec = 53
        g_dps = 15
        g_rounding = round_nearest

    def set_prec(self, n):
        global g_prec, g_dps
        g_prec = max(1, int(n))
        g_dps = prec_to_dps(n)

    def set_dps(self, n):
        global g_prec, g_dps
        g_prec = dps_to_prec(n)
        g_dps = max(1, int(n))

    def set_rounding(self, s):
        global g_rounding
        try:
            g_rounding = rounding_table[s]
        except KeyError:
            raise ValueError(("unknown rounding mode: %s.\n" % s) + 
                "Value must be one of: %r" % rounding_table.keys())

    prec = property(lambda self: g_prec, set_prec)
    dps = property(lambda self: g_dps, set_dps)
    rounding = property(lambda self: reverse_rounding_table[g_rounding], set_rounding)

mp = Context()
mp.default()


class PrecisionManager:

    def __init__(self, precfun, dpsfun, normalize_output=False):
        self.precfun = precfun
        self.dpsfun = dpsfun
        self.normalize_output = normalize_output

    def __call__(self, f):
        def g(*args, **kwargs):
            orig = mp.prec
            try:
                if self.precfun:
                    mp.prec = self.precfun(mp.prec)
                else:
                    mp.dps = self.dpsfun(mp.dps)
                if self.normalize_output:
                    return +f(*args, **kwargs)
                else:
                    return f(*args, **kwargs)
            finally:
                mp.prec = orig
        g.__name__ = f.__name__
        g.__doc__ = f.__doc__
        return g

    def __enter__(self):
        self.orig_prec = mp.prec
        if self.precfun:
            mp.prec = self.precfun(mp.prec)
        else:
            mp.dps = self.dpsfun(mp.dps)

    def __exit__(self, exc_type, exc_val, exc_tb):
        mp.prec = self.orig_prec
        return False

def extraprec(n, normalize_output=False):
    """
    The block

        with extraprec(n):
            <code>

    increases the precision n bits, executes <code>, and then
    restores the precision.

    extraprec(n)(f) returns a decorated version of the function f
    that increases the working precision by n bits before execution,
    and restores the parent precision afterwards. With
    normalize_output=True, it rounds the return value to the parent
    precision.
    """
    return PrecisionManager(lambda p: p + n, None, normalize_output)

def extradps(n, normalize_output=False):
    """
    This function is analogous to extraprec (see documentation)
    but changes the decimal precision instead of the number of bits.
    """
    return PrecisionManager(None, lambda d: d + n, normalize_output)

def workprec(n, normalize_output=False):
    """
    The block

        with workprec(n):
            <code>

    sets the precision to n bits, executes <code>, and then restores
    the precision.

    workprec(n)(f) returns a decorated version of the function f
    that sets the precision to n bits before execution,
    and restores the precision afterwards. With normalize_output=True,
    it rounds the return value to the parent precision.
    """
    return PrecisionManager(lambda p: n, None, normalize_output)

def workdps(n, normalize_output=False):
    """
    This function is analogous to workprec (see documentation)
    but changes the decimal precision instead of the number of bits.
    """
    return PrecisionManager(None, lambda d: n, normalize_output)

class mpnumeric(object):
    """Base class for mpf and mpc. Calling mpnumeric(x) returns an mpf
    if x can be converted to an mpf (if it is a float, int, mpf, ...),
    and an mpc if x is complex."""
    __slots__ = []
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
    sufficient precision to represent it exactly. If x is a str, just
    convert it to an mpf with the current working precision (perhaps
    this should be done differently...)"""
    if isinstance(x, mpnumeric):
        return x
    if isinstance(x, int_types):
        return make_mpf(from_int(x))
    if isinstance(x, float):
        return make_mpf(from_float(x, 53, round_floor))
    if isinstance(x, complex):
        return mpc(x)
    if isinstance(x, basestring):
        return make_mpf(from_str(x, g_prec, g_rounding))
    if hasattr(x, '_mpf_'):
        return make_mpf(x._mpf_)
    if hasattr(x, '_mpc_'):
        return mpc(*x._mpc_)
    raise TypeError("cannot create mpf from " + repr(x))

def mpf_convert_rhs(x):
    if isinstance(x, int_types):
        return make_mpf(from_int(x, bitcount(x), round_floor))
    if isinstance(x, float):
        return make_mpf(from_float(x, 53, round_floor))
    if hasattr(x, '_mpf_'):
        return make_mpf(x._mpf_)
    return NotImplemented

def mpf_convert_lhs(x):
    if isinstance(x, complex_types):
        return mpc(x)
    if isinstance(x, int_types):
        return make_mpf(from_int(x, bitcount(x), round_floor))
    if isinstance(x, float):
        return make_mpf(from_float(x, 53, round_floor))
    if hasattr(x, '_mpf_'):
        return make_mpf(x._mpf_)
    return NotImplemented

new = object.__new__

ctx = [53, round_nearest]

class mpf(mpnumeric):
    """An mpf instance holds a real-valued floating-point number. mpf:s
    work analogously to Python floats, but support arbitrary-precision
    arithmetic."""

    #__metaclass__ = context
    __slots__ = ['_mpf_']

    def __new__(cls, val=fzero):
        """A new mpf can be created from a Python float, an int, a
        or a decimal string representing a number in floating-point
        format."""
        if isinstance(val, mpf):
            sign, man, exp, bc = val._mpf_
            if (not man) and exp:
                return val
            return make_mpf(normalize(sign, man, exp, bc, g_prec, g_rounding))
        elif isinstance(val, tuple):
            if len(val) == 2:
                return make_mpf(from_man_exp(val[0], val[1], g_prec, g_rounding))
            if len(val) == 4:
                sign, man, exp, bc = val
                return make_mpf(normalize(sign, man, exp, bc, g_prec, g_rounding))
            raise ValueError
        return +convert_lossless(val)

    man_exp = property(lambda self: self._mpf_[1:3])
    man = property(lambda self: self._mpf_[1])
    exp = property(lambda self: self._mpf_[2])
    bc = property(lambda self: self._mpf_[3])

    def __repr__(s): return "mpf('%s')" % to_str(s._mpf_, repr_dps(g_prec))
    def __str__(s): return to_str(s._mpf_, g_dps)
    def __hash__(s): return fhash(s._mpf_)
    def __int__(s): return to_int(s._mpf_)
    def __float__(s): return to_float(s._mpf_)
    def __complex__(s): return float(s) + 0j

    def __nonzero__(s):
        return s._mpf_ != fzero

    def __abs__(s): return make_mpf(fabs(s._mpf_, g_prec, g_rounding))
    def __pos__(s): return make_mpf(fpos(s._mpf_, g_prec, g_rounding))
    def __neg__(s): return make_mpf(fneg(s._mpf_, g_prec, g_rounding))

    def __eq__(s, t):
        if isinstance(t, mpf):
            return feq(s._mpf_, t._mpf_)
        if isinstance(t, complex_types):
            return mpc(s) == t
        t = mpf_convert_rhs(t)
        if t is NotImplemented:
            return False
        return feq(s._mpf_, t._mpf_)

    def __ne__(s, t):
        return not s.__eq__(t)

    def __cmp__(s, t):
        if not isinstance(t, mpf):
            t = mpf_convert_rhs(t)
            if t is NotImplemented:
                return t
        return fcmp(s._mpf_, t._mpf_)

    def __lt__(s, t):
        if s._mpf_ == fnan:
            return False
        r = s.__cmp__(t)
        if r is NotImplemented:
            return r
        return r == -1

    def __gt__(s, t):
        if s._mpf_ == fnan:
            return False
        r = s.__cmp__(t)
        if r is NotImplemented:
            return r
        return r == 1

    def __le__(s, t):
        if s._mpf_ == fnan:
            return False
        r = s.__cmp__(t)
        if r is NotImplemented:
            return r
        return r <= 0

    def __ge__(s, t):
        if s._mpf_ == fnan:
            return False
        r = s.__cmp__(t)
        if r is NotImplemented:
            return r
        return r >= 0

    def binop(s, t, f):
        if isinstance(t, complex_types):
            s = mpc(s)
            if f is feq: return s == t
            if f is fmul: return s * t
            if f is fadd: return s + t
            if f is fsub: return s - t
            if f is fdiv: return s / t
            raise ValueError("bad operation")
        t = mpf_convert_rhs(t)
        if t is NotImplemented:
            return t
        return make_mpf(f(s._mpf_, t._mpf_, g_prec, g_rounding))

    def __add__(s, t):
        r = new(mpf)
        sval = s._mpf_
        try:
            r._mpf_ = fadd(sval, t._mpf_, g_prec, g_rounding)
            return r
        except:
            if isinstance(t, int_types):
                r._mpf_ = fadd(sval, from_int(t), g_prec, g_rounding)
                return r
            return s.binop(t, fadd)

    def __sub__(s, t):
        r = new(mpf)
        sval = s._mpf_
        try:
            r._mpf_ = fsub(sval, t._mpf_, g_prec, g_rounding)
            return r
        except:
            if isinstance(t, int_types):
                r._mpf_ = fadd(sval, from_int(-t), g_prec, g_rounding)
                return r
            return s.binop(t, fsub)

    def __mul__(s, t):
        r = new(mpf)
        sval = s._mpf_
        try:
            r._mpf_ = fmul(sval, t._mpf_, g_prec, g_rounding)
            return r
        except:
            if isinstance(t, int_types):
                r._mpf_ = fmuli(sval, t, g_prec, g_rounding)
                return r
            return s.binop(t, fmul)

    def __rmul__(s, t):
        if isinstance(t, int_types):
            r = new(mpf)
            r._mpf_ = fmuli(s._mpf_, t, g_prec, g_rounding)
            return r
        return s.binop(t, fmul)

    def __div__(s, t):
        if isinstance(t, mpf):
            return make_mpf(fdiv(s._mpf_, t._mpf_, g_prec, g_rounding))
        return s.binop(t, fdiv)

    def __mod__(s, t):
        if isinstance(t, mpf):
            return make_mpf(fmod(s._mpf_, t._mpf_, g_prec, g_rounding))
        return s.binop(t, fmod)

    def __pow__(s, t):
        if isinstance(t, int_types):
            return make_mpf(fpowi(s._mpf_, t, g_prec, g_rounding))
        if not isinstance(t, mpf):
            if isinstance(t, complex_types):
                return power(s, t)
            t = mpf_convert_rhs(t)
            if t is NotImplemented:
                return t
        try:
            return make_mpf(fpow(s._mpf_, t._mpf_, g_prec, g_rounding))
        except ValueError:
            # Turns complex
            return power(s, t)

    __radd__ = __add__
    #__rmul__ = __mul__

    def __rsub__(s, t): return mpf_convert_lhs(t) - s

    def __rdiv__(s, t):
        if isinstance(t, int_types):
            r = new(mpf)
            r._mpf_ = fdivi(t, s._mpf_, g_prec, g_rounding)
            return r
        return mpf_convert_lhs(t) / s

    def __rpow__(s, t): return mpf_convert_lhs(t) ** s
    def __rmod__(s, t): return mpf_convert_lhs(t) % s

    def sqrt(s):
        return sqrt(s)

    def ae(s, t, rel_eps=None, abs_eps=None):
        return almosteq(s, t, rel_eps, abs_eps)


def make_mpf(tpl, construct=object.__new__, cls=mpf):
    """Create mpf verbatim from a given tuple of data."""
    a = construct(cls)
    a._mpf_ = tpl
    return a


inf = make_mpf(finf)
ninf = make_mpf(fninf)
nan = make_mpf(fnan)

def isnan(x):
    if not isinstance(x, mpf):
        return False
    return x._mpf_ == fnan

class mpc(mpnumeric):
    """An mpc represents a complex number using a pair of mpf:s (one
    for the real part and another for the imaginary part.) The mpc
    class behaves fairly similarly to Python's complex type."""

    def __new__(cls, real=0, imag=0):
        s = object.__new__(cls)
        if isinstance(real, complex_types):
            real, imag = real.real, real.imag
        elif hasattr(real, "_mpc_"):
            real, imag = real._mpc_
        s.real = mpf(real)
        s.imag = mpf(imag)
        return s

    def __repr__(s):
        r = repr(s.real)[4:-1]
        i = repr(s.imag)[4:-1]
        return "mpc(real=%s, imag=%s)" % (r, i)

    def __str__(s):
        return "(%s)" % complex_to_str(s.real._mpf_, s.imag._mpf_, g_dps)

    def __complex__(s): return complex(float(s.real), float(s.imag))
    def __pos__(s): return mpc(s.real, s.imag)
    def __abs__(s): return hypot(s.real, s.imag)
    def __neg__(s): return mpc(-s.real, -s.imag)
    def __nonzero__(s): return bool(s.real) or bool(s.imag)
    def conjugate(s): return mpc(s.real, -s.imag)

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

    def __add__(s, t):
        if not isinstance(t, mpc):
            t = mpc(t)
        return mpc(s.real+t.real, s.imag+t.imag)

    def __sub__(s, t):
        if not isinstance(t, mpc):
            t = mpc(t)
        return mpc(s.real-t.real, s.imag-t.imag)

    def __mul__(s, t):
        if not isinstance(t, mpc):
            t = mpc(t)
        return mpc(*fcmul(s.real._mpf_, s.imag._mpf_, t.real._mpf_, t.imag._mpf_,
            g_prec, g_rounding))

    def __div__(s, t):
        if not isinstance(t, mpc):
            t = mpc(t)
        a = s.real; b = s.imag; c = t.real; d = t.imag
        mag = c*c + d*d
        return mpc((a*c+b*d)/mag, (b*c-a*d)/mag)

    def __pow__(s, n):
        if n == 0: return mpc(1)
        if n == 1: return +s
        if n == -1: return 1/s
        if n == 2: return s*s
        if isinstance(n, int_types) and n > 0:
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

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(s, t): return (-s) + t
    def __rpow__(s, t): return convert_lossless(t) ** s
    def __rdiv__(s, t): return mpc(t) / s

    def ae(s, t, rel_eps=None, abs_eps=None):
        return almosteq(s, t, rel_eps, abs_eps)


complex_types = (complex, mpc)

def make_mpc((re, im), construct=object.__new__, cls=mpc):
    a = construct(cls)
    a.real, a.imag = make_mpf(re), make_mpf(im)
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
    def _mpf_(self):
        return self.func(g_prec, g_rounding)

    def __repr__(self):
        return "<%s: %s~>" % (self.name, nstr(self))


_180 = from_int(180, 10, round_floor)

pi = constant(fpi, "pi")
degree = constant(fdegree, "degree")
e = constant(fe, "e")
euler = constant(fgamma, "Euler's constant gamma")
ln2 = constant(flog2, "log 2")
ln10 = constant(flog10, "log 10")
catalan = constant(fcatalan, "Catalan's constant")
eps = constant(lambda p, r: (0, 1, -p+1, 1), "epsilon of working precision")

def fraction(p, q):
    """Given Python integers p, q, return a lazy mpf with value p/q.
    The value is updated with the precision.

        >>> mp.dps = 15
        >>> a = fraction(1,100)
        >>> b = mpf(1)/100
        >>> print a; print b
        0.01
        0.01
        >>> mp.dps = 30
        >>> print a; print b
        0.01
        0.0100000000000000002081668171172
        >>> mp.dps = 15
    """
    return constant(lambda prec, rnd: from_rational(p, q, prec, rnd),
        '%s/%s' % (p, q))

def sqrt(x):
    """For real x >= 0, returns the square root of x. For negative or
    complex x, returns the principal branch of the complex square root
    of x."""
    x = convert_lossless(x)
    if isinstance(x, mpf) and not x._mpf_[0]:
        return make_mpf(fsqrt(x._mpf_, g_prec, g_rounding))
    x = mpc(x)
    return make_mpc(fcsqrt(x.real._mpf_, x.imag._mpf_, g_prec, g_rounding))

def hypot(x, y):
    """Returns the Euclidean distance sqrt(x*x + y*y). Both x and y
    must be real."""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return make_mpf(fhypot(x._mpf_, y._mpf_, g_prec, g_rounding))

def floor(x):
    """Computes the floor function of x. Note: returns an mpf, not a
    Python int. If x is larger than the precision, it will be rounded,
    not necessarily in the floor direction."""
    x = convert_lossless(x)
    return make_mpf(ffloor(x._mpf_, g_prec, g_rounding))

def ceil(x):
    """Computes the ceiling function of x. Note: returns an mpf, not a
    Python int. If x is larger than the precision, it will be rounded,
    not necessarily in the ceiling direction."""
    x = convert_lossless(x)
    return make_mpf(fceil(x._mpf_, g_prec, g_rounding))

def ldexp(x, n):
    """Calculate mpf(x) * 2**n efficiently. No rounding is performed."""
    x = convert_lossless(x)
    return make_mpf(fshift(x._mpf_, n))

# Since E-functions simply map reals to reals and complexes to complexes, we
# can construct all of them the same way (unlike log, sqrt, etc)
def ef(name, real_f, complex_f, doc):
    def f(x):
        try:
            return make_mpf(real_f(x._mpf_, g_prec, g_rounding))
        except AttributeError:
            pass
        x = convert_lossless(x)
        if isinstance(x, mpf):
            return make_mpf(real_f(x._mpf_, g_prec, g_rounding))
        else:
            return make_mpc(complex_f(x.real._mpf_, x.imag._mpf_, g_prec, g_rounding))
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
        return make_mpf(ftan(x._mpf_, g_prec, g_rounding))
    # the complex division can cause enormous cancellation.
    # TODO: handle more robustly
    mp.prec += 20
    t = sin(x) / cos(x)
    mp.prec -= 20
    return +t

def tanh(x):
    """Returns the hyperbolic tangent of x."""
    x = convert_lossless(x)
    oldprec = mp.prec
    a = abs(x)
    mp.prec += 10
    high = a.exp + a.bc
    if high < -10:
        if high < (-(mp.prec-10) * 0.3):
            return x - (x**3)/3 + 2*(x**5)/15
        mp.prec += (-high)
    a = exp(2*x)
    t = (a-1)/(a+1)
    mp.prec = oldprec
    return +t

@extraprec(5)
def arg(x):
    """Returns the complex argument (phase) of x. The returned value is
    an mpf instance. The argument is here defined to satisfy
    -pi < arg(x) <= pi. On the negative real half-axis, it is taken to
    be +pi."""
    x = mpc(x)
    return atan2(x.imag, x.real)

def log(x, b=None):
    """Returns the base-b logarithm of x. If b is unspecified, return
    the natural (base-e) logarithm. log(x, b) is defined as
    log(x)/log(b). log(0) raises ValueError.

    The natural logarithm is real if x > 0 and complex if x < 0 or if x
    is complex. The principal branch of the complex logarithm is chosen,
    for which Im(log(x)) = -pi < arg(x) <= pi. """
    if b is not None:
        mp.prec += 3
        a = log(x) / log(b)
        mp.prec -= 3
        return +a
    x = convert_lossless(x)
    if not x:
        raise ValueError, "logarithm of 0"
    if isinstance(x, mpf) and x > 0:
        return make_mpf(flog(x._mpf_, g_prec, g_rounding))
    else:
        return mpc(log(abs(x)), arg(x))

@extraprec(10, normalize_output=True)
def power(x, y):
    """Returns x**y = exp(y*log(x)) for real or complex x and y."""
    # TODO: better estimate for extra precision needed
    return exp(y * log(x))

def atan(x):
    """Returns the inverse tangent of x."""
    x = convert_lossless(x)
    if isinstance(x, mpf):
        return make_mpf(fatan(x._mpf_, g_prec, g_rounding))
    # TODO: maybe get this to agree with Python's cmath atan about the
    # branch to choose on the imaginary axis
    # TODO: handle cancellation robustly
    mp.prec += 10
    t = (0.5j)*(log(1-1j*x) - log(1+1j*x))
    mp.prec -= 10
    return +t

def atan2(y,x):
    """atan2(y, x) has the same magnitude as atan(y/x) but accounts for
    the signs of y and x. (Defined for real x and y only.)"""
    x = convert_lossless(x)
    y = convert_lossless(y)
    if y < 0:
        return -atan2(-y, x)
    if not x and not y:
        return mpf(0)
    if y > 0 and x == 0:
        mp.prec += 2
        t = pi/2
        mp.prec -= 2
        return t
    mp.prec += 2
    if x > 0:
        a = atan(y/x)
    else:
        a = pi - atan(-y/x)
    mp.prec -= 2
    return +a

# TODO: robustly deal with cancellation in all of the following functions

@extraprec(10, normalize_output=True)
def _asin_complex(z):
    return -1j * log(1j * z + sqrt(1 - z*z))

def asin(x):
    """Returns the inverse sine of x. Outside the range [-1, 1], the
    result is complex and defined as the principal branch value of
    -i * log(i * x + sqrt(1 - x**2))."""
    x = convert_lossless(x)
    if isinstance(x, mpf) and abs(x) <= 1:
        return _asin_complex(x).real
    return _asin_complex(x)

def _acos_complex(z):
    mp.prec += 10
    t = pi/2 + 1j * log(1j * z + sqrt(1 - z*z))
    mp.prec -= 10
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
    oldprec = mp.prec
    a = abs(x)
    mp.prec += 10
    high = a.exp + a.bc
    if high < -10:
        if high < (-(mp.prec-10) * 0.3):
            return x - (x**3)/6 + 3*(x**5)/40
        mp.prec += (-high)
    t = log(x + sqrt(x**2 + 1))
    mp.prec = oldprec
    return +t

@extraprec(10, normalize_output=True)
def acosh(x):
    """Returns the inverse hyperbolic cosine of x. The value is
    given by log(x + sqrt(1 + x**2)), where the principal branch is
    used when the result is complex."""
    x = convert_lossless(x)
    return log(x + sqrt(x-1)*sqrt(x+1))

def atanh(x):
    """Returns the inverse hyperbolic tangent of x. Outside the range
    [-1, 1], the result is complex and defined as the principal branch
    value of (log(1+x) - log(1-x))/2."""
    x = convert_lossless(x)
    oldprec = mp.prec
    a = abs(x)
    mp.prec += 10
    high = a.exp + a.bc
    if high < -10:
        #print mp.prec, x, x-(x**3)/3+(x**5)/5
        if high < (-(mp.prec-10) * 0.3):
            return x - (x**3)/3 + (x**5)/5
        mp.prec += (-high)
    t = 0.5*(log(1+x)-log(1-x))
    mp.prec = oldprec
    return +t

def rand():
    """Return an mpf chosen randomly from [0, 1)."""
    return make_mpf(frand(mp.prec))

def almosteq(s, t, rel_eps=None, abs_eps=None):
    """
    Determine whether the difference between s and t is smaller
    than a given epsilon.

    Both a maximum relative difference and a maximum difference
    ('epsilons') may be specified. The absolute difference is
    defined as |s-t| and the relative difference is defined
    as |s-t|/max(|s|, |t|).

    If only one epsilon is given, both are set to the same value.
    If none is given, both epsilons are set to 2**(-prec+m) where
    prec is the current working precision and m is a small integer.
    """
    t = convert_lossless(t)
    if abs_eps is None and rel_eps is None:
        rel_eps = abs_eps = make_mpf((0, 1, -g_prec+4, 1))
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

def nstr(x, n=6):
    """Convert an mpf or mpc to a decimal string literal with n significant
    digits. The small default value for n is chosen to make this function
    useful for printing collections of numbers.

    If x is a list or tuple, the function is applied to each element.
    For unrecognized classes, this simply returns str(x).

    There is also a companion function nprint that prints the string
    instead of returning it.

        >>> nstr([+pi, ldexp(1,-500)])
        '[3.14159, 3.05494e-151]'
        >>> print([+pi, ldexp(1,-500)])
        [3.14159, 3.05494e-151]
    """
    if isinstance(x, list):
        return "[%s]" % (", ".join(nstr(c, n) for c in x))
    if isinstance(x, tuple):
        return "(%s)" % (", ".join(nstr(c, n) for c in x))
    if isinstance(x, mpf):
        return to_str(x._mpf_, n)
    if isinstance(x, mpc):
        return "(%s + %sj)" % (to_str(x.real._mpf_, n), to_str(x.imag._mpf_, n))
    if isinstance(x, basestring):
        return repr(x)
    return str(x)

def nprint(x, n=6):
    """Print the result of nstr(x, n)."""
    print nstr(x, n)

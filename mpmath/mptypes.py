"""
This module defines the mpf, mpc classes, and standard functions for
operating with them.
"""
__docformat__ = 'plaintext'

__all__ = ["mpnumeric", "mpf", "mpc", "pi", "e", "ln2", "ln10",
  "j", "sqrt", "hypot", "exp", "log", "cos", "sin", "tan", "atan", "atan2",
  "power", "asin", "acos", "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
  "arg", "degree", "rand", "inf", "nan", "floor", "ceil", "isnan", "almosteq",
  "ldexp", "fraction", "nstr", "nprint", "mp", "extraprec",
  "extradps", "workprec", "workdps", "eps", "convert_lossless", "make_mpf",
  "make_mpc", "sec", "csc", "cot", "sech", "csch", "coth",
  "asec", "acsc", "acot", "asech", "acsch", "acoth", "arange",
  "ln", "log10", "frexp", "radians", "degrees", "modf", "cbrt", "nthroot",
  "sign", "plot", "isinf", "mpi", "isint"]

from lib import *
from libmpc import *
from libmpi import *

class mpnumeric(object):
    """Base class for mpf and mpc. Calling mpnumeric(x) returns an mpf
    if x can be converted to an mpf (if it is a float, int, mpf, ...),
    and an mpc if x is complex."""
    __slots__ = []
    def __new__(cls, val):
        # TODO: should maybe normalize here
        if isinstance(val, cls): return val
        if isinstance(val, complex): return mpc(val)
        return mpf(val)

def convert_lossless(x, strings=True):
    """Attempt to convert x to an mpf or mpc losslessly. If x is an
    mpf or mpc, return it unchanged. If x is an int, create an mpf with
    sufficient precision to represent it exactly. If x is a str, just
    convert it to an mpf with the current working precision (perhaps
    this should be done differently...)"""
    if isinstance(x, mpnumeric): return x
    if isinstance(x, int_types): return make_mpf(from_int(x))
    if isinstance(x, float): return make_mpf(from_float(x))
    if isinstance(x, complex): return mpc(x)
    if strings and isinstance(x, basestring): return make_mpf(from_str(x, *prec_rounding))
    if hasattr(x, '_mpf_'): return make_mpf(x._mpf_)
    if hasattr(x, '_mpc_'): return make_mpc(x._mpc_)
    if hasattr(x, '_mpmath_'): return convert_lossless(x._mpmath_(*prec_rounding))
    raise TypeError("cannot create mpf from " + repr(x))

def try_convert_mpf_value(x, prec, rounding):
    if isinstance(x, float): return from_float(x)
    if hasattr(x, '_mpf_'): return x._mpf_
    if hasattr(x, '_mpmath_'):
        t = convert_lossless(x._mpmath_(prec, rounding))
        if isinstance(t, mpf):
            return t._mpf_
    return NotImplemented

def mpf_convert_arg(x, prec, rounding):
    if isinstance(x, int_types): return from_int(x)
    if isinstance(x, float): return from_float(x)
    if isinstance(x, basestring): return from_str(x, prec, rounding)
    if isinstance(x, constant): return x.func(prec, rounding)
    if hasattr(x, '_mpf_'): return x._mpf_
    if hasattr(x, '_mpmath_'):
        t = convert_lossless(x._mpmath_(prec, rounding))
        if isinstance(t, mpf):
            return t._mpf_
    raise TypeError("cannot create mpf from " + repr(x))

def mpf_convert_rhs(x):
    if isinstance(x, int_types): return from_int(x)
    if isinstance(x, float): return from_float(x)
    if isinstance(x, complex_types): return mpc(x)
    if hasattr(x, '_mpf_'): return x._mpf_
    if hasattr(x, '_mpmath_'):
        t = convert_lossless(x._mpmath_(*prec_rounding))
        if isinstance(t, mpf):
            return t._mpf_
        return t
    return NotImplemented

def mpf_convert_lhs(x):
    x = mpf_convert_rhs(x)
    if type(x) is tuple:
        return make_mpf(x)
    return x

def mpc_convert_lhs(x):
    try:
        return convert_lossless(x)
    except TypeError:
        return NotImplemented

new = object.__new__

class mpf(mpnumeric):
    """
    An mpf instance holds a real-valued floating-point number. mpf:s
    work analogously to Python floats, but support arbitrary-precision
    arithmetic.
    """
    __slots__ = ['_mpf_']

    def __new__(cls, val=fzero, **kwargs):
        """A new mpf can be created from a Python float, an int, a
        or a decimal string representing a number in floating-point
        format."""
        prec, rounding = prec_rounding
        if kwargs:
            prec = kwargs.get('prec', prec)
            if 'dps' in kwargs:
                prec = dps_to_prec(kwargs['dps'])
            rounding = kwargs.get('rounding', rounding)
        if type(val) is cls:
            sign, man, exp, bc = val._mpf_
            if (not man) and exp:
                return val
            return make_mpf(normalize(sign, man, exp, bc, prec, rounding))
        elif type(val) is tuple:
            if len(val) == 2:
                return make_mpf(from_man_exp(val[0], val[1], prec, rounding))
            if len(val) == 4:
                sign, man, exp, bc = val
                return make_mpf(normalize(sign, MP_BASE(man), exp, bc, prec, rounding))
            raise ValueError
        else:
            return make_mpf(fpos(mpf_convert_arg(val, prec, rounding), prec, rounding))

    man_exp = property(lambda self: self._mpf_[1:3])
    man = property(lambda self: self._mpf_[1])
    exp = property(lambda self: self._mpf_[2])
    bc = property(lambda self: self._mpf_[3])

    real = property(lambda self: self)
    imag = property(lambda self: zero)

    def __getstate__(self): return to_pickable(self._mpf_)
    def __setstate__(self, val): self._mpf_ = from_pickable(val)

    def __repr__(s): return "mpf('%s')" % to_str(s._mpf_, repr_dps(mp.prec))
    def __str__(s): return to_str(s._mpf_, mp.dps)
    def __hash__(s): return fhash(s._mpf_)
    def __int__(s): return int(to_int(s._mpf_))
    def __long__(s): return long(to_int(s._mpf_))
    def __float__(s): return to_float(s._mpf_)
    def __complex__(s): return complex(float(s))
    def __nonzero__(s): return s._mpf_ != fzero
    def __abs__(s): return make_mpf(fabs(s._mpf_, *prec_rounding))
    def __pos__(s): return make_mpf(fpos(s._mpf_, *prec_rounding))
    def __neg__(s): return make_mpf(fneg(s._mpf_, *prec_rounding))

    def _cmp(s, t, func):
        if hasattr(t, '_mpf_'):
            t = t._mpf_
        else:
            t = mpf_convert_rhs(t)
            if t is NotImplemented:
                return t
        return func(s._mpf_, t)

    def __cmp__(s, t): return s._cmp(t, fcmp)
    def __lt__(s, t): return s._cmp(t, flt)
    def __gt__(s, t): return s._cmp(t, fgt)
    def __le__(s, t): return s._cmp(t, fle)
    def __ge__(s, t): return s._cmp(t, fge)

    def __ne__(s, t):
        v = s.__eq__(t)
        if v is NotImplemented:
            return v
        return not v

    def __rsub__(s, t):
        prec, rounding = prec_rounding
        if type(t) in int_types:
            return make_mpf(fsub(from_int(t), s._mpf_, prec, rounding))
        t = mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t - s

    def __rdiv__(s, t):
        prec, rounding = prec_rounding
        if isinstance(t, int_types):
            return make_mpf(fdivi(t, s._mpf_, prec, rounding))
        t = mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t / s

    def __rpow__(s, t):
        t = mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t ** s

    def __rmod__(s, t):
        t = mpf_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t % s

    def sqrt(s):
        return sqrt(s)

    def ae(s, t, rel_eps=None, abs_eps=None):
        return almosteq(s, t, rel_eps, abs_eps)


mpf_binary_op = """
def %NAME%(self, other):
    prec, rounding = prec_rounding
    sval = self._mpf_
    if hasattr(other, '_mpf_'):
    #try:
        tval = other._mpf_
        %WITH_MPF%
    #except AttributeError:
    #    pass
    ttype = type(other)
    if ttype in int_types:
        %WITH_INT%
    elif ttype is float:
        tval = from_float(other)
        %WITH_MPF%
    elif ttype is mpc:
        tval = other._mpc_
        %WITH_MPC%
    elif ttype is complex:
        tval = from_float(other.real), from_float(other.imag)
        %WITH_MPC%
    if isinstance(other, mpnumeric):
        return NotImplemented
    try:
        other = convert_lossless(other, strings=False)
    except TypeError:
        return NotImplemented
    return self.%NAME%(other)
"""

return_mpf = "; obj = new(mpf); obj._mpf_ = val; return obj"
return_mpc = "; obj = new(mpc); obj._mpc_ = val; return obj"

mpf_pow_same = """
        try:
            val = fpow(sval, tval, prec, rounding) %s
        except ComplexResult:
            if mp.trap_complex:
                raise
            val = mpc_pow((sval, fzero), (tval, fzero), prec, rounding) %s
""" % (return_mpf, return_mpc)

def binary_op(name, with_mpf='', with_int='', with_mpc=''):
    code = mpf_binary_op
    code = code.replace("%WITH_INT%", with_int)
    code = code.replace("%WITH_MPC%", with_mpc)
    code = code.replace("%WITH_MPF%", with_mpf)
    code = code.replace("%NAME%", name)
    np = {}
    exec code in globals(), np
    return np[name]

mpf.__eq__ = binary_op('__eq__',
    'return feq(sval, tval)',
    'return feq(sval, from_int(other))',
    'return (tval[1] == fzero) and feq(tval[0], sval)')

mpf.__add__ = binary_op('__add__',
    'val = fadd(sval, tval, prec, rounding)' + return_mpf,
    'val = fadd(sval, from_int(other), prec, rounding)' + return_mpf,
    'val = mpc_add_mpf(tval, sval, prec, rounding)' + return_mpc)

mpf.__sub__ = binary_op('__sub__',
    'val = fsub(sval, tval, prec, rounding)' + return_mpf,
    'val = fsub(sval, from_int(other), prec, rounding)' + return_mpf,
    'val = mpc_sub((sval, fzero), tval, prec, rounding)' + return_mpc)

mpf.__mul__ = binary_op('__mul__',
    'val = fmul(sval, tval, prec, rounding)' + return_mpf,
    'val = fmuli(sval, other, prec, rounding)' + return_mpf,
    'val = mpc_mul_mpf(tval, sval, prec, rounding)' + return_mpc)

mpf.__div__ = binary_op('__div__',
    'val = fdiv(sval, tval, prec, rounding)' + return_mpf,
    'val = fdiv(sval, from_int(other), prec, rounding)' + return_mpf,
    'val = mpc_div((sval, fzero), tval, prec, rounding)' + return_mpc)

mpf.__mod__ = binary_op('__mod__',
    'val = fmod(sval, tval, prec, rounding)' + return_mpf,
    'val = fmod(sval, from_int(other), prec, rounding)' + return_mpf,
    'raise NotImplementedError("complex modulo")')

mpf.__pow__ = binary_op('__pow__',
    mpf_pow_same,
    'val = fpowi(sval, other, prec, rounding)' + return_mpf,
    'val = mpc_pow((sval, fzero), tval, prec, rounding)' + return_mpc)

mpf.__radd__ = mpf.__add__
mpf.__rmul__ = mpf.__mul__
mpf.__truediv__ = mpf.__div__
mpf.__rtruediv__ = mpf.__rdiv__


class mpc(mpnumeric):
    """
    An mpc represents a complex number using a pair of mpf:s (one
    for the real part and another for the imaginary part.) The mpc
    class behaves fairly similarly to Python's complex type.
    """

    __slots__ = ['_mpc_']

    def __new__(cls, real=0, imag=0):
        s = object.__new__(cls)
        if isinstance(real, complex_types):
            real, imag = real.real, real.imag
        elif hasattr(real, "_mpc_"):
            s._mpc_ = real._mpc_
            return s
        real = mpf(real)
        imag = mpf(imag)
        s._mpc_ = (real._mpf_, imag._mpf_)
        return s

    real = property(lambda self: make_mpf(self._mpc_[0]))
    imag = property(lambda self: make_mpf(self._mpc_[1]))

    def __getstate__(self):
        return to_pickable(self._mpc_[0]), to_pickable(self._mpc_[1])

    def __setstate__(self, val):
        self._mpc_ = from_pickable(val[0]), from_pickable(val[1])

    def __repr__(s):
        r = repr(s.real)[4:-1]
        i = repr(s.imag)[4:-1]
        return "mpc(real=%s, imag=%s)" % (r, i)

    def __str__(s):
        return "(%s)" % complex_to_str(s.real._mpf_, s.imag._mpf_, mp.dps)

    def __complex__(s): return complex(float(s.real), float(s.imag))
    def __pos__(s): return mpc(s.real, s.imag)
    def __abs__(s): return make_mpf(mpc_abs(s._mpc_, *prec_rounding))
    def __neg__(s): return mpc(-s.real, -s.imag)
    def __nonzero__(s): return bool(s.real) or bool(s.imag)
    def conjugate(s): return mpc(s.real, -s.imag)

    def __eq__(s, t):
        if not isinstance(t, mpc):
            if isinstance(t, str):
                return False
            t = mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
        return s.real == t.real and s.imag == t.imag

    def __ne__(s, t):
        b = s.__eq__(t)
        if b is NotImplemented:
            return b
        return not b

    def _compare(*args):
        raise TypeError("no ordering relation is defined for complex numbers")

    __gt__ = _compare
    __le__ = _compare
    __gt__ = _compare
    __ge__ = _compare

    def __add__(s, t):
        prec, rounding = prec_rounding
        if not isinstance(t, mpc):
            t = mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if isinstance(t, mpf):
                return make_mpc(mpc_add_mpf(s._mpc_, t._mpf_, prec, rounding))
        return make_mpc(mpc_add(s._mpc_, t._mpc_, prec, rounding))

    def __sub__(s, t):
        prec, rounding = prec_rounding
        if not isinstance(t, mpc):
            t = mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if isinstance(t, mpf):
                return make_mpc(mpc_sub_mpf(s._mpc_, t._mpf_, prec, rounding))
        return make_mpc(mpc_sub(s._mpc_, t._mpc_, prec, rounding))

    def __mul__(s, t):
        prec, rounding = prec_rounding
        if not isinstance(t, mpc):
            if isinstance(t, int_types):
                return make_mpc(mpc_mul_int(s._mpc_, t, prec, rounding))
            t = mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if isinstance(t, mpf):
                return make_mpc(mpc_mul_mpf(s._mpc_, t._mpf_, prec, rounding))
            t = mpc(t)
        return make_mpc(mpc_mul(s._mpc_, t._mpc_, prec, rounding))

    def __div__(s, t):
        prec, rounding = prec_rounding
        if not isinstance(t, mpc):
            t = mpc_convert_lhs(t)
            if t is NotImplemented:
                return t
            if isinstance(t, mpf):
                return make_mpc(mpc_div_mpf(s._mpc_, t._mpf_, prec, rounding))
        return make_mpc(mpc_div(s._mpc_, t._mpc_, prec, rounding))

    def __pow__(s, t):
        prec, rounding = prec_rounding
        if isinstance(t, int_types):
            return make_mpc(mpc_pow_int(s._mpc_, t, prec, rounding))
        t = mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        if isinstance(t, mpf):
            return make_mpc(mpc_pow_mpf(s._mpc_, t._mpf_, prec, rounding))
        return make_mpc(mpc_pow(s._mpc_, t._mpc_, prec, rounding))

    __radd__ = __add__

    def __rsub__(s, t):
        t = mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t - s

    def __rmul__(s, t):
        prec, rounding = prec_rounding
        if isinstance(t, int_types):
            return make_mpc(mpc_mul_int(s._mpc_, t, prec, rounding))
        t = mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t * s

    def __rdiv__(s, t):
        t = mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t / s

    def __rpow__(s, t):
        t = mpc_convert_lhs(t)
        if t is NotImplemented:
            return t
        return t ** s

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def ae(s, t, rel_eps=None, abs_eps=None):
        return almosteq(s, t, rel_eps, abs_eps)


complex_types = (complex, mpc)


class mpi(mpnumeric):
    """
    Interval arithmetic class. Precision is controlled by mp.prec.
    """

    def __new__(cls, a, b=None):
        if isinstance(a, mpi):
            return a
        if b is None:
            b = a
        a = mpf(a, rounding=round_floor)
        b = mpf(b, rounding=round_ceiling)
        if isnan(a) or isnan(b):
            a, b = -inf, inf
        assert a <= b, "endpoints must be properly ordered"
        return make_mpi((a._mpf_, b._mpf_))

    @property
    def a(self):
        return make_mpf(self._val[0])

    @property
    def b(self):
        return make_mpf(self._val[1])

    @property
    def mid(self):
        return make_mpf(mpi_mid(self._val, mp.prec))

    @property
    def delta(self):
        return make_mpf(mpi_delta(self._val, mp.prec))

    def __contains__(self, t):
        t = mpi(t)
        return (self.a <= t.a) and (t.b <= self.b)

    def __repr__(self):
        return mpi_str(self._val, mp.prec)

    __str__ = __repr__

    def __eq__(self, other):
        if not isinstance(other, mpi):
            try:
                other = mpi(other)
            except:
                return NotImplemented
        return (self.a == other.a) and (self.b == other.b)

    def __abs__(self):
        return make_mpi(mpi_abs(self._val, mp.prec))

    def __pos__(self):
        return make_mpi(mpi_pos(self._val, mp.prec))

    def __neg__(self):
        return make_mpi(mpi_neg(self._val, mp.prec))

    def __add__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_add(self._val, other._val, mp.prec))

    def __sub__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_sub(self._val, other._val, mp.prec))

    def __mul__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_mul(self._val, other._val, mp.prec))

    def __div__(self, other):
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_div(self._val, other._val, mp.prec))

    def __pow__(self, other):
        if isinstance(other, (int, long)):
            return make_mpi(mpi_pow_int(self._val, int(other), mp.prec))
        if not isinstance(other, mpi):
            other = mpi(other)
        return make_mpi(mpi_pow(self._val, other._val, mp.prec))

    def __rsub__(s, t):
        return mpi(t) - s

    def __rdiv__(s, t):
        return mpi(t) / s

    def __rpow__(s, t):
        return mpi(t) ** s

    __radd__ = __add__
    __rmul__ = __mul__
    __truediv__ = __div__
    __rtruediv__ = __rdiv__
    __floordiv__ = __div__
    __rfloordiv__ = __rdiv__

def make_mpi(val, cls=mpi):
    a = new(cls)
    a._val = val
    return a

def make_mpf(v, cls=mpf):
    a = new(cls)
    a._mpf_ = v
    return a

def make_mpc(v, cls=mpc):
    a = new(cls)
    a._mpc_ = v
    return a

one = make_mpf(fone)
zero = make_mpf(fzero)
inf = make_mpf(finf)
ninf = make_mpf(fninf)
nan = make_mpf(fnan)
j = mpc(0,1)

def isnan(x):
    if not isinstance(x, mpf):
        return False
    return x._mpf_ == fnan

def isinf(x):
    if not isinstance(x, mpf):
        return False
    return x._mpf_ in (finf, fninf)

def isint(x):
    if isinstance(x, int_types):
        return True
    try:
        x = convert_lossless(x)
    except:
        return False
    if isinstance(x, mpf):
        if isnan(x) or isinf(x):
            return False
        return x == int(x)
    return False

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
        prec, rounding = prec_rounding
        return self.func(prec, rounding)

    def __repr__(self):
        return "<%s: %s~>" % (self.name, nstr(self))

pi = constant(fpi, "pi")
degree = constant(fdegree, "degree")
e = constant(fe, "e")
ln2 = constant(flog2, "log 2")
ln10 = constant(flog10, "log 10")
eps = constant(lambda p, r: (0, MP_ONE, -p+1, 1), "epsilon of working precision")

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



def funcwrapper(f):
    def g(*args, **kwargs):
        orig = mp.prec
        try:
            args = [convert_lossless(z) for z in args]
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
            x = convert_lossless(x)
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
                return make_mpc(complex_f((x._mpf_, fzero), prec, rounding))
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

sqrt = mpfunc('sqrt', fsqrt, mpc_sqrt, "principal square root", mpi_sqrt)
cbrt = mpfunc('cbrt', fcbrt, mpc_cbrt, "principal cubic root")
exp = mpfunc('exp', fexp, mpc_exp, "exponential function", mpi_exp)
ln = mpfunc('ln', flog, mpc_log, "natural logarithm", mpi_log)

cos = mpfunc('cos', fcos, mpc_cos, "cosine")
sin = mpfunc('sin', fsin, mpc_sin, "sine")
tan = mpfunc('tan', ftan, mpc_tan, "tangent")
cosh = mpfunc('cosh', fcosh, mpc_cosh, "hyperbolic cosine")
sinh = mpfunc('sinh', fsinh, mpc_sinh, "hyperbolic sine")
tanh = mpfunc('tanh', ftanh, mpc_tanh, "hyperbolic tangent")

acos = mpfunc('acos', facos, mpc_acos, "inverse cosine")
asin = mpfunc('asin', fasin, mpc_asin, "inverse sine")
atan = mpfunc('atan', fatan, mpc_atan, "inverse tangent")
asinh = mpfunc('asinh', fasinh, mpc_asinh, "inverse hyperbolic sine")
acosh = mpfunc('acosh', facosh, mpc_acosh, "inverse hyperbolic cosine")
atanh = mpfunc('atanh', fatanh, mpc_atanh, "inverse hyperbolic tangent")

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

@funcwrapper
def nthroot(x, n):
    """principal n-th root"""
    n = int(n)
    if isinstance(x, mpf):
        try:
            return make_mpf(fnthroot(x._mpf_, n, *prec_rounding))
        except ComplexResult:
            if mp.trap_complex:
                raise
            x = (x._mpf_, fzero)
    else:
        x = x._mpc_
    return make_mpc(mpc_nthroot(x, n, *prec_rounding))

def hypot(x, y):
    """Returns the Euclidean distance sqrt(x*x + y*y). Both x and y
    must be real."""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return make_mpf(fhypot(x._mpf_, y._mpf_, *prec_rounding))

def floor(x):
    """Computes the floor function of x. Note: returns an mpf, not a
    Python int. If x is larger than the precision, it will be rounded,
    not necessarily in the floor direction."""
    x = convert_lossless(x)
    return make_mpf(ffloor(x._mpf_, *prec_rounding))

def ceil(x):
    """Computes the ceiling function of x. Note: returns an mpf, not a
    Python int. If x is larger than the precision, it will be rounded,
    not necessarily in the ceiling direction."""
    x = convert_lossless(x)
    return make_mpf(fceil(x._mpf_, *prec_rounding))

def ldexp(x, n):
    """Calculate mpf(x) * 2**n efficiently. No rounding is performed."""
    x = convert_lossless(x)
    return make_mpf(fshift(x._mpf_, n))

def frexp(x):
    """Convert x to a scaled number y in the range [0.5, 1). Returns
    (y, n) such that x = y * 2**n. No rounding is performed."""
    x = convert_lossless(x)
    y, n = mpf_frexp(x._mpf_)
    return make_mpf(y), n

def sign(x):
    """Return sign(x), defined as x/abs(x), or 0 for x = 0."""
    x = convert_lossless(x)
    if not x or isnan(x):
        return x
    if isinstance(x, mpf):
        return cmp(x, 0)
    return x / abs(x)

@extraprec(5)
def arg(x):
    """Returns the complex argument (phase) of x. The returned value is
    an mpf instance. The argument is here defined to satisfy
    -pi < arg(x) <= pi. On the negative real half-axis, it is taken to
    be +pi."""
    x = mpc(x)
    return atan2(x.imag, x.real)

@funcwrapper
def log(x, b=None):
    """Returns the base-b logarithm of x. If b is unspecified, return
    the natural (base-e) logarithm. log(x, b) is defined as
    log(x)/log(b). log(0) raises ValueError.

    The natural logarithm is real if x > 0 and complex if x < 0 or if x
    is complex. The principal branch of the complex logarithm is chosen,
    for which Im(log(x)) = -pi < arg(x) <= pi. """
    if b is None:
        return ln(x)
    return ln(x) / ln(b)

def log10(x):
    """Base-10 logarithm. Equivalent to log(x,10)."""
    return log(x, 10)

def power(x, y):
    """Converts x and y to mpf or mpc and returns x**y = exp(y*log(x))."""
    return convert_lossless(x) ** convert_lossless(y)

def modf(x,y):
    """Converts x and y to mpf or mpc and returns x % y"""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return x % y

def degrees(x):
    """Convert x given in radians to degrees"""
    return x / degree

def radians(x):
    """Convert x given in degrees to radians"""
    return x * degree

def atan2(y,x):
    """atan2(y, x) has the same magnitude as atan(y/x) but accounts for
    the signs of y and x. (Defined for real x and y only.)"""
    x = convert_lossless(x)
    y = convert_lossless(y)
    return make_mpf(fatan2(y._mpf_, x._mpf_, *prec_rounding))

def rand():
    """Return an mpf chosen randomly from [0, 1)."""
    return make_mpf(frand(mp.prec))

from operator import gt, lt

def arange(*args):
    """arange([a,] b[, dt]) -> list [a, a + dt, a + 2*dt, ..., b]"""
    if not len(args) <= 3:
        raise TypeError('arange expected at most 3 arguments, got %i'
                        % len(args))
    if not len(args) >= 1:
        raise TypeError('arange expected at least 1 argument, got %i'
                        % len(args))
    # set default
    a = 0
    dt = 1
    # interpret arguments
    if len(args) == 1:
        b = args[0]
    elif len(args) >= 2:
        a = args[0]
        b = args[1]
    if len(args) == 3:
        dt = args[2]
    a, b, dt = mpf(a), mpf(b), mpf(dt)
    assert a + dt != a, 'dt is too small and would cause an infinite loop'
    # adapt code for sign of dt
    if a > b:
        if dt > 0:
            return []
        op = gt
    else:
        if dt < 0:
            return []
        op = lt
    # create list
    result = []
    i = 0
    t = a
    while 1:
        t = a + dt*i
        i += 1
        ##print i, t, op(t, b)
        if op(t, b):
            result.append(t)
        else:
            break
    return result

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
        rel_eps = abs_eps = make_mpf((0, 1, -mp.prec+4, 1))
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
        return "(" + complex_to_str(x._mpc_[0], x._mpc_[1], n)  + ")"
    if isinstance(x, basestring):
        return repr(x)
    from matrices import matrix
    if isinstance(x, matrix):
        return x.__nstr__(n)
    return str(x)

def nprint(x, n=6):
    """Print the result of nstr(x, n)."""
    print nstr(x, n)


plot_ignore = (ValueError, ArithmeticError, ZeroDivisionError)

def plot(f, xlim=[-5,5], ylim=None, points=200, file=None):
    """
    Shows a simple 2D plot of a function or list of functions over
    a given interval. Some examples:

        plot(lambda x: exp(x)*li(x), [1, 4])
        plot([cos, sin], [-4, 4])
        plot([fresnels, fresnelc], [-4, 4])
        plot([sqrt, cbrt], [-4, 4])
        plot(lambda t: zeta(0.5+t*j), [-20, 20])
        plot([floor, ceil, abs, sign], [-5, 5])

    Points where the function raises a numerical exception or
    returns an infinite value are removed from the graph.

    For parts where the function assumes complex values, the
    real part is plotted with dashes and the imaginary part
    is plotted with dots.

    NOTE: This function requires matplotlib (pylab).
    """
    import pylab
    if not isinstance(f, (tuple, list)):
        f = [f]
    a, b = xlim
    colors = ['b', 'r', 'g', 'm', 'k']
    for n, func in enumerate(f):
        x = arange(a, b, (b-a)/float(points))
        segments = []
        segment = []
        in_complex = False
        for i in xrange(len(x)):
            try:
                v = func(x[i])
                if isnan(v) or abs(v) == inf:
                    raise ValueError
                if isinstance(v, complex_types):
                    re = float(v.real)
                    im = float(v.imag)
                    if not in_complex:
                        in_complex = True
                        segments.append(segment)
                        segment = []
                    segment.append((float(x[i]), re, im))
                else:
                    if in_complex:
                        in_complex = False
                        segments.append(segment)
                        segment = []
                    segment.append((float(x[i]), v))
            except plot_ignore:
                if segment:
                    segments.append(segment)
                segment = []
        if segment:
            segments.append(segment)
        for segment in segments:
            x = [s[0] for s in segment]
            y = [s[1] for s in segment]
            if not x:
                continue
            c = colors[n % len(colors)]
            if len(segment[0]) == 3:
                z = [s[2] for s in segment]
                pylab.plot(x, y, '--'+c, linewidth=1.5)
                pylab.plot(x, z, ':'+c, linewidth=1.5)
            else:
                pylab.plot(x, y, c, linewidth=1.5)
    pylab.xlim(xlim)
    if ylim:
        pylab.ylim(ylim)
    pylab.grid(True)
    if file:
        pylab.savefig(file)
    else:
        pylab.show()

from lib import *

from decimal import Context, getcontext, setcontext, Decimal

def _convert(x):
    if isinstance(x, float):
        return float_from_pyfloat(x, mpf._prec, mpf._rounding)
    if isinstance(x, (int, long)):
        return float_from_int(x, mpf._prec, mpf._rounding)
    if isinstance(x, (Decimal, str)):
        return mpf(decimal_to_binary(x, mpf._prec, mpf._rounding))

class context(type):
    _prec = 53
    _dps = 15
    _rounding = ROUND_HALF_EVEN

    def _setprec(self, n):
        self._prec = max(1, int(n))
        self._dps = max(1, int(round(int(n)/LOG2_10)-1))

    prec = property(lambda self: self._prec, _setprec)

    def _setdps(self, n):
        self._prec = max(1, int(round((int(n)+1)*LOG2_10)))
        self._dps = max(1, int(n))

    dps = property(lambda self: self._dps, _setdps)

    extraprec = lambda self, n: self._setprec(self._prec+n)
    extradps = lambda self, n: self._setdps(self._dps+n)


class mpf(tuple):

    __metaclass__ = context

    def __new__(cls, val=fzero):
        if val.__class__ is tuple:
            return tuple.__new__(cls, val)
        elif isinstance(val, mpf):
            return +val
        else:
            return tuple.__new__(cls, _convert(val))

    man = property(lambda self: self[0])
    exp = property(lambda self: self[1])
    bc = property(lambda self: self[2])

    def __repr__(s):
        st = "mpf('%s')"
        return st % binary_to_decimal(s, mpf._dps+2)

    def __str__(s):
        return binary_to_decimal(s, mpf._dps)

    def __cmp__(s, t):
        if not isinstance(t, mpf): t = mpf(t)
        return fcmp(s, t)

    def __eq__(s, t):
        if not isinstance(t, mpf): t = mpf(t)
        return s[:] == t[:]

    def __ne__(s, t):
        if not isinstance(t, mpf): t = mpf(t)
        return s[:] != t[:]

    def __lt__(s, t): return s.__cmp__(t) < 0
    def __le__(s, t): return s.__cmp__(t) <= 0
    def __gt__(s, t): return s.__cmp__(t) > 0
    def __ge__(s, t): return s.__cmp__(t) >= 0

    def __int__(s):
        return float_to_int(s)

    def __float__(s):
        return float_to_pyfloat(s)

    def __complex__(s):
        return float_to_pyfloat(s) + 0j

    def __abs__(s):
        return mpf(fabs(s, mpf._prec, mpf._rounding))

    def __pos__(s):
        return mpf(normalize(s[0], s[1], mpf._prec, mpf._rounding))

    def __neg__(s):
        return mpf(fneg(s, mpf._prec, mpf._rounding))

    def __add__(s, t):
        if not isinstance(t, mpf): t = mpf(t)
        return mpf(fadd(s, t, mpf._prec, mpf._rounding))

    __radd__ = __add__

    def __sub__(s, t):
        if not isinstance(t, mpf): t = mpf(t)
        return mpf(fsub(s, t, mpf._prec, mpf._rounding))

    def __rsub__(s, t):
        if not isinstance(t, mpf): t = mpf(t)
        return mpf(fsub(t, s, mpf._prec, mpf._rounding))

    def __mul__(s, t):
        if not isinstance(t, mpf): t = mpf(t)
        return mpf(fmul(s, t, mpf._prec, mpf._rounding))

    def __div__(s, t):
        if not isinstance(t, mpf): t = mpf(t)
        return mpf(fdiv(s, t, mpf._prec, mpf._rounding))

    def __rdiv__(s, t):
        if not isinstance(t, mpf): t = mpf(t)
        return mpf(fdiv(t, s, mpf._prec, mpf._rounding))

    def __pow__(s, t):
        if isinstance(t, (int, long)):
            return mpf(fpow(s, t, mpf._prec, mpf._rounding))
        if not isinstance(t, mpf): t = mpf(t)
        if t == 0.5:
            return s.sqrt()
        intt = int(t)
        if t == intt:
            return mpf(fpow(s, intt, mpf._prec, mpf._rounding))

    def sqrt(s): return mpf(fsqrt(s, mpf._prec, mpf._rounding))
    #def exp(s):  return mpf(fexp(s, mpf._prec, mpf._rounding))
    def log(s):  return mpf(flog(s, mpf._prec, mpf._rounding))
    def sin(s):  return mpf(fsin(s, mpf._prec, mpf._rounding))
    def cos(s):  return mpf(fcos(s, mpf._prec, mpf._rounding))
    def tan(s):  return mpf(ftan(s, mpf._prec, mpf._rounding))

    @classmethod
    def const_pi(cls):
        """Return pi as an mpf rounded to the current working precision"""
        return mpf(fpi(cls._prec, cls._rounding))

    @classmethod
    def const_gamma(cls):
        """Return Euler's constant gamma as an mpf rounded to the
        current working precision"""
        return mpf(fgamma(cls._prec, cls._rounding))

    @classmethod
    def const_log2(cls):
        """Return log(2) as an mpf rounded to the current working
        precision"""
        return mpf(flog2(cls._prec, cls._rounding))

    @classmethod
    def const_log10(cls):
        """Return log(10) as an mpf rounded to the current working
        precision"""
        return mpf(flog10(cls._prec, cls._rounding))

def hypot(x, y):
    x = mpf(x)
    y = mpf(y)
    return mpf(fhypot(x, y, mpf._prec, mpf._rounding))


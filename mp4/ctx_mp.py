"""
Defines the multiprecision multiprecision floating-point context
with multiprecision low-level implementations of core functions.

"""

# Used for the low-level implementation
from mpf_impl import (V_cache, V_from, V_to_str, V_to_float, V_to_complex, V_hash,
    V_add, V_sub, V_mul, V_sqr, V_sum, V_dot, V_div, V_sqrt, V_recip, V_pow, V_0,
    V_re, V_im, V_cmp, V_sign, V_exp, V_pi, V_nthroot, V_ln, V_arg, V_mag,
    V_nint_distance, V_from_rational, V_agm,
    V_gamma_upper_int, V_expint_int, V_agm, V_hypot, V_atan2, V_psi,
    V_bernoulli, V_besselj,
    V_from_str,
    V_e, V_ln2, V_ln10, V_phi, V_degree, V_euler, V_catalan, V_glaisher,
    V_apery, V_khinchin, V_twinprime, V_mertens,
    V_to_mpf, V_to_mpc, V_from_mpf, V_from_mpc,
    V_mod,
    V_inf, V_ninf, V_ninfj, V_infj, V_nan)

from mpf_impl import S_HAS_NAN, S_HAS_INF, MPZ
from mpf_impl import inttypes

import rational

from mpmath.settings import repr_dps
from mpmath.settings import prec_to_dps
from mpmath.settings import dps_to_prec
from mpmath.mptypes import ComplexResult
from mpmath.libmpf import fzero
import mpmath.libmpf as lmpf
import mpmath.libmpc as lmpc
import mpmath.libelefun as lelf
import mpmath.gammazeta as lgz
import mpmath.libhyper as lhyp
from mpf_impl import i_trim, i_add, i_sqrt, MPZ_0, MPZ_1, S_NORMAL, from_man_exp

from ctx_base import StandardBaseContext
import function_docs
from hyp import make_hyp_summator

hyp_summators = {}

def V_from_rounded(x, prec, strings=False):
    if isinstance(x, rational.mpq):
        p, q = x
        return V_from_rational(p, q, prec)
    if strings and isinstance(x, basestring):
        return V_from_str(x, prec)
    return NotImplemented


def binop(name, commutative=True):
    code = """\
def __FUN__(x, y):
    _new, _type, _prec = x._info
    xval = x._v
    if hasattr(y, "_v"):
        yval = y._v
    elif y in V_cache:
        yval = V_cache[y]
    else:
        try:
            yval = V_from(y)
        except:
            yval = V_from_rounded(y, _prec)
            if yval is NotImplemented:
                return yval
    z = _new(_type)
    z._v = V_OP(%s, %s, _prec)
    return z
"""
    namespace = {}
    code = code.replace("OP", name)
    code1 = code.replace("FUN", name) % ("xval", "yval")
    if commutative:
        exec code1 in globals(), namespace
        lfun = rfun = namespace["__"+name+"__"]
    else:
        code2 = code.replace("FUN", "r"+name) % ("yval", "xval")
        exec code1 in globals(), namespace
        exec code2 in globals(), namespace
        lfun = namespace["__"+name+"__"]
        rfun = namespace["__r"+name+"__"]
    return lfun, rfun

def cmpop(name, action):
    code = """\
def FUN(x, y):
    if hasattr(y, "_v"):
        yval = y._v
    else:
        try:
            yval = V_from(y)
        except:
            _prec = x._info[-1]
            yval = V_from_rounded(y, _prec)
            if yval is NotImplemented:
                return yval
    return ACTION
"""
    code = code.replace("FUN", name).replace("ACTION", action%("x._v", "yval"))
    namespace = {}
    exec code in globals(), namespace
    return namespace[name]

class mpf(object):
    """
    Represents a complex floating-point number `(a+bi)` where both `a` and
    `b` have the form `m 2^n` for integers `m` and `n`. Arithmetic operations
    are rounded to the user-set working precision. Can also represent one of
    several special values (including infinities and indeterminate values).

    Some operations that exactly preserve number parts (possibly changing
    signs) are performed exactly without rounding. These operations include
    -x, x+0, x-0, x.real, x.imag, and abs(x) for pure real or pure
    imaginary x. A consequence is that for example -(-x) == x holds, even
    if x was generated at a higher precision.
    """

    __slots__ = ["_v"]

    def __new__(cls, value, value2=None):
        if value2 is not None:
            return cls(value) + cls(1j)*cls(value2)
        # XXX: faster lookup, avoid isinstance
        if type(value) is cls:
            return value
        self = object.__new__(cls)
        try:
            self._v = V_from(value)
        except:
            if hasattr(value, "_v"):
                v = value._v
            else:
                v = V_from_rounded(value, self._info[-1], strings=True)
                if v is NotImplemented:
                    raise ValueError
            self._v = v
        return self

    # -----------------------------------------------------------------
    # Conversion methods
    # -----------------------------------------------------------------

    def __repr__(self):
        if self._ctx.pretty:
            return self.__str__()
        p = self._info[-1]
        return "mpf('%s')" % V_to_str(self._v, repr_dps(p))

    def __str__(self):
        p = self._info[-1]
        return V_to_str(self._v, prec_to_dps(p))

    def __int__(x):
        am, ae, bm, be, special = x._v
        if bm or special:
            raise ValueError
        if ae >= 0:
            return int(am << ae)
        if am < 0:
            return int(-((-am) >> (-ae)))
        return int(am >> (-ae))

    __long__ = __int__

    def __hash__(x): return V_hash(x._v)
    def __float__(x): return V_to_float(x._v)
    def __complex__(x): return V_to_complex(x._v)
    def __nonzero__(x): return x._v != V_0

    __eq__ = cmpop("__eq__", "%s == %s")
    __ne__ = cmpop("__ne__", "%s != %s")
    __le__ = cmpop("__le__", "V_cmp(%s, %s) <= 0")
    __lt__ = cmpop("__lt__", "V_cmp(%s, %s) < 0")
    __ge__ = cmpop("__ge__", "V_cmp(%s, %s) >= 0")
    __gt__ = cmpop("__gt__", "V_cmp(%s, %s) > 0")
    __cmp__ = cmpop("__cmp__", "V_cmp(%s, %s)")

    __add__, __radd__ = binop("add")
    __sub__, __rsub__ = binop("sub", False)
    __mul__, __rmul__ = binop("mul")
    __div__, __rdiv__ = binop("div", False)
    __pow__, __rpow__ = binop("pow", False)
    __mod__, __rmod__ = binop("mod", False)
    __truediv__, __rtruediv__ = __div__, __rdiv__

    def ae(x, y, *args, **kwargs):
        return x._ctx.almosteq(x, y, *args, **kwargs)

    def is_special(self):
        return bool(self._v[-1])

    def to_fixed(x, prec):
        am, ae, bm, be, special = x._v
        if special or bm:
            raise ValueError
        shift = ae + prec
        if shift >= 0:
            return am << shift
        else:
            return am >> (-shift)

    # neg, conjugate, abs and related methods inlined for speed
    # The main speedup does not come from avoiding a function call,
    # but from being able to simply return the same object when
    # the value is e.g. already positive real (a very common case)

    # XXX
    def __pos__(x):
        _new, _type, _prec = x._info
        am, ae, bm, be, special = x._v
        #if special:
        #    #raise NotImplementedError
        #    return x
        if am: am, ae = i_trim(am, ae, _prec, 'n')
        if bm: bm, be = i_trim(bm, be, _prec, 'n')
        z = _new(_type)
        z._v = am, ae, bm, be, special
        return z

    def __neg__(x):
        _new, _type, _prec = x._info
        am, ae, bm, be, special = x._v
        z = _new(_type)
        z._v = -am, ae, -bm, be, special
        return z

    @property
    def real(x):
        v = x._v
        am, ae, bm, be, special = v
        if special:
            if v == V_inf: return x._ctx.inf
            if v == V_ninf: return x._ctx.ninf
            if v == V_nan: return x._ctx.nan
        if bm:
            if am:
                _new, _type, _prec = x._info
                z = _new(_type)
                z._v = am, ae, MPZ_0, 0, special
                return z
            return x._zero
        return x

    @property
    def imag(x):
        am, ae, bm, be, special = v = x._v
        if special:
            if v == V_infj: return x._ctx.infj
            if v == V_ninfj: return x._ctx.ninfj
            if v == V_nan: return x._ctx.nan
        if bm:
            _new, _type, _prec = x._info
            z = _new(_type)
            z._v = bm, be, MPZ_0, 0, special
            return z
        return x._zero

    def conjugate(x):
        am, ae, bm, be, special = x._v
        if bm:
            _new, _type, _prec = x._info
            z = _new(_type)
            z._v = am, ae, -bm, be, special
            return z
        return x

    def __abs__(x):
        am, ae, bm, be, special = x._v
        if special:
            if special & S_HAS_NAN: return x._ctx.nan
            if special & S_HAS_INF: return x._ctx.inf
            raise NotImplementedError
        if bm:
            if am:
                _new, _type, _prec = x._info
                mm, me = i_add(am*am, ae+ae, bm*bm, be+be, _prec+20, 'u')
                mm, me = i_sqrt(mm, me, _prec, 'n')
                z = _new(_type)
                z._v = mm, me, MPZ_0, 0, S_NORMAL
                return z
            v = abs(bm), be, MPZ_0, 0, S_NORMAL
        else:
            if am >= 0:
                return x
            v = -am, ae, MPZ_0, 0, S_NORMAL
        _new, _type, _prec = x._info
        z = _new(_type)
        z._v = v
        return z

    def __invert__(x):
        _new, _type, _prec = x._info
        z = _new(_type)
        z._v = V_recip(x._v, _prec)
        return z

    def sqrt(x):
        _new, _type, _prec = x._info
        z = _new(_type)
        z._v = V_sqrt(x._v, _prec)
        return z

class constant(mpf):

    def __new__(cls, v_fun):
        x = object.__new__(cls)
        x._f = v_fun
        return x

    # XXX
    def __call__(x, prec=None, dps=None):
        _new, _type, _prec = x._info
        if prec is None:
            prec = _prec
        if dps is not None:
            prec = dps_to_prec(dps)
        z = _new(_type)
        z._v = x._f(prec)
        return z

    def get_v(x):
        return x._f(x._info[-1])

    _v = property(get_v)

class MPContext(StandardBaseContext):

    def __init__(ctx):

        # XXX
        ctx.mpf = type('mpf', (mpf,), {})
        ctx.constant = type('constant', (constant,), {}) # XXX
        ctx.convert = ctx.mpc = ctx.mpf  # XXX
        ctx.mpf._ctx = ctx
        ctx.constant._ctx = ctx
        # note: must all be linked
        ctx.constant._info = ctx.mpf._info = ctx._info = [object.__new__, ctx.mpf, 53]

        StandardBaseContext.__init__(ctx)
        ctx.pretty = False
        ctx.set_default_prec()
        ctx.mpq = rational.mpq

        #SpecialFunctions.__init__(ctx)
        #QuadratureMethods.__init__(ctx)
        #CalculusMethods.__init__(ctx)
        #MatrixMethods.__init__(ctx)

        ctx.init_builtins()


    def init_builtins(ctx):

        mpf = ctx.mpf

        # Exact constants
        ctx.j = mpf(1j)
        ctx.one = mpf(1)
        ctx.zero = mpf._zero = mpf(0)
        ctx.inf = mpf(1e300*1e300) # XXX
        ctx.ninf = mpf(-1e300*1e300)
        ctx.nan = mpf(1e300*1e300 - 1e300*1e300)

        # Approximate constants
        ctx.pi = ctx.constant(V_pi)
        ctx.ln2 = ctx.constant(V_ln2)
        ctx.ln10 = ctx.constant(V_ln10)
        ctx.phi = ctx.constant(V_phi)
        ctx.e = ctx.constant(V_e)
        ctx.euler = ctx.constant(V_euler)
        ctx.catalan = ctx.constant(V_catalan)
        ctx.khinchin = ctx.constant(V_khinchin)
        ctx.glaisher = ctx.constant(V_glaisher)
        ctx.apery = ctx.constant(V_apery)
        ctx.degree = ctx.constant(V_degree)
        ctx.twinprime = ctx.constant(V_twinprime)
        ctx.mertens = ctx.constant(V_mertens)

        # XXX
        ctx.nthroot = ctx.root
        ctx.digamma = ctx.psi0
        ctx.trigamma = ctx.psi1
        ctx.tetragamma = ctx.psi2
        ctx.pentagamma = ctx.psi3
        ctx.bernfrac = lgz.bernfrac

        # Standard functions with new implementations
        ctx.sqrt = ctx.def_V_function(V_sqrt)
        ctx.re = ctx.def_V_function(V_re)
        ctx.im = ctx.def_V_function(V_im)
        ctx.sign = ctx.def_V_function(V_sign)

        ctx._nthroot = ctx.def_V_function_n(V_nthroot)
        ctx._gamma_upper_int = ctx.def_V_n_function(V_gamma_upper_int)
        ctx._expint_int = ctx.def_V_n_function(V_expint_int)
        ctx._agm = ctx.def_V_function_2(V_agm)
        ctx.hypot = ctx.def_V_function_2(V_hypot)
        ctx.atan2 = ctx.def_V_function_2(V_atan2)
        ctx._besselj = ctx.def_V_n_function(V_besselj)
        ctx.psi = ctx.polygamma = ctx.def_V_n_function(V_psi)
        ctx.bernoulli = ctx.def_V_n(V_bernoulli)

        ctx.exp2 = ctx.def_V_function(V_exp)

        # Standard functions
        ctx.cbrt = ctx.def_mp_function(lelf.mpf_cbrt, lmpc.mpc_cbrt)
        ctx.ln = ctx.def_mp_function(lelf.mpf_log, lmpc.mpc_log)
        ctx.atan = ctx.def_mp_function(lelf.mpf_atan, lmpc.mpc_atan)
        ctx.exp = ctx.def_mp_function(lelf.mpf_exp, lmpc.mpc_exp)
        ctx.sin = ctx.def_mp_function(lelf.mpf_sin, lmpc.mpc_sin)
        ctx.cos = ctx.def_mp_function(lelf.mpf_cos, lmpc.mpc_cos)
        ctx.tan = ctx.def_mp_function(lelf.mpf_tan, lmpc.mpc_tan)
        ctx.sinh = ctx.def_mp_function(lelf.mpf_sinh, lmpc.mpc_sinh)
        ctx.cosh = ctx.def_mp_function(lelf.mpf_cosh, lmpc.mpc_cosh)
        ctx.tanh = ctx.def_mp_function(lelf.mpf_tanh, lmpc.mpc_tanh)
        ctx.asin = ctx.def_mp_function(lelf.mpf_asin, lmpc.mpc_asin)
        ctx.acos = ctx.def_mp_function(lelf.mpf_acos, lmpc.mpc_acos)
        ctx.atan = ctx.def_mp_function(lelf.mpf_atan, lmpc.mpc_atan)
        ctx.asinh = ctx.def_mp_function(lelf.mpf_asinh, lmpc.mpc_asinh)
        ctx.acosh = ctx.def_mp_function(lelf.mpf_acosh, lmpc.mpc_acosh)
        ctx.atanh = ctx.def_mp_function(lelf.mpf_atanh, lmpc.mpc_atanh)
        ctx.sinpi = ctx.def_mp_function(lelf.mpf_sin_pi, lmpc.mpc_sin_pi)
        ctx.cospi = ctx.def_mp_function(lelf.mpf_cos_pi, lmpc.mpc_cos_pi)
        ctx.floor = ctx.def_mp_function(lmpf.mpf_floor, lmpc.mpc_floor)
        ctx.ceil = ctx.def_mp_function(lmpf.mpf_ceil, lmpc.mpc_ceil)
        ctx.fib = ctx.fibonacci = ctx.def_mp_function(lelf.mpf_fibonacci, lmpc.mpc_fibonacci)
        ctx.zeta = ctx.def_mp_function(lgz.mpf_zeta, lgz.mpc_zeta)
        ctx.altzeta = ctx.def_mp_function(lgz.mpf_altzeta, lgz.mpc_altzeta)
        ctx.gamma = ctx.def_mp_function(lgz.mpf_gamma, lgz.mpc_gamma)
        ctx.fac = ctx.factorial = ctx.def_mp_function(lgz.mpf_factorial, lgz.mpc_factorial)
        ctx.harmonic = ctx.def_mp_function(lgz.mpf_harmonic, lgz.mpc_harmonic)
        ctx.ei = ctx.def_mp_function(lhyp.mpf_ei, lhyp.mpc_ei)
        ctx.e1 = ctx.def_mp_function(lhyp.mpf_e1, lhyp.mpc_e1)
        ctx.ci = ctx.def_mp_function(lhyp.mpf_ci, lhyp.mpc_ci)
        ctx.si = ctx.def_mp_function(lhyp.mpf_si, lhyp.mpc_si)
        ctx.ellipk = ctx.def_mp_function(lhyp.mpf_ellipk, lhyp.mpc_ellipk)
        ctx.ellipe = ctx.def_mp_function(lhyp.mpf_ellipe, lhyp.mpc_ellipe)
        ctx.agm1 = ctx.def_mp_function(lhyp.mpf_agm1, lhyp.mpc_agm1)
        ctx._erf = ctx.def_mp_function(lhyp.mpf_erf, None)
        ctx._erfc = ctx.def_mp_function(lhyp.mpf_erfc, None)

        ctx.lnn = ctx.def_V_function(V_ln)
        ctx.arg = ctx.def_V_function(V_arg)

        ctx.absmin = ctx.absmax = abs

    def to_fixed(ctx, x, prec):
        return x.to_fixed(prec)

    def rand(ctx):
        _new, _type, _prec = ctx._info
        z = _new(_type)
        z._v = V_from_mpf(lmpf.mpf_rand(_prec))
        return z

    def fraction(ctx, p, q):
        return ctx.constant(lambda prec: V_from_rational(p, q, prec))

    def mpf_or_rational(ctx, x):
        hasx = hasattr(x, "_v")
        if hasx:
            v = x._v
        else:
            if type(x) in inttypes:
                return int(x), 'Z'
            if isinstance(x, tuple):
                p, q = x
                return ctx.mpq((p,q)), 'Q'
            if isinstance(x, basestring) and '/' in x:
                p, q = x.split('/')
                return ctx.mpq((int(p), int(q))), 'Q'
            v = V_from(x)
        a, b, c, d, special = v
        if not (c or special):
            if b >= 0:
                return int(a) << b, 'Z'
            if b >= -4:
                p, q = int(a), (1<<(-b))
                return ctx.mpq((p,q)), 'Q'
        flag = 'R'
        if c:
            flag = 'C'
        if hasx:
            return x, flag
        _new, _type, _prec = ctx._info
        x = _new(_type)
        x._v = v
        return x, flag

    def hypsum(ctx, p, q, flags, coeffs, z, **kwargs):
        v = z._v
        if v[2]:
            key = p, q, flags, 'C'
        else:
            key = p, q, flags, 'R'
        if key not in hyp_summators:
            hyp_summators[key] = make_hyp_summator(key)[1]
        _new, _type, _prec = ctx._info
        zv = hyp_summators[key](coeffs, v, _prec, **kwargs)
        z = _new(_type)
        z._v = zv
        return z

    def mag(ctx, z):
        if hasattr(z, "_v"):
            v = z._v
        else:
            v = ctx.convert(z)._v
        return V_mag(v)

    def is_real_type(ctx, z):
        # XXX
        if hasattr(z, "_v"):
            return not z._v[2]
        return True

    def is_complex_type(ctx, z):
        # XXX
        if hasattr(z, "_v"):
            return bool(z._v[2])
        raise NotImplementedError

    def nint_distance(ctx, z):
        if hasattr(z, "_v"):
            v = z._v
        else:
            v = ctx.convert(z)._v
        return V_nint_distance(v)

    def _V_from(ctx, x):
        try:
            v = V_from(x)
        except:
            v = V_from_rounded(x, ctx._info[-1]+10, True)
            if v is NotImplemented:
                raise NotImplementedError
        return v

    def def_V_function(ctx, Vf):
        _info = ctx._info
        def f(x, **kwargs):
            _new, _type, _prec = _info
            rounding = None
            if kwargs:
                _prec = kwargs.get("prec", _prec)
                if "rounding" in kwargs:
                    rounding = kwargs["rounding"]
            if hasattr(x, "_v"):
                v = x._v
            else:
                v = ctx._V_from(x)
            z = _new(_type)
            if rounding:
                z._v = Vf(v, _prec, rounding=rounding)
            else:
                z._v = Vf(v, _prec)
            return z
        return f

    def def_V_n(ctx, Vf):
        _info = ctx._info
        def f(n):
            _new, _type, _prec = _info
            z = _new(_type)
            z._v = Vf(int(n), _prec)
            return z
        return f

    def def_V_function_n(ctx, Vf):
        _info = ctx._info
        def f(x, n):
            _new, _type, _prec = _info
            n = int(n)
            if hasattr(x, "_v"):
                v = x._v
            else:
                v = ctx._V_from(x)
            z = _new(_type)
            z._v = Vf(v, n, _prec)
            return z
        return f

    def def_V_n_function(ctx, Vf):
        _info = ctx._info
        def f(n, x):
            _new, _type, _prec = _info
            n = int(n)
            if hasattr(x, "_v"):
                v = x._v
            else:
                v = ctx._V_from(x)
            z = _new(_type)
            z._v = Vf(n, v, _prec)
            return z
        return f

    def def_V_function_2(ctx, Vf):
        _info = ctx._info
        def f(x, y):
            _new, _type, _prec = _info
            if hasattr(x, "_v"):
                v = x._v
            else:
                v = ctx._V_from(x)
            if hasattr(y, "_v"):
                w = y._v
            else:
                w = ctx._V_from(y)
            z = _new(_type)
            z._v = Vf(v, w, _prec)
            return z
        return f

    def def_mp_function(ctx, mpf_f, mpc_f):
        _info = ctx._info
        def f(x, **kwargs):
            _new, _type, _prec = _info
            if kwargs:
                _prec = kwargs.get("prec", _prec)
                if 'dps' in kwargs:
                    _prec = dps_to_prec(kwargs['dps'])
                rounding = kwargs.get("rounding", 'n')
            else:
                rounding = 'n'
            if hasattr(x, "_v"):
                v = x._v
            else:
                v = ctx._V_from(x)
            am, ae, bm, be, special = v
            if bm:
                if special:
                    fv = V_to_mpc(v)
                else:
                    fv = from_man_exp(am, ae), from_man_exp(bm, be)
                re, im = mpc_f(fv, _prec, rounding)
                rsign, rman, rexp, rbc = re
                isign, iman, iexp, ibc = im
                if ((not rman) and rexp) or ((not iman) and iexp):
                    raise NotImplementedError
                if rsign: rman = -rman
                if isign: iman = -iman
                v = rman, rexp, iman, iexp, S_NORMAL
            else:
                try:
                    if special:
                        fv = V_to_mpf(v)
                    else:
                        fv = from_man_exp(am, ae)
                    sign, man, exp, bc = mpf_f(fv, _prec, rounding)
                    if (not man) and exp:
                        v = V_from_mpf((sign,man,exp,bc))
                    elif sign:
                        v = -man, exp, MPZ_0, 0, S_NORMAL
                    else:
                        v = man, exp, MPZ_0, 0, S_NORMAL
                except ComplexResult:
                    re, im = mpc_f((from_man_exp(am, ae), fzero), _prec, rounding)
                    rsign, rman, rexp, rbc = re
                    isign, iman, iexp, ibc = im
                    if ((not rman) and rexp) or ((not iman) and iexp):
                        v = V_from_mpc((re, im))
                    else:
                        if rsign: rman = -rman
                        if isign: iman = -iman
                    v = rman, rexp, iman, iexp, S_NORMAL
            z = _new(_type)
            z._v = v
            return z
        return f

    # Called by SpecialFunctions.__init__()
    @classmethod
    def wrap_specfun(cls, name, f, wrap):
        if wrap:
            def f_wrapped(ctx, *args, **kwargs):
                convert = ctx.convert
                args = [convert(a) for a in args]
                prec = ctx.prec
                try:
                    ctx.prec += 10
                    retval = f(ctx, *args, **kwargs)
                finally:
                    ctx.prec = prec
                return +retval
        else:
            f_wrapped = f

        f_wrapped.__doc__ = function_docs.__dict__.get(name, "<no doc>")

        setattr(cls, name, f_wrapped)

    def bad_domain(ctx, msg):
        raise ValueError(msg)

    def AS_POINTS(ctx, x):
        return x

    def fprod(ctx, factors):
        orig = ctx.prec
        try:
            v = ctx.one
            for p in factors:
                v *= p
        finally:
            ctx.prec = orig
        return +v

    def fmul(ctx, x, y, prec=None, exact=False):
        _new, _type, _prec = ctx.mpf._info
        if exact:
            prec = 0
        elif prec is None:
            prec = _prec
        v = V_mul(ctx.convert(x)._v, ctx.convert(y)._v, prec)
        z = _new(_type)
        z._v = v
        return z

    def isnan(ctx, x):
        return bool(ctx.convert(x)._v[-1] & S_HAS_NAN)

    def isinf(ctx, x):
        return bool(ctx.convert(x)._v[-1] & S_HAS_INF)

    def isint(ctx, x):
        if hasattr(x, "_v"):
            a, b, c, d, special = x._v
            return (not (c or special)) and b >= 0
        if type(x) in inttypes:
            return True
        if isinstance(x, ctx.mpq):
            # XXX: WRONG
            p, q = x
            if not p:
                return True
            if p == 1 or p == -1:
                return q == 1
            return not (q % p)
        print "MOO", type(x), x
        raise NotImplementedError

    def isnpint(ctx, x):
        if hasattr(x, "_v"):
            a, b, c, d, special = x._v
            return (not (c or special)) and b >= 0 and a <= 0
        if type(x) in inttypes:
            return x <= 0
        if isinstance(x, ctx.mpq):
            # XXX: WRONG
            p, q = x
            if not p:
                return True
            return (not (q % p)) and p <= 0
        raise NotImplementedError

    def is_special(ctx, x):
        return x.is_special()

    def set_default_prec(ctx):
        ctx._prec = ctx.mpf._info[-1] = 53
        ctx._dps = 15

    def set_prec(ctx, n):
        ctx._prec = ctx.mpf._info[-1] = max(1, int(n))
        ctx._dps = prec_to_dps(n)

    def set_dps(ctx, n):
        ctx._prec = ctx.mpf._info[-1] = dps_to_prec(n)
        ctx._dps = max(1, int(n))

    prec = property(lambda ctx: ctx._prec, set_prec)
    dps = property(lambda ctx: ctx._dps, set_dps)

    def get_eps(ctx):
        _new, _type, _prec = ctx.mpf._info
        z = _new(_type)
        z._v = MPZ_1, 1-_prec, MPZ_0, 0, S_NORMAL
        return z

    eps = property(get_eps)

    def ldexp(ctx, x, n):
        x = ctx.convert(x)
        _new, _type, _prec = x._info
        am, ae, bm, be, special = x._v
        if special:
            return x
        z = _new(_type)
        if am: ae += n
        if bm: be += n
        z._v = am, ae, bm, be, special
        return z

    def frexp(ctx, x):
        x = ctx.convert(x)
        y, n = lmpf.mpf_frexp(V_to_mpf(x._v))
        _new, _type, _prec = x._info
        z = _new(_type)
        z._v = V_from_mpf(y)
        return z, n

    def fdot(ctx, xs, ys=None):
        _new, _type, _prec = ctx.mpf._info
        _hasattr = hasattr
        _from = V_from
        if ys is None:
            xs = [((x._v if _hasattr(x, "_v") else _from(x)), \
                   (y._v if _hasattr(y, "_v") else _from(y))) \
                for x, y in xs]
            val = V_dot(xs, _prec)
        else:
            xs = [x._v if _hasattr(x, "_v") else _from(x) for x in xs]
            ys = [y._v if _hasattr(y, "_v") else _from(y) for y in ys]
            val = V_dot(zip(xs, ys), _prec)
        z = _new(_type)
        z._v = val
        return z

    def fsum(ctx, xs, absolute=False, squared=False):
        _new, _type, _prec = ctx.mpf._info
        _hasattr = hasattr
        _from = V_from
        if squared:
            # XXX!
            val = V_sum((V_sqr(x._v) if _hasattr(x, "_v") else V_sqr(_from(x)) for x in xs), _prec, absolute)
        else:
            val = V_sum((x._v if _hasattr(x, "_v") else _from(x) for x in xs), _prec, absolute)
        z = _new(_type)
        z._v = val
        return z

    def nstr(ctx, x, n=6, **kwargs):
        if isinstance(x, list):
            return "[%s]" % (", ".join(ctx.nstr(c, n) for c in x))
        if isinstance(x, tuple):
            return "(%s)" % (", ".join(ctx.nstr(c, n) for c in x))
        if hasattr(x, '_v'):
            return V_to_str(x._v, n)
        if isinstance(x, basestring):
            return repr(x)
        if isinstance(x, ctx.matrix):
            return x.__nstr__(n)
        #if hasattr(x, '_mpi_'):
        #    return ctx.mpi_to_str(x, n, **kwargs)
        return str(x)

    def nprint(ctx, x, n=6, **kwargs):
        print ctx.nstr(x, n, **kwargs)



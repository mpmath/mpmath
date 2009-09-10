from ctx_base import StandardBaseContext

import math
import cmath
import math2

from gammazeta import mpf_bernoulli
from libmpf import to_float

# XXX
from hyp import NoConvergence

class FPContext(StandardBaseContext):
    """
    Context for fast low-precision arithmetic (53-bit precision, giving at most
    about 15-digit accuracy), using Python's builtin float and complex.
    """

    def __init__(ctx):
        StandardBaseContext.__init__(ctx)

        # Override SpecialFunctions implementation
        ctx.loggamma = math2.loggamma

        ctx._bernoulli_cache = {}

    def _get_prec(ctx): return 53
    def _set_prec(ctx, p): return
    def _get_dps(ctx): return 15
    def _set_dps(ctx, p): return

    _fixed_precision = True

    prec = property(_get_prec, _set_prec)
    dps = property(_get_dps, _set_dps)

    zero = 0.0
    one = 1.0
    eps = math2.EPS
    inf = math2.INF
    ninf = math2.NINF
    nan = math2.NAN
    j = 1j

    def bernoulli(ctx, n):
        cache = ctx._bernoulli_cache
        if n in cache:
            return cache[n]
        cache[n] = to_float(mpf_bernoulli(n, 53, 'n'), strict=True)
        return cache[n]

    pi = math2.pi
    e = math2.e
    euler = math2.euler

    absmin = absmax = abs

    def AS_POINTS(ctx, x):
        return x

    def fsum(ctx, args, absolute=False, squared=False):
        if absolute:
            if squared:
                return sum((abs(x)**2 for x in args), ctx.zero)
            return sum((abs(x) for x in args), ctx.zero)
        if squared:
            return sum((x**2 for x in args), ctx.zero)
        return sum((x for x in args), ctx.zero)

    def fdot(ctx, xs, ys):
        if ys:
            xs = zip(xs, ys)
        return sum((x*y for (x,y) in xs), ctx.zero)

    def is_special(ctx, x):
        return x - x != 0.0

    def isnan(ctx, x):
        return x != x

    def isinf(ctx, x):
        return abs(x) == math2.INF

    def isnpint(ctx, x):
        if type(x) is complex:
            if x.imag:
                return False
            x = x.real
        return x <= 0.0 and round(x) == x

    def mpf(ctx, x):
        try:
            return float(x)
        except:
            return complex(x)

    mpc = complex

    convert = mpf

    power = staticmethod(math2.pow)
    sqrt = staticmethod(math2.sqrt)
    exp = staticmethod(math2.exp)
    ln = log = staticmethod(math2.log)
    cos = staticmethod(math2.cos)
    sin = staticmethod(math2.sin)
    tan = staticmethod(math2.tan)
    cos_sin = staticmethod(math2.cos_sin)
    acos = staticmethod(math2.acos)
    asin = staticmethod(math2.asin)
    atan = staticmethod(math2.atan)
    cosh = staticmethod(math2.cosh)
    sinh = staticmethod(math2.sinh)
    tanh = staticmethod(math2.tanh)
    gamma = staticmethod(math2.gamma)

    # XXX: math2
    def arg(ctx, z):
        z = complex(z)
        return math.atan2(z.imag, z.real)

    def sinpi(ctx, x):
        return ctx.sin(ctx.pi*x)

    def re(ctx, x):
        if type(x) is float:
            return x
        return x.real

    def im(ctx, x):
        if type(x) is float:
            return 0.0
        return x.imag

    ldexp = math.ldexp
    frexp = math.frexp

    def mag(ctx, z):
        return ctx.frexp(abs(z))[1]

    def isint(ctx, z):
        if z.imag:
            return False
        z = z.real
        try:
            return z == int(z)
        except:
            return False

    def nint_distance(ctx, z):
        n = round(z)
        return n, ctx.mag(abs(z-n))

    def mpf_or_rational(ctx, z):
        if z == int(z):
            return int(z), 'Z'
        return z, 'R'

    def is_real_type(ctx, z):
        return True

    def hypsum(ctx, p, q, types, coeffs, z, maxterms=6000): 
        s = t = 1.0
        k = 0
        coeffs = list(coeffs)
        num = range(p)
        den = range(p,p+q)
        while 1:
            for i in num: t *= (coeffs[i]+k)
            for i in den: t /= (coeffs[i]+k)
            k += 1
            t /= k
            t *= z
            s += t
            if abs(t) < ctx.eps:
                break
            if k > maxterms:
                raise NoConvergence
        return s


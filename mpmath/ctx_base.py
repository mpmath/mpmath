from operator import gt, lt

from functions import SpecialFunctions
from quadrature import QuadratureMethods
from calculus import CalculusMethods
from matrices import MatrixMethods
from linalg import LinearAlgebraMethods
from identification import IdentificationMethods
from optimization import OptimizationMethods
from odes import ODEMethods
from visualization import VisualizationMethods

class Context(object):
    pass

class StandardBaseContext(Context,
    SpecialFunctions,
    QuadratureMethods,
    CalculusMethods,
    MatrixMethods,
    LinearAlgebraMethods,
    IdentificationMethods,
    OptimizationMethods,
    ODEMethods,
    VisualizationMethods):

    def __init__(ctx):
        # Call those that need preinitialization (e.g. for wrappers)
        SpecialFunctions.__init__(ctx)
        QuadratureMethods.__init__(ctx)
        CalculusMethods.__init__(ctx)
        MatrixMethods.__init__(ctx)

    _fixed_precision = False

    # XXX
    verbose = False

    def warn(ctx, msg):
        print "Warning:", msg

    def bad_domain(ctx, msg):
        raise ValueError(msg)

    def chop(ctx, x, tol=None):
        if tol is None:
            tol = 100*ctx.eps
        try:
            x = ctx.convert(x)
            absx = abs(x)
            if abs(x) < tol:
                return ctx.zero
            if ctx.im(x):
                if abs(x.imag) < min(tol, absx*tol):
                    return x.real
                if abs(x.real) < min(tol, absx*tol):
                    return ctx.mpc(0, x.imag)
        except (TypeError, ValueError):
            if hasattr(x, "__iter__"):
                return [ctx.chop(a, tol) for a in x]
            raise
        return x

    def almosteq(ctx, s, t, rel_eps=None, abs_eps=None):
        t = ctx.convert(t)
        if abs_eps is None and rel_eps is None:
            rel_eps = abs_eps = ctx.ldexp(1, -ctx.prec+4)
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

    def arange(ctx, *args):
        r"""
        This is a generalized version of Python's :func:`range` function
        that accepts fractional endpoints and step sizes and
        returns a list of ``mpf`` instances. Like :func:`range`,
        :func:`arange` can be called with 1, 2 or 3 arguments:

        ``arange(b)``
            `[0, 1, 2, \ldots, x]`
        ``arange(a, b)``
            `[a, a+1, a+2, \ldots, x]`
        ``arange(a, b, h)``
            `[a, a+h, a+h, \ldots, x]`

        where `b-1 \le x < b` (in the third case, `b-h \le x < b`).

        Like Python's :func:`range`, the endpoint is not included. To
        produce ranges where the endpoint is included, :func:`linspace`
        is more convenient.

        **Examples**

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> arange(4)
            [mpf('0.0'), mpf('1.0'), mpf('2.0'), mpf('3.0')]
            >>> arange(1, 2, 0.25)
            [mpf('1.0'), mpf('1.25'), mpf('1.5'), mpf('1.75')]
            >>> arange(1, -1, -0.75)
            [mpf('1.0'), mpf('0.25'), mpf('-0.5')]

        """
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
        a, b, dt = ctx.mpf(a), ctx.mpf(b), ctx.mpf(dt)
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
            if op(t, b):
                result.append(t)
            else:
                break
        return result

    def linspace(ctx, *args, **kwargs):
        """
        ``linspace(a, b, n)`` returns a list of `n` evenly spaced
        samples from `a` to `b`. The syntax ``linspace(mpi(a,b), n)``
        is also valid.

        This function is often more convenient than :func:`arange`
        for partitioning an interval into subintervals, since
        the endpoint is included::

            >>> from mpmath import *
            >>> mp.dps = 15; mp.pretty = False
            >>> linspace(1, 4, 4)
            [mpf('1.0'), mpf('2.0'), mpf('3.0'), mpf('4.0')]
            >>> linspace(mpi(1,4), 4)
            [mpf('1.0'), mpf('2.0'), mpf('3.0'), mpf('4.0')]

        You may also provide the keyword argument ``endpoint=False``::

            >>> linspace(1, 4, 4, endpoint=False)
            [mpf('1.0'), mpf('1.75'), mpf('2.5'), mpf('3.25')]

        """
        if len(args) == 3:
            a = ctx.mpf(args[0])
            b = ctx.mpf(args[1])
            n = int(args[2])
        elif len(args) == 2:
            assert hasattr(args[0], '_mpi_')
            a = args[0].a
            b = args[0].b
            n = int(args[1])
        else:
            raise TypeError('linspace expected 2 or 3 arguments, got %i' \
                            % len(args))
        if n < 1:
            raise ValueError('n must be greater than 0')
        if not 'endpoint' in kwargs or kwargs['endpoint']:
            if n == 1:
                return [ctx.mpf(a)]
            step = (b - a) / ctx.mpf(n - 1)
            y = [i*step + a for i in xrange(n)]
            y[-1] = b
        else:
            step = (b - a) / ctx.mpf(n)
            y = [i*step + a for i in xrange(n)]
        return y

    def cos_sin(ctx, z, **kwargs):
        return ctx.cos(z, **kwargs), ctx.sin(z, **kwargs)

    def _default_hyper_maxprec(ctx, p):
        return int(1000 * p**0.25 + 4*p)

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

    def warn(ctx, msg):
        print "Warning:", msg

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
        import operator
        # adapt code for sign of dt
        if a > b:
            if dt > 0:
                return []
            op = operator.gt
        else:
            if dt < 0:
                return []
            op = operator.lt
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


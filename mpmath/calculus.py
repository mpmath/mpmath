"""
High-level calculus-oriented functions.

* Numerical differentiation
* Numerical polynomial operations
* Numerical root-finding
* Numerical integration

"""

__docformat__ = 'plaintext'

from mptypes import *
from specfun import factorial, bernoulli2n

#----------------------------------------------------------------------------#
#                                Differentiation                             #
#----------------------------------------------------------------------------#

def diff(f, x, direction=0):
    """
    Compute f'(x) using a simple finite difference approximation.

    With direction = 0, use the central difference f(x-h), f(x+h)
    With direction = 1, use the forward difference f(x), f(x+h)
    With direction = -1, use the backward difference f(x-h), f(x)

        >>> print diff(cos, 1)
        -0.841470984807897
        >>> print diff(abs, 0, 0)
        0.0
        >>> print diff(abs, 0, 1)
        1.0
        >>> print diff(abs, 0, -1)
        -1.0

    The step size is taken similar to the epsilon of the precision.
    To eliminate cancellation errors, diff temporarily doubles the
    working precision while calculating the function values.
    """
    prec = mp.prec
    extra = 5
    h = ldexp(1, -prec-extra)
    try:
        mp.prec = 2*(prec+extra)
        if   direction == 0:  return (f(x+h) - f(x-h)) * ldexp(1, prec+extra-1)
        elif direction == 1:  return (f(x+h) - f(x)) * ldexp(1, prec+extra)
        elif direction == -1: return (f(x) - f(x-h)) * ldexp(1, prec+extra)
        else:
            raise ValueError("invalid difference direction: %r" % direction)
    finally:
        mp.prec = prec

def diffc(f, x, n=1, radius=mpf(0.5)):
    """
    Compute an approximation of the nth derivative of f at the point x
    using the Cauchy integral formula. This only works for analytic
    functions. A circular path with the given radius is used.

    diffc increases the working precision slightly to avoid simple
    rounding errors. Note that, especially for large n, differentiation
    is extremely ill-conditioned, so this precaution does not
    guarantee a correct result. (Provided there are no singularities
    in the way, increasing the radius may help.)

    The returned value will be a complex number; a large imaginary part
    for a derivative that should be real may indicate a large numerical
    error.
    """
    prec = mp.prec
    try:
        mp.prec += 10
        def g(t):
            rei = radius*exp(j*t)
            z = x + rei
            return f(z) / rei**n
        d = quadts(g, 0, 2*pi)
        return d * factorial(n) / (2*pi)
    finally:
        mp.prec = prec


#----------------------------------------------------------------------------#
#                           Generic root-finding                             #
#----------------------------------------------------------------------------#

msg1 = "Cannot perform a step with the secant method because the " \
  "function values are equal at the two chosen start points. Try " \
  "different start points."

msg2 = "The derivative cannot be computed. The previous value " \
  "will be reused for the next iteration."

def secant(f, x0, x1=None, maxsteps=20, verbose=False):
    """Solve the equation f(x) = 0 using the secant method, starting
    at the given initial point x0 and performing up to `maxsteps`
    steps or quitting when the difference between successive x values
    is smaller than the epsilon of the current working precision.

    The secant method requires a second starting point x1 with both
    x0 and x1 located close to the root. If only x0 is provided, x1
    is automatically generated as x0 + 1/4."""
    weps = 2*eps
    x = x0 * mpf(1)
    if x1 is None:
        xprev = x0 + mpf(0.25)
    else:
        xprev = x1 * mpf(1)
    deriv_prev = None
    fxprev = f(xprev)
    for i in xrange(maxsteps):
        if verbose:
            print "Step", i
            print "x =", x
        fx = f(x)
        ydiff = fx - fxprev
        xdiff = x - xprev
        if verbose:
            print "f(x) =", fx
            print "xdiff = ", xdiff
            print "ydiff = ", ydiff
        try:
            deriv = xdiff / ydiff
            deriv_prev = deriv
        except ZeroDivisionError:
            if deriv_prev is None:
                raise ZeroDivisionError(msg1)
            if verbose and abs(xdiff) > weps:
                print msg2
            deriv = deriv_prev
        x, xprev = x - fx*deriv, x
        fxprev = fx
        if verbose:
            print
        if abs(xdiff) <= weps:
            break
    return x


#----------------------------------------------------------------------------#
#                                Polynomials                                 #
#----------------------------------------------------------------------------#

def polyval(coeffs, x, derivative=False):
    """
    Given coefficients [c0, c1, c2, ..., cn], evaluate
    P(x) = c0 + c1*x + c2*x**2 + ... + cn*x**n.

    If derivative=True is set, a tuple (P(x), P'(x)) is returned.
    """
    p = mpnumeric(coeffs[-1])
    q = mpf(0)
    for c in coeffs[-2::-1]:
        if derivative:
            q = p + x*q
        p = c + x*p
    if derivative:
        return p, q
    else:
        return p

def polyroots(coeffs, maxsteps=20):
    """
    Numerically locate all (complex) roots of a polynomial using the
    Durand-Kerner method.

    This function returns a tuple (roots, err) where roots is a list of
    complex numbers sorted by absolute value, and err is an estimate of
    the maximum error. The polynomial should be given as a list of
    coefficients.

        >>> nprint(polyroots([24,-14,-1,1]), 4)
        ([(2.0 + 8.968e-44j), (3.0 + 1.156e-33j), (-4.0 + 0.0j)], 5.921e-16)
        >>> nprint(polyroots([2,3,4]))
        ([(-0.375 + -0.599479j), (-0.375 + 0.599479j)], 2.22045e-16)

    """
    deg = len(coeffs) - 1
    # Must be monic
    lead = mpnumeric(coeffs[-1])
    if lead == 1:
        coeffs = map(mpnumeric, coeffs)
    else:
        coeffs = [c/lead for c in coeffs]
    f = lambda x: polyval(coeffs, x)
    roots = [mpc((0.4+0.9j)**n) for n in range(deg)]
    error = [mpf(1) for n in range(deg)]
    for step in range(maxsteps):
        if max(error).ae(0):
            break
        for i in range(deg):
            if not error[i].ae(0):
                p = roots[i]
                x = f(p)
                for j in range(deg):
                    if i != j:
                        try:
                            x /= (p-roots[j])
                        except ZeroDivisionError:
                            continue
                roots[i] = p - x
                error[i] = abs(x)
    roots.sort(key=abs)
    err = max(error)
    err = max(err, ldexp(1, -mp.prec+1))
    return roots, err


##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#       Implementation of tanh-sinh (doubly exponential) quadrature          #
#----------------------------------------------------------------------------#

"""
The implementation of the tanh-sinh algorithm is based on the
description given in Borwein, Bailey & Girgensohn, "Experimentation
in Mathematics - Computational Paths to Discovery", A K Peters,
2003, pages 312-313.

Various documents are available online, e.g.
http://crd.lbl.gov/~dhbailey/dhbpapers/dhb-tanh-sinh.pdf
http://users.cs.dal.ca/~jborwein/tanh-sinh.pdf
"""

def transform(f, a, b):
    """Given an integrand f defined over the interval [a, b], return an
    equivalent integrand g defined on the standard interval [-1, 1].

    If a and b are finite, this is achived by means of a linear change
    of variables. If at least one point is infinite, the substitution
    t = 1/x is used."""
    a = mpf(a)
    b = mpf(b)
    if (a, b) == (-1, 1):
        return f
    one = mpf(1)
    half = mpf(0.5)
    # The transformation 1/x sends [1, inf] to [0, 1], which in turn
    # can be transformed to [-1, 1] the usual way. For a double
    # infinite interval, we simply evaluate the function symmetrically
    if (a, b) == (-inf, inf):
        # return transform(lambda x: (f(-1/x+1)+f(1/x-1))/x**2, 0, 1)
        def y(x):
            u = 2/(x+one)
            w = one - u
            return half * (f(w)+f(-w)) * u**2
        return y
    if a == -inf:
        # return transform(lambda x: f(-1/x+b+1)/x**2, 0, 1)
        b1 = b+1
        def y(x):
            u = 2/(x+one)
            return half * f(b1-u) * u**2
        return y
    if b == inf:
        # return transform(lambda x: f(1/x+a-1)/x**2, 0, 1)
        a1 = a-1
        def y(x):
            u = 2/(x+one)
            return half * f(a1+u) * u**2
        return y
    # Simple linear change of variables
    C = (b-a)/2
    D = (b+a)/2
    def g(x):
        return C * f(D + C*x)
    return g

def TS_estimate_error(res, prec, eps):
    """Estimate error of the calculation at the present level by
    comparing it to the results from two previous levels. The
    algorithm is that described by D. H. Bailey."""
    try:
        D1 = log(abs(res[-1]-res[-2]), 10)
        D2 = log(abs(res[-1]-res[-3]), 10)
    except ValueError:
        return eps
    D3 = -prec
    D4 = min(0, max(D1**2/D2, 2*D1, D3))
    return mpf('0.1') ** -int(D4)

def TS_guess_level(prec):
    """Guess a reasonable first level of tanh-sinh quadrature for a
    given precision. The level should not be too high, or time will
    be wasted on computing unneeded nodes, and not too low or
    the integrator may fail or have to restart unnecessarily. This
    function gives e.g.
        50 bits -> 4
        100 bits -> 5
        500 bits -> 8
        3000 bits -> 10
    These numbers are based purely on a limited amount of
    experimentation and will sometimes be wrong."""
    return int(4 + max(0, log(prec/30.0, 2)))


def TS_node(k, hn, prec, a, ar):
    """Calculate an (abscissa, weight) pair for tanh-sinh quadrature.

        x[k] = tanh(pi/2 * sinh(k*h))
        w[k] = pi/2 * cosh(k*h) / cosh(pi/2 sinh(k*h))**2
    """
    oldprec = mp.prec
    mp.prec = prec
    mp.rounding = 'up'
    t = ldexp(mpf(k), -hn)
    # We only need to calculate one exponential
    #sinht, cosht = ldexp(a-ar, -1), ldexp(a+ar, -1)
    b = exp(a - ar); br = 1/b
    sinhb, coshb = ldexp(b-br, -1), ldexp(b+br, -1)
    x, w = sinhb/coshb, (a + ar)/coshb**2
    mp.rounding = 'default'
    mp.prec = oldprec
    return x, w


TS_cache = {}

def TS_nodes(prec, m, verbose=False):
    """
    Return a list of abscissas and a list of corresponding weights for
    tanh-sinh quadrature at level m with precision prec.
    """
    if (prec, m) in TS_cache:
        return TS_cache[(prec, m)]
    eps = ldexp(1, -prec)
    t = ldexp(1, -m)
    a0 = exp(t); a0m = 1/a0
    a1 = a1m = ldexp(pi, -2)
    xs = [0]
    ws = [ldexp(pi, -1)]
    for k in xrange(1, 20 * 2**m + 1):
        a1 = a1*a0
        a1m = a1m*a0m
        x, w = TS_node(k, m, prec, a1, a1m)
        diff = abs(x-1)
        if diff <= eps:
            break
        if verbose and m > 6 and k % 300 == 150:
            # note: the number displayed is rather arbitrary. should
            # figure out how to print something that looks more like a
            # percentage
            print "Calculating nodes:", nstr(-log(diff, 10) / prec)
        xs.append(x)
        ws.append(w)
    TS_cache[(prec, m)] = (xs, ws)
    return xs, ws

def TS_eval(f, nodes, target_prec, working_prec, m, verbose=False):
    """Evaluate f at the given set of tanh-sinh nodes."""
    eps = ldexp(1, -target_prec)
    S = mpf(0)
    h = mpf(1)
    xs, ws = nodes
    res = []
    for k in xrange(1, m+1):
        if verbose:
            print "Evaluating integral (level %s of %s)" % (k, m)
        h = h / 2
        for i in xrange(0, len(xs), 2**(m-k)):
            if i % (2**(m-k+1)) != 0 or k == 1:
                if i == 0:
                    S = S + ws[0]*f(mpf(0))
                else:
                    S = S + (ws[i])*(f(-xs[i]) + f(xs[i]))
        res.append(h*S)
        if k > 2:
            err = TS_estimate_error(res, target_prec, eps)
            if verbose:
                print "Estimated error:", nstr(err)
            if err <= eps:
                break
    return +res[-1], TS_estimate_error(res, target_prec, eps)


def TS_adaptive(f, target_prec, working_prec, min_level, max_level, verbose):
    """Repeatedly attempt to integrate f, trying different levels. Quit
    as soon as the estimated error is small enough, or if that doesn't
    happen, when the max level has been tried."""
    eps = ldexp(1, -target_prec)
    for m in xrange(min_level, max_level+1):
        if verbose:
            print "Using tanh-sinh algorithm with level ", m
        nodes = TS_nodes(working_prec, m, verbose=verbose)
        s, err = TS_eval(f, nodes, target_prec, working_prec, m,
            verbose=verbose)
        steps = 2*len(nodes[0])
        if err <= eps:
            return s, err
    if verbose:
        print "Failed to reach full accuracy. Estimated error:", \
            nstr(err)
    return s, err


def quadts(f, a, b, **options):
    """
    Integrate f(x) dx over the interval [a, b], using tanh-sinh
    quadrature. Use quadts(f, (a, b), (c, d)) to calculate the
    two-dimensional integral over [a, b] x [c, d].

        >>> print quadts(lambda x: x**2, -2, 4)
        24.0
        >>> print quadts(lambda x, y: x+y, (0, 1), (0, 2))
        3.0

    Options
    =======

    target_prec
        The number of accurate bits to aim for in the result. If not
        specified, the current working precision mp.prec is used.

    working_prec
        Precision to use when evaluating the function. This should
        be set slightly higher than the target precision to eliminate
        the effects of rounding and cancellation errors.

    min_level
    max_level
        The quadts function first attempts to perform tanh-sinh
        quadrature at min_level; if that fails, at min_level+1, etc, up
        to max_level. One 'level' corresponds to roughly 2**n
        evaluation points. The levels should be integers roughly
        of the size 5-10. If not specified, reasonable values are
        inferred from the target precision.

    error
        Set to True to obtain an error estimate along with the result.

    verbose
        Set to True to display progress messages.
    """

    verbose = options.get('verbose', False)
    target_prec = options.get('target_prec', mp.prec)
    working_prec = options.get('working_prec', target_prec + 20)
    min_level = options.get('min_level', TS_guess_level(target_prec))
    max_level = options.get('max_level', min_level + 2)

    oldprec = mp.prec
    mp.prec = working_prec

    # Handle double integrals
    if isinstance(a, tuple):
        (a, b), (c, d) = a, b
        if c == d:
            return mpf(0)
        g = f
        # Define the inner integral to recursively call quadts. We must
        # be careful to pass along the right settings.
        def f(x):
            return quadts(lambda y: g(x,y), c, d,
                min_level=min_level, max_level=max_level,
                target_prec=target_prec, working_prec=working_prec)
    if a == b:
        return mpf(0)

    # Based on experience, integrals on (semi)infinite intervals require
    # a little extra work
    if inf in (abs(a), abs(b)):
        min_level += 1; max_level += 1

    # Standardize to [-1, 1] and evaluate
    f = transform(f, a, b)
    val, err = TS_adaptive(f, target_prec, working_prec,
        min_level, max_level, verbose=verbose)

    mp.prec = oldprec
    if options.get('error', False):
        return val, err
    return val

##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                               Numerical summation                          #
#----------------------------------------------------------------------------#

@extraprec(15, normalize_output=True)
def sumem(f, a=0, b=inf, N=None, integral=None, fderiv=None, verbose=False):
    """
    Calculate the sum of f(n) for n = a..b using Euler-Maclaurin
    summation. This algorithm is efficient for slowly convergent
    nonoscillatory sums; the essential condition is that f must be
    analytic. The method relies on approximating the sum by an
    integral, so f must be smooth and well-behaved enough to be
    integrated numerically.

    A tuple (s, err) is returned where s is the calculated sum and err
    is the estimated magnitude of the error. With verbose=True,
    detailed information about progress and errors is printed.

        >>> mp.dps = 15
        >>> s, err = sumem(lambda n: 1/n**2, 1, inf)
        >>> print s
        1.64493406684823
        >>> print pi**2 / 6
        1.64493406684823
        >>> nprint(err)
        2.22045e-16

    N is the number of terms to compute directly before using the
    Euler-Maclaurin formula to approximate the tail. It must be set
    high enough; often roughly N ~ dps is the right size.

    High-order derivatives of f are also needed. By default, these
    are computed using numerical integration, which is the most
    expensive part of the calculation. The default method assumes
    that all poles of f are located close to the origin. A custom
    nth derivative function fderiv(x, n) can be provided as a
    keyword parameter.

    This is much more efficient:

        >>> f = lambda n: 1/n**2
        >>> fp = lambda x, n: (-1)**n * factorial(n+1) * x**(-2-n)
        >>> mp.dps = 50
        >>> s, err = sumem(lambda n: 1/n**2, 1, inf, fderiv=fp)
        >>> print s
        1.6449340668482264364724151666460251892189499012068
        >>> print pi**2 / 6
        1.6449340668482264364724151666460251892189499012068

    If b = inf, f and its derivatives are all assumed to vanish
    at infinity. It is assumed that a is finite, so doubly
    infinite sums cannot be evaluated directly.
    """
    if N is None:
        N = 3*mp.dps + 20
    a, b, N = mpf(a), mpf(b), mpf(N)
    infinite = (b == inf)
    weps = eps * 2**8
    if verbose:
        print "Summing f(k) from k = %i to %i" % (a, a+N-1)
    S = sum(f(mpf(k)) for k in xrange(a, a+N))
    if integral is None:
        if verbose:
            print "Integrating f(x) from x = %i to %s" % (a+N, nstr(b))
        I, ierr = quadts(f, a+N, b, error=1)
    else:
        I, ierr = integral(a+N, b), mpf(0)
    # There is little hope if the tail cannot be integrated
    # accurately. Estimate magnitude of tail as the error.
    if ierr > weps:
        if verbose:
            print "Failed to converge to target accuracy (integration failed)"
        return S+I, abs(I) + ierr
    if infinite:
        C = f(a+N) / 2
    else:
        C = (f(a+N) + f(b)) / 2
    # Default (inefficient) approach for derivatives
    if not fderiv:
        fderiv = lambda x, n: diffc(f, x, n, radius=N*0.75)
    k = 1
    prev = 0
    if verbose:
        print "Summing tail"
    B = bernoulli2n()
    fac = 2
    while 1:
        if infinite:
            D = fderiv(a+N, 2*k-1)
        else:
            D = fderiv(a+N, 2*k-1) - fderiv(b, 2*k-1)
        # B(2*k) / fac(2*k)
        term = B.next() / fac * D
        mag = abs(term)
        if verbose:
            print "term", k, "magnitude =", nstr(mag)
        # Error can be estimated as the magnitude of the smallest term
        if k >= 2:
            if mag < weps:
                if verbose:
                    print "Converged to target accuracy"
                res, err = I + C + S, eps * 2**15
                break
            if mag > abs(prev):
                if verbose:
                    print "Failed to converge to target accuracy (N too low)"
                res, err = I + C + S, abs(term)
                break
        S -= term
        k += 1
        fac *= (2*k) * (2*k-1)
        prev = term
    if isinstance(res, mpc) and not isinstance(I, mpc):
        return res.real, err
    return res, err

def ODE_integrate(t_list, x0, derivs, step):
    """
    Given the list t_list of values, returns the solution at these points.
    """
    x = x0
    result = [x]
    for i in range(len(t_list)-1):
        dt = t_list[i+1] - t_list[i]
        x = step(t_list[i], x, dt, derivs)
        result.append(x)
    return result

def smul(a, x):
    """Multiplies the vector "x" by the scalar "a"."""
    R = []
    for i in range(len(x)):
        R.append(a*x[i])
    return R

def vadd(*args):
    """Adds vectors "x", "y", ... together."""
    assert len(args) >= 2
    n = len(args[0])
    for x in args:
        assert len(x) == n
    R = []
    for i in range(n):
        s = 0.
        for x in args:
            s += x[i]
        R.append(s)
    return R

def ODE_step_euler(x, y, h, derivs):
    """
    Advances the solution y(x) from x to x+h using the Euler method.

    derivs .... a python function f(x, (y1, y2, y3, ...)) returning
    a tuple (y1', y2', y3', ...) where y1' is the derivative of y1 at x.
    """
    X = derivs(x,y)
    return vadd(y, smul(h, X))

half = mpf(0.5)

def ODE_step_rk4(x, y, h, derivs):
    """
    Advances the solution y(x) from x to x+h using the 4th-order Runge-Kutta
    method.

    derivs .... a python function f(x, (y1, y2, y3, ...)) returning
    a tuple (y1', y2', y3', ...) where y1' is the derivative of y1 at x.
    """
    h2 = h/2
    third = mpf(1)/3
    sixth = mpf(1)/6
    k1 = smul(h, derivs(x, y))
    k2 = smul(h, derivs(x+h2, vadd(y, smul(half, k1))))
    k3 = smul(h, derivs(x+h2, vadd(y, smul(half, k2))))
    k4 = smul(h, derivs(x+h, vadd(y, k3)))
    return vadd(y, smul(sixth, k1), smul(third, k2), smul(third, k3), 
            smul(sixth, k4))

def arange(a, b, dt):
    """
    Returns a list [a, a + dt, a + 2*dt, ..., b]
    """
    a, b, dt = mpf(a), mpf(b), mpf(dt)
    t = a
    result = [t]
    while t <= b:
        t += dt
        result.append(t)
    return result

def hypser(a, b, c, z):
    """
    Calculates 2F1 using a series expansion.
    """
    a, b, c = mpf(a), mpf(b), mpf(c)
    z = mpc(z)
    fac = mpf("1.0")
    temp = fac
    deriv = mpf("0.0")
    # The upper bound (1000) should be either passed as a parameter, or
    # automatically determined for some given precision:
    for n in range(1,1000):
        fac *= a*b/c
        deriv += fac
        fac *= z/n
        series = temp + fac
        temp = series
        a += 1
        b += 1
        c += 1
    return series, deriv

def hypgeo(a, b, c, z):
    """
    Calculates 2F1 by solving a differential equation.
    """

    def derivs(s, (yy1, yy2, yy3, yy4)):
        """
        F .... y1 = yy1 + I*yy2
        F' ... y2 = yy3 + I*yy4

        """
        z1 = z
        y1 = yy1 + 1j*yy2
        y2 = yy3 + 1j*yy4
        Z = z0 + s * (z1-z0)
        A = y2 * (z1-z0)
        B = (z1-z0)*(a*b*y1-(c-(a+b+1)*Z)*y2)/(Z*(1-Z))
        return A.real, A.imag, B.real, B.imag

    a, b, c = mpf(a), mpf(b), mpf(c)
    z = mpc(z)

    if z.real**2 + z.imag**2 <= 0.25:
        series, deriv = hypser(a, b, c, z)
        return series
    if z.real < 0.0: 
        z0 = -mpf("0.5")
    elif z.real <= 1.0:
        z0 = mpf("0.5")
    else:
        if z.imag >= 0.:
            z0 = mpf("0.5")*1j
        else:
            z0 = -mpf("0.5")*1j

    solver = ODE_step_rk4
    s = arange(0, 1, 0.0001)
    series, deriv = hypser(a, b, c, z0)
    yy1, yy2 = series.real, series.imag
    yy3, yy4 = deriv.real, deriv.imag
    sol = ODE_integrate(s, (yy1, yy2, yy3, yy4), derivs, solver)
    yy1 = sol[-1][0]
    yy2 = sol[-1][1]
    y1 = yy1 + 1j*yy2
    F = y1
    return F

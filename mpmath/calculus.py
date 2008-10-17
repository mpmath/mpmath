"""
High-level, mostly calculus-oriented functions.

* Numerical differentiation
* Numerical polynomial operations
* Numerical root-finding
* Numerical integration

etc
"""

__docformat__ = 'plaintext'

from settings import (mp, extraprec)
from mptypes import (mpnumeric, convert_lossless, mpf, mpc, j, inf, eps,
    AS_POINTS, arange, nstr, nprint)
from functions import (ldexp, factorial, exp, cos, pi, bernoulli, sign)

from quadrature import quadgl, quadts

#----------------------------------------------------------------------------#
#                       General extrapolation methods                        #
#----------------------------------------------------------------------------#

def richardson_extrapolation(f, n, N):
    if not callable(f):
        g = f; f = lambda k: g.__getitem__(int(k))
    orig = mp.prec
    try:
        mp.prec = 2*orig
        s = mpf(0)
        for j in range(0, N+1):
            c = (n+j)**N * (-1)**(j+N) / (factorial(j) * factorial(N-j))
            s += c * f(mpf(n+j))
    finally:
        mp.prec = orig
    return +s

def shanks_extrapolation(f, n, m=1):
    if not callable(f):
        g = f; f = lambda k: g.__getitem__(int(k))
    orig = mp.prec
    try:
        mp.prec = 2*orig
        table = [f(mpf(j)) for j in range(n+m+2)]
        table2 = table[:]
        for i in range(1, m+1):
            for j in range(i, n+m+1):
                x, y, z = table[j-1], table[j], table[j+1]
                try:
                    table2[j] = (z*x - y**2) / (z + x - 2*y)
                except ZeroDivisionError:
                    return table2[j-1]
            table = table2[:]
    finally:
        mp.prec = orig
    return +table[n]

def limit(f, x, direction=-1, n=None, N=None):
    """Compute lim of f(t) as t -> x using Richardson extrapolation.
    For infinite x, the function values [f(n), ... f(n+N)] are used.
    For finite x, [f(x-direction/n)), ... f(x-direction/(n+N))] are
    used. If x is inf, f can be also be a precomputed sequence
    with a __getitem__ method."""
    if callable(f):
        if not n: n = 3 + int(mp.dps * 0.5)
        if not N: N = 2*n
    else:
        # If a sequence, take as many terms are are available
        g = f; f = lambda k: g.__getitem__(int(k))
        if not N: N = len(g)-1
        if not n: n = 0
    if   x == inf:  return richardson_extrapolation(lambda k: f(mpf(k)), n, N)
    elif x == -inf: return richardson_extrapolation(lambda k: f(mpf(-k)), n, N)
    direction *= mpf(1)
    def g(k):
        return f(x - direction/(1+k))
    return richardson_extrapolation(g, n, N)


#----------------------------------------------------------------------------#
#                                Differentiation                             #
#----------------------------------------------------------------------------#

def diff(f, x, n=1, method='step', scale=1, direction=0):
    """
    Numerically computes the derivative of f'(x). Optionally, computes
    the nth derivative f^(n)(x), for any order n. Basic examples:

        >>> print diff(lambda x: x**2 + x, 1.0)
        3.0
        >>> print diff(lambda x: x**2 + x, 1.0, 2)
        2.0
        >>> print diff(lambda x: x**2 + x, 1.0, 3)
        0.0

    The exponential function is invariant under differentiation:

        >>> nprint([diff(exp, 3, n) for n in range(5)])
        [20.0855, 20.0855, 20.0855, 20.0855, 20.0855]

    Method
    ------

    method = 'step':

    By default, the derivative is computed using a finite difference
    approximation, with a small step h. This requires n+1 function
    evaluations and must be performed at (n+1) times the target
    precison. Accordingly, f must support fast evaluation at high
    precision.

    method = 'quad':

    Alternatively, the derivative can be computed using complex
    numerical integration. This requires a larger number of function
    evaluations, but the advantage is that not much extra precision
    is required. For high order derivatives, this method may thus
    be faster if f is very expensive to evaluate at high precision.

    With method='quad', n may be any complex number, not necessarily
    an integer. This allows fractional derivatives to be computed.

    Scale
    -----

    The scale option specifies the scale of variation of f. The step
    size in the finite difference is taken to be approximately
    eps*scale. Thus, for example if f(x) = cos(1000*x), the scale
    should be set to 1/1000 and if f(x) = cos(x/1000), the scale
    should be 1000. By default, scale = 1.

    (In practice, the default scale will work even for cos(1000*x) or
    cos(x/1000). Changing this parameter is a good idea if the scale
    is something *preposterous*.)

    If numerical integration is used, the radius of integration is
    taken to be equal to scale/2. Note that f must not have any
    singularities within the circle of radius scale/2 centered around
    x. If possible, a larger scale value is preferable because it
    typically makes the integration faster and more accurate.

    Direction
    ---------

    By default, diff uses a central difference approximation.
    This corresponds to direction=0. Alternatively, it can compute a
    left difference (direction=-1) or right difference (direction=1).
    This is useful for computing left- or right-sided derivatives
    of nonsmooth functions:

        >>> print diff(abs, 0, direction=0)
        0.0
        >>> print diff(abs, 0, direction=1)
        1.0
        >>> print diff(abs, 0, direction=-1)
        -1.0

    More generally, if the direction is nonzero, a right difference
    is computed where the step size is multiplied by sign(direction).
    For example, with direction=+j, the derivative from the positive
    imaginary direction will be computed.

    This option only makes sense with method='step'. If integration
    is used, it is assumed that f is analytic, implying that the
    derivative is the same in all directions.

    """
    if n == 0:
        return f(x)
    orig = mp.prec
    try:
        if method == 'step':
            mp.prec = (orig+20) * (n+1)
            h = ldexp(scale, -orig-10)
            v = mpf(0)
            # Applying the finite difference formula recursively n times,
            # we get a step sum weighted by a row of binomial coefficients
            # Directed: steps x, x+h, ... x+n*h
            if direction:
                h *= sign(direction)
                steps = xrange(0, n+1)
                norm = h**n
            # Central: steps x-n*h, x-(n-2)*h ..., x, ..., x+(n-2)*h, x+n*h
            else:
                steps = xrange(-n, n+1, 2)
                norm = (2*h)**n
            # Initial binomial coefficient is 1 or -1
            b = (-1) ** (n & 1)
            k = 1
            for i in steps:
                v += b * f(x+i*h)
                # Binomial recurrence
                b = (b * (n-k+1)) // (-k)
                k += 1
            v = v / norm
        elif method == 'quad':
            mp.prec += 10
            radius = mpf(scale)/2
            def g(t):
                rei = radius*exp(j*t)
                z = x + rei
                return f(z) / rei**n
            d = quadts(g, [0, 2*pi])
            v = d * factorial(n) / (2*pi)
        else:
            raise ValueError("unknown method: %r" % method)
    finally:
        mp.prec = orig
    return +v

def diffun(f, n=1, **options):
    """
    Given a function f, returns a function g(x) that evaluates the nth
    derivative f^(n)(x).

        >>> cos2 = diffun(sin)
        >>> sin2 = diffun(sin, 4)
        >>> print cos(1.3), cos2(1.3)
        0.267498828624587 0.267498828624587
        >>> print sin(1.3), sin2(1.3)
        0.963558185417193 0.963558185417193

    The function f must support arbitrary precision evaluation.
    See diff() for additional details and supported keyword options.
    """
    if n == 0:
        return f
    def g(x):
        return diff(f, x, n, **options)
    return g

def taylor(f, x, n, **options):
    """
    Produce a degree-n Taylor polynomial around the point x of the
    given function f. The coefficients are returned as a list.

        >>> nprint(taylor(sin, 0, 5))
        [0.0, 1.0, 0.0, -0.166667, 1.08755e-55, 8.33333e-3]

    The coefficients are computed using high-order numerical
    differentiation. The function must be possible to evaluate
    to arbitrary precision. See diff() for additional details and
    supported keyword options.

    Note that to evaluate the Taylor polynomial as an approximation
    of f, e.g. with polyval, the coefficients must be reversed, and
    the point of the Taylor expansion must be subtracted from
    the argument:

        >>> p = taylor(exp, 2.0, 10)
        >>> print polyval(p[::-1], 2.5 - 2.0)
        12.1824939606092
        >>> print exp(2.5)
        12.1824939607035

    """
    return [diff(f, x, i, **options) / factorial(i) for i in xrange(n+1)]


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
    Given coefficients [cn, ..., c2, c1, c0], evaluate
    P(x) = cn*x**n + ... + c2*x**2 + c1*x + c0.

    If derivative=True is set, a tuple (P(x), P'(x)) is returned.
    """
    if not coeffs:
        return mpf(0)
    p = mpnumeric(coeffs[0])
    q = mpf(0)
    for c in coeffs[1:]:
        if derivative:
            q = p + x*q
        p = c + x*p
    if derivative:
        return p, q
    else:
        return p

def polyroots(coeffs, maxsteps=50, cleanup=True, extraprec=10, error=False):
    """
    Numerically locate all (complex) roots of a polynomial using the
    Durand-Kerner method.

    With error=True, this function returns a tuple (roots, err) where roots
    is a list of complex numbers sorted by absolute value, and err is an
    estimate of the maximum error. The polynomial should be given as a list
    of coefficients, in the same format as accepted by polyval(). The
    leading coefficient must be nonzero.

    These are the roots of x^3 - x^2 - 14*x + 24 and 4x^2 + 3x + 2:

        >>> nprint(polyroots([1,-1,-14,24]), 4)
        [-4.0, 2.0, 3.0]
        >>> nprint(polyroots([4,3,2], error=True))
        ([(-0.375 - 0.599479j), (-0.375 + 0.599479j)], 2.22045e-16)

    """
    if len(coeffs) <= 1:
        if not coeffs or not coeffs[0]:
            raise ValueError("Input to polyroots must not be the zero polynomial")
        # Constant polynomial with no roots
        return []

    orig = mp.prec
    weps = +eps
    try:
        mp.prec += 10
        deg = len(coeffs) - 1
        # Must be monic
        lead = convert_lossless(coeffs[0])
        if lead == 1:
            coeffs = map(convert_lossless, coeffs)
        else:
            coeffs = [c/lead for c in coeffs]
        f = lambda x: polyval(coeffs, x)
        roots = [mpc((0.4+0.9j)**n) for n in xrange(deg)]
        err = [mpf(1) for n in xrange(deg)]
        for step in xrange(maxsteps):
            if max(err).ae(0):
                break
            for i in xrange(deg):
                if not err[i].ae(0):
                    p = roots[i]
                    x = f(p)
                    for j in range(deg):
                        if i != j:
                            try:
                                x /= (p-roots[j])
                            except ZeroDivisionError:
                                continue
                    roots[i] = p - x
                    err[i] = abs(x)
        if cleanup:
            for i in xrange(deg):
                if abs(roots[i].imag) < weps:
                    roots[i] = roots[i].real
                elif abs(roots[i].real) < weps:
                    roots[i] = roots[i].imag * 1j
        roots.sort(key=lambda x: (abs(x.imag), x.real))
    finally:
        mp.prec = orig
    if error:
        err = max(err)
        err = max(err, ldexp(1, -orig+1))
        return [+r for r in roots], +err
    else:
        return [+r for r in roots]

def quadosc(f, interval, period=None, zeros=None, alt=1):
    """
    Integrates f(x) over interval = [a, b] where at least one of a and
    b is infinite and f is a slowly decaying oscillatory function. The
    zeros of f must be provided, either by specifying a period (suitable
    when f contains a pure sine or cosine factor) or providing a function
    that returns the nth zero (suitable when the oscillation is not
    strictly periodic).
    """
    a, b = AS_POINTS(interval)
    a = convert_lossless(a)
    b = convert_lossless(b)
    if period is None and zeros is None:
        raise ValueError( \
            "either the period or zeros keyword parameter must be specified")
    if a == -inf and b == inf:
        s1 = quadosc(f, [a, 0], zeros=zeros, period=period, alt=alt)
        s2 = quadosc(f, [0, b], zeros=zeros, period=period, alt=alt)
        return s1 + s2
    if a == -inf:
        if zeros:
            return quadosc(lambda x:f(-x), [-b,-a], lambda n: zeros(-n), alt=alt)
        else:
            return quadosc(lambda x:f(-x), [-b,-a], period=period, alt=alt)
    if b != inf:
        raise ValueError("quadosc requires an infinite integration interval")
    if not zeros:
        zeros = lambda n: n*period/2
    for n in range(1,10):
        p = zeros(n)
        if p > a:
            break
    if n >= 9:
        raise ValueError("zeros do not appear to be correctly indexed")
    if alt == 0:
        s = quadgl(f, [a, zeros(n+1)])
        s += sumrich(lambda k: quadgl(f, [zeros(2*k), zeros(2*k+2)]), [n, inf])
    else:
        s = quadgl(f, [a, zeros(n)])
        s += sumsh(lambda k: quadgl(f, [zeros(k), zeros(k+1)]), [n, inf])
    return s


##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                               Numerical summation                          #
#----------------------------------------------------------------------------#

def sumrich(f, interval, n=None, N=None):
    """
    Sum f(k) for k = a, a+1, ..., b where [a, b] = interval,
    using Richardson extrapolation. This function is essentially
    equivalent to using limit() on the sequence of partial sums.
    """
    a, b = AS_POINTS(interval)
    assert b == inf
    if not n: n = 3 + int(mp.dps * 0.5)
    if not N: N = 2*n
    orig = mp.prec
    try:
        mp.prec = 2*orig
        s = mpf(0)
        tbl = []
        for k in range(a, a+n+N+1):
            s += f(mpf(k))
            tbl.append(s)
        s = richardson_extrapolation(lambda k: tbl[int(k)], n, N)
    finally:
        mp.prec = orig
    return +s

def sumsh(f, interval, n=None, m=None):
    """
    Sum f(k) for k = a, a+1, ..., b where [a, b] = interval,
    using an n-term Shanks transformation. With m > 1, the Shanks
    transformation is applied recursively m times.

    Shanks summation often works well for slowly convergent and/or
    alternating Taylor series.
    """
    a, b = AS_POINTS(interval)
    assert b == inf
    if not n: n = 5 + int(mp.dps * 1.2)
    if not m: m = 2 + n//3
    orig = mp.prec
    try:
        mp.prec = 2*orig
        s = mpf(0)
        tbl = []
        for k in range(a, a+n+m+2):
            s += f(mpf(k))
            tbl.append(s)
        s = shanks_extrapolation(tbl, n, m)
    finally:
        mp.prec = orig
    return +s

@extraprec(15, normalize_output=True)
def sumem(f, interval, N=None, integral=None, fderiv=None, error=False,
    verbose=False):
    """
    Sum f(k) for k = a, a+1, ..., b where [a, b] = interval,
    using Euler-Maclaurin summation. This algorithm is efficient
    for slowly convergent nonoscillatory sums; the essential condition
    is that f must be analytic. The method relies on approximating the
    sum by an integral, so f must be smooth and well-behaved enough
    to be integrated numerically.

    With error=True, a tuple (s, err) is returned where s is the
    calculated sum and err is the estimated magnitude of the error.
    With verbose=True, detailed information about progress and errors
    is printed.

        >>> mp.dps = 15
        >>> s, err = sumem(lambda n: 1/n**2, 1, inf, error=True)
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
        >>> print sumem(lambda n: 1/n**2, 1, inf, fderiv=fp)
        1.6449340668482264364724151666460251892189499012068
        >>> print pi**2 / 6
        1.6449340668482264364724151666460251892189499012068

    If b = inf, f and its derivatives are all assumed to vanish
    at infinity. It is assumed that a is finite, so doubly
    infinite sums cannot be evaluated directly.
    """
    a, b = AS_POINTS(interval)
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
        I, ierr = quadts(f, [a+N, b], error=1)
        # XXX: hack for relative error
        ierr /= abs(I)
    else:
        I, ierr = integral(a+N, b), mpf(0)
    # There is little hope if the tail cannot be integrated
    # accurately. Estimate magnitude of tail as the error.
    if ierr > weps:
        if verbose:
            print "Failed to converge to target accuracy (integration failed)"
        if error:
            return S+I, abs(I) + ierr
        else:
            return S+I
    if infinite:
        C = f(a+N) / 2
    else:
        C = (f(a+N) + f(b)) / 2
    # Default (inefficient) approach for derivatives
    if not fderiv:
        fderiv = lambda x, n: diff(f, x, n)
    k = 1
    prev = 0
    if verbose:
        print "Summing tail"
    fac = 2
    while 1:
        if infinite:
            D = fderiv(a+N, 2*k-1)
        else:
            D = fderiv(a+N, 2*k-1) - fderiv(b, 2*k-1)
        # B(2*k) / fac(2*k)
        term = bernoulli(2*k) / fac * D
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
        res, err = res.real, err
    if error:
        return res, err
    else:
        return res

#----------------------------------------------------------------------------#
#                                  ODE solvers                               #
#----------------------------------------------------------------------------#

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
    rest = args[1:]
    for x in args:
        assert len(x) == n
    R = []
    for i in range(n):
        s = args[0][i]
        for x in rest:
            s += x[i]
        R.append(s)
    return R

def ODE_step_euler(x, y, h, derivs):
    """
    Advances the solution y(x) from x to x+h using the Euler method.

    derivs .... a python function f(x, (y1, y2, y3, ...)) returning
    a tuple (y1', y2', y3', ...) where y1' is the derivative of y1 at x.
    """
    X = derivs(y,x)
    return vadd(y, smul(h, X))

half = mpf(0.5)

def ODE_step_rk4(x, y, h, derivs):
    """
    Advances the solution y(x) from x to x+h using the 4th-order Runge-Kutta
    method.

    derivs .... a python function f(x, (y1, y2, y3, ...)) returning
    a tuple (y1', y2', y3', ...) where y1' is the derivative of y1 at x.
    """
    h2 = ldexp(h, -1)
    third = mpf(1)/3
    k1 = smul(h, derivs(y, x))
    k2 = smul(h, derivs(vadd(y, smul(half, k1)), x+h2))
    k3 = smul(h, derivs(vadd(y, smul(half, k2)), x+h2))
    k4 = smul(h, derivs(vadd(y, k3), x+h))
    v = []
    for i in range(len(y)):
        v.append(y[i] + third*(k2[i]+k3[i] + half*(k1[i]+k4[i])))
    return v

def odeint(derivs, x0, t_list, step=ODE_step_rk4):
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

#----------------------------------------------------------------------------#
#                              Approximation methods                         #
#----------------------------------------------------------------------------#

# The Chebyshev approximation formula is given at:
# http://mathworld.wolfram.com/ChebyshevApproximationFormula.html

# The only major changes in the following code is that we return the
# expanded polynomial coefficients instead of Chebyshev coefficients,
# and that we automatically transform [a,b] -> [-1,1] and back
# for convenience.

# Coefficient in Chebyshev approximation
def chebcoeff(f,a,b,j,N):
    s = mpf(0)
    h = mpf(0.5)
    for k in range(1, N+1):
        t = cos(pi*(k-h)/N)
        s += f(t*(b-a)*h + (b+a)*h) * cos(pi*j*(k-h)/N)
    return 2*s/N

# Generate Chebyshev polynomials T_n(ax+b) in expanded form
def chebT(a=1, b=0):
    Tb = [1]
    yield Tb
    Ta = [b, a]
    while 1:
        yield Ta
        # Recurrence: T[n+1](ax+b) = 2*(ax+b)*T[n](ax+b) - T[n-1](ax+b)
        Tmp = [0] + [2*a*t for t in Ta]
        for i, c in enumerate(Ta): Tmp[i] += 2*b*c
        for i, c in enumerate(Tb): Tmp[i] -= c
        Ta, Tb = Tmp, Ta

def chebyfit(f, interval, N, error=False):
    """
    Chebyshev approximation: returns coefficients of a degree N-1
    polynomial that approximates f on the interval [a, b]. With error=True,
    also returns an estimate of the maximum error.
    """
    a, b = AS_POINTS(interval)
    orig = mp.prec
    try:
        mp.prec = orig + int(N**0.5) + 20
        c = [chebcoeff(f,a,b,k,N) for k in range(N)]
        d = [mpf(0)] * N
        d[0] = -c[0]/2
        h = mpf(0.5)
        T = chebT(mpf(2)/(b-a), mpf(-1)*(b+a)/(b-a))
        for k in range(N):
            Tk = T.next()
            for i in range(len(Tk)):
                d[i] += c[k]*Tk[i]
        d = d[::-1]
        # Estimate maximum error
        err = mpf(0)
        for k in range(N):
            x = cos(pi*k/N) * (b-a)*h + (b+a)*h
            err = max(err, abs(f(x) - polyval(d, x)))
    finally:
        mp.prec = orig
        if error:
            return d, +err
        else:
            return d

"""
This module implements high-level calculus-oriented functions:

* Numerical differentiation
* Numerical polynomial operations
* Numerical root-finding
* Numerical integration

"""

__docformat__ = 'plaintext'

from mptypes import *
from apps.extrafun import factorial

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
            return rei * f(z) / (z-x)**(n+1)
        d = quadts(g, 0, 2*pi)
        return d * factorial(n) / (2*pi)
    finally:
        mp.prec = prec


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
    err = max(err, ldexp(1, -getprec()+1))
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
    # The transformation 1/x sends [1, inf] to [0, 1], which in turn
    # can be transformed to [-1, 1] the usual way. For a double
    # infinite interval, we simply evaluate the function symmetrically
    if (a, b) == (-inf, inf):
        return transform(lambda x: (f(-1/x+1)+f(1/x-1))/x**2, 0, 1)
    if a == -inf:
        return transform(lambda x: f(-1/x+b+1)/x**2, 0, 1)
    if b == inf:
        return transform(lambda x: f(1/x+a-1)/x**2, 0, 1)
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


def TS_node(k, h, prec):
    """Calculate an (abscissa, weight) pair for tanh-sinh quadrature.

        x[k] = tanh(pi/2 * sinh(k*h))
        w[k] = pi/2 * cosh(k*h) / cosh(pi/2 sinh(k*h))**2
    """
    oldprec = mp.prec
    mp.prec = prec
    mp.rounding = 'up'
    t = mpf(k) * h
    # We only need to calculate one exponential
    a = exp(t); ar = 1/a
    sinht, cosht = (a-ar)/2, (a+ar)/2
    b = exp((pi * sinht) / 2); br = 1/b
    sinhb, coshb = (b-br)/2, (b+br)/2
    x, w = sinhb/coshb, (pi/2)*cosht/coshb**2
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
    h = ldexp(1, -m)
    xs = []
    ws = []
    for k in xrange(20 * 2**m + 1):
        x, w = TS_node(k, h, prec)
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

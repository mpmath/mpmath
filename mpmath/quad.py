"""
The implementation of the tanh-sinh algorithm is based on the
description given in Borwein, Bailey & Girgensohn, "Experimentation
in Mathematics - Computational Paths to Discovery", A K Peters,
2003, pages 312-313.

Various documents are available online, e.g.
http://crd.lbl.gov/~dhbailey/dhbpapers/dhb-tanh-sinh.pdf
http://users.cs.dal.ca/~jborwein/tanh-sinh.pdf
"""

from mpmath import *


#----------------------------------------------------------------------
# Utilities
#

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
    # The substitution 1/x sends [0, inf] to [0, 1], which in turn
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

def smallstr(f):
    """Represent x as a string with just a few digits"""
    prec = mpf.prec
    mpf.dps = 5
    s = str(f)
    mpf.prec = prec
    return s


#----------------------------------------------------------------------
# Tanh-sinh (doubly exponential) quadrature
#

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
    oldprec = mpf.prec; mpf.prec = prec
    mpf.round_up()
    t = mpf(k) * h
    # We only need to calculate one exponential
    a = exp(t); ar = 1/a
    sinht, cosht = (a-ar)/2, (a+ar)/2
    b = exp((pi * sinht) / 2); br = 1/b
    sinhb, coshb = (b-br)/2, (b+br)/2
    x, w = sinhb/coshb, (pi/2)*cosht/coshb**2
    mpf.round_default()
    mpf.prec = oldprec
    return x, w

TS_cache = {}

def TS_nodes(prec, m, verbose=False):
    """
    Return a list of abscissas and a list of corresponding weights for
    tanh-sinh quadrature at level m with precision prec.
    """
    if (prec, m) in TS_cache:
        return TS_cache[(prec, m)]
    eps = mpf((1, -prec))
    h = mpf((1, -m))
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
            print "calculating nodes:", smallstr(-log(diff, 10) / prec)
        xs.append(x)
        ws.append(w)
    TS_cache[(prec, m)] = (xs, ws)
    return xs, ws

def TS_eval(f, nodes, target_prec, working_prec, m, verbose=False):
    """Evaluate f at the given set of tanh-sinh nodes."""
    eps = mpf((1, -target_prec))
    S = mpf(0)
    h = mpf(1)
    xs, ws = nodes
    res = []
    for k in xrange(1, m+1):
        if verbose:
            print "evaluating integral (level %s of %s)" % (k, m)
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
                print "the estimated error is", smallstr(err)
            if err <= eps:
                break
    return +res[-1], TS_estimate_error(res, target_prec, eps)


def TS_adaptive(f, target_prec, working_prec, min_level, max_level, verbose):
    """Repeatedly attempt to integrate f, trying different levels. Quit
    as soon as the estimated error is small enough, or if that doesn't
    happen, when the max level has been tried."""
    eps = mpf((1, -target_prec))
    for m in xrange(min_level, max_level+1):
        if verbose:
            print "using tanh-sinh algorithm with m =", m
        nodes = TS_nodes(working_prec, m, verbose=verbose)
        s, err = TS_eval(f, nodes, target_prec, working_prec, m,
            verbose=verbose)
        steps = 2*len(nodes[0])
        if err <= eps:
            return s, err
    if verbose:
        print "Warning: failed to reach full accuracy. Estimated error:", \
            smallstr(err)
    return s, err


def quadts(f, a, b, **options):
    """
    Integrate f(x) dx over the interval [a, b], using tanh-sinh
    quadrature. Use quadts(f, (a, b), (c, d)) to calculate the
    two-dimensional integral over [a, b] x [c, d].

    Options
    =======

    target_prec
        The number of accurate bits to aim for in the result. If not
        specified, the current working precision mpf.prec is used.

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
        Return error estimate along with the result.

    verbose
        Set to True to display progress messages.
    """

    verbose = options.get('verbose', False)
    target_prec = options.get('target_prec', mpf.prec)
    working_prec = options.get('working_prec', target_prec + 30)
    min_level = options.get('min_level', TS_guess_level(target_prec))
    max_level = options.get('max_level', min_level + 2)

    oldprec = mpf.prec
    mpf.prec = working_prec

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

    mpf.prec = oldprec
    if options.get('error', False):
        return val, err
    return val

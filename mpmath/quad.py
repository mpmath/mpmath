from mpmath import *

# Utilities

def transform(f, a, b):
    a = mpf(a)
    b = mpf(b)
    if (a, b) == (-1, 1):
        return f
    if (a, b) == (-inf, inf):
        return transform(lambda x: (f(-1/x+1)+f(1/x-1))/x**2, 0, 1)
    if a == -inf: return transform(lambda x: f(-1/x+b+1)/x**2, 0, 1)
    if b == inf: return transform(lambda x: f(1/x+a-1)/x**2, 0, 1)
    C = (b-a)/2
    D = (b+a)/2
    def g(x):
        return C * f(D + C*x)
    return g

def smallstr(f):
    prec = mpf.prec
    mpf.dps = 5
    s = str(f)
    mpf.prec = prec
    return s

# Tanh-sinh (doubly exponential) quadrature

def TS_node(k, h, prec):
    # x[k] = tanh(pi/2 * sinh(k*h))
    # w[k] = pi/2 * cosh(k*h) / cosh(pi/2 sinh(k*h))**2
    oldprec = mpf.prec; mpf.prec = prec
    mpf.round_up()
    t = mpf(k) * h
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
    if (prec, m) in TS_cache:
        return TS_cache[(prec, m)]
    eps = mpf((1, -prec))
    h = mpf((1, -m))
    prec2 = prec + 10
    xs = []
    ws = []
    for k in xrange(20 * 2**m + 1):
        x, w = TS_node(k, h, prec2)
        xs.append(x)
        ws.append(w)
        diff = abs(xs[-1]-1)
        if diff <= eps:
            break
        if verbose and m > 6 and k % 300 == 150:
            print "calculating nodes:", smallstr(-log(diff, 10) / prec2)
    TS_cache[(prec, m)] = (xs, ws)
    return xs, ws

def TS_eval(f, nodes, m, verbose=False):
    prec = mpf.prec
    eps = mpf((1, -prec))
    mpf.prec += 30
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
            err = TS_estimate_error(res, prec, eps)
            if verbose:
                print "the estimated error is", smallstr(err)
            if err <= eps:
                break
    mpf.prec -= 30
    return +res[-1], TS_estimate_error(res, prec, eps)



def TS_estimate_error(res, prec, eps):
    try:
        D1 = log(abs(res[-1]-res[-2]), 10)
        D2 = log(abs(res[-1]-res[-3]), 10)
    except ValueError:
        return eps
    D3 = -prec
    D4 = min(0, max(D1**2/D2, 2*D1, D3))
    return mpf('0.1') ** -int(D4)



def TS_guess_level(prec):
    return int(4 + max(0, log(prec/30.0, 2)))


def quadts(f, a, b, **options):
    """
    Integrate f(x) dx over the interval [a, b], using "doubly
    exponential" or "tanh-sinh" quadrature. This algorithm tends to
    work exceptionally well at high precision levels when the
    integrand is singular at the endpoints or when the interval
    [a, b] is infinite.

    The algorithm is nonadaptive, however, and does not perform well if
    there are singularities between the endpoints or if the integrand
    is bumpy or oscillatory.
    """
    verbose=options.get('verbose', False)
    extra_level = 0

    if isinstance(a, tuple):
        (a, b), (c, d) = a, b
        if a == b or c == d:
            return mpf(0)
        g = f
        def f(y):
            return quadts(lambda x: g(x,y), a, b)
        if inf in (abs(c), abs(d)):
            extra_level = 1
        f = transform(f, c, d)
    else:
        if a == b:
            return mpf(0)
        if inf in (abs(a), abs(b)):
            extra_level = 1
        f = transform(f, a, b)
    prec = mpf.prec
    m = TS_guess_level(prec) + extra_level
    nodes = TS_nodes(prec, m, verbose)
    val, err = TS_eval(f, nodes, m, verbose)
    return val

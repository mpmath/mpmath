# Implementation overview
# -----------------------
#
# heung evaluates Hl (DLMF 31.3.1), the solution of the general Heun equation
# (DLMF 31.2.1) with H(0) = 1, by analytic continuation along the straight
# segment from 0 to z. The first step uses the Maclaurin series of Hl (radius
# min(1, |a|)) over at most a third of that radius. Every later step is a
# Taylor series about the current center w, started from (H(w), H'(w)), with
# a displacement of at most one third of the distance from w to the nearest
# finite singular point. The plane cut along the rays
# {t : t >= 1} and {t*a : t >= 1} is star-shaped about 0, so the segment never
# crosses a cut.
#
# Cancellation can lose many more bits than any fixed number of guard bits
# (a large delta, or a path that passes close to a singular point), so the
# working precision adapts. The whole path is evaluated at two precisions that
# differ by CHECK_BITS. The result is accepted when both values agree to the
# caller's precision. It is returned as 0 when it lies within an estimate of
# the rounding noise (2**-wp times the largest series sum on the path) and
# that estimate is below 2**-(2*prec + 40) (see heung). The estimate does not
# include error growth along the path. Otherwise the guard bits double, up to
# maxprec.
#
# The origin recurrence is DLMF 31.3.3 with n = j + 1. The recurrence at an
# ordinary point is the plain Taylor recurrence of the equation multiplied by
# z(z-1)(z-a).

from .functions import defun

CHECK_BITS = 20
# Rounding errors add up over the steps of a path, so the noise estimate
# 2**-wp * scale is widened by this many bits. The estimate still does not
# include error growth along the path.
NOISE_BITS = 12
INITIAL_GUARD_BITS = 30


def _min_terms(a, q, alpha, beta, gamma, delta, epsilon):
    # Terms can grow for about as many indices as the parameters are large,
    # so a series may not stop on small terms before this index.
    size = (abs(q) + abs(a*delta) + abs(epsilon) + abs(gamma) + abs(alpha) +
            abs(beta))
    return int(size) + 2


def _origin_terms(ctx, h, a, q, alpha, beta, gamma, c):
    # Yields t[1], t[2], ... where t[n] = b[n]*h**n and the Maclaurin
    # coefficients b[n] satisfy (b[-1] = 0, b[0] = 1)
    #   a*n*(gamma+n-1)*b[n] = A_n*b[n-1] - B_n*b[n-2],
    #   A_n = q + (n-1)*((a+1)*(gamma+n-2) + c),  c = epsilon + a*delta,
    #   B_n = (n-2+alpha)*(n-2+beta).
    a1 = a + 1
    h2 = h*h
    prev2 = ctx.zero
    prev1 = ctx.one
    n = 0
    while True:
        n += 1
        A = q + (n-1)*(a1*(gamma+n-2) + c)
        B = (n-2+alpha)*(n-2+beta)
        t = (A*h*prev1 - B*h2*prev2) / (a*n*(gamma+n-1))
        yield t
        prev2 = prev1
        prev1 = t


def _ordinary_terms(ctx, w, h, tn, tn1, a, q, alpha, beta, gamma, delta,
                    epsilon):
    # Starts from tn = t[0] and tn1 = t[1]. Yields t[2], t[3], ... where
    # t[n] = c[n]*h**n and H(w+u) = sum c[n] u**n.
    # With P = z(z-1)(z-a), Q = gamma(z-1)(z-a) + delta*z(z-a) + epsilon*z(z-1)
    # and T = alpha*beta*z - q, the equation is P H'' + Q H' + T H = 0. Write
    # P(w+u) = p0 + p1 u + p2 u^2 + u^3, Q(w+u) = Q0 + Q1 u + Q2 u^2 and
    # T(w+u) = r0 + r1 u. Then, with c[-1] = 0, for n >= 0:
    #   p0 (n+2)(n+1) c[n+2] = - (n+1) (n p1 + Q0) c[n+1]
    #       - (n(n-1) p2 + n Q1 + r0) c[n]
    #       - ((n-1)(n-2) + (n-1) Q2 + r1) c[n-1]
    p0 = w*(w-1)*(w-a)
    p1 = 3*w*w - 2*(1+a)*w + a
    p2 = 3*w - (1+a)
    Q2 = gamma + delta + epsilon
    k = gamma*(1+a) + delta*a + epsilon
    Q1 = 2*Q2*w - k
    Q0 = (Q2*w - k)*w + gamma*a
    ab = alpha*beta
    r0 = ab*w - q
    r1 = ab
    h2 = h*h
    h3 = h2*h
    tnm1 = ctx.zero
    n = 0
    while True:
        t = (-(n+1)*(n*p1 + Q0)*h*tn1
             - (n*(n-1)*p2 + n*Q1 + r0)*h2*tn
             - ((n-1)*(n-2) + (n-1)*Q2 + r1)*h3*tnm1) / (p0*(n+2)*(n+1))
        yield t
        tnm1 = tn
        tn = tn1
        tn1 = t
        n += 1


def _series(ctx, terms, n, S0, S1, eps, run, nmin, maxterms):
    # Sums terms whose first yielded term has index n + 1. S0 = sum t[m] and
    # S1 = sum m*t[m] for m <= n. Stops after `run` consecutive terms that are
    # exactly zero (the recurrence is then dead), or, once n >= nmin, after
    # `run` consecutive terms with m*|t[m]| <= eps*max(|S0|, |S1|). The caller
    # sets `run` above the number of earlier terms the recurrence uses (2 for
    # the origin series, 3 for an ordinary center), so run zeros end a series.
    # Raises NoConvergence before summing a term with index above maxterms.
    small = zeros = 0
    for t in terms:
        n += 1
        if n > maxterms:
            raise ctx.NoConvergence(f"heung: more than maxterms={maxterms} "
                                    "terms needed in one series")
        S0 += t
        S1 += n*t
        zeros = zeros + 1 if not t else 0
        if zeros >= run:
            break
        if n*abs(t) <= eps*max(abs(S0), abs(S1)):
            small += 1
            if small >= run and n >= nmin:
                break
        else:
            small = 0
    return S0, S1


def _next_center(ctx, w, z, lim):
    # Returns (center, h) for the step from w toward z, at most lim long. The
    # exact differences keep w + h equal to the next center, so the series is
    # summed at the point the next step starts from.
    rem = ctx.fsub(z, w, exact=True)
    rr = abs(rem)
    if rr <= lim:
        return z, rem
    wn = w + rem*(lim/rr)
    return wn, ctx.fsub(wn, w, exact=True)


def _on_ray(ctx, z, b):
    # True when z = t*b for a real t >= 1, tested exactly (no rounding).
    # That holds exactly when z*conj(b) = pr + i*pim has pim == 0 and
    # pr >= |b|**2.
    # .real and .imag return the parts unrounded, unlike ctx.re and ctx.im.
    zr, zi, br, bi = z.real, z.imag, b.real, b.imag
    pr = ctx.fadd(ctx.fmul(zr, br, exact=True), ctx.fmul(zi, bi, exact=True),
                  exact=True)
    pim = ctx.fsub(ctx.fmul(zi, br, exact=True), ctx.fmul(zr, bi, exact=True),
                  exact=True)
    b2 = ctx.fadd(ctx.fmul(br, br, exact=True), ctx.fmul(bi, bi, exact=True),
                  exact=True)
    return pim == 0 and pr >= b2


def _evaluate(ctx, wp, a, q, alpha, beta, gamma, delta, z, maxterms):
    # Evaluates H(z) along the segment 0 -> z at working precision wp and
    # returns (H(z), scale), where scale is the largest |H| or |h*H'| of any
    # series sum on the path. 2**-wp * scale estimates the absolute rounding
    # noise, without the growth of errors along the path. Returns None when a
    # step can no longer be resolved at this precision, so that the caller
    # retries with more bits. Sets ctx.prec to wp; the caller restores it.
    # Raises NoConvergence when a series needs more than maxterms terms.
    ctx.prec = wp
    eps = ctx.ldexp(ctx.one, -wp)
    epsilon = 1 + alpha + beta - gamma - delta
    nmin = _min_terms(a, q, alpha, beta, gamma, delta, epsilon)

    # First step: Maclaurin series.
    w, h = _next_center(ctx, ctx.zero, z, min(ctx.one, abs(a))/3)
    terms = _origin_terms(ctx, h, a, q, alpha, beta, gamma,
                          epsilon + a*delta)
    S0, S1 = _series(ctx, terms, 0, ctx.one, ctx.zero, eps, 3, nmin,
                     maxterms)
    scale = max(abs(S0), abs(S1))
    H = S0
    Hp = S1/h

    # Later steps: Taylor series about the current center.
    while w != z:
        d = min(abs(w), abs(w-1), abs(w-a))
        wn, h = _next_center(ctx, w, z, d/3)
        if wn == w:
            return None
        t1 = Hp*h
        terms = _ordinary_terms(ctx, w, h, H, t1, a, q, alpha, beta, gamma,
                                delta, epsilon)
        S0, S1 = _series(ctx, terms, 1, H + t1, t1, eps, 4, nmin, maxterms)
        scale = max(scale, abs(S0), abs(S1))
        w = wn
        H = S0
        Hp = S1/h
    return H, scale


@defun
def heung(ctx, a, q, alpha, beta, gamma, delta, z, *, maxprec=None,
          maxterms=None):
    r"""
    Evaluates the general Heun function `H_\ell(a, q; \alpha, \beta, \gamma,
    \delta; z)` (DLMF 31.3.1), the solution of the Heun equation

    .. math ::

        \frac{d^2 H}{dz^2} + \left(\frac{\gamma}{z} + \frac{\delta}{z-1} +
        \frac{\epsilon}{z-a}\right) \frac{dH}{dz} +
        \frac{\alpha\beta z - q}{z(z-1)(z-a)} H = 0,
        \qquad \epsilon = \alpha + \beta + 1 - \gamma - \delta,

    that is analytic at `z = 0` and normalized by `H(0) = 1`. The parameters
    are `a` (the finite singular point other than 0 and 1, `a \ne 0, 1`), the
    accessory parameter `q`, and the exponent parameters `\alpha`, `\beta`,
    `\gamma`, `\delta`. The function is continued analytically from the
    origin to `z` along the straight segment, in the plane cut along the two
    rays
    `\{t : t \in \mathbb{R}, t \ge 1\}` and
    `\{t a : t \in \mathbb{R}, t \ge 1\}`.

    The power series at the origin is used for `|z| \le \min(1, |a|)/3`.
    Farther out, the solution is continued with Taylor series about points on
    the segment, each step at most one third of the distance to the nearest
    singular point. The working precision is increased automatically until
    two evaluations at different precisions agree, or until the value is
    within an estimate of the rounding noise and that estimate is below
    `2^{-(2p+40)}` (see below).

    A ``ValueError`` is raised if `a` is 0 or 1, if `\gamma` is a
    nonpositive integer (the solution at the origin is then not defined by
    this normalization), if any argument is infinite or NaN, or if `z` lies
    exactly on one of the two cuts, including their endpoints `1` and `a`.
    The test for the cuts is exact: a point one unit in the last place away
    from a cut is accepted. Evaluation very close to a cut can be slow.

    A true value whose absolute size is below about `2^{-(2p+40)}`, where
    `p` is the precision in bits, may be returned as exactly zero, because it
    cannot be told apart from rounding noise. This applies when the true
    value is that small, for example near a root of a polynomial solution.
    Larger values are not returned as zero. The noise estimate does not
    include the growth of rounding errors along the path, so an exact zero
    is returned as 0 only when the estimate is small enough. An exact zero
    reached along a path that passes close to a singular point can raise
    ``NoConvergence`` instead.

    The keyword arguments *maxprec* (largest working precision in bits) and
    *maxterms* (largest term index in one series) limit the work, and
    ``NoConvergence`` is raised when either is exceeded. By default *maxprec*
    is the limit :func:`~mpmath.hyp2f1` uses and *maxterms* is
    ``2*maxprec + 1000``. The contexts ``fp`` and ``iv`` raise
    ``NotImplementedError``.

    **Examples**

    A value inside the disk of convergence of the power series::

        >>> from mpmath import *
        >>> mp.dps = 25; mp.pretty = True
        >>> heung(2, 1/mpf(3), 0.5, 0.75, 2/mpf(3), 4/mpf(5), 0.25)
        1.071417042933064090803712

    A value outside that disk, at a complex argument::

        >>> heung(2, 1/mpf(3), 0.5, 0.75, 2/mpf(3), 4/mpf(5),
        ...       mpc(1.25, 0.25))
        (1.382209817660756516579802 + 0.5411476473179165717001981j)

    With `a = 4`, `q = 9/4`, `\alpha = \beta = 3/2`, `\gamma = 1/2`,
    `\delta = 2` the function is elementary,
    `2 / ((1-z) \sqrt{4-z})`::

        >>> z = mpf(1)/5
        >>> heung(4, 9/mpf(4), 1.5, 1.5, 0.5, 2, z)
        1.282472940106442619320844
        >>> 2/((1-z)*sqrt(4-z))
        1.282472940106442619320844

    If `q = a \alpha \beta` and `\delta = \alpha + \beta + 1 - \gamma`, then
    `\epsilon = 0` and the function reduces to the Gauss hypergeometric
    function, `H_\ell = {}_2F_1(\alpha, \beta; \gamma; z)` (DLMF 31.7.1)::

        >>> mp.dps = 15
        >>> a = mpf(3); A = mpf(1)/3; B = mpf(5)/7; G = mpf(3)/4
        >>> heung(a, a*A*B, A, B, G, A+B+1-G, 0.3)
        1.11963268820494
        >>> hyp2f1(A, B, G, 0.3)
        1.11963268820494

    The second solution at the origin, `z^{1-\gamma} H_\ell(a, q +
    (1-\gamma)(\epsilon + a\delta); \alpha+1-\gamma, \beta+1-\gamma, 2-\gamma,
    \delta; z)` (DLMF 31.3.5), is also available through ``heung``::

        >>> a, q = 2, mpf(1)/3
        >>> alpha, beta, gamma, delta = 0.5, 0.75, mpf(2)/3, 4/mpf(5)
        >>> epsilon = 1 + alpha + beta - gamma - delta
        >>> s = 1 - gamma
        >>> z = mpf(1)/4 + 1j/4
        >>> z**s * heung(a, q + s*(epsilon + a*delta), alpha + s, beta + s,
        ...              2 - gamma, delta, z)
        (0.724567559984042 + 0.297155917497087j)

    **References**

    1. [DLMF]_ chapter 31, Heun functions.

    2. O. V. Motygin, "On evaluation of the Heun functions",
       arXiv:1506.03848 (2015).
    """
    # A top-level import is circular: ctx_mp imports ctx_base, which imports
    # this package.
    from ..ctx_mp import MPContext
    if not isinstance(ctx, MPContext):
        raise NotImplementedError("heung is only implemented for mp contexts")
    args = [ctx.convert(x) for x in (a, q, alpha, beta, gamma, delta, z)]
    for x in args:
        if not ctx.isfinite(x):
            raise ValueError("heung: arguments must be finite")
    a, q, alpha, beta, gamma, delta, z = args
    if a == 0 or a == 1:
        raise ValueError("heung: a must differ from 0 and 1")
    if ctx.isnpint(gamma):
        raise ValueError("heung: gamma must not be a nonpositive integer")

    if _on_ray(ctx, z, ctx.one):
        raise ValueError("heung: z lies on the cut [1, inf)")
    if _on_ray(ctx, z, a):
        raise ValueError("heung: z lies on the cut {t*a : t >= 1}")

    # Computation stays real when every argument is real. The result type is
    # mpc whenever any argument is an mpc, even with zero imaginary part.
    cplx = any(isinstance(x, ctx.mpc) for x in args)
    if not cplx or all(x.imag == 0 for x in args):
        a, q, alpha, beta, gamma, delta, z = [x.real for x in args]
    if not z:
        return ctx.mpc(1) if cplx else ctx.one

    prec = ctx.prec
    if maxprec is None:
        maxprec = ctx._default_hyper_maxprec(prec)
    if maxterms is None:
        maxterms = 2*maxprec + 1000
    try:
        g = INITIAL_GUARD_BITS
        while True:
            wp1 = prec + g
            wp2 = wp1 + CHECK_BITS
            if wp2 > maxprec:
                raise ctx.NoConvergence("heung: working precision above "
                                        f"maxprec={maxprec} needed")
            r1 = _evaluate(ctx, wp1, a, q, alpha, beta, gamma, delta, z,
                           maxterms)
            r2 = None
            if r1 is not None:
                r2 = _evaluate(ctx, wp2, a, q, alpha, beta, gamma, delta, z,
                               maxterms)
            if r2 is not None:
                (v1, _), (v2, s2) = r1, r2
                # Both values come from the same path at different
                # precisions, so agreement means the guard bits sufficed.
                if abs(v1 - v2) <= abs(v2)*ctx.ldexp(ctx.one, -(prec+2)):
                    break
                # At a zero of Hl only rounding noise remains, and two noisy
                # values rarely agree. Return 0 only when v2 is within the
                # noise estimate and the estimate is below 2**-(2*prec + 40),
                # far below the scale H(0) = 1. A higher estimate cannot tell
                # a zero from too few guard bits or from a small true value,
                # so the guard bits keep doubling. The estimate misses error
                # growth along the path, so for an exact zero reached near a
                # singular point it may never get small enough, and the loop
                # ends in NoConvergence.
                noise = s2*ctx.ldexp(ctx.one, NOISE_BITS - wp2)
                if abs(v2) <= noise <= ctx.ldexp(ctx.one, -(2*prec + 40)):
                    v2 = ctx.zero
                    break
            g *= 2
    finally:
        ctx.prec = prec
    result = +v2
    return ctx.mpc(result) if cplx else result

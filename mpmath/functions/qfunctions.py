from functions import defun, defun_wrapped

@defun
def qp(ctx, a, q=None, n=None, **kwargs):
    r"""
    Evaluates the q-Pochhammer symbol (or q-rising factorial)

    .. math ::

        (a; q)_n = \prod_{k=0}^{n-1} (1-a q^k)

    where `n = \infty` is permitted if `|q| < 1`. Called with two arguments,
    ``qp(a,q)`` computes `(a;q)_{\infty}`; with a single argument, ``qp(q)``
    computes `(q;q)_{\infty}`. The special case

    .. math ::

        \phi(q) = (q; q)_{\infty} = \prod_{k=1}^{\infty} (1-q^k) =
            \sum_{k=-\infty}^{\infty} (-1)^k q^{(3k^2-k)/2}

    is also known as the Euler function, or (up to a factor `q^{-1/24}`)
    the Dirichlet eta function.

    **Examples**

    If `n` is a positive integer, the function amounts to a finite product::

        >>> from mpmath import *
        >>> mp.dps = 25; mp.pretty = True
        >>> qp(2,3,5)
        -725305.0
        >>> fprod(1-2*3**k for k in range(5))
        -725305.0

    Complex arguments are allowed::

        >>> qp(2-1j, 0.75j)
        (0.4628842231660149089976379 + 4.481821753552703090628793j)

    The regular Pochhammer symbol `(a)_n` is obtained in the
    following limit as `q \to 1`::

        >>> a, n = 4, 7
        >>> limit(lambda q: qp(q**a,q,n) / (1-q)**n, 1)
        604800.0
        >>> rf(a,n)
        604800.0

    The Taylor series of the reciprocal Euler function gives
    the partition function `P(n)`, i.e. the number of ways of writing
    `n` as a sum of positive integers::

        >>> taylor(lambda q: 1/qp(q), 0, 10)
        [1.0, 1.0, 2.0, 3.0, 5.0, 7.0, 11.0, 15.0, 22.0, 30.0, 42.0]

    Special values include::

        >>> qp(0)
        1.0
        >>> findroot(diffun(qp), -0.4)   # location of maximum
        -0.4112484791779547734440257
        >>> qp(_)
        1.228348867038575112586878

    The q-Pochhammer symbol is related to the Jacobi theta functions.
    For example, the following identity holds::

        >>> q = mpf(0.5)    # arbitrary
        >>> qp(q)
        0.2887880950866024212788997
        >>> root(3,-2)*root(q,-24)*jtheta(2,pi/6,root(q,6))
        0.2887880950866024212788997

    """
    a = ctx.convert(a)
    if n is None:
        n = ctx.inf
    else:
        n = ctx.convert(n)
    assert n >= 0
    if q is None:
        q = a
    else:
        q = ctx.convert(q)
    infinite = (n == ctx.inf)
    same = (a == q)
    if infinite:
        if abs(q) >= 1:
            if same and (q == -1 or q == 1):
                return ctx.zero * q
            raise ValueError("q-function only defined for |q| < 1")
        elif q == 0:
            return ctx.one - a
    maxterms = kwargs.get('maxterms', 50*ctx.prec)
    if infinite and same:
        # Euler's pentagonal theorem
        def terms():
            t = 1
            yield t
            k = 1
            x1 = q
            x2 = q**2
            while 1:
                yield (-1)**k * x1
                yield (-1)**k * x2
                x1 *= q**(3*k+1)
                x2 *= q**(3*k+2)
                k += 1
                if k > maxterms:
                    raise ctx.NoConvergence
        return ctx.sum_accurately(terms)
    # return ctx.nprod(lambda k: 1-a*q**k, [0,n-1])
    def factors():
        k = 0
        r = ctx.one
        while 1:
            yield 1 - a*r
            r *= q
            k += 1
            if k >= n:
                raise StopIteration
            if k > maxterms:
                raise ctx.NoConvergence
    return ctx.mul_accurately(factors)

@defun_wrapped
def qgamma(ctx, z, q, **kwargs):
    return ctx.qp(q, q, None, **kwargs) / \
        ctx.qp(q**z, q, None, **kwargs) * (1-q)**(1-z)

@defun_wrapped
def qfac(ctx, z, q, **kwargs):
    if ctx.isint(z) and ctx._re(z) > 0:
        n = int(ctx._re(z))
        print "yass", q, q, n
        return ctx.qp(q, q, n, **kwargs) / (1-q)**n
    return ctx.qgamma(z+1, q, **kwargs)

@defun
def qhyper(ctx, a_s, b_s, q, z, **kwargs):
    r"""
    Evaluates the basic hypergeometric series or hypergeometric q-series

    .. math ::

        \,_r\phi_s \left[\begin{matrix} 
            a_1 & a_2 & \ldots & a_r \\ 
            b_1 & b_2 & \ldots & b_s
        \end{matrix} ; q,z \right] =
        \sum_{n=0}^\infty
        \frac{(a_1;q)_n, \ldots, (a_r;q)_n}
             {(b_1;q)_n, \ldots, (b_s;q)_n}
        \left((-1)^n q^{n\choose 2}\right)^{1+s-r}
        \frac{z^n}{(q;q)_n}

    where `(a;q)_n` denotes the q-Pochhammer symbol (see :func:`qp`).
    """
    #a_s = [ctx._convert_param(a)[0] for a in a_s]
    #b_s = [ctx._convert_param(b)[0] for b in b_s]
    #q = ctx._convert_param(q)[0]
    a_s = map(ctx.convert, a_s)
    b_s = map(ctx.convert, b_s)
    q = ctx.convert(q)
    z = ctx.convert(z)
    r = len(a_s)
    s = len(b_s)
    d = 1+s-r
    maxterms = kwargs.get('maxterms', 50*ctx.prec)
    def terms():
        t = ctx.one
        yield t
        qk = 1
        k = 0
        x = 1
        while 1:
            for a in a_s:
                p = 1 - a*qk
                t *= p
            for b in b_s:
                p = 1 - b*qk
                if not p:
                    raise ValueError
                t /= p
            t *= z
            x *= (-1)**d * qk ** d
            qk *= q
            t /= (1 - qk)
            k += 1
            yield t * x
            if k > maxterms:
                raise ctx.NoConvergence
    return ctx.sum_accurately(terms)

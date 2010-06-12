from functions import defun, defun_wrapped

@defun_wrapped
def hermite(ctx, n, z, **kwargs):
    if not z:
        try:
            return 2**n * ctx.sqrt(ctx.pi) / ctx.gamma(0.5*(1-n))
        except ValueError:
            return 0.0*(n+z)
    if ctx.re(z) > 0 or (ctx.re(z) == 0 and ctx.im(z) > 0) or ctx.isnpint(-n):
        prec = ctx.prec
        ctx.prec = ctx.prec*4+20
        z2 = -z**(-2)
        ctx.prec = prec
        return (2*z)**n * ctx.hyp2f0(-0.5*n, -0.5*(n-1), z2, **kwargs)
    else:
        prec = ctx.prec
        ctx.prec = ctx.prec*4+20
        z2 = z**2
        ctx.prec = prec
        return ctx.hermite(n,-z) + 2**(n+2)*ctx.sqrt(ctx.pi) * (-z) / \
            ctx.gamma(-0.5*n) * ctx.hyp1f1((1-n)*0.5, 1.5, z2, **kwargs)

@defun_wrapped
def gegenbauer(ctx, n, a, z, **kwargs):
    # Special cases: a+0.5, a*2 poles
    if ctx.isnpint(a):
        return 0*(z+n)
    if ctx.isnpint(a+0.5):
        # TODO: something else is required here
        # E.g.: gegenbauer(-2, -0.5, 3) == -12
        if ctx.isnpint(n+1):
            raise NotImplementedError("Gegenbauer function with two limits")
        def h(a):
            a2 = 2*a
            T = [], [], [n+a2], [n+1, a2], [-n, n+a2], [a+0.5], 0.5*(1-z)
            return [T]
        return ctx.hypercomb(h, [a], **kwargs)
    def h(n):
        a2 = 2*a
        T = [], [], [n+a2], [n+1, a2], [-n, n+a2], [a+0.5], 0.5*(1-z)
        return [T]
    return ctx.hypercomb(h, [n], **kwargs)

@defun_wrapped
def jacobi(ctx, n, a, b, x, **kwargs):
    if not ctx.isnpint(a):
        def h(n):
            return (([], [], [a+n+1], [n+1, a+1], [-n, a+b+n+1], [a+1], (1-x)*0.5),)
        return ctx.hypercomb(h, [n], **kwargs)
    if not ctx.isint(b):
        def h(n, a):
            return (([], [], [-b], [n+1, -b-n], [-n, a+b+n+1], [b+1], (x+1)*0.5),)
        return ctx.hypercomb(h, [n, a], **kwargs)
    # XXX: determine appropriate limit
    return ctx.binomial(n+a,n) * ctx.hyp2f1(-n,1+n+a+b,a+1,(1-x)/2, **kwargs)

@defun_wrapped
def laguerre(ctx, n, a, z, **kwargs):
    # XXX: limits, poles
    #if ctx.isnpint(n):
    #    return 0*(a+z)
    def h(a):
        return (([], [], [a+n+1], [a+1, n+1], [-n], [a+1], z),)
    return ctx.hypercomb(h, [a], **kwargs)

@defun_wrapped
def legendre(ctx, n, x, **kwargs):
    if ctx.isint(n):
        n = int(n)
        # Accuracy near zeros
        if (n + (n < 0)) & 1:
            if not x:
                return x
            mag = ctx.mag(x)
            if mag < -2*ctx.prec-10:
                return x
            if mag < -5:
                ctx.prec += -mag
    return ctx.hyp2f1(-n,n+1,1,(1-x)/2, **kwargs)

@defun
def legenp(ctx, n, m, z, type=2, **kwargs):
    # Legendre function, 1st kind
    n = ctx.convert(n)
    m = ctx.convert(m)
    # Faster
    if not m:
        return ctx.legendre(n, z, **kwargs)
    # TODO: correct evaluation at singularities
    if type == 2:
        def h(n,m):
            g = m*0.5
            T = [1+z, 1-z], [g, -g], [], [1-m], [-n, n+1], [1-m], 0.5*(1-z)
            return (T,)
        return ctx.hypercomb(h, [n,m], **kwargs)
    if type == 3:
        def h(n,m):
            g = m*0.5
            T = [z+1, z-1], [g, -g], [], [1-m], [-n, n+1], [1-m], 0.5*(1-z)
            return (T,)
        return ctx.hypercomb(h, [n,m], **kwargs)
    raise ValueError("requires type=2 or type=3")

@defun
def legenq(ctx, n, m, z, type=2, **kwargs):
    # Legendre function, 2nd kind
    n = ctx.convert(n)
    m = ctx.convert(m)
    z = ctx.convert(z)
    if z in (1, -1):
        #if ctx.isint(m):
        #    return ctx.nan
        #return ctx.inf  # unsigned
        return ctx.nan
    if type == 2:
        def h(n, m):
            cos, sin = ctx.cospi_sinpi(m)
            s = 2 * sin / ctx.pi
            c = cos
            a = 1+z
            b = 1-z
            u = m/2
            w = (1-z)/2
            T1 = [s, c, a, b], [-1, 1, u, -u], [], [1-m], \
                [-n, n+1], [1-m], w
            T2 = [-s, a, b], [-1, -u, u], [n+m+1], [n-m+1, m+1], \
                [-n, n+1], [m+1], w
            return T1, T2
        return ctx.hypercomb(h, [n, m], **kwargs)
    if type == 3:
        # The following is faster when there only is a single series
        # Note: not valid for -1 < z < 0 (?)
        if abs(z) > 1:
            def h(n, m):
                T1 = [ctx.expjpi(m), 2, ctx.pi, z, z-1, z+1], \
                     [1, -n-1, 0.5, -n-m-1, 0.5*m, 0.5*m], \
                     [n+m+1], [n+1.5], \
                     [0.5*(2+n+m), 0.5*(1+n+m)], [n+1.5], z**(-2)
                return [T1]
            return ctx.hypercomb(h, [n, m], **kwargs)
        else:
            # not valid for 1 < z < inf ?
            def h(n, m):
                s = 2 * ctx.sinpi(m) / ctx.pi
                c = ctx.expjpi(m)
                a = 1+z
                b = z-1
                u = m/2
                w = (1-z)/2
                T1 = [s, c, a, b], [-1, 1, u, -u], [], [1-m], \
                    [-n, n+1], [1-m], w
                T2 = [-s, c, a, b], [-1, 1, -u, u], [n+m+1], [n-m+1, m+1], \
                    [-n, n+1], [m+1], w
                return T1, T2
            return ctx.hypercomb(h, [n, m], **kwargs)
    raise ValueError("requires type=2 or type=3")

@defun_wrapped
def chebyt(ctx, n, x, **kwargs):
    return ctx.hyp2f1(-n,n,(1,2),(1-x)/2, **kwargs)

@defun_wrapped
def chebyu(ctx, n, x, **kwargs):
    return (n+1) * ctx.hyp2f1(-n, n+2, (3,2), (1-x)/2, **kwargs)

@defun
def spherharm(ctx, l, m, theta, phi, **kwargs):
    l = ctx.convert(l)
    m = ctx.convert(m)
    theta = ctx.convert(theta)
    phi = ctx.convert(phi)
    l_isint = ctx.isint(l)
    l_natural = l_isint and l >= 0
    m_isint = ctx.isint(m)
    if l_isint and l < 0 and m_isint:
        return ctx.spherharm(-(l+1), m, theta, phi, **kwargs)
    if theta == 0 and m_isint and m < 0:
        return ctx.zero * 1j
    if l_natural and m_isint:
        if abs(m) > l:
            return ctx.zero * 1j
        # http://functions.wolfram.com/Polynomials/
        #     SphericalHarmonicY/26/01/02/0004/
        def h(l,m):
            absm = abs(m)
            C = [-1, ctx.expj(m*phi),
                 (2*l+1)*ctx.fac(l+absm)/ctx.pi/ctx.fac(l-absm),
                 ctx.sin(theta)**2,
                 ctx.fac(absm), 2]
            P = [0.5*m*(ctx.sign(m)+1), 1, 0.5, 0.5*absm, -1, -absm-1]
            return ((C, P, [], [], [absm-l, l+absm+1], [absm+1],
                ctx.sin(0.5*theta)**2),)
    else:
        # http://functions.wolfram.com/HypergeometricFunctions/
        #     SphericalHarmonicYGeneral/26/01/02/0001/
        def h(l,m):
            if ctx.isnpint(l-m+1) or ctx.isnpint(l+m+1) or ctx.isnpint(1-m):
                return (([0], [-1], [], [], [], [], 0),)
            cos, sin = ctx.cos_sin(0.5*theta)
            C = [0.5*ctx.expj(m*phi), (2*l+1)/ctx.pi,
                 ctx.gamma(l-m+1), ctx.gamma(l+m+1),
                 cos**2, sin**2]
            P = [1, 0.5, 0.5, -0.5, 0.5*m, -0.5*m]
            return ((C, P, [], [1-m], [-l,l+1], [1-m], sin**2),)
    return ctx.hypercomb(h, [l,m], **kwargs)

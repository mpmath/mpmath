from functions import defun, defun_wrapped

@defun
def j0(ctx, x):
    """Computes the Bessel function `J_0(x)`. See :func:`~mpmath.besselj`."""
    return ctx.besselj(0, x)

@defun
def j1(ctx, x):
    """Computes the Bessel function `J_1(x)`.  See :func:`~mpmath.besselj`."""
    return ctx.besselj(1, x)

@defun
def besselj(ctx, n, z, derivative=0, **kwargs):
    if type(n) is int:
        n_isint = True
    else:
        n = ctx.convert(n)
        n_isint = ctx.isint(n)
        if n_isint:
            n = int(n)
    if n_isint and n < 0:
        return (-1)**n * ctx.besselj(-n, z, derivative, **kwargs)
    z = ctx.convert(z)
    M = ctx.mag(z)
    if derivative:
        d = ctx.convert(derivative)
        # TODO: the integer special-casing shouldn't be necessary.
        # However, the hypergeometric series gets inaccurate for large d
        # because of inaccurate pole cancellation at a pole far from
        # zero (needs to be fixed in hypercomb or hypsum)
        if ctx.isint(d) and d >= 0:
            d = int(d)
            orig = ctx.prec
            try:
                ctx.prec += 15
                v = ctx.fsum((-1)**k * ctx.binomial(d,k) * ctx.besselj(2*k+n-d,z)
                    for k in range(d+1))
            finally:
                ctx.prec = orig
            v *= ctx.mpf(2)**(-d)
        else:
            def h(n,d):
                r = ctx.fmul(ctx.fmul(z, z, prec=ctx.prec+M), -0.25, exact=True)
                B = [0.5*(n-d+1), 0.5*(n-d+2)]
                T = [([2,ctx.pi,z],[d-2*n,0.5,n-d],[],B,[(n+1)*0.5,(n+2)*0.5],B+[n+1],r)]
                return T
            v = ctx.hypercomb(h, [n,d], **kwargs)
    else:
    # Fast case: J_n(x), n int, appropriate magnitude for fixed-point calculation
        if (not derivative) and n_isint and abs(M) < 10 and abs(n) < 20:
            try:
                return ctx._besselj(n, z)
            except NotImplementedError:
                pass
        if not z:
            if not n:
                v = ctx.one + n+z
            elif ctx.re(n) > 0:
                v = n*z
            else:
                v = ctx.inf + z + n
        else:
            #v = 0
            orig = ctx.prec
            try:
                # XXX: workaround for accuracy in low level hypergeometric series
                # when alternating, large arguments
                ctx.prec += min(3*abs(M), ctx.prec)
                w = ctx.fmul(z, 0.5, exact=True)
                def h(n):
                    r = ctx.fneg(ctx.fmul(w, w, prec=max(0,ctx.prec+M)), exact=True)
                    return [([w], [n], [], [n+1], [], [n+1], r)]
                v = ctx.hypercomb(h, [n], **kwargs)
            finally:
                ctx.prec = orig
        v = +v
    return v

@defun
def besseli(ctx, n, z, derivative=0, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    if not z:
        if derivative:
            raise ValueError
        if not n:
            # I(0,0) = 1
            return 1+n+z
        if ctx.isint(n):
            return 0*(n+z)
        r = ctx.re(n)
        if r == 0:
            return ctx.nan*(n+z)
        elif r > 0:
            return 0*(n+z)
        else:
            return ctx.inf+(n+z)
    M = ctx.mag(z)
    if derivative:
        d = ctx.convert(derivative)
        def h(n,d):
            r = ctx.fmul(ctx.fmul(z, z, prec=ctx.prec+M), 0.25, exact=True)
            B = [0.5*(n-d+1), 0.5*(n-d+2), n+1]
            T = [([2,ctx.pi,z],[d-2*n,0.5,n-d],[n+1],B,[(n+1)*0.5,(n+2)*0.5],B,r)]
            return T
        v = ctx.hypercomb(h, [n,d], **kwargs)
    else:
        def h(n):
            w = ctx.fmul(z, 0.5, exact=True)
            r = ctx.fmul(w, w, prec=max(0,ctx.prec+M))
            return [([w], [n], [], [n+1], [], [n+1], r)]
        v = ctx.hypercomb(h, [n], **kwargs)
    return v

@defun_wrapped
def bessely(ctx, n, z, derivative=0, **kwargs):
    if not z:
        if derivative:
            # Not implemented
            raise ValueError
        if not n:
            # ~ log(z/2)
            return -ctx.inf + (n+z)
        if ctx.im(n):
            return nan * (n+z)
        r = ctx.re(n)
        q = n+0.5
        if ctx.isint(q):
            if n > 0:
                return -ctx.inf + (n+z)
            else:
                return 0 * (n+z)
        if r < 0 and int(ctx.floor(q)) % 2:
            return ctx.inf + (n+z)
        else:
            return ctx.ninf + (n+z)
    # XXX: use hypercomb
    ctx.prec += 10
    m, d = ctx.nint_distance(n)
    if d < -ctx.prec:
        h = +ctx.eps
        ctx.prec *= 2
        n += h
    elif d < 0:
        ctx.prec -= d
    # TODO: avoid cancellation for imaginary arguments
    cos, sin = ctx.cospi_sinpi(n)
    return (ctx.besselj(n,z,derivative,**kwargs)*cos - \
        ctx.besselj(-n,z,derivative,**kwargs))/sin

@defun_wrapped
def besselk(ctx, n, z, **kwargs):
    if not z:
        return ctx.inf
    M = ctx.mag(z)
    if M < 1:
        # Represent as limit definition
        def h(n):
            r = (z/2)**2
            T1 = [z, 2], [-n, n-1], [n], [], [], [1-n], r
            T2 = [z, 2], [n, -n-1], [-n], [], [], [1+n], r
            return T1, T2
    # We could use the limit definition always, but it leads
    # to very bad cancellation (of exponentially large terms)
    # for large real z
    # Instead represent in terms of 2F0
    else:
        ctx.prec += M
        def h(n):
            return [([ctx.pi/2, z, ctx.exp(-z)], [0.5,-0.5,1], [], [], \
                [n+0.5, 0.5-n], [], -1/(2*z))]
    return ctx.hypercomb(h, [n], **kwargs)

@defun_wrapped
def hankel1(ctx,n,x,**kwargs):
    return ctx.besselj(n,x,**kwargs) + ctx.j*ctx.bessely(n,x,**kwargs)

@defun_wrapped
def hankel2(ctx,n,x,**kwargs):
    return ctx.besselj(n,x,**kwargs) - ctx.j*ctx.bessely(n,x,**kwargs)

@defun_wrapped
def whitm(ctx,k,m,z,**kwargs):
    if z == 0:
        # M(k,m,z) = 0^(1/2+m)
        if ctx.re(m) > -0.5:
            return z
        elif ctx.re(m) < -0.5:
            return ctx.inf + z
        else:
            return ctx.nan * z
    x = ctx.fmul(-0.5, z, exact=True)
    y = 0.5+m
    return ctx.exp(x) * z**y * ctx.hyp1f1(y-k, 1+2*m, z, **kwargs)

@defun_wrapped
def whitw(ctx,k,m,z,**kwargs):
    if z == 0:
        g = abs(ctx.re(m))
        if g < 0.5:
            return z
        elif g > 0.5:
            return ctx.inf + z
        else:
            return ctx.nan * z
    x = ctx.fmul(-0.5, z, exact=True)
    y = 0.5+m
    return ctx.exp(x) * z**y * ctx.hyperu(y-k, 1+2*m, z, **kwargs)

@defun
def hyperu(ctx, a, b, z, **kwargs):
    a, atype = ctx._convert_param(a)
    b, btype = ctx._convert_param(b)
    z = ctx.convert(z)
    if not z:
        if ctx.re(b) <= 1:
            return ctx.gammaprod([1-b],[a-b+1])
        else:
            return ctx.inf + z
    bb = 1+a-b
    bb, bbtype = ctx._convert_param(bb)
    try:
        orig = ctx.prec
        try:
            ctx.prec += 10
            v = ctx.hypsum(2, 0, (atype, bbtype), [a, bb], -1/z, maxterms=ctx.prec)
            return v / z**a
        finally:
            ctx.prec = orig
    except ctx.NoConvergence:
        pass
    def h(a,b):
        w = ctx.sinpi(b)
        T1 = ([ctx.pi,w],[1,-1],[],[a-b+1,b],[a],[b],z)
        T2 = ([-ctx.pi,w,z],[1,-1,1-b],[],[a,2-b],[a-b+1],[2-b],z)
        return T1, T2
    return ctx.hypercomb(h, [a,b], **kwargs)

@defun
def struveh(ctx,n,z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/StruveH/26/01/02/
    def h(n):
        return [([z/2, 0.5*ctx.sqrt(ctx.pi)], [n+1, -1], [], [n+1.5], [1], [1.5, n+1.5], -(z/2)**2)]
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def struvel(ctx,n,z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/StruveL/26/01/02/
    def h(n):
        return [([z/2, 0.5*ctx.sqrt(ctx.pi)], [n+1, -1], [], [n+1.5], [1], [1.5, n+1.5], (z/2)**2)]
    return ctx.hypercomb(h, [n], **kwargs)

def _anger(ctx,which,v,z,**kwargs):
    v = ctx._convert_param(v)[0]
    z = ctx.convert(z)
    def h(v):
        b = ctx.mpq_1_2
        u = v*b
        m = b*3
        a1,a2,b1,b2 = m-u, m+u, 1-u, 1+u
        c, s = ctx.cospi_sinpi(u)
        if which == 0:
            A, B = [b*z, s], [c]
        if which == 1:
            A, B = [b*z, -c], [s]
        w = ctx.square_exp_arg(z, mult=-0.25)
        T1 = A, [1, 1], [], [a1,a2], [1], [a1,a2], w
        T2 = B, [1], [], [b1,b2], [1], [b1,b2], w
        return T1, T2
    return ctx.hypercomb(h, [v], **kwargs)

@defun
def angerj(ctx, v, z, **kwargs):
    return _anger(ctx, 0, v, z, **kwargs)

@defun
def webere(ctx, v, z, **kwargs):
    return _anger(ctx, 1, v, z, **kwargs)

@defun
def lommels1(ctx, u, v, z, **kwargs):
    u = ctx._convert_param(u)[0]
    v = ctx._convert_param(v)[0]
    z = ctx.convert(z)
    def h(u,v):
        b = ctx.mpq_1_2
        w = ctx.square_exp_arg(z, mult=-0.25)
        return ([u-v+1, u+v+1, z], [-1, -1, u+1], [], [], [1], \
            [b*(u-v+3),b*(u+v+3)], w),
    return ctx.hypercomb(h, [u,v], **kwargs)

@defun
def lommels2(ctx, u, v, z, **kwargs):
    u = ctx._convert_param(u)[0]
    v = ctx._convert_param(v)[0]
    z = ctx.convert(z)
    # Asymptotic expansion (GR p. 947) -- need to be careful
    # not to use for small arguments
    # def h(u,v):
    #    b = ctx.mpq_1_2
    #    w = -(z/2)**(-2)
    #    return ([z], [u-1], [], [], [b*(1-u+v)], [b*(1-u-v)], w),
    def h(u,v):
        b = ctx.mpq_1_2
        w = ctx.square_exp_arg(z, mult=-0.25)
        T1 = [u-v+1, u+v+1, z], [-1, -1, u+1], [], [], [1], [b*(u-v+3),b*(u+v+3)], w
        T2 = [2, z], [u+v-1, -v], [v, b*(u+v+1)], [b*(v-u+1)], [], [1-v], w
        T3 = [2, z], [u-v-1, v], [-v, b*(u-v+1)], [b*(1-u-v)], [], [1+v], w
        #c1 = ctx.cospi((u-v)*b)
        #c2 = ctx.cospi((u+v)*b)
        #s = ctx.sinpi(v)
        #r1 = (u-v+1)*b
        #r2 = (u+v+1)*b
        #T2 = [c1, s, z, 2], [1, -1, -v, v], [], [-v+1], [], [-v+1], w
        #T3 = [-c2, s, z, 2], [1, -1, v, -v], [], [v+1], [], [v+1], w
        #T2 = [c1, s, z, 2], [1, -1, -v, v+u-1], [r1, r2], [-v+1], [], [-v+1], w
        #T3 = [-c2, s, z, 2], [1, -1, v, -v+u-1], [r1, r2], [v+1], [], [v+1], w
        return T1, T2, T3
    return ctx.hypercomb(h, [u,v], **kwargs)

@defun
def ber(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinBer2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        cos, sin = ctx.cospi_sinpi(-0.75*n)
        T1 = [cos, z/2], [1, n], [], [n+1], [], [0.5, 0.5*(n+1), 0.5*n+1], r
        T2 = [sin, z/2], [1, n+2], [], [n+2], [], [1.5, 0.5*(n+3), 0.5*n+1], r
        return T1, T2
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def bei(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinBei2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        cos, sin = ctx.cospi_sinpi(0.75*n)
        T1 = [cos, z/2], [1, n+2], [], [n+2], [], [1.5, 0.5*(n+3), 0.5*n+1], r
        T2 = [sin, z/2], [1, n], [], [n+1], [], [0.5, 0.5*(n+1), 0.5*n+1], r
        return T1, T2
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def ker(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinKer2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        cos1, sin1 = ctx.cospi_sinpi(0.25*n)
        cos2, sin2 = ctx.cospi_sinpi(0.75*n)
        T1 = [2, z, 4*cos1], [-n-3, n, 1], [-n], [], [], [0.5, 0.5*(1+n), 0.5*(n+2)], r
        T2 = [2, z, -sin1], [-n-3, 2+n, 1], [-n-1], [], [], [1.5, 0.5*(3+n), 0.5*(n+2)], r
        T3 = [2, z, 4*cos2], [n-3, -n, 1], [n], [], [], [0.5, 0.5*(1-n), 1-0.5*n], r
        T4 = [2, z, -sin2], [n-3, 2-n, 1], [n-1], [], [], [1.5, 0.5*(3-n), 1-0.5*n], r
        return T1, T2, T3, T4
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def kei(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinKei2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        cos1, sin1 = ctx.cospi_sinpi(0.75*n)
        cos2, sin2 = ctx.cospi_sinpi(0.25*n)
        T1 = [-cos1, 2, z], [1, n-3, 2-n], [n-1], [], [], [1.5, 0.5*(3-n), 1-0.5*n], r
        T2 = [-sin1, 2, z], [1, n-1, -n], [n], [], [], [0.5, 0.5*(1-n), 1-0.5*n], r
        T3 = [-sin2, 2, z], [1, -n-1, n], [-n], [], [], [0.5, 0.5*(n+1), 0.5*(n+2)], r
        T4 = [-cos2, 2, z], [1, -n-3, n+2], [-n-1], [], [], [1.5, 0.5*(n+3), 0.5*(n+2)], r
        return T1, T2, T3, T4
    return ctx.hypercomb(h, [n], **kwargs)

# TODO: do this more generically?
def c_memo(f):
    name = f.__name__
    def f_wrapped(ctx):
        cache = ctx._misc_const_cache
        prec = ctx.prec
        p,v = cache.get(name, (-1,0))
        if p >= prec:
            return +v
        else:
            cache[name] = (prec, f(ctx))
            return cache[name][1]
    return f_wrapped

@c_memo
def _airyai_C1(ctx):
    return 1 / (ctx.cbrt(9) * ctx.gamma(ctx.mpf(2)/3))

@c_memo
def _airyai_C2(ctx):
    return -1 / (ctx.cbrt(3) * ctx.gamma(ctx.mpf(1)/3))

@c_memo
def _airybi_C1(ctx):
    return 1 / (ctx.nthroot(3,6) * ctx.gamma(ctx.mpf(2)/3))

@c_memo
def _airybi_C2(ctx):
    return ctx.nthroot(3,6) / ctx.gamma(ctx.mpf(1)/3)

def _airybi_n2_inf(ctx):
    prec = ctx.prec
    try:
        v = ctx.power(3,'2/3')*ctx.gamma('2/3')/(2*ctx.pi)
    finally:
        ctx.prec = prec
    return +v

# Derivatives at z = 0
# TODO: could be expressed more elegantly using triple factorials
def _airyderiv_0(ctx, z, n, ntype, which):
    if ntype == 'Z':
        if n < 0:
            return z
        r = ctx.mpq_1_3
        prec = ctx.prec
        try:
            ctx.prec += 10
            v = ctx.gamma((n+1)*r) * ctx.power(3,n*r) / ctx.pi
            if which == 0:
                v *= ctx.sinpi(2*(n+1)*r)
                v /= ctx.power(3,'2/3')
            else:
                v *= abs(ctx.sinpi(2*(n+1)*r))
                v /= ctx.power(3,'1/6')
        finally:
            ctx.prec = prec
        return +v + z
    else:
        # singular (does the limit exist?)
        raise NotImplementedError

@defun
def airyai(ctx, z, derivative=0, **kwargs):
    z = ctx.convert(z)
    if derivative:
        n, ntype = ctx._convert_param(derivative)
    else:
        n = 0
    # Values at infinities
    if not ctx.isnormal(z) and z:
        if n and ntype == 'Z':
            if n == -1:
                if z == ctx.inf:
                    return ctx.mpf(1)/3 + 1/z
                if z == ctx.ninf:
                    return ctx.mpf(-2)/3 + 1/z
            if n < -1:
                if z == ctx.inf:
                    return z
                if z == ctx.ninf:
                    return (-1)**n * (-z)
        if (not n) and z == ctx.inf or z == ctx.ninf:
            return 1/z
        # TODO: limits
        raise ValueError("essential singularity of Ai(z)")
    # Account for exponential scaling
    if z:
        extraprec = max(0, int(1.5*ctx.mag(z)))
    else:
        extraprec = 0
    if n:
        if n == 1:
            def h():
                # http://functions.wolfram.com/03.07.06.0005.01
                if ctx._re(z) > 4:
                    ctx.prec += extraprec
                    w = z**1.5; r = -0.75/w; u = -2*w/3
                    ctx.prec -= extraprec
                    C = -ctx.exp(u)/(2*ctx.sqrt(ctx.pi))*ctx.nthroot(z,4)
                    return ([C],[1],[],[],[(-1,6),(7,6)],[],r),
                # http://functions.wolfram.com/03.07.26.0001.01
                else:
                    ctx.prec += extraprec
                    w = z**3 / 9
                    ctx.prec -= extraprec
                    C1 = _airyai_C1(ctx) * 0.5
                    C2 = _airyai_C2(ctx)
                    T1 = [C1,z],[1,2],[],[],[],[ctx.mpq_5_3],w
                    T2 = [C2],[1],[],[],[],[ctx.mpq_1_3],w
                    return T1, T2
            return ctx.hypercomb(h, [], **kwargs)
        else:
            if z == 0:
                return _airyderiv_0(ctx, z, n, ntype, 0)
            # http://functions.wolfram.com/03.05.20.0004.01
            def h(n):
                ctx.prec += extraprec
                w = z**3/9
                ctx.prec -= extraprec
                q13,q23,q43 = ctx.mpq_1_3, ctx.mpq_2_3, ctx.mpq_4_3
                a1=q13; a2=1; b1=(1-n)*q13; b2=(2-n)*q13; b3=1-n*q13
                T1 = [3, z], [n-q23, -n], [a1], [b1,b2,b3], \
                    [a1,a2], [b1,b2,b3], w
                a1=q23; b1=(2-n)*q13; b2=1-n*q13; b3=(4-n)*q13
                T2 = [3, z, -z], [n-q43, -n, 1], [a1], [b1,b2,b3], \
                    [a1,a2], [b1,b2,b3], w
                return T1, T2
            v = ctx.hypercomb(h, [n], **kwargs)
            if ctx._is_real_type(z) and ctx.isint(n):
                v = ctx._re(v)
            return v
    else:
        def h():
            if ctx._re(z) > 4:
                # We could use 1F1, but it results in huge cancellation;
                # the following expansion is better.
                # TODO: asymptotic series for derivatives
                ctx.prec += extraprec
                w = z**1.5; r = -0.75/w; u = -2*w/3
                ctx.prec -= extraprec
                C = ctx.exp(u)/(2*ctx.sqrt(ctx.pi)*ctx.nthroot(z,4))
                return ([C],[1],[],[],[(1,6),(5,6)],[],r),
            else:
                ctx.prec += extraprec
                w = z**3 / 9
                ctx.prec -= extraprec
                C1 = _airyai_C1(ctx)
                C2 = _airyai_C2(ctx)
                T1 = [C1],[1],[],[],[],[ctx.mpq_2_3],w
                T2 = [z*C2],[1],[],[],[],[ctx.mpq_4_3],w
                return T1, T2
        return ctx.hypercomb(h, [], **kwargs)

@defun
def airybi(ctx, z, derivative=0, **kwargs):
    z = ctx.convert(z)
    if derivative:
        n, ntype = ctx._convert_param(derivative)
    else:
        n = 0
    # Values at infinities
    if not ctx.isnormal(z) and z:
        if n and ntype == 'Z':
            if z == ctx.inf:
                return z
            if z == ctx.ninf:
                if n == -1:
                    return 1/z
                if n == -2:
                    return _airybi_n2_inf(ctx)
                if n < -2:
                    return (-1)**n * (-z)
        if not n:
            if z == ctx.inf:
                return z
            if z == ctx.ninf:
                return 1/z
        # TODO: limits
        raise ValueError("essential singularity of Bi(z)")
    if z:
        extraprec = max(0, int(1.5*ctx.mag(z)))
    else:
        extraprec = 0
    if n:
        if n == 1:
            # http://functions.wolfram.com/03.08.26.0001.01
            def h():
                ctx.prec += extraprec
                w = z**3 / 9
                ctx.prec -= extraprec
                C1 = _airybi_C1(ctx)*0.5
                C2 = _airybi_C2(ctx)
                T1 = [C1,z],[1,2],[],[],[],[ctx.mpq_5_3],w
                T2 = [C2],[1],[],[],[],[ctx.mpq_1_3],w
                return T1, T2
            return ctx.hypercomb(h, [], **kwargs)
        else:
            if z == 0:
                return _airyderiv_0(ctx, z, n, ntype, 1)
            def h(n):
                ctx.prec += extraprec
                w = z**3/9
                ctx.prec -= extraprec
                q13,q23,q43 = ctx.mpq_1_3, ctx.mpq_2_3, ctx.mpq_4_3
                q16 = ctx.mpq_1_6
                q56 = ctx.mpq_5_6
                a1=q13; a2=1; b1=(1-n)*q13; b2=(2-n)*q13; b3=1-n*q13
                T1 = [3, z], [n-q16, -n], [a1], [b1,b2,b3], \
                    [a1,a2], [b1,b2,b3], w
                a1=q23; b1=(2-n)*q13; b2=1-n*q13; b3=(4-n)*q13
                T2 = [3, z], [n-q56, 1-n], [a1], [b1,b2,b3], \
                    [a1,a2], [b1,b2,b3], w
                return T1, T2
            v = ctx.hypercomb(h, [n], **kwargs)
            if ctx._is_real_type(z) and ctx.isint(n):
                v = ctx._re(v)
            return v
    else:
        def h():
            ctx.prec += extraprec
            w = z**3 / 9
            ctx.prec -= extraprec
            C1 = _airybi_C1(ctx)
            C2 = _airybi_C2(ctx)
            T1 = [C1],[1],[],[],[],[ctx.mpq_2_3],w
            T2 = [z*C2],[1],[],[],[],[ctx.mpq_4_3],w
            return T1, T2
        return ctx.hypercomb(h, [], **kwargs)

def _airy_zero(ctx, which, k, derivative, complex=False):
    # Asymptotic formulas are given in DLMF section 9.9
    def U(t): return t**(2/3.)*(1-7/(t**2*48))
    def T(t): return t**(2/3.)*(1+5/(t**2*48))
    k = int(k)
    assert k >= 1
    assert derivative in (0,1)
    if which == 0:
        if derivative:
            return ctx.findroot(lambda z: ctx.airyai(z,1),
                -U(3*ctx.pi*(4*k-3)/8))
        return ctx.findroot(ctx.airyai, -T(3*ctx.pi*(4*k-1)/8))
    if which == 1 and complex == False:
        if derivative:
            return ctx.findroot(lambda z: ctx.airybi(z,1),
                -U(3*ctx.pi*(4*k-1)/8))
        return ctx.findroot(ctx.airybi, -T(3*ctx.pi*(4*k-3)/8))
    if which == 1 and complex == True:
        if derivative:
            t = 3*ctx.pi*(4*k-3)/8 + 0.75j*ctx.ln2
            s = ctx.expjpi(ctx.mpf(1)/3) * T(t)
            return ctx.findroot(lambda z: ctx.airybi(z,1), s)
        t = 3*ctx.pi*(4*k-1)/8 + 0.75j*ctx.ln2
        s = ctx.expjpi(ctx.mpf(1)/3) * U(t)
        return ctx.findroot(ctx.airybi, s)

@defun
def airyaizero(ctx, k, derivative=0):
    return _airy_zero(ctx, 0, k, derivative, False)

@defun
def airybizero(ctx, k, derivative=0, complex=False):
    return _airy_zero(ctx, 1, k, derivative, complex)


@defun_wrapped
def coulombc(ctx, l, eta, _cache={}):
    if (l, eta) in _cache and _cache[l,eta][0] >= ctx.prec:
        return +_cache[l,eta][1]
    G3 = ctx.loggamma(2*l+2)
    G1 = ctx.loggamma(1+l+ctx.j*eta)
    G2 = ctx.loggamma(1+l-ctx.j*eta)
    v = 2**l * ctx.exp((-ctx.pi*eta+G1+G2)/2 - G3)
    if not (ctx.im(l) or ctx.im(eta)):
        v = ctx.re(v)
    _cache[l,eta] = (ctx.prec, v)
    return v

@defun_wrapped
def coulombf(ctx, l, eta, z, w=1, chop=True, **kwargs):
    # Regular Coulomb wave function
    # Note: w can be either 1 or -1; the other may be better in some cases
    # TODO: check that chop=True chops when and only when it should
    #ctx.prec += 10
    def h(l, eta):
        try:
            jw = ctx.j*w
            jwz = ctx.fmul(jw, z, exact=True)
            jwz2 = ctx.fmul(jwz, -2, exact=True)
            C = ctx.coulombc(l, eta)
            T1 = [C, z, ctx.exp(jwz)], [1, l+1, 1], [], [], [1+l+jw*eta], \
                [2*l+2], jwz2
        except ValueError:
            T1 = [0], [-1], [], [], [], [], 0
        return (T1,)
    v = ctx.hypercomb(h, [l,eta], **kwargs)
    if chop and (not ctx.im(l)) and (not ctx.im(eta)) and (not ctx.im(z)) and \
        (ctx.re(z) >= 0):
        v = ctx.re(v)
    return v

@defun_wrapped
def _coulomb_chi(ctx, l, eta, _cache={}):
    if (l, eta) in _cache and _cache[l,eta][0] >= ctx.prec:
        return _cache[l,eta][1]
    def terms():
        l2 = -l-1
        jeta = ctx.j*eta
        return [ctx.loggamma(1+l+jeta) * (-0.5j),
            ctx.loggamma(1+l-jeta) * (0.5j),
            ctx.loggamma(1+l2+jeta) * (0.5j),
            ctx.loggamma(1+l2-jeta) * (-0.5j),
            -(l+0.5)*ctx.pi]
    v = ctx.sum_accurately(terms, 1)
    _cache[l,eta] = (ctx.prec, v)
    return v

@defun_wrapped
def coulombg(ctx, l, eta, z, w=1, chop=True, **kwargs):
    # Irregular Coulomb wave function
    # Note: w can be either 1 or -1; the other may be better in some cases
    # TODO: check that chop=True chops when and only when it should
    if not ctx._im(l):
        l = ctx._re(l)  # XXX: for isint
    def h(l, eta):
        # Force perturbation for integers and half-integers
        if ctx.isint(l*2):
            T1 = [0], [-1], [], [], [], [], 0
            return (T1,)
        l2 = -l-1
        try:
            chi = ctx._coulomb_chi(l, eta)
            jw = ctx.j*w
            s = ctx.sin(chi); c = ctx.cos(chi)
            C1 = ctx.coulombc(l,eta)
            C2 = ctx.coulombc(l2,eta)
            u = ctx.exp(jw*z)
            x = -2*jw*z
            T1 = [s, C1, z, u, c], [-1, 1, l+1, 1, 1], [], [], \
                [1+l+jw*eta], [2*l+2], x
            T2 = [-s, C2, z, u],   [-1, 1, l2+1, 1],    [], [], \
                [1+l2+jw*eta], [2*l2+2], x
            return T1, T2
        except ValueError:
            T1 = [0], [-1], [], [], [], [], 0
            return (T1,)
    v = ctx.hypercomb(h, [l,eta], **kwargs)
    if chop and (not ctx._im(l)) and (not ctx._im(eta)) and (not ctx._im(z)) and \
        (ctx._re(z) >= 0):
        v = ctx._re(v)
    return v

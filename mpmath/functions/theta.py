from mpmath.libmp.libintmath import jacobi_symbol

from .functions import defun, defun_wrapped


@defun
def _djacobi_theta2(ctx, z, q, nd):
    # the loops below break when the fixed precision quantities
    # a and b go to zero;
    # right shifting small negative numbers by wp one obtains -1, not zero,
    # so the condition a**2 + b**2 > MIN is used to break the loops.
    MIN = 2
    extra1 = 10
    extra2 = 20
    if not ctx._im(q) and not ctx._im(z):
        wp = ctx.prec + extra1
        x = ctx.to_fixed(ctx._re(q), wp)
        x2 = (x*x) >> wp
        a = b = x2
        c1, s1 = ctx.cos_sin(ctx._re(z), prec=wp)
        cn = c1 = ctx.to_fixed(c1, wp)
        sn = s1 = ctx.to_fixed(s1, wp)
        c2 = (c1*c1 - s1*s1) >> wp
        s2 = (c1 * s1) >> (wp - 1)
        cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        if nd&1:
            s = s1 + ((a * sn * 3**nd) >> wp)
        else:
            s = c1 + ((a * cn * 3**nd) >> wp)
        n = 2
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            if nd&1:
                s += (a * sn * (2*n+1)**nd) >> wp
            else:
                s += (a * cn * (2*n+1)**nd) >> wp
            n += 1
        s = -(s << 1)
        s = ctx.ldexp(s, -wp)
    # case z real, q complex
    elif not ctx._im(z):
        wp = ctx.prec + extra2
        xre = ctx.to_fixed(ctx._re(q), wp)
        xim = ctx.to_fixed(ctx._im(q), wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = x2re
        aim = bim = x2im
        c1, s1 = ctx.cos_sin(ctx._re(z), prec=wp)
        cn = c1 = ctx.to_fixed(c1, wp)
        sn = s1 = ctx.to_fixed(s1, wp)
        c2 = (c1*c1 - s1*s1) >> wp
        s2 = (c1 * s1) >> (wp - 1)
        cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        if nd&1:
            sre = s1 + ((are * sn * 3**nd) >> wp)
            sim = ((aim * sn * 3**nd) >> wp)
        else:
            sre = c1 + ((are * cn * 3**nd) >> wp)
            sim = ((aim * cn * 3**nd) >> wp)
        n = 5
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp

            if nd&1:
                sre += ((are * sn * n**nd) >> wp)
                sim += ((aim * sn * n**nd) >> wp)
            else:
                sre += ((are * cn * n**nd) >> wp)
                sim += ((aim * cn * n**nd) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)
    # case z complex, q real
    elif not ctx._im(q):
        wp = ctx.prec + extra2
        x = ctx.to_fixed(ctx._re(q), wp)
        x2 = (x*x) >> wp
        a = b = x2
        c1, s1 = ctx.cos_sin(z, prec=wp)
        cnre = c1re = ctx.to_fixed(ctx._re(c1), wp)
        cnim = c1im = ctx.to_fixed(ctx._im(c1), wp)
        snre = s1re = ctx.to_fixed(ctx._re(s1), wp)
        snim = s1im = ctx.to_fixed(ctx._im(s1), wp)
        c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
        c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
        s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
        s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
        t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
        t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
        t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
        t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
        cnre = t1
        cnim = t2
        snre = t3
        snim = t4
        if nd&1:
            sre = s1re + ((a * snre * 3**nd) >> wp)
            sim = s1im + ((a * snim * 3**nd) >> wp)
        else:
            sre = c1re + ((a * cnre * 3**nd) >> wp)
            sim = c1im + ((a * cnim * 3**nd) >> wp)
        n = 5
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if nd&1:
                sre += ((a * snre * n**nd) >> wp)
                sim += ((a * snim * n**nd) >> wp)
            else:
                sre += ((a * cnre * n**nd) >> wp)
                sim += ((a * cnim * n**nd) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)
    # case z and q complex
    else:
        wp = ctx.prec + extra2
        xre = ctx.to_fixed(ctx._re(q), wp)
        xim = ctx.to_fixed(ctx._im(q), wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = x2re
        aim = bim = x2im
        c1, s1 = ctx.cos_sin(z, prec=wp)
        cnre = c1re = ctx.to_fixed(ctx._re(c1), wp)
        cnim = c1im = ctx.to_fixed(ctx._im(c1), wp)
        snre = s1re = ctx.to_fixed(ctx._re(s1), wp)
        snim = s1im = ctx.to_fixed(ctx._im(s1), wp)
        c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
        c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
        s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
        s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
        t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
        t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
        t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
        t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
        cnre = t1
        cnim = t2
        snre = t3
        snim = t4
        if nd&1:
            sre = s1re + (((are * snre - aim * snim) * 3**nd) >> wp)
            sim = s1im + (((are * snim + aim * snre)* 3**nd) >> wp)
        else:
            sre = c1re + (((are * cnre - aim * cnim) * 3**nd) >> wp)
            sim = c1im + (((are * cnim + aim * cnre)* 3**nd) >> wp)
        n = 5
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if nd&1:
                sre += (((are * snre - aim * snim) * n**nd) >> wp)
                sim += (((aim * snre + are * snim) * n**nd) >> wp)
            else:
                sre += (((are * cnre - aim * cnim) * n**nd) >> wp)
                sim += (((aim * cnre + are * cnim) * n**nd) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)
    s *= ctx.nthroot(q, 4)
    return (-1)**(1 - (nd&1) + nd//2) * s

@defun
def _djacobi_theta3(ctx, z, q, nd):
    MIN = 2
    extra1 = 10
    extra2 = 20
    if not ctx._im(q) and not ctx._im(z):
        s = 0
        wp = ctx.prec + extra1
        x = ctx.to_fixed(ctx._re(q), wp)
        a = (1 << wp)
        b = x
        x2 = (x*x) >> wp
        c1, s1 = ctx.cos_sin(ctx._re(z)*2, prec=wp)
        c1 = ctx.to_fixed(c1, wp)
        s1 = ctx.to_fixed(s1, wp)
        cn = c1
        sn = s1
        if nd&1:
            s += (a * sn) >> wp
        else:
            s += (a * cn) >> wp
        n = 2
        while True:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            if abs(a) <= MIN:
                break
            cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            if nd&1:
                s += (a * sn * n**nd) >> wp
            else:
                s += (a * cn * n**nd) >> wp
            n += 1
        s = -(s << (nd+1))
        s = ctx.ldexp(s, -wp)*q
    # case z real, q complex
    elif not ctx._im(z):
        wp = ctx.prec + extra2
        xre = ctx.to_fixed(ctx._re(q), wp)
        xim = ctx.to_fixed(ctx._im(q), wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = (1 << wp)
        aim = 0
        bre = xre
        bim = xim
        c1, s1 = ctx.cos_sin(ctx._re(z)*2, prec=wp)
        c1 = ctx.to_fixed(c1, wp)
        s1 = ctx.to_fixed(s1, wp)
        cn = c1
        sn = s1
        if nd&1:
            sre = (are * sn) >> wp
            sim = (aim * sn) >> wp
        else:
            sre = (are * cn) >> wp
            sim = (aim * cn) >> wp
        n = 2
        while True:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            if are**2 + aim**2 <= MIN:
                break
            cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            if nd&1:
                sre += (are * sn * n**nd) >> wp
                sim += (aim * sn * n**nd) >> wp
            else:
                sre += (are * cn * n**nd) >> wp
                sim += (aim * cn * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)*q
    # case z complex, q real
    elif not ctx._im(q):
        wp = ctx.prec + extra2
        x = ctx.to_fixed(ctx._re(q), wp)
        a = (1 << wp)
        b = x
        x2 = (x*x) >> wp
        c1, s1 = ctx.cos_sin(2*z, prec=wp)
        cnre = c1re = ctx.to_fixed(ctx._re(c1), wp)
        cnim = c1im = ctx.to_fixed(ctx._im(c1), wp)
        snre = s1re = ctx.to_fixed(ctx._re(s1), wp)
        snim = s1im = ctx.to_fixed(ctx._im(s1), wp)
        if nd&1:
            sre = (a * snre) >> wp
            sim = (a * snim) >> wp
        else:
            sre = (a * cnre) >> wp
            sim = (a * cnim) >> wp
        n = 2
        while True:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            if abs(a) <= MIN:
                break
            t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
            t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
            t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
            t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if nd&1:
                sre += (a * snre * n**nd) >> wp
                sim += (a * snim * n**nd) >> wp
            else:
                sre += (a * cnre * n**nd) >> wp
                sim += (a * cnim * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)*q
    # case z and q complex
    else:
        wp = ctx.prec + extra2
        xre = ctx.to_fixed(ctx._re(q), wp)
        xim = ctx.to_fixed(ctx._im(q), wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = (1 << wp)
        aim = 0
        bre = xre
        bim = xim
        c1, s1 = ctx.cos_sin(2*z, prec=wp)
        cnre = c1re = ctx.to_fixed(ctx._re(c1), wp)
        cnim = c1im = ctx.to_fixed(ctx._im(c1), wp)
        snre = s1re = ctx.to_fixed(ctx._re(s1), wp)
        snim = s1im = ctx.to_fixed(ctx._im(s1), wp)
        if nd&1:
            sre = (are * snre - aim * snim) >> wp
            sim = (aim * snre + are * snim) >> wp
        else:
            sre = (are * cnre - aim * cnim) >> wp
            sim = (aim * cnre + are * cnim) >> wp
        n = 2
        while True:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            if are**2 + aim**2 <= MIN:
                break
            t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
            t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
            t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
            t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if nd&1:
                sre += ((are * snre - aim * snim) * n**nd) >> wp
                sim += ((aim * snre + are * snim) * n**nd) >> wp
            else:
                sre += ((are * cnre - aim * cnim) * n**nd) >> wp
                sim += ((aim * cnre + are * cnim) * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)*q
    if nd&1:
        return (-1)**(nd//2) * s
    else:
        return (-1)**(1 + nd//2) * s + (ctx.zero if nd else ctx.one)

@defun
def _reduce_psl2z(ctx, z):
    """
    Returns the cumulative transformation matrix, that reduces a complex
    number z to the fundamental domain of PSL(2, Z), chosen to be
    |Re(z)| ≤ 0.5 and |z| ≥ 1.
    """
    z = ctx.convert(z)
    assert z.imag > 0, f"Expected point from upper half-plane, got {ctx.mpc(z)}"

    a = d = 1
    b = c = 0

    z_orig = z
    with ctx.extraprec(30):
        while True:
            # Translate to center in |Re(z)| ≤ 1/2
            n = round(z.real)
            if n:
                z -= n
                a -= n*c
                b -= n*d

            # Maybe apply an inversion
            if z.real**2 + z.imag**2 < 1:
                z = -1/z
                a, c = -c, a
                b, d = -d, b
                if abs(z.real) <= 0.5:
                    break
            else:
                break

        # Canonicalize matrix
        if c < 0 or (c == 0 and d < 0):
            a, b, c, d = -a, -b, -c, -d

    return a, b, c, d

#
# General modular transformations for jtheta()
#
# References:
# * Hans Rademacher (1973), "Topics in Analytic Number Theory",
#   Springer. Section 81.
# * [DLMF]_, §20.7(viii).
#

_T_map = {(0, 0): 1, (0, 1): 2, (1, 0): 4, (1, 1): 3}

def _jtheta_permutation(n, a, b, c, d):
    if n == 2:
        return _T_map[(c%2, d%2)]
    if n == 3:
        return _T_map[((a + c)%2, (b + d)%2)]
    if n == 4:
        return _T_map[(a%2, b%2)]
    return 1

@defun
def _jtheta_eps(ctx, n, a, b, c, d):
    if n != 1:
        if n == 2:
            phi = (c - 2)*d - 2 + 2*(1 - c)*((d + 1)%2)
        elif n == 3:
            phi = (a + c - 2)*(b + d) - 3 + 2*(1 - a - c)*((b + d + 1)%2)
        else:
            phi = (a - 2)*b - 4 + 2*(1 - a)*((b + 1)%2)
        k = ctx._jtheta_eps(1, -d, b, c, -a)
    else:
        if c % 2 == 0:
            phi = d*(b - c - 1) + 2
            k = jacobi_symbol(c, d)
        else:
            phi = c*(a + d + 1) - 3
            k = jacobi_symbol(d, c)
    return ctx.expjpi(ctx.convert(phi)/4)/k

@defun
def _jtheta_needs_modular(ctx, z, q):
    if not z.imag:
        return False
    tau = ctx.taufrom(q=q)
    assert abs(q) < 1 and tau.imag > 0
    return abs(tau.real) > 0.5 or tau.real**2 + tau.imag**2 < 1

@defun
def _jtheta_modular(ctx, g, n, z, q, nd):
    a, b, c, d = g
    tau = ctx.taufrom(q=q)
    v = -1/(c*tau + d)
    alpha = 1j*v*c/ctx.pi

    assert abs(q) < 1 and tau.imag > 0

    new_n = _jtheta_permutation(n, -d, b, c, -a)
    new_z = z*v
    new_tau = (a*tau + b)/(c*tau + d)
    new_q = ctx.qfrom(tau=new_tau)

    assert abs(new_tau.real) <= 0.5 and new_tau.real**2 + new_tau.imag**2 >= 1

    def terms():
        Him1, Hi = ctx.zero, ctx.one
        a2 = alpha*2
        a2z = a2*z
        for i in range(nd + 1):
            yield (ctx.binomial(nd, i) * Hi * v**(nd - i)
                   * ctx.jtheta(new_n, new_z, new_q, nd - i))
            Him1, Hi = Hi, a2z*Hi + a2*i*Him1

    C = ctx._jtheta_eps(n, -d, b, c, -a)*ctx.sqrt(v/1j)
    X = alpha*z**2
    return C*ctx.exp(X)*sum(terms())

@defun
def jtheta(ctx, n, z, q, derivative=0):
    n = int(n)
    z = ctx.convert(z)
    q = ctx.convert(q)
    nd = int(derivative)

    if n not in range(1, 5):
        raise ValueError("First argument expected to be 1, 2, 3 or 4")
    if abs(q) >= 1:
        raise ValueError(f"abs(q) >= 1")

    # We use Fourier series (DLMF, §20.2(i)) to compute functions, when
    # |q| is not near 1.  Else, transform τ to the fundamental
    # domain (|Re(τ)| ≤ 0.5 and |τ| ≥ 1), applying transformations
    # of lattice parameter (DLMF, §20.7(viii)).

    if ctx._jtheta_needs_modular(z, q):
        tau = ctx.taufrom(q=q)
        g = ctx._reduce_psl2z(tau)

        # Estimate exponential factor
        c, d = g[2:]
        extra = 10*(nd + 1) + max(0, ctx.mag(c/(c*tau + d)*z**2))

        return ctx.extraprec(extra, True)(ctx._jtheta_modular)(g, n, z, q, nd)

    # At that point, τ is in the fundamental domain and thus Im(τ) ≥ √3π/2.
    # Using quasi-periodicity property (see DLMF, §20.2(ii)) brings
    # z to the domain |Im(z)| ≤ π |Im(τ)|/2.

    if abs(z.imag) > abs(ctx.log(q).real)/2:
        with ctx.extraprec(10):
            tau = ctx.taufrom(q=q)
            tau_pi = tau*ctx.pi
            k = round(z.imag/tau_pi.imag)
            assert k != 0
            beta = -ctx.j*2*k
            C = q**(k**2)*ctx.exp(beta*z)
            if n in (1, 4) and k & 1:
                C = -C
            new_z = z - k*tau_pi

            def terms():
                for i in range(nd + 1):
                    yield (ctx.binomial(nd, i) * beta**i
                           * ctx.jtheta(n, new_z, q, nd - i))

            res = C*sum(terms())
        return +res

    extra = 10 + ctx.prec * nd // 10
    if z:
        M = ctx.mag(z)
        if M > 5 or ((n != 1 if nd else n == 1) and M < -5):
            extra += 2*abs(M)
    with ctx.extraprec(extra):
        if n < 3:
            z_inner = z - ctx.pi/2 if n == 1 else z
            res = ctx._djacobi_theta2(z_inner, q, nd)
        else:
            q_inner = -q if n == 4 else q
            res = ctx._djacobi_theta3(z, q_inner, nd)
    return +res

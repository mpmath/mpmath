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
def _djacobi_theta2a(ctx, z, q, nd):
    """
    case ctx._im(z) != 0
    dtheta(2, z, q, nd) =
    j*nd q**1/4 * Sum(q**(n*n + n) * (2*n+1)*nd * exp(j*(2*n + 1)*z), n=-inf, inf)
    max term for (2*n0+1)*log(q).real - 2* ctx._im(z) ~= 0
    n0 = int(ctx._im(z)/log(q).real - 1/2)
    """
    n = n0 = int(z.imag/ctx.log(q).real - 1/2)
    e2 = ctx.expj(2*z)
    e = e0 = ctx.expj((2*n + 1)*z)
    a = q**(n*n + n)
    # leading term
    term = (2*n+1)**nd * a * e
    s = term
    eps1 = ctx.eps*abs(term)
    while 1:
        n += 1
        e = e * e2
        term = (2*n+1)**nd * q**(n*n + n) * e
        if abs(term) < eps1:
            break
        s += term
    e = e0
    e2 = ctx.expj(-2*z)
    n = n0
    while 1:
        n -= 1
        e = e * e2
        term = (2*n+1)**nd * q**(n*n + n) * e
        if abs(term) < eps1:
            break
        s += term
    return ctx.j**nd * s * ctx.nthroot(q, 4)

@defun
def _djacobi_theta3a(ctx, z, q, nd):
    """
    case ctx._im(z) != 0
    djtheta3(z, q, nd) = (2*j)**nd *
      Sum(q**(n*n) * n**nd * exp(j*2*n*z), n, -inf, inf)
    max term for minimum n*abs(log(q).real) + ctx._im(z)
    """
    n = n0 = int(-z.imag/abs(ctx.log(q).real))
    e2 = ctx.expj(2*z)
    e = e0 = ctx.expj(2*n*z)
    a = q**(n*n) * e
    s = term = n**nd * a
    eps1 = ctx.eps*abs(term if term else a)
    while 1:
        n += 1
        e = e * e2
        a = q**(n*n) * e
        term = n**nd * a
        aterm = abs(term if term else a)
        if aterm < eps1:
            break
        s += term
    e = e0
    e2 = ctx.expj(-2*z)
    n = n0
    while 1:
        n -= 1
        e = e * e2
        a = q**(n*n) * e
        term = n**nd * a
        aterm = abs(term if term else a)
        if aterm < eps1:
            break
        s += term
    return (2*ctx.j)**nd * s


def _theta_modular_reduce(ctx, n, z, q):
    """Apply PSL(2,Z) reduction so that tau = log(q)/(i*pi) lies in the
    fundamental domain |Re tau| <= 1/2, |tau| >= 1, and reduce z modulo the
    quasi-period pi*tau, so that |Im(z_inner)| <= Im(pi*tau)/2.

    Returns (n_new, z_inner, q_new, C, alpha, beta, T) such that

        theta_n(z, q) = C * exp(alpha*z^2 + beta*z)
                        * theta_{n_new}(z_inner, q_new)

    where z_inner is an affine function of z with dz_inner/dz = 1/T.
    """
    pi = ctx.pi
    j = ctx.j

    n_cur = int(n)
    tau = ctx.log(q) / (j * pi)

    # Constant scalar accumulator
    C = ctx.one
    # coefficient of z^2 in the prefix exponent
    alpha = ctx.zero
    # coefficient of z in the prefix exponent
    beta = ctx.zero
    # cumulative scale from the inversions; z_inner = z/T
    T = ctx.one

    # Classical reduction to the fundamental domain.  It terminates: each
    # inversion strictly increases Im(tau), and the SL(2,Z) orbit of tau
    # has only finitely many points with Im above any given bound.
    while True:
        # Step 1: bring Re(tau) into [-1/2, 1/2] (DLMF 20.7.26-20.7.29).
        # theta_1, theta_2 gain e^(i*pi/4) per unit shift, theta_3 <-> theta_4.
        re_tau = ctx._re(tau)
        k = ctx.nint(re_tau)
        if k:
            tau -= k
            re_tau = ctx._re(tau)
            if n_cur in (3, 4) and (int(k) & 1):
                n_cur = 7 - n_cur
            if n_cur in (1, 2):
                C *= ctx.expjpi(k / 4)
        # Step 2: if |tau| < 1, invert via tau' = -1/tau (DLMF 20.7.30-20.7.33):
        #   theta_n(z|tau) = (-i*tau)^(-1/2) * exp(-i*z^2/(pi*tau))
        #                    * theta_m(z*tau' | tau'),   z*tau' = -z/tau,
        # index map 1->1 (extra factor -i), 2<->4, 3->3.  So C gains
        # (-i*tau)^(-1/2) (and -i for theta_1), alpha gains -i/(pi*tau) (scaled
        # by 1/T^2), and T accumulates -tau (z_inner = z/T = -z/tau = z*tau').
        im_tau = ctx._im(tau)
        if re_tau * re_tau + im_tau * im_tau < 1:
            C /= ctx.sqrt(-j * tau)
            if n_cur == 1:
                C *= -j
            elif n_cur in (2, 4):
                n_cur = 6 - n_cur
            alpha -= j / (pi * tau * T * T)
            T *= -tau
            tau = -1 / tau
        else:
            break

    q_new = ctx.expjpi(tau)

    # Step 3: subtract one quasi-period pi*tau in z (DLMF 20.2.6-20.2.9).
    # In the fundamental domain Im(pi*tau) >= pi*sqrt(3)/2, so one step
    # brings |Im(z_inner)| within Im(pi*tau)/2.  The e^(-2i*k_z*z_inner)
    # factor (z_inner = z/T) adds -2i*k_z/T to beta.
    pi_tau = pi * tau
    z_inner = z / T
    k_z = ctx.nint(ctx._im(z_inner) / ctx._im(pi_tau))
    if k_z:
        if n_cur in (1, 4) and (int(k_z) & 1):
            C = -C
        C *= q_new ** (k_z * k_z)
        beta -= 2 * j * k_z / T
        z_inner -= k_z * pi_tau

    return n_cur, z_inner, q_new, C, alpha, beta, T


def _theta_need_modular(ctx, z, q):
    # jtheta by direct series loses precision when tau = log(q)/(i*pi) lies
    # outside the PSL(2,Z) fundamental domain (|Re tau| <= 1/2 and
    # |tau| >= 1).  Only route through modular reduction in that regime;
    # otherwise fall through to the standard dispatcher.
    if not ctx._im(z):
        return False
    tau = ctx.log(q) / (ctx.j * ctx.pi)
    re = ctx._re(tau)
    im = ctx._im(tau)
    if abs(re) > ctx.mpf(0.5):
        return True
    if re*re + im*im < 1:
        return True
    return False


def _theta_inner_a(ctx, n, z, q, nd=0):
    """Dispatch to the shifted-series helpers, used inside the modular
    path after reduction (Im(z) may be non-negligible)."""
    if n == 1:
        return ctx._djacobi_theta2a(z - ctx.pi/2, q, nd)
    if n == 2:
        return ctx._djacobi_theta2a(z, q, nd)
    if n == 3:
        return ctx._djacobi_theta3a(z, q, nd)
    return ctx._djacobi_theta3a(z, -q, nd)


@defun
def jtheta(ctx, n, z, q, derivative=0):
    n = int(n)
    z = ctx.convert(z)
    q = ctx.convert(q)
    nd = int(derivative)

    if n not in range(1, 5):
        raise ValueError("First argument expected to be 1, 2, 3 or 4")
    if abs(q) > ctx.THETA_Q_LIM:
        raise ValueError(f"abs(q) > THETA_Q_LIM = {ctx.THETA_Q_LIM}")

    # Implementation note
    # If ctx._im(z) is close to zero, _jacobi_theta2 and _jacobi_theta3
    # are used,
    # which compute the series starting from n=0 using fixed precision
    # numbers;
    # otherwise  _jacobi_theta2a and _jacobi_theta3a are used, which compute
    # the series starting from n=n0, which is the largest term.

    # TODO: write _jacobi_theta2a and _jacobi_theta3a using fixed-point

    extra = 10 + ctx.prec * nd // 10
    if z:
        M = ctx.mag(z)
        if M > 5 or ((n != 1 if nd else n == 1) and M < -5):
            extra += 2*abs(M)
    cz = 0.5
    extra2 = 50
    with ctx.extraprec(extra):
        if _theta_need_modular(ctx, z, q):
            # tau outside the PSL(2,Z) fundamental domain: reduce to the
            # domain (where the nome series converges fast and accurately),
            # then undo the transform.  For nd > 0 differentiate by Leibniz
            # over prefix(z) = exp(alpha*z^2 + beta*z) and theta_inner(z/T):
            # prefix^(i) = prefix * H_i, with the Hermite-type recurrence
            # H_{i+1} = (2*alpha*z + beta)*H_i + 2*alpha*i*H_{i-1} (H_0 = 1).
            # nd == 0 collapses the sum to the single bare term.
            ctx.prec += 30 + 10 * nd
            n_new, z_inner, q_new, C, alpha, beta, T = \
                _theta_modular_reduce(ctx, n, z, q)
            # exp(X) loses ~mag(X) bits when |q| -> 1; add that many and redo.
            X = alpha * z * z + beta * z
            mag_X = ctx.mag(X)
            if mag_X > 0:
                ctx.prec += mag_X + 10
                n_new, z_inner, q_new, C, alpha, beta, T = \
                    _theta_modular_reduce(ctx, n, z, q)
                X = alpha * z * z + beta * z
            prefix = C * ctx.exp(X)
            gp = 2 * alpha * z + beta
            T_inv = 1 / T
            def terms():
                Hm1, Hi = ctx.zero, ctx.one     # H_{-1} (unused), H_0
                for i in range(nd + 1):
                    yield (ctx.binomial(nd, i) * Hi * T_inv**(nd - i)
                           * _theta_inner_a(ctx, n_new, z_inner, q_new, nd - i))
                    Hm1, Hi = Hi, gp * Hi + 2 * alpha * i * Hm1
            res = prefix * ctx.fsum(terms())
        elif n in [1, 2]:
            z_inner = z - ctx.pi/2 if n == 1 else z
            if z.imag:
                if abs(z.imag) < cz * abs(ctx.log(q).real):
                    ctx.dps += extra2
                    res = ctx._djacobi_theta2(z_inner, q, nd)
                else:
                    ctx.dps += 10
                    res = ctx._djacobi_theta2a(z_inner, q, nd)
            else:
                res = ctx._djacobi_theta2(z_inner, q, nd)
        else:
            q_inner = -q if n == 4 else q
            if z.imag:
                if abs(z.imag) < cz * abs(ctx.log(q).real):
                    ctx.dps += extra2
                    res = ctx._djacobi_theta3(z, q_inner, nd)
                else:
                    ctx.dps += 10
                    res = ctx._djacobi_theta3a(z, q_inner, nd)
            else:
                res = ctx._djacobi_theta3(z, q_inner, nd)
    return +res

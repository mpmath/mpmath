#!/usr/bin/env python
"""
    elliptic.py

    Implements the Jacobi theta and Jacobi elliptic functions, using
    arbitrary precision math library

    Author of the first version: M.T. Taschuk

    References:

    [1] Abramowitz & Stegun. 'Handbook of Mathematical Functions, 9th Ed.',
        (Dover duplicate of 1972 edition)
    [2] Whittaker 'A Course of Modern Analysis, 4th Ed.', 1946,
        Cambridge Univeristy Press

"""
import sys

from mptypes import (mpf, mpc, mp, convert_lossless, eps, one, zero)
from functions import (pi, sqrt, cos, sin, exp, tanh, ellipk, sech, nthroot)
from libmpf import to_fixed, MP_ZERO, mpf_shift, from_man_exp
from mpmath.libelefun import cos_sin

def calculate_nome(k):
    """
    Calculate the nome, q, from the value for k.

    Useful factoids:

    k**2 = m;   m is used in Abramowitz
    """
    k = convert_lossless(k)

    if abs(k) > one:             # range error
        raise ValueError

    if k == zero:
        return zero
    elif k == one:
        return one
    else:
        kprimesquared = one - k**2
        kprime = sqrt(kprimesquared)
        top = ellipk(kprimesquared)
        bottom = ellipk(k**2)

        argument = -pi*top/bottom

        nome = exp(argument)
        return nome

def calculate_k(q):
    """
    Calculates the value of k for a particular nome, q,
    using jacobi theta functions.
    """
    #zero = mpf('0')
    #one = mpf('1')

    q = convert_lossless(q)
    #if q > one or q < zero:
    #    raise ValueError

    v2 = jacobi_theta(2, 0, q)
    v3 = jacobi_theta(3, 0, q)
    m = v2**2/v3**2
    return m

def _jacobi_theta2(z, q):
    if abs(q) >= 1:
        raise ValueError
    extra1 = 10
    extra2 = 20
    if z == zero:
        if isinstance(q, mpf):
            wp = mp.prec + extra1
            x = to_fixed(q._mpf_, wp)
            x2 = (x*x) >> wp
            a = b = x2
            s = x2
            while(a):
                b = (b*x2) >> wp
                a = (a*b) >> wp
                s += a
            s = (1 << (wp+1)) + (s << 1)
            s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
        else:
            wp = mp.prec + extra1
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = x2re
            aim = bim = x2im
            sre = (1<<wp) + are
            sim = aim
            while (are or aim):
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                sre += are
                sim += aim
            sre = (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
    else:
        if isinstance(q, mpf) and isinstance(z, mpf):
            wp = mp.prec + extra1
            x = to_fixed(q._mpf_, wp)
            x2 = (x*x) >> wp
            a = b = x2
            c1, s1 = cos_sin(z._mpf_, wp)
            cn = c1 = to_fixed(c1, wp)
            sn = s1 = to_fixed(s1, wp)
            c2 = (c1*c1 - s1*s1) >> wp
            s2 = (c1 * s1) >> (wp - 1)
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            s = c1 + ((a * cn) >> wp)
            while a:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
                s += (a * cn) >> wp
            s = (s << 1)
            s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
            s *= nthroot(q, 4)
            return s
        # case z real, q complex
        elif isinstance(z, mpf):
            wp = mp.prec + extra2
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = x2re
            aim = bim = x2im
            c1, s1 = cos_sin(z._mpf_, wp)
            cn = c1 = to_fixed(c1, wp)
            sn = s1 = to_fixed(s1, wp)
            c2 = (c1*c1 - s1*s1) >> wp
            s2 = (c1 * s1) >> (wp - 1)
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            sre = c1 + ((are * cn) >> wp)
            sim = ((aim * cn) >> wp)
            while (are or aim):
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp

                sre += ((are * cn) >> wp)
                sim += ((aim * cn) >> wp)
            sre = (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
        #case z complex, q real
        elif isinstance(q, mpf):
            wp = mp.prec + extra2
            x = to_fixed(q._mpf_, wp)
            x2 = (x*x) >> wp
            a = b = x2
            prec0 = mp.prec
            mp.prec = wp
            c1 = cos(z)
            s1 = sin(z)
            mp.prec = prec0
            cnre = c1re = to_fixed(c1.real._mpf_, wp)
            cnim = c1im = to_fixed(c1.imag._mpf_, wp)
            snre = s1re = to_fixed(s1.real._mpf_, wp)
            snim = s1im = to_fixed(s1.imag._mpf_, wp)
            #c2 = (c1*c1 - s1*s1) >> wp
            c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
            c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
            #s2 = (c1 * s1) >> (wp - 1)
            s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
            s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
            #cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            
            sre = c1re + ((a * cnre) >> wp)
            sim = c1im + ((a * cnim) >> wp)
            while (a):
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
                sre += ((a * cnre) >> wp)
                sim += ((a * cnim) >> wp)
            sre = (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
        # case z and q complex
        else:
            wp = mp.prec + extra2
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = x2re
            aim = bim = x2im
            prec0 = mp.prec
            mp.prec = wp
            # cos(2*z), siz(2*z) with z complex
            c1 = cos(z)
            s1 = sin(z)
            mp.prec = prec0
            cnre = c1re = to_fixed(c1.real._mpf_, wp)
            cnim = c1im = to_fixed(c1.imag._mpf_, wp)
            snre = s1re = to_fixed(s1.real._mpf_, wp)
            snim = s1im = to_fixed(s1.imag._mpf_, wp)
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
            sre = c1re + ((are * cnre - aim * cnim) >> wp)
            sim = c1im + ((are * cnim + aim * cnre) >> wp) 
            while (are or aim):
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                #cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
                t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
                t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
                t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
                t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                sre += ((are * cnre - aim * cnim) >> wp)
                sim += ((aim * cnre + are * cnim) >> wp)
            sre = (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
    s *= nthroot(q, 4)
    return s

def _djacobi_theta2(z, q):
    if abs(q) >= 1:
        raise ValueError
    extra1 = 10
    extra2 = 20
    if isinstance(q, mpf) and isinstance(z, mpf):
        wp = mp.prec + extra1
        x = to_fixed(q._mpf_, wp)
        x2 = (x*x) >> wp
        a = b = x2
        c1, s1 = cos_sin(z._mpf_, wp)
        cn = c1 = to_fixed(c1, wp)
        sn = s1 = to_fixed(s1, wp)
        c2 = (c1*c1 - s1*s1) >> wp
        s2 = (c1 * s1) >> (wp - 1)
        cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        s = s1 + ((a * sn * 3) >> wp)
        n = 5
        while a:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            s += (a * sn * n) >> wp
	    n += 2
        s = -(s << 1)
        s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
        s *= nthroot(q, 4)
        return s
        # case z real, q complex
    elif isinstance(z, mpf):
        wp = mp.prec + extra2
        xre, xim = q._mpc_
        xre = to_fixed(xre, wp)
        xim = to_fixed(xim, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = x2re
        aim = bim = x2im
        c1, s1 = cos_sin(z._mpf_, wp)
        cn = c1 = to_fixed(c1, wp)
        sn = s1 = to_fixed(s1, wp)
        c2 = (c1*c1 - s1*s1) >> wp
        s2 = (c1 * s1) >> (wp - 1)
        cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        sre = s1 + ((are * sn * 3) >> wp)
        sim = ((aim * sn * 3) >> wp)
        n = 5
        while (are or aim):
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp

            sre += ((are * sn * n) >> wp)
            sim += ((aim * sn * n) >> wp)
	    n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    #case z complex, q real
    elif isinstance(q, mpf):
        wp = mp.prec + extra2
        x = to_fixed(q._mpf_, wp)
        x2 = (x*x) >> wp
        a = b = x2
        prec0 = mp.prec
        mp.prec = wp
        c1 = cos(z)
        s1 = sin(z)
        mp.prec = prec0
        cnre = c1re = to_fixed(c1.real._mpf_, wp)
        cnim = c1im = to_fixed(c1.imag._mpf_, wp)
        snre = s1re = to_fixed(s1.real._mpf_, wp)
        snim = s1im = to_fixed(s1.imag._mpf_, wp)
        #c2 = (c1*c1 - s1*s1) >> wp
        c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
        c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
        #s2 = (c1 * s1) >> (wp - 1)
        s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
        s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
        #cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
        t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
        t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
        t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
        cnre = t1
        cnim = t2
        snre = t3
        snim = t4
        
        sre = s1re + ((a * snre * 3) >> wp)
        sim = s1im + ((a * snim * 3) >> wp)
        n = 5
        while (a):
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
            sre += ((a * snre * n) >> wp)
            sim += ((a * snim * n) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    # case z and q complex
    else:
        wp = mp.prec + extra2
        xre, xim = q._mpc_
        xre = to_fixed(xre, wp)
        xim = to_fixed(xim, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = x2re
        aim = bim = x2im
        prec0 = mp.prec
        mp.prec = wp
        # cos(2*z), siz(2*z) with z complex
        c1 = cos(z)
        s1 = sin(z)
        mp.prec = prec0
        cnre = c1re = to_fixed(c1.real._mpf_, wp)
        cnim = c1im = to_fixed(c1.imag._mpf_, wp)
        snre = s1re = to_fixed(s1.real._mpf_, wp)
        snim = s1im = to_fixed(s1.imag._mpf_, wp)
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
        sre = s1re + (((are * snre - aim * snim) * 3) >> wp)
        sim = s1im + (((are * snim + aim * snre)* 3) >> wp) 
        n = 5
        while (are or aim):
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            #cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            sre += (((are * snre - aim * snim) * n) >> wp)
            sim += (((aim * snre + are * snim) * n) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    s *= nthroot(q, 4)
    return s

def _jacobi_theta3(z, q):
    if abs(q) >= 1:
        raise ValueError
    extra1 = 10
    extra2 = 20
    if z == zero:
        if isinstance(q, mpf):
            wp = mp.prec + extra1
            x = to_fixed(q._mpf_, wp)
            s = x
            a = b = x
            x2 = (x*x) >> wp
            while(a):
                b = (b*x2) >> wp
                a = (a*b) >> wp
                s += a
            s = (1 << wp) + (s << 1)
            s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
            return s
        else:
            wp = mp.prec + extra1
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            sre = are = bre = xre
            sim = aim = bim = xim
            while (are or aim):
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                sre += are
                sim += aim
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
            return s
    else:
        if isinstance(q, mpf) and isinstance(z, mpf):
            s = MP_ZERO
            wp = mp.prec + extra1
            x = to_fixed(q._mpf_, wp)
            a = b = x
            x2 = (x*x) >> wp
            c1, s1 = cos_sin(mpf_shift(z._mpf_, 1), wp)
            c1 = to_fixed(c1, wp)
            s1 = to_fixed(s1, wp)
            cn = c1
            sn = s1
            s += (a * cn) >> wp
            while a:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
                s += (a * cn) >> wp
            s = (1 << wp) + (s << 1)
            s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
            return s
        # case z real, q complex
        elif isinstance(z, mpf):
            wp = mp.prec + extra2
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = xre
            aim = bim = xim
            c1, s1 = cos_sin(mpf_shift(z._mpf_, 1), wp)
            c1 = to_fixed(c1, wp)
            s1 = to_fixed(s1, wp)
            cn = c1
            sn = s1
            sre = (are * cn) >> wp
            sim = (aim * cn) >> wp
            while (are or aim):
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp

                sre += (are * cn) >> wp
                sim += (aim * cn) >> wp
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
            return s
        #case z complex, q real
        elif isinstance(q, mpf):
            wp = mp.prec + extra2
            x = to_fixed(q._mpf_, wp)
            a = b = x
            x2 = (x*x) >> wp
            prec0 = mp.prec
            mp.prec = wp
            c1 = cos(2*z)
            s1 = sin(2*z)
            mp.prec = prec0
            cnre = c1re = to_fixed(c1.real._mpf_, wp)
            cnim = c1im = to_fixed(c1.imag._mpf_, wp)
            snre = s1re = to_fixed(s1.real._mpf_, wp)
            snim = s1im = to_fixed(s1.imag._mpf_, wp)
            sre = (a * cnre) >> wp
            sim = (a * cnim) >> wp
            while (a):
                b = (b*x2) >> wp
                a = (a*b) >> wp
                t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
                t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
                t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
                t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                sre += (a * cnre) >> wp
                sim += (a * cnim) >> wp
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
            return s
        # case z and q complex
        else:
            wp = mp.prec + extra2
            xre, xim = q._mpc_
            xre = to_fixed(xre, wp)
            xim = to_fixed(xim, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = xre
            aim = bim = xim
            prec0 = mp.prec
            mp.prec = wp
            # cos(2*z), sin(2*z) with z complex
            c1 = cos(2*z)
            s1 = sin(2*z)
            mp.prec = prec0
            cnre = c1re = to_fixed(c1.real._mpf_, wp)
            cnim = c1im = to_fixed(c1.imag._mpf_, wp)
            snre = s1re = to_fixed(s1.real._mpf_, wp)
            snim = s1im = to_fixed(s1.imag._mpf_, wp)
            sre = (are * cnre - aim * cnim) >> wp
            sim = (aim * cnre + are * cnim) >> wp
            while (are or aim):
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
                t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
                t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
                t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                sre += (are * cnre - aim * cnim) >> wp
                sim += (aim * cnre + are * cnim) >> wp
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
            return s

def _djacobi_theta3(z, q):
    if abs(q) >= 1:
        raise ValueError
    extra1 = 10
    extra2 = 20
    if isinstance(q, mpf) and isinstance(z, mpf):
        s = MP_ZERO
        wp = mp.prec + extra1
        x = to_fixed(q._mpf_, wp)
        a = b = x
        x2 = (x*x) >> wp
        c1, s1 = cos_sin(mpf_shift(z._mpf_, 1), wp)
        c1 = to_fixed(c1, wp)
        s1 = to_fixed(s1, wp)
        cn = c1
        sn = s1
        s += (a * sn) >> wp
	n = 2
        while a:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            s += (a * sn * n) >> wp
	    n += 1
        s = -(s << 2)
        s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
        return s
    # case z real, q complex
    elif isinstance(z, mpf):
        wp = mp.prec + extra2
        xre, xim = q._mpc_
        xre = to_fixed(xre, wp)
        xim = to_fixed(xim, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = xre
        aim = bim = xim
        c1, s1 = cos_sin(mpf_shift(z._mpf_, 1), wp)
        c1 = to_fixed(c1, wp)
        s1 = to_fixed(s1, wp)
        cn = c1
        sn = s1
        sre = (are * sn) >> wp
        sim = (aim * sn) >> wp
	n = 2
        while (are or aim):
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp

            sre += (are * sn * n) >> wp
            sim += (aim * sn * n) >> wp
	    n += 1
        sre = -(sre << 2)
        sim = -(sim << 2)
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
        return s
    #case z complex, q real
    elif isinstance(q, mpf):
        wp = mp.prec + extra2
        x = to_fixed(q._mpf_, wp)
        a = b = x
        x2 = (x*x) >> wp
        prec0 = mp.prec
        mp.prec = wp
        c1 = cos(2*z)
        s1 = sin(2*z)
        mp.prec = prec0
        cnre = c1re = to_fixed(c1.real._mpf_, wp)
        cnim = c1im = to_fixed(c1.imag._mpf_, wp)
        snre = s1re = to_fixed(s1.real._mpf_, wp)
        snim = s1im = to_fixed(s1.imag._mpf_, wp)
        sre = (a * snre) >> wp
        sim = (a * snim) >> wp
	n = 2
        while (a):
            b = (b*x2) >> wp
            a = (a*b) >> wp
            t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
            t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
            t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
            t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            sre += (a * snre * n) >> wp
            sim += (a * snim * n) >> wp
	    n += 1
        sre = -(sre << 2)
        sim = -(sim << 2)
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
        return s
    # case z and q complex
    else:
        wp = mp.prec + extra2
        xre, xim = q._mpc_
        xre = to_fixed(xre, wp)
        xim = to_fixed(xim, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = xre
        aim = bim = xim
        prec0 = mp.prec
        mp.prec = wp
        # cos(2*z), sin(2*z) with z complex
        c1 = cos(2*z)
        s1 = sin(2*z)
        mp.prec = prec0
        cnre = c1re = to_fixed(c1.real._mpf_, wp)
        cnim = c1im = to_fixed(c1.imag._mpf_, wp)
        snre = s1re = to_fixed(s1.real._mpf_, wp)
        snim = s1im = to_fixed(s1.imag._mpf_, wp)
        sre = (are * snre - aim * snim) >> wp
        sim = (aim * snre + are * snim) >> wp
	n = 2
        while (are or aim):
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
            t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
            t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
            t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            sre += ((are * snre - aim * snim) * n) >> wp
            sim += ((aim * snre + are * snim) * n) >> wp
	    n += 1
        sre = -(sre << 2)
        sim = -(sim << 2)
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
        return s

def jtheta(n, z, q):
    """
    Jacobi theta functions as functions of the nome q
    n = 1,2,3,4 
    z complex number
    q complex number in the unit disk
    theta(1, z, q) = 
      2 * q**1/4 * Sum((-)**n * q**(n*n + n) * sin((2*n + 1)*z), n=0, inf)
    theta(2, z, q) = 
      2 * q**1/4 * Sum(q**(n*n + n) * cos((2*n + 1)*z), n=0, inf)
    theta(3, z, q) = 1 + 2 * Sum(q**(n**2) * cos(2*n*z), n=1, inf)
    theta(4, z, q) = 1 + 2 * Sum((-q)**(n**2) * cos(2*n*z), n=1, inf)
    """
    z = convert_lossless(z)
    q = convert_lossless(q)

    extra = 10
    prec0 = mp.prec
    try:
        mp.prec += extra
        if n == 1:
            res = _jacobi_theta2(z - pi/2, q)
        elif n == 2:
            res = _jacobi_theta2(z, q)
        elif n == 3:
            res = _jacobi_theta3(z, q)
        elif n == 4:
            res = _jacobi_theta3(z, -q)
        else:
            raise ValueError
    finally:
        mp.prec = prec0
    return res

def djtheta(n, z, q):
    """
    derivative of the Jacobi theta functions as functions of the nome q
    n = 1,2,3,4 
    z complex number
    q complex number in the unit disk
    djtheta(1, z, q) = 
      2 * q**1/4 * Sum((-)**n * q**(n*n + n) * 
      (2*n + 1) * cos((2*n + 1)*z), n=0, inf)
    djtheta(2, z, q) = 
      -2 * q**1/4 * Sum(q**(n*n + n) * 
      (2*n + 1) * cos((2*n + 1)*z), n=0, inf)
    djtheta(3, z, q) = -2 * Sum(q**(n**2) * 2 * n * sin(2*n*z), n=1, inf)
    djtheta(4, z, q) = -2 * Sum((-q)**(n**2) * 2 * n * sin(2*n*z), n=1, inf)
    """
    z = convert_lossless(z)
    q = convert_lossless(q)

    extra = 10
    prec0 = mp.prec
    try:
        mp.prec += extra
        if n == 1:
            res = _djacobi_theta2(z - pi/2, q)
        elif n == 2:
            res = _djacobi_theta2(z, q)
        elif n == 3:
            res = _djacobi_theta3(z, q)
        elif n == 4:
            res = _djacobi_theta3(z, -q)
        else:
            raise ValueError
    finally:
        mp.prec = prec0
    return res

def jacobi_theta_1(z, m):
    """
    The jacobi theta function 1 is defined by the series 
    expansion found in Abramowitz & Stegun [4]
    theta1(z, q) = 
    2 * q**1/4 * Sum((-)**n * q**(n*n + n) * sin((2*n + 1)*z), n=0, inf)
    The nome q is computed from the parameter m
    z is any complex number, q is a complex number in the unit circle
    """
    m = convert_lossless(m)
    z = convert_lossless(z)

    k = sqrt(m)
    q = calculate_nome(k)

    if abs(q) >= 1:
        raise ValueError

    res = _jacobi_theta2(z - pi/2, q)
    return res

def jacobi_theta_2(z, m):
    """
    The jacobi theta function 2 is defined by the series
    expansion found in Abramowitz & Stegun [4].
    theta2(z, q) = 
    2 * q**1/4 * Sum(q**(n*n + n) * cos((2*n + 1)*z), n=0, inf)
    The nome q is computed from the parameter m
    z is any complex number, q is a complex number in the unit circle
    """
    m = convert_lossless(m)
    z = convert_lossless(z)

    k = sqrt(m)
    q = calculate_nome(k)

    if abs(q) >= 1:
        raise ValueError

    return _jacobi_theta2(z, q)

def jacobi_theta_3(z, m):
    """
    The jacobi theta function 3 is defined by the series expansion
    found in Abramowitz & Stegun [4]
    theta3(z, q) = 1 + 2 * Sum(q**(n**2) * cos(2*n*z), n=1, inf)
    The nome q is computed from the parameter m
    z is any complex number, q is a complex number in the unit circle
    """
    m = convert_lossless(m)
    z = convert_lossless(z)

    k = sqrt(m)
    q = calculate_nome(k)

    if abs(q) >= 1:
        raise ValueError

    return _jacobi_theta3(z, q)

def jacobi_theta_4(z, m):
    """
    The jacobi theta function 4 is defined by the series expansion
    found in Abramowitz & Stegun [4]
    theta4(z, q) = 1 + 2 * Sum((-q)**(n**2) * cos(2*n*z), n=1, inf)
    The nome q is computed from the parameter m
    z is any complex number, q is a complex number in the unit circle
    """
    m = convert_lossless(m)
    z = convert_lossless(z)

    k = sqrt(m)
    q = calculate_nome(k)

    if abs(q) >= 1:
        raise ValueError

    res = _jacobi_theta3(z, -q)
    return res


def jsn(u, m):
    """
    Implementation of the jacobi elliptic sn function in term
    of jacoby theta functions.
    u is any complex number, m must be in the unit disk
    """
    if abs(m) < eps:
        return sin(u)
    elif m == one:
        return tanh(u)
    else:
        extra = 10
	try:
            mp.prec += extra
            q = calculate_nome(sqrt(m))
    
            v3 = jtheta(3, zero, q)
            v2 = jtheta(2, zero, q)        # mathworld says v4
            arg1 = u / (v3*v3)
            v1 = jtheta(1, arg1, q)
            v4 = jtheta(4, arg1, q)
    
            sn = (v3/v2)*(v1/v4)
	finally:
            mp.prec -= extra

        return sn


def jcn(u, m):
    """
    Implementation of the jacobi elliptic cn function in term
    of theta functions.
    u is any complex number, m must be in the unit disk
    """
    if abs(m) < eps:
        return cos(u)
    elif m == one:
        return sech(u)
    else:
        extra = 10
	try:
            mp.prec += extra
            q = calculate_nome(sqrt(m))
    
            v3 = jtheta(3, zero, q)
            v2 = jtheta(2, zero, q)
            v04 = jtheta(4, zero, q)
    
            arg1 = u / (v3*v3)
    
            v1 = jtheta(2, arg1, q)
            v4 = jtheta(4, arg1, q)
    
            cn = (v04/v2)*(v1/v4)
        finally:
	    mp.prec -= extra
        return cn


def jdn(u, m):
    """
    Implementation of the jacobi elliptic dn function in term
    of theta functions.
    u is any complex number, m must be in the unit disk
    """
    if m == zero:
        return one
    elif m == one:
        return sech(u)
    else:
        extra = 10
	try:
            mp.prec += extra
            q = calculate_nome(sqrt(m))
    
            v3 = jtheta(3, zero, q)
            v2 = jtheta(2, zero, q)
            v04 = jtheta(4, zero, q)
    
            arg1 = u / (v3*v3)
    
            v1 = jtheta(3, arg1, q)
            v4 = jtheta(4, arg1, q)
    
            cn = (v04/v3)*(v1/v4)
        finally:
	    mp.prec -= extra
        return cn

jacobi_elliptic_sn = jsn
jacobi_elliptic_cn = jcn
jacobi_elliptic_dn = jdn

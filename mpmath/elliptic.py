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

from mptypes import (mpf, mpc, mp, mpmathify, eps, one, zero)
from functions import (pi, sqrt, cos, sin, exp, tanh, ellipk, sech, nthroot)
from libmpf import to_fixed, MP_ZERO, mpf_shift, from_man_exp
from mpmath.libelefun import cos_sin

# The series for the Jacobi theta functions converge for |q| < 1;
# in the current implementation they throw a ValueError for
# abs(q) > Q_LIM
Q_LIM = 1 - 10**-7

def calculate_nome(k):
    """
    Calculate the nome, q, from the value for k.

    Useful factoids:

    k**2 = m;   m is used in Abramowitz
    """
    k = mpmathify(k)

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

    q = mpmathify(q)

    v2 = jacobi_theta(2, 0, q)
    v3 = jacobi_theta(3, 0, q)
    m = v2**2/v3**2
    return m


def _jacobi_theta2(z, q):
    if abs(q) > Q_LIM:
        raise ValueError('abs(q) > Q_LIM = %f' % Q_LIM)
    extra1 = 10
    extra2 = 20
    # the loops below break when the fixed precision quantities
    # a and b go to zero;
    # right shifting small negative numbers by wp one obtains -1, not zero,
    # so the condition a**2 + b**2 > MIN is used to break the loops.
    MIN = 2
    if z == zero:
        if isinstance(q, mpf):
            wp = mp.prec + extra1
            x = to_fixed(q._mpf_, wp)
            x2 = (x*x) >> wp
            a = b = x2
            s = x2
            while abs(a) > MIN:
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
            while are**2 + aim**2 > MIN:
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
            while abs(a) > MIN:
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
            while are**2 + aim**2 > MIN:
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
            # cos(z), siz(z) with z complex
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
            n = 1
            termre = c1re
            termim = c1im
            sre = c1re + ((are * cnre - aim * cnim) >> wp)
            sim = c1im + ((are * cnim + aim * cnre) >> wp)

            n = 3
            termre = ((are * cnre - aim * cnim) >> wp)
            termim = ((are * cnim + aim * cnre) >> wp)
            sre = c1re + ((are * cnre - aim * cnim) >> wp)
            sim = c1im + ((are * cnim + aim * cnre) >> wp)

            n = 5
            while are**2 + aim**2 > MIN:
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
                termre = ((are * cnre - aim * cnim) >> wp)
                termim = ((aim * cnre + are * cnim) >> wp)
                sre += ((are * cnre - aim * cnim) >> wp)
                sim += ((aim * cnre + are * cnim) >> wp)
                n += 2
            sre = (sre << 1)
            sim = (sim << 1)
            sre = from_man_exp(sre, -wp, mp.prec, 'n')
            sim = from_man_exp(sim, -wp, mp.prec, 'n')
            s = mpc(sre, sim)
    s *= nthroot(q, 4)
    return s

def _djacobi_theta2(z, q, nd):
    if abs(q) > Q_LIM:
        raise ValueError('abs(q) > Q_LIM = %f' % Q_LIM)
    MIN = 2
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
        if (nd&1):
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
        s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
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
        if (nd&1):
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

            if (nd&1):
                sre += ((are * sn * n**nd) >> wp)
                sim += ((aim * sn * n**nd) >> wp)
            else:
                sre += ((are * cn * n**nd) >> wp)
                sim += ((aim * cn * n**nd) >> wp)
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

        if (nd&1):
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
            if (nd&1):
                sre += ((a * snre * n**nd) >> wp)
                sim += ((a * snim * n**nd) >> wp)
            else:
                sre += ((a * cnre * n**nd) >> wp)
                sim += ((a * cnim * n**nd) >> wp)
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
        if (nd&1):
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
            #cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if (nd&1):
                sre += (((are * snre - aim * snim) * n**nd) >> wp)
                sim += (((aim * snre + are * snim) * n**nd) >> wp)
            else:
                sre += (((are * cnre - aim * cnim) * n**nd) >> wp)
                sim += (((aim * cnre + are * cnim) * n**nd) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    s *= nthroot(q, 4)
    if (nd&1):
        return (-1)**(nd//2) * s
    else:
        return (-1)**(1 + nd//2) * s

def _jacobi_theta3(z, q):
    if abs(q) > Q_LIM:
        raise ValueError('abs(q) > Q_LIM = %f' % Q_LIM)
    extra1 = 10
    extra2 = 20
    MIN = 2
    if z == zero:
        if isinstance(q, mpf):
            wp = mp.prec + extra1
            x = to_fixed(q._mpf_, wp)
            s = x
            a = b = x
            x2 = (x*x) >> wp
            while abs(a) > MIN:
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
            while are**2 + aim**2 > MIN:
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
            while abs(a) > MIN:
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
            while are**2 + aim**2 > MIN:
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
            while abs(a) > MIN:
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
            while are**2 + aim**2 > MIN:
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

def _djacobi_theta3(z, q, nd):
    """nd=1,2,3 order of the derivative with respect to z"""
    if abs(q) > Q_LIM:
        raise ValueError('abs(q) > Q_LIM = %f' % Q_LIM)
    MIN = 2
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
        if (nd&1):
            s += (a * sn) >> wp
        else:
            s += (a * cn) >> wp
        n = 2
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            if nd&1:
                s += (a * sn * n**nd) >> wp
            else:
                s += (a * cn * n**nd) >> wp
            n += 1
        s = -(s << (nd+1))
        s = mpf(from_man_exp(s, -wp, mp.prec, 'n'))
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
        if (nd&1):
            sre = (are * sn) >> wp
            sim = (aim * sn) >> wp
        else:
            sre = (are * cn) >> wp
            sim = (aim * cn) >> wp
        n = 2
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
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
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
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
        if (nd&1):
            sre = (a * snre) >> wp
            sim = (a * snim) >> wp
        else:
            sre = (a * cnre) >> wp
            sim = (a * cnim) >> wp
        n = 2
        while abs(a) > MIN:
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
            if (nd&1):
                sre += (a * snre * n**nd) >> wp
                sim += (a * snim * n**nd) >> wp
            else:
                sre += (a * cnre * n**nd) >> wp
                sim += (a * cnim * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
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
        if (nd&1):
            sre = (are * snre - aim * snim) >> wp
            sim = (aim * snre + are * snim) >> wp
        else:
            sre = (are * cnre - aim * cnim) >> wp
            sim = (aim * cnre + are * cnim) >> wp
        n = 2
        while are**2 + aim**2 > MIN:
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
            if(nd&1):
                sre += ((are * snre - aim * snim) * n**nd) >> wp
                sim += ((aim * snre + are * snim) * n**nd) >> wp
            else:
                sre += ((are * cnre - aim * cnim) * n**nd) >> wp
                sim += ((aim * cnre + are * cnim) * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
        sre = from_man_exp(sre, -wp, mp.prec, 'n')
        sim = from_man_exp(sim, -wp, mp.prec, 'n')
        s = mpc(sre, sim)
    if (nd&1):
        return (-1)**(nd//2) * s
    else:
        return (-1)**(1 + nd//2) * s

def jtheta(n, z, q):
    r"""
    Computes the Jacobi theta function theta(n, z, q), where
    n = 1, 2, 3 or 4. The theta functions are functions of two
    variables:

    - z is the *argument*, an arbitrary real or complex number

    - q is the *nome*, which must be a real or complex number
      in the unit disk (i.e. abs(q) < 1)

    One also commonly encounters the notation theta(n, z, tau) in the
    literature. The variable tau is called the *parameter* and can
    be converted to a nome using the formula ``q = exp(j*pi*tau)``.
    Note the condition abs(q) < 1 requires imag(tau) > 0; i.e.
    Jacobi theta functions are defined for tau in the upper half
    plane.

    Other notations are also in use. For example, some authors use
    the single-argument form theta(n, x). Depending on context, this
    can mean ``jtheta(n, 0, x)``, ``jtheta(n, x, q)``, or possibly
    something else. Needless to say, it is a good idea to cross-check
    the definitions when working with theta functions.

    **Definition**

    The four Jacobi theta functions as implemented by :func:`jtheta`
    are defined by the following infinite series::

                               oo
                               ___           2
                          1/4 \         n  (n  + n)
        theta(1,z,q) = 2 q     )    (-1)  q         sin((2n+1) z)
                              /___
                              n = 0

                               oo
                               ___     2
                          1/4 \      (n  + n)
        theta(2,z,q) = 2 q     )    q         cos((2n+1) z)
                              /___
                              n = 0

                               oo
                               ___     2
                              \      (n )
        theta(3,z,q) = 1 +  2  )    q     cos(2n z)
                              /___
                              n = 1

                               oo
                               ___        2
                              \         (n )
        theta(4,z,q) = 1 +  2  )    (-q)     cos(2n z)
                              /___
                              n = 1

    For abs(q) << 1, these series converge very quickly, so the
    Jacobi theta functions can efficiently be evaluated to high
    precision.

    **Examples and basic properties**

    Considered as functions of z, the Jacobi theta functions may be
    viewed as generalizations of the ordinary trigonometric functions
    cos and sin. They are periodic functions::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print jtheta(1, 0.1, 1/5.)
        0.117756191842059
        >>> print jtheta(1, 0.1 + 2*pi, 1/5.)
        0.117756191842059

    Indeed, the series defining the theta functions are essentially
    trigonometric Fourier series. The coefficients can be retrieved
    using :func:`fourier`::

        >>> nprint(fourier(lambda x: jtheta(2, x, 0.5), [-pi, pi], 4))
        ([0.0, 1.68179, 0.0, 0.420448, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0])

    The Jacobi theta functions are also so-called quasiperiodic
    functions of z and tau, meaning that for fixed tau, theta(z, q)
    and theta(z+tau*pi, q) are the same except for an exponential
    factor::

        >>> tau = 0.3*j
        >>> q = exp(pi*j*tau)
        >>> z = 10
        >>> print jtheta(4, z+tau*pi, q)
        (-0.682420280786035 + 1.5266839997214j)
        >>> print -exp(-2*j*z)/q * jtheta(4, z, q)
        (-0.682420280786035 + 1.5266839997214j)

    The Jacobi theta functions satisfy a huge number of other
    functional equations, such as the following identity (valid for
    any q)::

        >>> q = 0.3
        >>> print jtheta(3,0,q)**4
        6.82374408935276
        >>> print jtheta(2,0,q)**4 + jtheta(4,0,q)**4
        6.82374408935276

    Extensive listings of identities satisfied by the Jacobi theta
    functions can be found in standard reference works.

    The Jacobi theta functions are related to the gamma function
    for special arguments::

        >>> print jtheta(3, 0, exp(-pi))
        1.08643481121331
        >>> print pi**(1/4.) / gamma(3/4.)
        1.08643481121331

    :func:`jtheta` supports arbitrary precision evaluation and complex
    arguments::

        >>> mp.dps = 50
        >>> print jtheta(4, sqrt(2), 0.5)
        2.0549510717571539127004115835148878097035750653737
        >>> mp.dps = 25
        >>> print jtheta(4, 1+2j, (1+j)/5)
        (7.180331760146805926357227 - 1.634292858119162417301859j)

    **Possible issues**

    For abs(q) >= 1 or imag(tau) <= 0, :func:`jtheta` raises
    ``ValueError``. This exception is also raised abs(q) extremely
    close to 1 (or equivalently tau very close to 0), since the
    series would converge too slowly::

        >>> jtheta(1, 10, 0.99999999 * exp(0.5*j))
        Traceback (most recent call last):
          ...
        ValueError: abs(q) > Q_LIM = 1.000000

    """
    z = mpmathify(z)
    q = mpmathify(q)

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

def djtheta(n, z, q, nd=1):
    """
    For an integer nd >= 1, computes the nd:th derivative with respect
    to z of the Jacobi theta function ``jtheta(n, z, q)``. For
    example::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> print djtheta(3, 7, 0.2)
        -0.795947847483158
        >>> print diff(lambda x: jtheta(3, x, 0.2), 7)
        -0.795947847483158

    For additional details, see :func:`jtheta`.
    """

    z = mpmathify(z)
    q = mpmathify(q)

    extra = 10
    prec0 = mp.prec
    try:
        mp.prec += extra
        if n == 1:
            res = _djacobi_theta2(z - pi/2, q, nd)
        elif n == 2:
            res = _djacobi_theta2(z, q, nd)
        elif n == 3:
            res = _djacobi_theta3(z, q, nd)
        elif n == 4:
            res = _djacobi_theta3(z, -q, nd)
        else:
            raise ValueError
    finally:
        mp.prec = prec0
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

if __name__ == '__main__':
    import doctest
    doctest.testmod()

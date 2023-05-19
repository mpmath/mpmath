"""
Torture tests for asymptotics and high precision evaluation of
special functions.

(Other torture tests may also be placed here.)

Running this file (gmpy recommended!) takes several CPU minutes.
The multiprocessing module is used automatically to run tests
in parallel if many cores are available. (A single test may take between
a second and several minutes; possibly more.)

The idea:

* We evaluate functions at positive, negative, imaginary, 45- and 135-degree
  complex values with magnitudes between 10^-20 to 10^20, at precisions between
  5 and 150 digits (we can go even higher for fast functions).

* Comparing the result from two different precision levels provides
  a strong consistency check (particularly for functions that use
  different algorithms at different precision levels).

* That the computation finishes at all (without failure), within reasonable
  time, provides a check that evaluation works at all: that the code runs,
  that it doesn't get stuck in an infinite loop, and that it doesn't use
  some extremely slowly algorithm where it could use a faster one.

TODO:

* Speed up those functions that take long to finish!
* Generalize to test more cases; more options.
* Implement a timeout mechanism.
* Some functions are notably absent, including the following:
  * inverse trigonometric functions (some become inaccurate for complex arguments)
  * ci, si (not implemented properly for large complex arguments)
  * zeta functions (need to modify test not to try too large imaginary values)
  * and others...

"""

import pytest

from mpmath import (agm, airyai, airybi, apery, barnesg, bernfrac, bernoulli,
                    besseli, besselj, besselk, bessely, catalan, cbrt, chi, ci,
                    cos, cosh, coulombf, coulombg, e, e1, ei, ellipe, ellipk,
                    erf, erfc, erfi, euler, exp, expint, expm1, gamma,
                    gammainc, glaisher, hermite, hyp0f1, hyp1f1, hyp1f2,
                    hyp2f0, hyp2f1, hyp2f2, hyp2f3, hyperu, j, jtheta,
                    khinchin, lambertw, legendre, legenp, legenq, li, ln, ln2,
                    ln10, loggamma, mertens, mp, mpf, phi, pi, polylog, power,
                    root, shi, si, sin, sinh, sqrt, stieltjes, tan, tanh,
                    twinprime, workprec)


a1, a2, a3, a4, a5 = 1.5, -2.25, 3.125, 4, 2

@pytest.mark.parametrize('f,maxdps,huge_range',
                         [(lambda z: +pi, 10000, False),
                          (lambda z: +e, 10000, False),
                          (lambda z: +ln2, 10000, False),
                          (lambda z: +ln10, 10000, False),
                          (lambda z: +phi, 10000, False),
                          (lambda z: +catalan, 5000, False),
                          (lambda z: +euler, 5000, False),
                          (lambda z: +glaisher, 1000, False),
                          (lambda z: +khinchin, 1000, False),
                          (lambda z: +twinprime, 150, False),
                          (lambda z: stieltjes(2), 150, False),
                          (lambda z: +mertens, 150, False),
                          (lambda z: +apery, 5000, False),
                          (sqrt, 10000, True),
                          (cbrt, 5000, True),
                          (lambda z: root(z,4), 5000, True),
                          (lambda z: root(z,-5), 5000, True),
                          (exp, 5000, True),
                          (expm1, 1500, False),
                          (ln, 5000, True),
                          (cosh, 5000, False),
                          (sinh, 5000, False),
                          (tanh, 1500, False),
                          (sin, 5000, True),
                          (cos, 5000, True),
                          (tan, 1500, False),
                          (agm, 1500, True),
                          (ellipk, 1500, False),
                          (ellipe, 1500, False),
                          (lambertw, 150, True),
                          (lambda z: lambertw(z,-1), 150, False),
                          (lambda z: lambertw(z,1), 150, False),
                          (lambda z: lambertw(z,4), 150, False),
                          (gamma, 150, False),
                          (loggamma, 150, False),  # True ?
                          (ei, 150, False),
                          (e1, 150, False),
                          (li, 150, True),
                          (ci, 150, False),
                          (si, 150, False),
                          (chi, 150, False),
                          (shi, 150, False),
                          (erf, 150, False),
                          (erfc, 150, False),
                          (erfi, 150, False),
                          (lambda z: besselj(2, z), 150, False),
                          (lambda z: bessely(2, z), 150, False),
                          (lambda z: besseli(2, z), 150, False),
                          (lambda z: besselk(2, z), 150, False),
                          (lambda z: besselj(-2.25, z), 150, False),
                          (lambda z: bessely(-2.25, z), 150, False),
                          (lambda z: besseli(-2.25, z), 150, False),
                          (lambda z: besselk(-2.25, z), 150, False),
                          (airyai, 150, False),
                          (airybi, 150, False),
                          (lambda z: hyp0f1(a1, z), 150, False),
                          (lambda z: hyp1f1(a1, a2, z), 150, False),
                          (lambda z: hyp1f2(a1, a2, a3, z), 150, False),
                          (lambda z: hyp2f0(a1, a2, z), 150, False),
                          (lambda z: hyperu(a1, a2, z), 150, False),
                          (lambda z: hyp2f1(a1, a2, a3, z), 150, False),
                          (lambda z: hyp2f2(a1, a2, a3, a4, z), 150, False),
                          (lambda z: hyp2f3(a1, a2, a3, a4, a5, z), 150, False),
                          (lambda z: coulombf(a1, a2, z), 150, False),
                          (lambda z: coulombg(a1, a2, z), 150, False),
                          (lambda z: polylog(2,z), 150, False),
                          (lambda z: polylog(3,z), 150, False),
                          (lambda z: polylog(-2,z), 150, False),
                          (lambda z: expint(4, z), 150, False),
                          (lambda z: expint(-4, z), 150, False),
                          (lambda z: expint(2.25, z), 150, False),
                          (lambda z: gammainc(2.5, z, 5), 150, False),
                          (lambda z: gammainc(2.5, 5, z), 150, False),
                          (lambda z: hermite(3, z), 150, False),
                          (lambda z: hermite(2.5, z), 150, False),
                          (lambda z: legendre(3, z), 150, False),
                          (lambda z: legendre(4, z), 150, False),
                          (lambda z: legendre(2.5, z), 150, False),
                          (lambda z: legenp(a1, a2, z), 150, False),
                          (lambda z: legenq(a1, a2, z), 90, False),  # abnormally slow
                          (lambda z: jtheta(1, z, 0.5), 150, False),
                          (lambda z: jtheta(2, z, 0.5), 150, False),
                          (lambda z: jtheta(3, z, 0.5), 150, False),
                          (lambda z: jtheta(4, z, 0.5), 150, False),
                          (lambda z: jtheta(1, z, 0.5, 1), 150, False),
                          (lambda z: jtheta(2, z, 0.5, 1), 150, False),
                          (lambda z: jtheta(3, z, 0.5, 1), 150, False),
                          (lambda z: jtheta(4, z, 0.5, 1), 150, False),
                          (barnesg, 90, False)])
def test_asymp(f, maxdps, huge_range):
    dps = [5,15,25,50,90,150,500,1500,5000,10000]
    dps = [p for p in dps if p <= maxdps]
    def check(x,y,p,inpt):
        assert abs(x-y)/abs(y) < workprec(20)(power)(10, -p+1)
    exponents = list(range(-20,20))
    if huge_range:
        exponents += [-1000, -100, -50, 50, 100, 1000]
    for n in exponents:
        mp.dps = 25
        xpos = mpf(10)**n / 1.1287
        xneg = -xpos
        ximag = xpos*j
        xcomplex1 = xpos*(1+j)
        xcomplex2 = xpos*(-1+j)
        for i in range(len(dps)):
            mp.dps = dps[i]
            new = f(xpos), f(xneg), f(ximag), f(xcomplex1), f(xcomplex2)
            if i != 0:
                p = dps[i-1]
                check(prev[0], new[0], p, xpos)
                check(prev[1], new[1], p, xneg)
                check(prev[2], new[2], p, ximag)
                check(prev[3], new[3], p, xcomplex1)
                check(prev[4], new[4], p, xcomplex2)
            prev = new


def test_bernoulli_huge():
    p, q = bernfrac(9000)
    assert p % 10**10 == 9636701091
    assert q == 4091851784687571609141381951327092757255270
    mp.dps = 15
    assert str(bernoulli(10**100)) == '-2.58183325604736e+987675256497386331227838638980680030172857347883537824464410652557820800494271520411283004120790908623'
    mp.dps = 50
    assert str(bernoulli(10**100)) == '-2.5818332560473632073252488656039475548106223822913e+987675256497386331227838638980680030172857347883537824464410652557820800494271520411283004120790908623'

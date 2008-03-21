"""
Miscellaneous useful (or just interesting :) mathematical constants.
"""

from mpmath import *
from mpmath.lib import *
from mpmath.specfun.misc import bernoulli2n, logk
from mpmath.mptypes import constant


#----------------------------------------------------------------------#
# The golden ratio is given by phi = (1 + sqrt(5))/2

@constant_memo
def phi_fixed(prec):
    prec += 10
    sqrt = [sqrt_fixed2, sqrt_fixed][prec < 20000]
    a = sqrt(5<<prec, prec) + (1 << prec)
    return a >> 11


#----------------------------------------------------------------------#

# Catalan's constant is computed using Lupas's rapidly convergent series
# (listed on http://mathworld.wolfram.com/CatalansConstant.html)
#            oo
#            ___       n-1  8n     2                   3    2
#        1  \      (-1)    2   (40n  - 24n + 3) [(2n)!] (n!)
#  K =  ---  )     -----------------------------------------
#       64  /___               3               2
#                             n  (2n-1) [(4n)!]
#           n = 1

@constant_memo
def catalan_fixed(prec):
    prec = prec + 20
    a = one = 1 << prec
    s, t, n = 0, 1, 1
    while t:
        a *= 32 * n**3 * (2*n-1)
        a //= (3-16*n+16*n**2)**2
        t = a * (-1)**(n-1) * (40*n**2-24*n+3) // (n**3 * (2*n-1))
        s += t
        n += 1
    return s >> (20 + 6)

#----------------------------------------------------------------------#

# Euler's constant (gamma) is computed using the Brent-McMillan formula,
# gamma ~= A(n)/B(n) - log(n), where

#   A(n) = sum_{k=0,1,2,...} (n**k / k!)**2 * H(k)
#   B(n) = sum_{k=0,1,2,...} (n**k / k!)**2
#   H(k) = 1 + 1/2 + 1/3 + ... + 1/k

# The error is bounded by O(exp(-4n)). Choosing n to be a power
# of two, 2**p, the logarithm becomes particularly easy to calculate.

# Reference:
# Xavier Gourdon & Pascal Sebah, The Euler constant: gamma
# http://numbers.computation.free.fr/Constants/Gamma/gamma.pdf

@constant_memo
def euler_fixed(prec):
    prec += 30
    # choose p such that exp(-4*(2**p)) < 2**-n
    p = int(math.log((prec/4) * math.log(2), 2)) + 1
    n = 1<<p
    r = one = 1<<prec
    H, A, B, npow, k, d = 0, 0, 0, 1, 1, 1
    while r:
        A += (r * H) >> prec
        B += r
        r = r * (n*n) // (k*k)
        H += one // k
        k += 1
    S = ((A<<prec) // B) - p*log2_fixed(prec)
    return S >> 30


#----------------------------------------------------------------------#

# Khinchin's constant is relatively difficult to compute. Here
# we use the rational zeta series

#                    oo                2*n-1
#                   ___                ___
#                   \   ` zeta(2*n)-1  \   ` (-1)^(k+1)
#  log(K)*log(2) =   )    ------------  )    ----------
#                   /___.      n       /___.      k
#                   n = 1              k = 1

# which adds half a digit per term. The essential trick for achieving
# reasonable efficiency is to recycle both the values of the zeta
# function (essentially Bernoulli numbers) and the partial terms of
# the inner sum.

# An alternative might be to use K = 2*exp[1/log(2) X] where

#      / 1     1       [ pi*x*(1-x^2) ]
#  X = |    ------ log [ ------------ ].
#      / 0  x(1+x)     [  sin(pi*x)   ]

# and integrate numerically. In practice, this seems to be slightly
# slower than the zeta series at high precision.

@constant_memo
def khinchin_fixed(prec):
    orig = mp.prec
    try:
        mp.prec = int(prec + prec**0.5 + 5)
        s = mpf(0)
        t = one = mpf(1)
        B = bernoulli2n()
        fac = mpf(4)
        pipow = twopi2 = (2*pi)**2
        n = 1
        while 1:
            zeta2n = (-1)**(n+1) * B.next() * pipow / fac
            term = ((zeta2n - 1) * t) / n
            # print n, nstr(term)
            if term < eps:
                break
            s += term
            t += (one/(2*n+1) - one/(2*n))
            n += 1
            fac *= (2*n)*(2*n-1)
            pipow *= twopi2
        return to_fixed(exp(s/ln2)._mpf_, prec)
    finally:
        mp.prec = orig

#----------------------------------------------------------------------#

# Glaisher's constant is defined as A = exp(1/2 - zeta'(-1)).

# One way to compute it would be to perform direct numerical
# differentiation, but computing arbitrary Riemann zeta function
# values at high precision is expensive. We instead use the formula

#     A = exp((6 (-zeta'(2))/pi^2 + log 2 pi + gamma)/12)

# and compute zeta'(2) from the series representation

#              oo
#              ___
#             \     log k
#  -zeta'(2) = )    -----
#             /___     2
#                    k
#            k = 2

# This series converges exceptionally slowly, but can be accelerated
# using Euler-Maclaurin formula. The important insight is that the
# E-M integral can be done in closed form and that the high order
# are given by

#    n  /       \
#   d   | log x |   a + b log x
#   --- | ----- | = -----------
#     n |   2   |      2 + n
#   dx  \  x    /     x

# where a and b are integers given by a simple recurrence. Note
# that just one logarithm is needed. However, lots of integer
# logarithms are required for the initial summation.

# This algorithm could possibly be turned into a faster algorithm
# for general evaluation of zeta(s) or zeta'(s); this should be
# looked into.

@constant_memo
def glaisher_fixed(prec):
    orig = mp.prec
    try:
        dps = mp.dps
        mp.prec = prec + 20
        N = int(1.0*dps + 5)
        logs = logk()
        s = mpf(0)
        # E-M step 1: sum log(k)/k**2 for k = 2..N-1
        for n in range(2, N):
            # print n, N
            logn = logs.next()
            s += logn / n**2
        logN = logs.next()
        # E-M step 2: integral of log(x)/x**2 from N to inf
        s += (1+logN)/N
        # E-M step 3: endpoint correction term f(N)/2
        s += logN/(N**2 * 2)
        # E-M step 4: the series of derivatives
        pN, a, b, j, B2k, fac, k = N**3, 1, -2, 3, bernoulli2n(), 2, 1
        while 1:
            # D(2*k-1) * B(2*k) / fac(2*k) [D(n) = nth derivative]
            D = (a+b*logN)/pN
            B = B2k.next()
            term = B * D / fac
            if abs(term) < eps:
                break
            # print k, nstr(term)
            s -= term
            # Advance derivative twice
            a, b, pN, j = b-a*j, -j*b, pN*N, j+1
            a, b, pN, j = b-a*j, -j*b, pN*N, j+1
            k += 1
            fac *= (2*k) * (2*k-1)
        A = exp((6*s/pi**2 + log(2*pi) + euler)/12)
        return to_fixed(A._mpf_, prec)
    finally:
        mp.prec = orig

#----------------------------------------------------------------------#

# Apery's constant can be computed using the very rapidly convergent
# series
#              oo
#              ___              2                      10
#             \         n  205 n  + 250 n + 77     (n!)
#  zeta(3) =   )    (-1)   -------------------  ----------
#             /___               64                      5
#             n = 0                             ((2n+1)!)

@constant_memo
def apery_fixed(prec):
    prec += 20
    d = 1 << prec
    term = 77 << prec
    n = 1
    s = 0
    while term:
        s += term
        d *= (n**10)
        d //= (((2*n+1)**5) * (2*n)**5)
        term = (-1)**n * (205*(n**2) + 250*n + 77) * d
        n += 1
    return s >> (20 + 6)


#----------------------------------------------------------------------#

fme = from_man_exp

def mpf_phi(p, r): return fme(phi_fixed(p+5), -p-5, p, r)
def mpf_khinchin(p, r): return fme(khinchin_fixed(p+5), -p-5, p, r)
def mpf_glaisher(p, r): return fme(glaisher_fixed(p+5), -p-5, p, r)
def mpf_apery(p, r): return fme(apery_fixed(p+5), -p-5, p, r)
def mpf_euler(p, r): return fme(euler_fixed(p+5), -p-5, p, r)
def mpf_catalan(p, r): return fme(catalan_fixed(p+5), -p-5, p, r)

phi = constant(mpf_phi, "Golden ratio (phi)")
catalan = constant(mpf_catalan, "Catalan's constant")
euler = constant(mpf_euler, "Euler's constant (gamma)")
khinchin = constant(mpf_khinchin, "Khinchin's constant")
glaisher = constant(mpf_glaisher, "Glaisher's constant")
apery = constant(mpf_apery, "Apery's constant")

"""
This script calculates solutions to some of the problems from the
"Many Digits" competition:
http://www.cs.ru.nl/~milad/manydigits/problems.php

Run with:

    python manydigits.py

"""

from mpmath import *
from mpmath.libmp import to_fixed, bin_to_radix

dps = 100
mp.dps = dps + 10

def pr(x):
    """Return the first dps digits after the decimal point"""
    x = x._mpf_
    p = int(dps*3.33 + 10)
    t = to_fixed(x, p)
    d = bin_to_radix(t, p, 10, dps)
    s = str(d).zfill(dps)[-dps:]
    return s[:dps//2] + "\n" + s[dps//2:]

print """
This script prints answers to a selection of the "Many Digits"
competition problems: http://www.cs.ru.nl/~milad/manydigits/problems.php

The output for each problem is the first 100 digits after the
decimal point in the result.
"""

print "C01: sin(tan(cos(1)))"
print pr(sin(tan(cos(1))))
print

print "C02: sqrt(e/pi)"
print pr(sqrt(e/pi))
print

print "C03: sin((e+1)^3)"
print pr(sin((e+1)**3))
print

print "C04: exp(pi*sqrt(2011))"
mp.dps += 65
print pr(exp(pi*sqrt(2011)))
mp.dps -= 65
print

print "C05: exp(exp(exp(1/2)))"
print pr(exp(exp(exp(0.5))))
print

print "C06: arctanh(1-arctanh(1-arctanh(1-arctanh(1/pi))))"
print pr(atanh(1-atanh(1-atanh(1-atanh(1/pi)))))
print

print "C07: pi^1000"
mp.dps += 505
print pr(pi**1000)
mp.dps -= 505
print

print "C08: sin(6^(6^6))"
print pr(sin(6**(6**6)))
print

print "C09: sin(10*arctan(tanh(pi*(2011^(1/2))/3)))"
mp.dps += 150
print pr(sin(10*atan(tanh(pi*sqrt(2011)/3))))
mp.dps -= 150
print

print "C10: (7+2^(1/5)-5*(8^(1/5)))^(1/3) + 4^(1/5)-2^(1/5)"
a = mpf(1)/5
print pr(((7 + 2**a - 5*(8**a))**(mpf(1)/3) + 4**a - 2**a))
print

print "C11: tan(2^(1/2))+arctanh(sin(1))"
print pr((tan(sqrt(2)) + atanh(sin(1))))
print

print "C12: arcsin(1/e^2) + arcsinh(e^2)"
print pr(asin(1/exp(2)) + asinh(exp(2)))
print

print "C17: S= -4*Zeta(2) - 2*Zeta(3) + 4*Zeta(2)*Zeta(3) + 2*Zeta(5)"
print pr(-4*zeta(2) - 2*zeta(3) + 4*zeta(2)*zeta(3) + 2*zeta(5))
print

print "C18: Catalan G = Sum{i=0}{\infty}(-1)^i/(2i+1)^2"
print pr(catalan)
print

print "C21: Equation exp(cos(x)) = x"
print pr(findroot(lambda x: exp(cos(x))-x, 1))
print

print "C22: J = integral(sin(sin(sin(x)))), x=0..1"
print pr(quadts(lambda x: sin(sin(sin(x))), [0, 1]))
print

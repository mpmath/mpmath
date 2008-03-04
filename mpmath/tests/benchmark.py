# Note: the code in this file is a big pile of ugly

from mpmath import *
from mpmath.lib import *
from decimal import getcontext, Decimal
import dmath
from random import seed, random, randint

def getrandom(type):
    if type == 'mpf':
        return mpf(random() * 2.0**randint(-10, 10)) ** 0.5
    if type == 'mpfval':
        return (mpf(random() * 2.0**randint(-10, 10)) ** 0.5).val
    if type == 'Decimal':
        return Decimal(repr(random() * 2.0**randint(-10, 10))).sqrt()
    raise TypeError

def rndnums(type, N):
    seed(1234)
    xs = [getrandom(type) for i in xrange(N)]
    ys = xs[::-1]
    xys = zip(xs, ys)
    return xs, ys, xys

def setprec(type, prec):
    if type == 'Decimal':
        getcontext().prec = prec
    else:
        mpf.dps = prec
        # change prec value to bits for mpfval use
        prec = mpf.prec
    return prec

testcode = \
"""
def testit(prec, N):
    from time import clock
    RF = round_half_even
    prec = setprec('TYPE', prec)
    xs, ys, xys = rndnums('TYPE', N)
    t = 1e100
    for i in range(3):
        t1 = clock()
        for x, y in xys:
            OP; OP; OP; OP; OP; OP; OP; OP; OP; OP;
        t2 = clock()
        t = min(t, (t2-t1)/10)
    return t
"""

tests = []
atests = [
  ('Convert to integer (int(x))', 'int(x)', 'int(x)', 'to_int(x)'),
  ('Convert to string (str(x))', 'str(x)', 'str(x)', 'to_str(x, int(prec/3.321))'),
  ('Convert to float (float(x))', 'float(x)', 'float(x)', 'to_float(x)'),
  ('Equality (x==y)', 'x==y', 'x==y', 'feq(x, y)'),
  ('Comparison (x<y)', 'x<y', 'x<y', 'fcmp(x, y) < 0'),
  ('Addition (x+y)', 'x+y', 'x+y', 'fadd(x, y, prec, RF)'),
  ('Subtraction (x-y)', 'x+y', 'x+y', 'fsub(x, y, prec, RF)'),
  ('Multiplication (x*y)', 'x*y', 'x*y', 'fmul(x, y, prec, RF)'),
  ('Division (x/y)', 'x/y', 'x/y', 'fdiv(x, y, prec, RF)'),
  ('Square root (x^0.5)', 'x.sqrt()', 'sqrt(x)', 'fsqrt(x, prec, RF)'),
  ('Integer power (x^42)', 'x**42', 'x**42', 'fpowi(x, 42, prec, RF)'),
#  ('Exponential function (exp(x))', 'dmath.exp(x)', 'exp(x)', 'fexp(x, prec, RF)'),
#  ('Natural logarithm (log(x))', 'dmath.log(x+1)', 'log(x)', 'flog(x, prec, RF)'),
#  ('Sine (sin(x))', 'dmath.sin(x)', 'sin(x)', 'fsin(x, prec, RF)'),
#  ('Tangent (tan(x))', 'dmath.tan(x)', 'tan(x)', 'ftan(x, prec, RF)'),
#  ('Inverse tangent(atan(x))', 'dmath.atan(x)', 'atan(x)', 'fatan(x, prec, RF)'),
#  ('Hyperbolic cosine (cosh(x))', 'dmath.cosh(x)', 'cosh(x)', 'fcosh(x, prec, RF)')
]

slow = ["power", "exp", "log", "sin", "tan", "cos"]

for op in atests:
    cases = [op[0]]
    if op[1]:
        exec testcode.replace("OP", op[1]).replace("TYPE", "Decimal")
        cases += [testit]
    else:
        cases += [None]
    exec testcode.replace("OP", op[2]).replace("TYPE", "mpf")
    cases += [testit]
    exec testcode.replace("OP", op[3]).replace("TYPE", "mpfval")
    cases += [testit]
    tests.append(cases)

def rd(x):
    if x > 100000: return int(x // 10000) * 10000
    if x > 10000: return int(x // 1000) * 1000
    if x > 1000: return int(x // 100) * 100
    if x > 100: return int(x // 10) * 10
    return int(x)

def runtests():
    results = []
    for test in tests:
        name, dectest, mpftest, mpfvaltest = test
        if any(s in name for s in slow):
            N = 1
            precs = [15, 30, 100, 300]
        else:
            N = 10
            precs = [15, 30, 100, 300, 1000]
        header_name = "*" + name + "*"
        rows = []
        for prec in precs:
            print name, prec
            if dectest is None:
                t1 = 1e1000-1e1000
            else:
                t1 = dectest(prec, N)
            t2 = mpftest(prec, N)
            t3 = mpfvaltest(prec, N)
            s = []
            s += ["%i" % prec]
            s += [str(rd(N/t1))]
            s += [str(rd(N/t2)) + (" (%.1fx)" % (t1/t2))]
            s += [str(rd(N/t3)) + (" (%.1fx)" % (t1/t3))]
            rows.append(s)
        results.append((header_name, rows))
    return results

import gc
gc.disable()
results1 = runtests()
import psyco
psyco.full()
results2 = runtests()

import sys
sys.stdout = open("results.txt", "w")

for r1, r2 in zip(results1, results2):
    name = r1[0]
    print name
    print "|| *digits* || *Decimal* || *mpf* ||  *raw mpf* || *Decimal+psyco* || *mpf+psyco* || *raw mpf+psyco* ||"
    for a, b in zip(r1[1], r2[1]):
        cols = a + b[1:]
        print "|| " + (" || ".join(cols)) + " ||"
    print


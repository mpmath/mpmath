from time import clock
from random import *
from mpmath import *

w = """
def f(N):
    xs = []
    seed(1234)
    for i in range(N):
        x = mpf(random() * 2.0**randint(-10, 10)) ** 0.5
        y = mpf(random() * 2.0**randint(-10, 10)) ** 0.5
        #x *= choice([1, -1])
        #y *= choice([1, -1])
        xs.append((x, y))
    x, y = xs[0]
    OP
    t = 100.0
    for i in range(5):
        t1 = clock()
        for x, y in xs:
            OP
        t2 = clock()
        t = min(t, t2-t1)
    times = N/t
    if times > 100000: return int(times // 10000) * 10000
    if times > 10000: return int(times // 1000) * 1000
    if times > 1000: return int(times // 100) * 100
    if times > 100: return int(times // 10) * 10
    return int(times)

tests.append(("OP", f))
"""

tests = []
def deftest(op):
    exec w.replace("OP", op)

atests = ["int(x)", "float(x)", "str(x)", "+x", "-x",
  "x * 10", "x + 10", "x / 10",
  "x + y", "x - y", "x * y", "x / y", "x == y",
  "x < y", "abs(x)", "sqrt(x)", "x**0.5",
  "exp(x)", "log(x)", "x**42", "x**y",
  "sin(x)", "tan(x)", "atan(x)", "cosh(x)", "asin(x/32)",
  "erf(x)", "gamma(x+1)", "zeta(x)"]

for test in atests:
    deftest(test)

precs = [15, 30, 100, 500, 1000]

def runtests():
    print "\n    op / dps :",
    for prec in precs:
        print "%7s" % prec,
    print
    print "-" * 75
    for name, test in tests:
        print ("%12s" % name), ":",
        for prec in precs:
            mpf.dps = prec
            print "%7s" % test(10),
            mpf.dps = 15
        print

print "\nmpf timings (operations / second)"

#mpf.round_floor()
#import psyco
#psyco.full()

runtests()

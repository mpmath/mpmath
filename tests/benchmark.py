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
        x *= choice([1, -1])
        y *= choice([1, -1])
        xs.append((x, y))
    x, y = xs[0]
    OP
    t1 = clock()
    for x, y in xs:
        OP
    t2 = clock()
    return int(N/(t2-t1))

tests.append(("OP", f))
"""

tests = []
def deftest(op):
    exec w.replace("OP", op)

atests = ["x + y", "x - y", "x * y", "x / y", "x == y",
  "x < y", "abs(x)", "sqrt(abs(x))", "abs(x)**0.5",
  "exp(x)", "log(abs(x))", "x**42", "x**y",
  "sin(x)", "tan(x)", "atan(x)", "cosh(x)", "asin(x)",
  "gamma(x)", "erf(x)", "zeta(x)"]

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
            print "%7s" % test(max(50, 5000//prec)),
            mpf.dps = 15
        print

print "\nmpf timings (operations / second)"

#mpf.round_floor()
#import psyco
#psyco.full()

runtests()

"""
Interval arithmetic demo: estimating error of numerical Taylor series.

This module can be run interactively with

    python taylor.py

"""

from mpmath import *

def taylor(x, n):
    print "-"*75
    t = x = mpi(x)
    s = 1
    print "adding 1"
    print s, "\n"
    s += t
    print "adding x"
    print s, "\n"
    for k in range(2, n+1):
        t = (t * x) / k
        s += t
        print "adding x^%i / %i! ~= %s" % (k, k, t.mid)
        print s, "\n"
    print "-"*75
    return s

# Note: this should really be computed using interval arithmetic too!
def remainder(x, n):
    xi = max(0, x)
    r = exp(xi) / factorial(n+1)
    r = r * x**(n+1)
    return abs(r)

def exponential(x, n):
    """
    Compute exp(x) using n terms of the Taylor series for exp using
    intervals, and print detailed error analysis.
    """
    t = taylor(x, n)
    r = remainder(x, n)
    expx = exp(x)
    print "Correct value of exp(x):    ", expx
    print
    print "Computed interval:          "
    print t
    print
    print "Computed value (midpoint):  ", t.mid
    print
    print "Estimated rounding error:   ", t.delta
    print "Estimated truncation error: ", r
    print "Estimated total error:      ", t.delta + r
    print "Actual error                ", abs(expx - t.mid)
    print
    u = t + mpi(-r, r)
    print "Interval with est. truncation error added:"
    print u
    print
    print "Correct value contained in computed interval:", t.a <= expx <= t.b
    print "When accounting for truncation error:", u.a <= expx <= u.b

if __name__ == "__main__":
    print "Interval arithmetic demo"
    print
    print "This script sums the Taylor series for exp(x) using interval arithmetic,"
    print "and then compares the numerical errors due to rounding and truncation."
    print
    x = mpf(raw_input("Enter the value of x (e.g. 3.5): "))
    n = int(raw_input("Enter the number of terms n (e.g. 10): "))
    print
    exponential(x, n)

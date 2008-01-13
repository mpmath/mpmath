"""
Interval arithmetic demo: estimating error of numerical Taylor series.

"""

from mpmath import *

def taylor(x, n):
    t = x = mpi(x)
    s = 1 + t
    for k in range(2, n+1):
        t = (t * x) / k
        s += t
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
    print "Correct value:       ", expx
    print "Computed midpoint:   ", t.mid
    print "Computed interval:   ", t
    print "Est rounding error:  ", t.delta
    print "Est truncation error:", r
    print "Est error:           ", t.delta + r
    print "Actual error         ", abs(expx - t.mid)

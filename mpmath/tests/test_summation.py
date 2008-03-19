from mpmath import *

def test_sumem():
    f = lambda x: x**-3
    fp = lambda x, k: (-1)**k * factorial(k+2)/2 * x**(-k-3)
    for prec in [15, 50]:
        s, err = sumem(f, 1, inf, fderiv=fp)
        assert s.ae(zeta(3))
    mp.dps = 15
    assert sumem(lambda k: k**4 + 3*k + 1, 10, 100, N=5)[0].ae(2050333103)

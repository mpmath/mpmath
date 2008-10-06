from mpmath import *
from mpmath.optimization import *

def test_findroot():
    # old tests, assuming secant
    mp.dps = 15
    assert findroot(lambda x: 4*x-3, mpf(5)).ae(0.75)
    assert findroot(sin, mpf(3)).ae(pi)
    assert findroot(sin, (mpf(3), mpf(3.14))).ae(pi)
    assert findroot(lambda x: x*x+1, mpc(2+2j)).ae(1j)
    # test all solvers with 1 starting point
    f = lambda x: cos(x)
    for solver in [Secant, MNewton, Muller]:
        x = findroot(f, 2., solver=solver)
        assert abs(f(x)) < eps
    # test all solvers with interval of 2 points
    # FIXME: 
    for solver in [Secant, Muller, Bisection, Illinois, Pegasus, Anderson,
                   Ridder]:
        x = findroot(f, (1., 2.), solver=solver)
        assert abs(f(x)) < eps

# TODO: test MNewton, Muller

def test_multiplicity():
    for i in xrange(1, 5):
        assert multiplicity(lambda x: (x - 1)**i, 1) == i
    assert multiplicity(lambda x: x**2, 1) == 0

import cmath
import math

from mpmath import (chop, cos, diff, diffs, diffun, e, exp, fp, j, log, sqrt,
                    sin, taylor)


def test_diff():
    assert diff(log, 2.0, n=0).ae(log(2))
    assert diff(cos, 1.0).ae(-sin(1))
    assert diff(abs, 0.0) == 0
    assert diff(abs, 0.0, direction=1) == 1
    assert diff(abs, 0.0, direction=-1) == -1
    assert diff(exp, 1.0).ae(e)
    assert diff(exp, 1.0, n=5).ae(e)
    assert diff(exp, 2.0, n=5, direction=3*j).ae(e**2)
    assert diff(lambda x: x**2, 3.0, method='quad').ae(6)
    assert diff(lambda x: 3+x**5, 3.0, n=2, method='quad').ae(540)
    assert diff(lambda x: 3+x**5, 3.0, n=2, method='step').ae(540)
    assert diffun(sin)(2).ae(cos(2))
    assert diffun(sin, n=2)(2).ae(-sin(2))

def test_diffs():
    assert [chop(d) for d in diffs(sin, 0, 1)] == [0, 1]
    assert [chop(d) for d in diffs(sin, 0, 1, method='quad')] == [0, 1]
    assert [chop(d) for d in diffs(sin, 0, 2)] == [0, 1, 0]
    assert [chop(d) for d in diffs(sin, 0, 2, method='quad')] == [0, 1, 0]

def test_diff_relative():
    # relative=True must make the step proportional to mag(x)
    for x in [1e-30, 1e-5, 1, 1e10, 1e30]:
        assert diff(log, x, relative=True).ae(1/x)
    # and must not fail at x = 0
    assert diff(exp, 0, relative=True).ae(1)

def test_taylor():
    # Easy to test since the coefficients are exact in floating-point
    assert taylor(sqrt, 1, 4) == [1, 0.5, -0.125, 0.0625, -0.0390625]

def test_diff_partial():
    x,y,z = xyz = 2,3,7
    f = lambda x,y,z: 3*x**2 * (y+2)**3 * z**5
    assert diff(f, xyz, (0,0,0)).ae(25210500)
    assert diff(f, xyz, (0,0,1)).ae(18007500)
    assert diff(f, xyz, (0,0,2)).ae(10290000)
    assert diff(f, xyz, (0,1,0)).ae(15126300)
    assert diff(f, xyz, (0,1,1)).ae(10804500)
    assert diff(f, xyz, (0,1,2)).ae(6174000)
    assert diff(f, xyz, (0,2,0)).ae(6050520)
    assert diff(f, xyz, (0,2,1)).ae(4321800)
    assert diff(f, xyz, (0,2,2)).ae(2469600)
    assert diff(f, xyz, (1,0,0)).ae(25210500)
    assert diff(f, xyz, (1,0,1)).ae(18007500)
    assert diff(f, xyz, (1,0,2)).ae(10290000)
    assert diff(f, xyz, (1,1,0)).ae(15126300)
    assert diff(f, xyz, (1,1,1)).ae(10804500)
    assert diff(f, xyz, (1,1,2)).ae(6174000)
    assert diff(f, xyz, (1,2,0)).ae(6050520)
    assert diff(f, xyz, (1,2,1)).ae(4321800)
    assert diff(f, xyz, (1,2,2)).ae(2469600)
    assert diff(f, xyz, (2,0,0)).ae(12605250)
    assert diff(f, xyz, (2,0,1)).ae(9003750)
    assert diff(f, xyz, (2,0,2)).ae(5145000)
    assert diff(f, xyz, (2,1,0)).ae(7563150)
    assert diff(f, xyz, (2,1,1)).ae(5402250)
    assert diff(f, xyz, (2,1,2)).ae(3087000)
    assert diff(f, xyz, (2,2,0)).ae(3025260)
    assert diff(f, xyz, (2,2,1)).ae(2160900)
    assert diff(f, xyz, (2,2,2)).ae(1234800)

def test_issue_474():
    # gh-474: the default method returned 0.0 in the fp context
    assert fp.almosteq(fp.diff(fp.exp, 0), 1, rel_eps=1e-14)

def test_diff_fp():
    tol = {1: 1e-13, 2: 1e-11, 3: 1e-9, 4: 1e-7}
    for x in [-1.0, 0.5, 3.0]:
        for n, eps in tol.items():
            assert fp.almosteq(fp.diff(fp.exp, x, n), math.exp(x), rel_eps=eps)
            assert fp.almosteq(fp.diff(fp.sin, x, n),
                               math.sin(x + n*math.pi/2), rel_eps=eps)
    assert fp.almosteq(fp.diff(lambda x: x**3 - 2*x, 2.0, 2), 12, rel_eps=1e-12)
    assert fp.almosteq(fp.diff(fp.exp, 1+2j), cmath.exp(1+2j), rel_eps=1e-13)
    # one-sided differences, also in a complex direction
    assert fp.almosteq(fp.diff(abs, 0, direction=1), 1, rel_eps=1e-13)
    assert fp.almosteq(fp.diff(abs, 0, direction=-1), -1, rel_eps=1e-13)
    assert fp.almosteq(fp.diff(fp.exp, 2.0, direction=3j), math.exp(2),
                       rel_eps=1e-11)
    # relative=True scales the step with |x|; x = 0 falls back to 1
    assert fp.almosteq(fp.diff(fp.log, 1e6, relative=True), 1e-6, rel_eps=1e-12)
    assert fp.almosteq(fp.diff(fp.log, 1e-6, relative=True), 1e6, rel_eps=1e-12)
    assert fp.almosteq(fp.diff(fp.exp, 0, relative=True), 1, rel_eps=1e-14)
    # removable singularity at 0
    sinc = lambda x: fp.sin(x)/x
    assert fp.almosteq(fp.diff(sinc, 0, 0, singular=True), 1, rel_eps=1e-14)
    assert abs(fp.diff(sinc, 0, singular=True)) < 1e-13
    assert fp.almosteq(fp.diff(sinc, 0, 2, singular=True), -1/3, rel_eps=1e-11)
    # a first derivative takes 8 evaluations, none of them at x
    calls = []
    fp.diff(lambda t: calls.append(t) or fp.exp(t), 0)
    assert len(calls) == 8 and 0.0 not in calls
    # an explicit h keeps the plain (n+1)-point difference
    assert fp.diff(fp.exp, 0, h=1e-5) == (fp.exp(1e-5) - fp.exp(-1e-5))/2e-5
    # partial derivatives, diffs() and taylor()
    assert fp.almosteq(fp.diff(lambda x, y: x*y**2, (1.0, 2.0), (1, 1)), 4,
                       rel_eps=1e-11)
    for d in fp.diffs(fp.exp, 0, 4):
        assert fp.almosteq(d, 1, rel_eps=1e-7)
    expected = [0, 1, 0, -1/6, 0, 1/120]
    for c, c0 in zip(fp.taylor(fp.sin, 0, 5), expected):
        assert abs(c - c0) < 1e-8
    # functions built on diff()
    assert fp.almosteq(fp.findroot(lambda x: x**2 - 2, 1.0, solver='newton'),
                       math.sqrt(2))
    assert fp.multiplicity(lambda x: fp.sin(x) - 1, fp.pi/2) == 2

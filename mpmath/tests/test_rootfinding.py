import pytest

from mpmath import (cos, eps, findroot, fp, inf, iv, jacobian, matrix, mnorm,
                    mp, mpc, mpf, multiplicity, norm, pi, polyval, sin, sqrt,
                    workprec)
from mpmath.calculus.optimization import (Anderson, ANewton, Bisection,
                                          Illinois, MDNewton, MNewton, Muller,
                                          Newton, Pegasus, Ridder, Secant, ModAB, Brent)


def test_findroot():
    # old tests, assuming secant
    assert findroot(lambda x: 4*x-3, mpf(5)).ae(0.75)
    assert findroot(sin, mpf(3)).ae(pi)
    assert findroot(sin, (mpf(3), mpf(3.14))).ae(pi)
    assert findroot(lambda x: x*x+1, mpc(2+2j)).ae(1j)
    # test all solvers with 1 starting point
    f = lambda x: cos(x)
    for solver in [Newton, Secant, MNewton, Muller, ANewton]:
        x = findroot(f, 2., solver=solver)
        assert abs(f(x)) < eps
    # test all solvers with interval of 2 points
    for solver in [Secant, Muller, Bisection, Illinois, Pegasus, Anderson,
                   Ridder, ModAB, Brent]:
        x = findroot(f, (1., 2.), solver=solver)
        assert abs(f(x)) < eps
    # test types
    f = lambda x: (x - 2)**2

    assert isinstance(findroot(f, 1, tol=1e-10), mpf)
    assert isinstance(iv.findroot(f, 1., tol=1e-10), iv.mpf)
    assert isinstance(fp.findroot(f, 1, tol=1e-10), float)
    assert isinstance(fp.findroot(f, 1+0j, tol=1e-10), complex)

    # issue 401
    with pytest.raises(ValueError):
        with workprec(2):
            findroot(lambda x: x**2 - 4456178*x + 60372201703370,
                     mpc(real='5.278e+13', imag='-5.278e+13'))

    # issue 192
    with pytest.raises(ValueError):
        findroot(lambda x: -1, 0)

    # issue 387
    with pytest.raises(ValueError):
        findroot(lambda p: (1 - p)**30 - 1, 0.9)

def test_bisection():
    # issue 273
    assert findroot(lambda x: x**2-1,(0,2),solver='bisect') == 1

    with pytest.raises(ValueError):
        findroot(lambda x: x**2-1, (4, 2), solver='bisect') == 1

    # issue 285
    mp.dps = 240
    sol = -mp.ceil(mp.log(abs(findroot(lambda x: mp.sign(x - 3), (1, 4),
                                       solver='bisect', verify=False,
                                       tol=1e-200) - 3))/mp.log(10))
    assert sol.ae(200)

    # issue 339
    mp.dps = 15
    res = mpf('0.73908513321516064')
    for dps in [100, 200, 300, 1000]:
        with mp.workdps(dps):
            sol = findroot(lambda x: cos(x) - x, [0, 1], solver='bisect')
        assert (+sol).ae(res)

def test_mnewton():
    f = lambda x: polyval([1, 3, 3, 1], x)
    x = findroot(f, -0.9, solver='mnewton')
    assert abs(f(x)) < eps

def test_anewton():
    f = lambda x: (x - 2)**100
    x = findroot(f, 1., solver=ANewton)
    assert abs(f(x)) < eps

def test_muller():
    f = lambda x: (2 + x)**3 + 2
    x = findroot(f, 1., solver=Muller)
    assert abs(f(x)) < eps

def test_ridder():
    f = lambda x: cos(x)/x
    x = findroot(f, (1, 2), solver='ridder')
    assert abs(f(x)) < eps

def test_brent():
    f = lambda x: cos(x)/x
    x = findroot(f, (1, 2), solver='brent')
    assert abs(f(x)) < eps

    with pytest.raises(ValueError, match="expected interval of 2 points"):
        findroot(lambda x: x**2 - 1, (0,), solver='brent')

    with pytest.raises(ValueError, match="Function must have opposite signs"):
        findroot(lambda x: x**2 - 1, (2, 4), solver='brent')

    assert findroot(lambda x: x, (-1, 2), solver='brent') == 0.0

    assert findroot(lambda x: x, (-1, 1), solver='brent') == 0.0

def test_modAB():
    assert findroot(lambda x: x**2 - 1, (0, 2), solver='modAB') == 1

    # test ordering
    assert findroot(lambda x: x**2 - 1, (2, 0), solver='modAB') == 1

    with pytest.raises(ValueError, match="expected interval of 2 points"):
        findroot(lambda x: x**2 - 1, (0,), solver='modAB')

    with pytest.raises(ValueError, match="Function must have opposite signs"):
        findroot(lambda x: x**2 - 1, (2, 4), solver='modAB')

    # test exact zero hit
    assert findroot(lambda x: x, (-1, 1), solver='modAB') == 0.0

    # test bisection to secant switch for a purely linear function
    f_linear = lambda x: 2*x - 4
    assert mp.almosteq(findroot(f_linear, (0, 5), solver='modAB'), 2.0)

    f_convex = lambda x: x**10 - 1
    assert mp.almosteq(findroot(f_convex, (0.1, 2.0), solver='modAB'), 1.0)

    f_concave = lambda x: 1 - x**10
    assert mp.almosteq(findroot(f_concave, (2.0, 0.1), solver='modAB'), 1.0)

    f_cubic_inflection = lambda x: x**3 - 3*x + 3
    root = findroot(f_cubic_inflection, (-3, 2), solver='modAB')
    assert abs(f_cubic_inflection(root)) < eps

    # test reset to Bisection if the interval width exceeds the threshold
    f_step = lambda x: mp.sin(x) if x > 1 else x - 1
    assert mp.almosteq(findroot(f_step, (0.4, 3.0), solver='modAB'), 1.0)

def test_multiplicity():
    for i in range(1, 5):
        assert multiplicity(lambda x: (x - 1)**i, 1) == i
    assert multiplicity(lambda x: x**2, 1) == 0

def test_multidimensional(capsys):
    def f(*x):
        return [3*x[0]**2-2*x[1]**2-1, x[0]**2-2*x[0]+x[1]**2+2*x[1]-8]
    assert mnorm(jacobian(f, (1,-2)) - matrix([[6,8],[0,-2]]),1) < 1.e-7
    for x, error in MDNewton(mp, f, (1,-2), verbose=0,
                             norm=lambda x: norm(x, inf)):
        pass
    assert norm(f(*x), 2) < 1e-14
    for x, error in MDNewton(mp, f, (1,-2), verbose=1,
                             norm=lambda x: norm(x, inf)):
        pass
    assert norm(f(*x), 2) < 1e-14
    captured = capsys.readouterr()
    assert captured.out.find("canceled, won't get more exact") >= 0
    # The Chinese mathematician Zhu Shijie was the very first to solve this
    # nonlinear system 700 years ago
    f1 = lambda x, y: -x + 2*y
    f2 = lambda x, y: (x**2 + x*(y**2 - 2) - 4*y)  /  (x + 4)
    f3 = lambda x, y: sqrt(x**2 + y**2)
    def f(x, y):
        f1x = f1(x, y)
        return (f2(x, y) - f1x, f3(x, y) - f1x)
    x = findroot(f, (10, 10))
    assert [round(i) for i in x] == [3, 4]

def test_trivial():
    assert findroot(lambda x: 0, 1) == 1
    assert findroot(lambda x: x, 0) == 0
    #assert findroot(lambda x, y: x + y, (1, -1)) == (1, -1)

def test_issue_869():
    f = [lambda x: sqrt(x) + 1]
    pytest.raises(mp.ComplexResult, lambda: findroot(f, [-1]))

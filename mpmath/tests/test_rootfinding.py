import pytest

from mpmath import (cos, eps, findroot, fp, inf, iv, jacobian, matrix, mnorm,
                    mp, mpc, mpf, multiplicity, norm, pi, polyval, sin, sqrt,
                    workprec)
from mpmath.calculus.optimization import (Anderson, ANewton, Bisection,
                                          Illinois, MDNewton, MNewton, Muller,
                                          Newton, Pegasus, Ridder, Secant, Brent)


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
                   Ridder,Brent]:
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

def test_mnewton():
    f = lambda x: polyval([1,3,3,1],x,asc=True)
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

def test_multiplicity():
    for i in range(1, 5):
        assert multiplicity(lambda x: (x - 1)**i, 1) == i
    assert multiplicity(lambda x: x**2, 1) == 0

def test_brent():
    f = lambda x: x**2-2
    x = findroot(f,(1,2),solver='brent')
    assert abs(f(x))<eps

def test_brent_invalid_interval():
    f= lambda x: x**2+1
    with pytest.raises(ValueError):
        findroot(f,(0,1),solver="brent")

def test_brent_swapping_branch():
    f= lambda x: (x-1.5)*(x-2.5)
    x= findroot(f,(1.4,2),solver="brent")
    assert abs(f(x))<1e-12


def test_brent_interpolation_branch():
    f= lambda x: x**2 - 2
    x= findroot(f,(1,2),solver="brent")
    assert abs(x-sqrt(2))<1e-12


def test_brent_fallback_to_bisection():
    f=lambda x:x-3
    x=findroot(f,(0,5),solver="brent")
    assert abs(f(x))<1e-12


def test_brent_bracket_update():
    f= lambda x:x**3-x-2
    x= findroot(f,(1,2),solver="brent")
    assert abs(f(x))<1e-12

def test_brent_bisection_fallback():
    f= lambda x:(x-1)*(x-2)
    x= findroot(f,(1.5,3.0),solver="brent")
    assert abs(f(x))<1e-12


def test_brent_small_step_safeguard():
    f= lambda x:(x-1e-14)
    x= findroot(f,(0.0,1.0),solver="brent")
    assert abs(f(x))<1e-12


def test_brent_maxsteps_reached():
    f= lambda x:x-0.5
    solver= Brent(mp,f,(0.0,1.0),tol=mp.eps**5)
    last_x= None
    for i,(x,err) in enumerate(solver):
        last_x=x
        if i>solver.maxsteps:
            break
    assert abs(f(last_x))<1e-2

def test_brent_raises_without_interval():
    f = lambda x: x**2 - 2
    # Providing only one point instead of an interval
    with pytest.raises(ValueError, match="Brent's method requires an interval"):
        findroot(f, 1.0, solver="brent")

def test_brent_triggers_b_adjustment():
    # Root very close to left endpoint
    f = lambda x: x - 1e-12
    # Bracket [0, 1]: f(0) = -1e-12 (negative), f(1) = 1 - 1e-12 (positive)
    root = findroot(f, (0, 1), solver="brent")
    assert abs(root - 1e-12) < 1e-8

def test_brent_requires_interval():
    f = lambda x: x**2 - 2
    with pytest.raises(ValueError, match="requires an interval"):
        findroot(f, 1, solver="brent")   # single start, not (a, b)

def test_brent_requires_interval():
    f = lambda x: x**2 - 2
    # Brent must be given an interval (a, b), so this raises
    with pytest.raises(ValueError, match="requires an interval"):
        findroot(f, 1, solver="brent")   # wrong usage: only one starting point

def test_brent_triggers_b_adjustment():
    # Root very close to the right endpoint, forces b-adjustment line
    eps = 1e-15
    f = lambda x: x - (1 - eps)
    root = findroot(f, (0, 1), solver="brent")
    assert abs(root - (1 - eps)) < 1e-12
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
    assert [int(round(i)) for i in x] == [3, 4]

def test_trivial():
    assert findroot(lambda x: 0, 1) == 1
    assert findroot(lambda x: x, 0) == 0
    #assert findroot(lambda x, y: x + y, (1, -1)) == (1, -1)

def test_issue_869():
    f = [lambda x: sqrt(x) + 1]
    pytest.raises(mp.ComplexResult, lambda: findroot(f, [-1]))

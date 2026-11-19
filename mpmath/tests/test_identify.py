from mpmath import e, exp, findpoly, identify, log, mp, pi, pslq, sqrt, zeta


def test_pslq():
    assert pslq([3*pi+4*e/7, pi, e, log(2)]) == [7, -21, -4, 0]
    assert pslq([4.9999999999999991, 1]) == [1, -5]
    assert pslq([2,1]) == [1, -2]

def test_pslq_verbose(capsys):
    assert pslq([-1, pi], tol=0.01, verbose=True) == [22, 7]
    captured = capsys.readouterr()
    assert captured.out.find("PSLQ using prec") >= 0
    assert captured.out.find("FOUND relation at iter") >= 0
    assert captured.out.find("0/100:  Error:") >= 0
    assert pslq([-1, pi], verbose=True) is None
    captured = capsys.readouterr()
    assert captured.out.find("CANCELLING after step") >= 0

def test_identify():
    mp.dps = 20
    assert identify(zeta(4), ['log(2)', 'pi**4']) == '((1/90)*pi**4)'
    mp.dps = 15
    assert identify(exp(5)) == 'exp(5)'
    assert identify(exp(4)) == 'exp(4)'
    assert identify(log(5)) == 'log(5)'
    assert identify(exp(3*pi), ['pi']) == 'exp((3*pi))'
    assert identify(3, full=True) == ['3', '3', '1/(1/3)', 'sqrt(9)',
        '1/sqrt((1/9))', '(sqrt(12)/2)**2', '1/(sqrt(12)/6)**2']
    assert '1/10' in identify(0.1, full=True)
    assert '((2+sqrt(8))/2)' in identify(1+sqrt(2), full=True)
    assert '-(sqrt((8/9)))' in identify(-sqrt(8)/3, full=True)
    assert identify(pi+1, {'a':+pi}) == '(1 + 1*a)'

def test_findpoly_deprecated():
    assert findpoly(1+sqrt(2), 2, asc=False) == [1, -2, -1]

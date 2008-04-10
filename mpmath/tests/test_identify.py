from mpmath import *

def test_pslq():
    mp.dps = 15
    assert pslq([3*pi+4*e/7, pi, e, log(2)]) == [7, -21, -4, 0]
    assert pslq([4.9999999999999991, 1]) == [1, -5]
    assert pslq([2,1]) == [1, -2]

def test_identify():
    mp.dps = 20
    assert identify(zeta(4), ['log(2)', 'pi**4']) == ['((1/90)*pi**4)', '(1/90)*pi**4', '1/90*pi**4']
    mp.dps = 15
    assert identify(exp(3*pi), ['pi']) == ['exp((3*pi))']
    assert identify(exp(5)) == ['exp(5)']
    assert identify(exp(4)) == ['exp(4)']
    assert identify(log(5)) == ['log(5)']

from mpmath import *

def test_special():
    assert inf == inf
    assert inf != -inf
    assert -inf == -inf
    assert inf != nan
    assert nan != nan
    assert nan is nan
    assert --inf == inf
    assert abs(inf) == inf
    assert abs(-inf) == inf
    assert abs(nan) != abs(nan)

    assert inf - inf is nan
    assert inf + (-inf) is nan
    assert -inf - (-inf) is nan

    assert inf + nan is nan
    assert -inf + nan is nan

    assert mpf(2) + inf == inf
    assert 2 + inf == inf
    assert mpf(2) - inf == -inf
    assert 2 - inf == -inf

    assert inf > 3
    assert 3 < inf

    assert 3 > -inf
    assert -inf < 3

    assert inf * 0 is nan
    assert -inf * 0 is nan
    assert inf * 3 == inf
    assert inf * -3 == -inf
    assert -inf * 3 == -inf
    assert -inf * -3 == inf
    assert inf * inf == inf
    assert -inf * -inf == inf

    assert nan / 3 is nan
    assert inf / -3 == -inf
    assert inf / 3 == inf
    assert 3 / inf == 0
    assert -3 / inf == 0
    assert 0 / inf == 0
    assert inf / inf is nan
    assert inf / -inf is nan
    assert inf / nan is nan

    assert mpf('inf') == mpf('+inf') == inf
    assert mpf('-inf') == -inf
    assert mpf('nan') is nan

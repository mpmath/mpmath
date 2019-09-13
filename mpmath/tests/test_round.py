import mpmath
from mpmath import *
import sys

def test_basics():
    assert round(mpf('3.5')) == mpf('4')
    assert round(mpf('3.5'), 1) == mpf('3.5')

def test_default():
    a = mpf('20') / mpf('3')
    assert round(a, 0) == round(a)

def test_lowprec():
    with(workdps(64)):
        a = mpf('2000') / mpf('3')
        assert str(round(a, -1)) in ('670.0', '670')
        assert str(round(a)) in ('667.0', '667')
        assert str(round(a, 1)) == '666.7'
        assert str(round(a, 2)) == '666.67'

def test_highprec():
    if sys.version_info[0] >= 3:
        with(workdps(64)):
            a = mpf('2000') / mpf('3')
            assert str(round(a, 32)) == '666.66666666666666666666666666666667'

def test_types():
    if sys.version_info[0] < 3:
        a = mpf('2000') / mpf('3')
        assert isinstance(round(a, 32), float)
    else:
        with(workdps(64)):
            a = mpf('2000') / mpf('3')
            assert isinstance(round(a, 32), mpf)
            assert isinstance(round(a, 15), mpf)
            assert isinstance(round(a), (int, mpf))  # TODO?
            assert isinstance(round(a, -1), (int, mpf))

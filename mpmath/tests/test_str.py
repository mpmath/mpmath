import math
import random

import hypothesis.strategies as st
from hypothesis import example, given

import mpmath
from mpmath import inf, matrix, mp, mpc, mpf, nstr


A1 = matrix([])
A2 = matrix([[]])
A3 = matrix(2)
A4 = matrix([1, 2, 3])


def test_nstr():
    m = matrix([[0.75, 0.190940654, -0.0299195971],
                [0.190940654, 0.65625, 0.205663228],
                [-0.0299195971, 0.205663228, 0.64453125e-20]])
    assert nstr(m, 4, min_fixed=-inf) == \
    '''[    0.75  0.1909                    -0.02992]
[  0.1909  0.6562                      0.2057]
[-0.02992  0.2057  0.000000000000000000006445]'''
    assert nstr(m, 4) == \
    '''[    0.75  0.1909   -0.02992]
[  0.1909  0.6562     0.2057]
[-0.02992  0.2057  6.445e-21]'''
    # Check that kwargs works properly for mpc
    assert nstr(mpc(1.23e-4+4.56e-4j)) == '(0.000123 + 0.000456j)'
    assert nstr(mpc(1.23e-4+4.56e-4j), min_fixed=-4) == '(1.23e-4 + 4.56e-4j)'

def test_matrix_repr():
    assert repr(A1) == \
    '''matrix(
[])'''
    assert repr(A2) == \
    '''matrix(
[[]])'''
    assert repr(A3) == \
    '''matrix(
[['0.0', '0.0'],
 ['0.0', '0.0']])'''
    assert repr(A4) == \
    '''matrix(
[['1.0'],
 ['2.0'],
 ['3.0']])'''

def test_matrix_str():
    assert str(A1) == ''
    assert str(A2) == '[]'
    assert str(A3) == \
    '''[0.0  0.0]
[0.0  0.0]'''
    assert str(A4) == \
'''[1.0]
[2.0]
[3.0]'''


@given(st.floats(allow_subnormal=False,
                 allow_nan=False,
                 allow_infinity=False))
@example(1.0)
@example(-10.0)
@example(3.411330784663857e+16)
@example(5.960464477539063e-08)
def test_float_short_repr(f):
    if not f and math.copysign(1, f) == -1:
        return
    s = str(f)
    mp._legacy = False
    m = mpf(f)
    sm = str(m)
    assert s == sm
    assert m == mpf(sm)


def test_short_repr_specials():
    mp._legacy = False
    assert str(mpf(0)) == '0.0'
    assert str(mpf('inf')) == 'inf'
    assert str(mpf('-inf')) == '-inf'
    assert str(mpf('nan')) == 'nan'


def test_short_repr_roundtrip():
    mp._legacy = False
    for dps in [15, 20, 30, 50, 100, 300]:
        with mp.workdps(dps):
            for _ in range(10000):
                f = random.choice([(mpmath.rand()-0.5)*2]*10
                                  + [(mpmath.rand()-0.5)*2*10**5]*5
                                  + [(mpmath.rand()-0.5)*2*10**100]*2)
                s = str(f)
                b = mpf(s)
                assert f == b  # round-trip
                if '.' not in s or len(s) < 2 or s[-2:] == '.0':
                    continue
                # test that short repr is really minimal
                integer, frac = s.split('.')
                frac, *exponent = frac.split('e')
                exponent = 'e' + exponent[0] if exponent else ''
                assert f == mpf(str(integer + '.' + frac + exponent))
                frac = frac[:-random.randint(1,len(frac))]
                assert f != mpf(str(integer + '.' + frac + exponent))

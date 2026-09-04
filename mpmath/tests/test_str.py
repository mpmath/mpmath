import math
import random

import hypothesis.strategies as st
from hypothesis import example, given

from mpmath import inf, matrix, mp, mpc, mpf, nstr, rand


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


@given(st.floats(allow_subnormal=True,
                 allow_nan=False,
                 allow_infinity=False),
       st.sampled_from(list('nfcud')))
@example(x=6.170920920537087e+17, rnd='f')
def test_eval_repr_roundtrip(x, rnd):
    mp.rounding = rnd
    mp.pretty = True
    mx = mp.mpf(x)
    smx = repr(mx)
    assert mx == mp.mpf(smx)


@given(st.floats(allow_subnormal=False,
                 allow_nan=False,
                 allow_infinity=False))
@example(1.0)
@example(-10.0)
@example(3.411330784663857e+16)
@example(5.960464477539063e-08)
@example(562949953421312.2)
def test_float_short_repr(f):
    if not f and math.copysign(1, f) == -1:
        return
    s = str(f)
    m = mpf(f)
    sm = str(m)
    assert s == sm
    assert f"mpf('{s}')" == repr(m)
    assert m == mpf(sm)


@given(st.complex_numbers(allow_subnormal=False,
                          allow_nan=False,
                          allow_infinity=False))
@example(1+0.1j)
def test_complex_short_repr(z):
    mp.pretty = False
    if ((not z.real and math.copysign(1, z.real) == -1)
            or (not z.imag and math.copysign(1, z.imag) == -1)):
        return  # skip negative zero
    s = str(z)
    mz = mpc(z)
    smz = str(mz)
    assert s == smz
    assert f"mpc(real='{mz.real!s}', imag='{mz.imag!s}')" == repr(mz)
    assert mz == mpc(smz)
    mp.pretty = True
    assert smz == repr(mz)


def test_short_repr_specials():
    assert str(mpf(0)) == '0.0'
    assert str(mpf('inf')) == 'inf'
    assert str(mpf('-inf')) == '-inf'
    assert str(mpf('nan')) == 'nan'


def test_short_repr_roundtrip():
    for dps in [15, 20, 30, 50, 100, 300]:
        with mp.workdps(dps):
            for _ in range(1000):
                f = random.choice([(rand()-0.5)*2 for _ in range(10)]
                                  + [(rand()-0.5)*2*10**5 for _ in range(5)]
                                  + [(rand()-0.5)*2/10**5 for _ in range(5)]
                                  + [(rand()-0.5)*2*10**100 for _ in range(2)]
                                  + [(rand()-0.5)*2*10**10000 for _ in range(2)]
                                  + [(rand()-0.5)*2/10**10000 for _ in range(2)])
                s = str(f)
                b = mpf(s)
                assert f == b  # round-trip

                integer, *frac = s.split('.')
                if not frac:
                    continue
                frac = frac[0]
                if len(frac) < 2:
                    continue
                frac, *exponent = frac.split('e')
                exponent = 'e' + exponent[0] if exponent else ''

                # round-trip:
                assert f == mpf(str(integer + '.' + frac + exponent))

                # test that short repr is really minimal
                frac = frac[:-1]
                for d in range(10):
                    frac = frac[:-1] + str(d)
                    assert f != mpf(str(integer + '.' + frac + exponent))

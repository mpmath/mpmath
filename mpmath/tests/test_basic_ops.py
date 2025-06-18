import decimal
import math
import operator
import random

import pytest
from hypothesis import example, given, settings
from hypothesis import strategies as st

import mpmath
from mpmath import (ceil, fadd, fdiv, floor, fmul, fneg, fp, frac, fsub, inf,
                    isinf, isint, isnan, isnormal, isspecial, iv, monitor, mp,
                    mpc, mpf, mpi, nan, ninf, nint, nint_distance, nstr, pi,
                    workprec)
from mpmath.libmp import (MPQ, MPZ, finf, fnan, fninf, fnone, fone, from_float,
                          from_int, from_pickable, from_str, isprime, mpf_add,
                          mpf_mul, mpf_sub, round_down, round_nearest,
                          round_up, to_int, to_man_exp, to_pickable)


def test_type_compare():
    assert mpf(2) == mpc(2,0)
    assert mpf(0) == mpc(0)
    assert mpf(2) != mpc(2, 0.00001)
    assert mpf(2) == 2.0
    assert mpf(2) != 3.0
    assert mpf(2) == 2
    assert mpf(2) != '2.0'
    assert mpc(2) != '2.0'

def test_add():
    assert mpf(2.5) + mpf(3) == 5.5
    assert mpf(2.5) + 3 == 5.5
    assert mpf(2.5) + 3.0 == 5.5
    assert 3 + mpf(2.5) == 5.5
    assert 3.0 + mpf(2.5) == 5.5
    assert (3+0j) + mpf(2.5) == 5.5
    assert mpc(2.5) + mpf(3) == 5.5
    assert mpc(2.5) + 3 == 5.5
    assert mpc(2.5) + 3.0 == 5.5
    assert mpc(2.5) + (3+0j) == 5.5
    assert 3 + mpc(2.5) == 5.5
    assert 3.0 + mpc(2.5) == 5.5
    assert (3+0j) + mpc(2.5) == 5.5

def test_sub():
    assert mpf(2.5) - mpf(3) == -0.5
    assert mpf(2.5) - 3 == -0.5
    assert mpf(2.5) - 3.0 == -0.5
    assert 3 - mpf(2.5) == 0.5
    assert 3.0 - mpf(2.5) == 0.5
    assert (3+0j) - mpf(2.5) == 0.5
    assert mpc(2.5) - mpf(3) == -0.5
    assert mpc(2.5) - 3 == -0.5
    assert mpc(2.5) - 3.0 == -0.5
    assert mpc(2.5) - (3+0j) == -0.5
    assert 3 - mpc(2.5) == 0.5
    assert 3.0 - mpc(2.5) == 0.5
    assert (3+0j) - mpc(2.5) == 0.5

def test_mul():
    assert mpf(2.5) * mpf(3) == 7.5
    assert mpf(2.5) * 3 == 7.5
    assert mpf(2.5) * 3.0 == 7.5
    assert 3 * mpf(2.5) == 7.5
    assert 3.0 * mpf(2.5) == 7.5
    assert (3+0j) * mpf(2.5) == 7.5
    assert mpc(2.5) * mpf(3) == 7.5
    assert mpc(2.5) * 3 == 7.5
    assert mpc(2.5) * 3.0 == 7.5
    assert mpc(2.5) * (3+0j) == 7.5
    assert 3 * mpc(2.5) == 7.5
    assert 3.0 * mpc(2.5) == 7.5
    assert (3+0j) * mpc(2.5) == 7.5

def test_div():
    assert mpf(6) / mpf(3) == 2.0
    assert mpf(6) / 3 == 2.0
    assert mpf(6) / 3.0 == 2.0
    assert 6 / mpf(3) == 2.0
    assert 6.0 / mpf(3) == 2.0
    assert (6+0j) / mpf(3.0) == 2.0
    assert mpc(6) / mpf(3) == 2.0
    assert mpc(6) / 3 == 2.0
    assert mpc(6) / 3.0 == 2.0
    assert mpc(6) / (3+0j) == 2.0
    assert 6 / mpc(3) == 2.0
    assert 6.0 / mpc(3) == 2.0
    assert (6+0j) / mpc(3) == 2.0
    assert 1/mpc(inf, 1) == 0.0
    assert (1+1j)/mpc(2, inf) == 0.0
    assert mpc(inf, 1)**-1 == 0.0

def test_mod():
    assert mpf(3.1) % decimal.Decimal(5.3) == mpf('3.1000000000000001')
    assert mpf(2.53) % inf == mpf(2.53)
    assert mpf(2.53) % ninf == mpf(2.53)

def test_floordiv():
    assert mpf(30.21) // mpf(2.53) == mpf(11)

def test_divmod():
    assert divmod(mpf(30.21), mpf(2.53)) == (mpf(11), mpf(2.380000000000003))

def test_pow():
    assert mpf(6) ** mpf(3) == 216.0
    assert mpf(6) ** 3 == 216.0
    assert mpf(6) ** 3.0 == 216.0
    assert 6 ** mpf(3) == 216.0
    assert 6.0 ** mpf(3) == 216.0
    assert (6+0j) ** mpf(3.0) == 216.0
    assert mpc(6) ** mpf(3) == 216.0
    assert mpc(6) ** 3 == 216.0
    assert mpc(6) ** 3.0 == 216.0
    assert mpc(6) ** (3+0j) == 216.0
    assert 6 ** mpc(3) == 216.0
    assert 6.0 ** mpc(3) == 216.0
    assert (6+0j) ** mpc(3) == 216.0
    assert inf ** mpf(0) == mpf(1)
    assert ninf ** mpf(0) == mpf(1)
    assert nan ** mpf(0) == mpf(1)

def test_mixed_misc():
    assert 1 + mpf(3) == mpf(3) + 1 == 4
    assert 1 - mpf(3) == -(mpf(3) - 1) == -2
    assert 3 * mpf(2) == mpf(2) * 3 == 6
    assert 6 / mpf(2) == mpf(6) / 2 == 3
    assert 1.0 + mpf(3) == mpf(3) + 1.0 == 4
    assert 1.0 - mpf(3) == -(mpf(3) - 1.0) == -2
    assert 3.0 * mpf(2) == mpf(2) * 3.0 == 6
    assert 6.0 / mpf(2) == mpf(6) / 2.0 == 3

def test_add_misc():
    assert mpf(4) + mpf(-70) == -66
    assert mpf(1) + mpf(1.1)/80 == 1 + 1.1/80
    assert mpf((1, 10000000000)) + mpf(3) == mpf((1, 10000000000))
    assert mpf(3) + mpf((1, 10000000000)) == mpf((1, 10000000000))
    assert mpf((1, -10000000000)) + mpf(3) == mpf(3)
    assert mpf(3) + mpf((1, -10000000000)) == mpf(3)
    assert mpf(1) + 1e-15 != 1
    assert mpf(1) + 1e-20 == 1
    assert mpf(1.07e-22) + 0 == mpf(1.07e-22)
    assert mpf(0) + mpf(1.07e-22) == mpf(1.07e-22)

def test_mpf_init():
    a1 = mpf(0.3, prec=20)
    a2 = mpf(0.3, dps=5)
    a3 = mpf(0.3)
    assert a1 == a2
    assert a1 != a3
    assert str(a1) == '0.300000190734863'
    assert str(a3) == '0.3'
    pytest.raises(ValueError, lambda: mpf((1, 2, 3)))
    pytest.raises(ValueError, lambda: mpf(mpi(1, 2)))
    pytest.raises(TypeError, lambda: mpf(object()))
    pytest.raises(TypeError, lambda: mpf(1 + 1j))
    class SomethingReal:
        def _mpmath_(self, prec, rounding):
            return mp.make_mpf(from_str('1.3', prec, rounding))
    class SomethingComplex:
        def _mpmath_(self, prec, rounding):
            return mp.make_mpc((from_str('1.3', prec, rounding), \
                from_str('1.7', prec, rounding)))
    class mympf:
        @property
        def _mpf_(self):
            return mpf(3.5)._mpf_
    assert mpf(SomethingReal(), prec=20) == mpf('1.3', prec=20)
    pytest.raises(TypeError, lambda: mpf(SomethingComplex()))
    assert mpf(mympf()) == mpf(3.5)
    assert mympf() - mpf(0.5) == mpf(3.0)
    assert mpf(decimal.Decimal('1.5')) == mpf('1.5')
    assert mpf(decimal.Decimal('+inf')) == +inf
    assert mpf(decimal.Decimal('-inf')) == -inf
    assert isnan(mpf(decimal.Decimal('nan')))
    assert mpf(decimal.Decimal(1).exp(), dps=5) == mpf('2.7182807922363281', dps=5)
    assert mpf(decimal.Decimal(1).exp(), prec=0) == mpf('2.718281828459045235360287471', prec=93)
    assert mpf('0x1.4ace478p+33') == mpf(11100000000.0)
    assert mpf('0x1.4ace478p+33', base=0) == mpf(11100000000.0)
    assert mpf('1.4ace478p+33', base=16) == mpf(11100000000.0)

    assert mpf(float('+inf')) == +inf
    assert mpf(float('-inf')) == -inf
    assert isnan(mpf(float('nan')))

def test_mpc_init():
    class mympc:
        @property
        def _mpc_(self):
            return (mpf(7)._mpf_, mpf(-1)._mpf_)
    assert mpc(3+1j, 7-1j) == mpc(real='4.0', imag='8.0')
    assert mpc(3+1j, mympc()) == mpc(real='4.0', imag='8.0')
    assert mpc('(1+2j)') == mpc(real='1.0', imag='2.0')

def test_mpf_props():
    a = mpf(0.5)
    assert a.man_exp == (1, -1)
    pytest.raises(ValueError, lambda: inf.man_exp)
    pytest.raises(ValueError, lambda: nan.man_exp)
    assert a.man == 1
    assert a.exp == -1
    assert a.bc == 1

def test_mpf_methods():
    assert mpf(0.5).as_integer_ratio() == (1, 2)
    assert mpf('0.3').as_integer_ratio() == (5404319552844595,
                                             18014398509481984)

def test_mpf_magic():
    assert complex(mpf(0.5)) == complex(0.5)

def test_complex_misc():
    # many more tests needed
    assert 1 + mpc(2) == 3
    assert not mpc(2).ae(2 + 1e-13)
    assert mpc(2+1e-15j).ae(2)

def test_complex_zeros():
    for a in [0,2]:
        for b in [0,3]:
            for c in [0,4]:
                for d in [0,5]:
                    assert mpc(a,b)*mpc(c,d) == complex(a,b)*complex(c,d)

def test_hash():
    for i in range(-256, 256):
        assert hash(mpf(i)) == hash(i)
    assert hash(mpf(0.5)) == hash(0.5)
    assert hash(mpc(2,3)) == hash(2+3j)
    # Check that this doesn't fail
    assert hash(inf)
    hash(nan)
    # Check that overflow doesn't assign equal hashes to large numbers
    assert hash(mpf('1e1000')) != hash('1e10000')
    assert hash(mpc(100,'1e1000')) != hash(mpc(200,'1e1000'))
    assert hash(MPQ(1,3))
    assert hash(MPQ(0,1)) == 0
    assert hash(MPQ(-1,1)) == hash(-1)
    assert hash(MPQ(1,1)) == hash(1)
    assert hash(MPQ(5,1)) == hash(5)
    assert hash(MPQ(1,2)) == hash(0.5)
    assert hash(mpf(1)*2**2000) == hash(2**2000)
    assert hash(mpf(1)/2**2000) == hash(MPQ(1,2**2000))

# Advanced rounding test
def test_add_rounding():
    a = from_float(1e-50)
    assert mpf_sub(mpf_add(fone, a, 53, round_up), fone, 53, round_up) == from_float(2.2204460492503131e-16)
    assert mpf_sub(fone, a, 53, round_up) == fone
    assert mpf_sub(fone, mpf_sub(fone, a, 53, round_down), 53, round_down) == from_float(1.1102230246251565e-16)
    assert mpf_add(fone, a, 53, round_down) == fone

def test_almost_equal():
    assert mpf(1.2).ae(mpf(1.20000001), 1e-7)
    assert not mpf(1.2).ae(mpf(1.20000001), 1e-9)
    assert not mpf(-0.7818314824680298).ae(mpf(-0.774695868667929))
    assert inf.ae(inf)
    assert not inf.ae(-inf)
    assert not mpf(1.2).ae(nan)
    assert not mpf(1.2).ae(inf)
    assert not nan.ae(nan)
    assert not nan.ae(inf)

def test_arithmetic_functions():
    ops = [(operator.add, fadd), (operator.sub, fsub), (operator.mul, fmul),
        (operator.truediv, fdiv)]
    a = mpf(0.27)
    b = mpf(1.13)
    c = mpc(0.51+2.16j)
    d = mpc(1.08-0.99j)
    for x in [a,b,c,d]:
        for y in [a,b,c,d]:
            for op, fop in ops:
                if fop is not fdiv:
                    mp.prec = 200
                    z0 = op(x,y)
                mp.prec = 60
                z1 = op(x,y)
                mp.prec = 53
                z2 = op(x,y)
                assert fop(x, y, prec=60) == z1
                assert fop(x, y) == z2
                if fop is not fdiv:
                    assert fop(x, y, prec=inf) == z0
                    assert fop(x, y, dps=inf) == z0
                    assert fop(x, y, exact=True) == z0
                assert fneg(fneg(z1, exact=True), prec=inf) == z1
                assert fneg(z1) == -(+z1)

def test_exact_integer_arithmetic():
    random.seed(0)
    for prec in [6, 10, 25, 40, 100, 250, 725]:
        for rounding in ['d', 'u', 'f', 'c', 'n']:
            mp.dps = prec
            mp.rounding = rounding
            M = 10**(prec-2)
            M2 = 10**(prec//2-2)
            for i in range(10):
                a = random.randint(-M, M)
                b = random.randint(-M, M)
                assert mpf(a, rounding=rounding) == a
                assert int(mpf(a, rounding=rounding)) == a
                assert int(mpf(str(a), rounding=rounding)) == a
                assert mpf(a) + mpf(b) == a + b
                assert mpf(a) - mpf(b) == a - b
                assert -mpf(a) == -a
                a = random.randint(-M2, M2)
                b = random.randint(-M2, M2)
                assert mpf(a) * mpf(b) == a*b
                assert mpf_mul(from_int(a), from_int(b), mp.prec, rounding) == from_int(a*b)

def test_odd_int_bug():
    assert to_int(from_int(3), round_nearest) == 3

def test_str_1000_digits():
    mp.dps = 1001
    # last digit may be wrong
    assert str(mpf(2)**0.5)[-10:-1] == '9518488472'[:9]
    assert str(pi)[-10:-1] == '2164201989'[:9]

def test_str_10000_digits():
    mp.dps = 10001
    # last digit may be wrong
    assert str(mpf(2)**0.5)[-10:-1] == '5873258351'[:9]
    assert str(pi)[-10:-1] == '5256375678'[:9]

def test_monitor():
    f = lambda x: x**2
    a = []
    b = []
    g = monitor(f, a.append, b.append)
    assert g(3) == 9
    assert g(4) == 16
    assert a[0] == ((3,), {})
    assert b[0] == 9

def test_nint_distance():
    assert nint_distance(mpf(-3)) == (-3, -inf)
    assert nint_distance(mpc(-3)) == (-3, -inf)
    assert nint_distance(mpf(-3.1)) == (-3, -3)
    assert nint_distance(mpf(-3.01)) == (-3, -6)
    assert nint_distance(mpf(-3.001)) == (-3, -9)
    assert nint_distance(mpf(-3.0001)) == (-3, -13)
    assert nint_distance(mpf(-2.9)) == (-3, -3)
    assert nint_distance(mpf(-2.99)) == (-3, -6)
    assert nint_distance(mpf(-2.999)) == (-3, -9)
    assert nint_distance(mpf(-2.9999)) == (-3, -13)
    assert nint_distance(mpc(-3+0.1j)) == (-3, -3)
    assert nint_distance(mpc(-3+0.01j)) == (-3, -6)
    assert nint_distance(mpc(-3.1+0.1j)) == (-3, -3)
    assert nint_distance(mpc(-3.01+0.01j)) == (-3, -6)
    assert nint_distance(mpc(-3.001+0.001j)) == (-3, -9)
    assert nint_distance(mpf(0)) == (0, -inf)
    assert nint_distance(mpf(0.01)) == (0, -6)
    assert nint_distance(mpf('1e-100')) == (0, -332)
    pytest.raises(ValueError, lambda: nint_distance(mpc(1, inf)))
    pytest.raises(ValueError, lambda: nint_distance(mpc(inf, 1)))

def test_floor_ceil_nint_frac():
    for n in range(-10,10):
        assert floor(n) == n
        assert floor(n+0.5) == n
        assert ceil(n) == n
        assert ceil(n+0.5) == n+1
        assert nint(n) == n
        # nint rounds to even
        if n % 2 == 1:
            assert nint(n+0.5) == n+1
        else:
            assert nint(n+0.5) == n
    assert floor(inf) == inf
    assert floor(ninf) == ninf
    assert isnan(floor(nan))
    assert ceil(inf) == inf
    assert ceil(ninf) == ninf
    assert isnan(ceil(nan))
    assert nint(inf) == inf
    assert nint(ninf) == ninf
    assert isnan(nint(nan))
    assert floor(0.1) == 0
    assert floor(0.9) == 0
    assert floor(-0.1) == -1
    assert floor(-0.9) == -1
    assert floor(10000000000.1) == 10000000000
    assert floor(10000000000.9) == 10000000000
    assert floor(-10000000000.1) == -10000000000-1
    assert floor(-10000000000.9) == -10000000000-1
    assert floor(1e-100) == 0
    assert floor(-1e-100) == -1
    assert floor(1e100) == 1e100
    assert floor(-1e100) == -1e100
    assert ceil(0.1) == 1
    assert ceil(0.9) == 1
    assert ceil(-0.1) == 0
    assert ceil(-0.9) == 0
    assert ceil(10000000000.1) == 10000000000+1
    assert ceil(10000000000.9) == 10000000000+1
    assert ceil(-10000000000.1) == -10000000000
    assert ceil(-10000000000.9) == -10000000000
    assert ceil(1e-100) == 1
    assert ceil(-1e-100) == 0
    assert ceil(1e100) == 1e100
    assert ceil(-1e100) == -1e100
    assert nint(0.1) == 0
    assert nint(0.9) == 1
    assert nint(-0.1) == 0
    assert nint(-0.9) == -1
    assert nint(10000000000.1) == 10000000000
    assert nint(10000000000.9) == 10000000000+1
    assert nint(-10000000000.1) == -10000000000
    assert nint(-10000000000.9) == -10000000000-1
    assert nint(1e-100) == 0
    assert nint(-1e-100) == 0
    assert nint(1e100) == 1e100
    assert nint(-1e100) == -1e100
    assert floor(3.2+4.6j) == 3+4j
    assert ceil(3.2+4.6j) == 4+5j
    assert nint(3.2+4.6j) == 3+5j
    for n in range(-10,10):
        assert frac(n) == 0
    assert frac(0.25) == 0.25
    assert frac(1.25) == 0.25
    assert frac(2.25) == 0.25
    assert frac(-0.25) == 0.75
    assert frac(-1.25) == 0.75
    assert frac(-2.25) == 0.75
    assert frac('1e100000000000000') == 0
    u = mpf('1e-100000000000000')
    assert frac(u) == u
    assert frac(-u) == 1  # rounding!
    u = mpf('1e-400')
    assert frac(-u, prec=0) == fsub(1, u, exact=True)
    assert frac(3.25+4.75j) == 0.25+0.75j

def test_isnan_etc():
    assert isnan(nan) is True
    assert isnan(3) is False
    assert isnan(mpf(3)) is False
    assert isnan(inf) is False
    assert isnan(mpc(2, nan)) is True
    assert isnan(mpc(2, nan)) is True
    assert isnan(mpc(nan, nan)) is True
    assert isnan(mpc(2, 2)) is False
    assert isnan(mpc(nan, inf)) is True
    assert isnan(mpc(inf, inf)) is False
    assert isnan(MPQ(3, 2)) is False
    assert isnan(MPQ(0, 1)) is False
    assert isinf(inf) is True
    assert isinf(-inf) is True
    assert isinf(3) is False
    assert isinf(nan) is False
    assert isinf(3 + 4j) is False
    assert isinf(mpc(inf)) is True
    assert isinf(mpc(3, inf)) is True
    assert isinf(mpc(inf, 3)) is True
    assert isinf(mpc(inf, inf)) is True
    assert isinf(mpc(nan, inf)) is True
    assert isinf(mpc(inf, nan)) is True
    assert isinf(mpc(nan, nan)) is False
    assert isinf(MPQ(3, 2)) is False
    assert isinf(MPQ(0, 1)) is False
    pytest.raises(TypeError, lambda: isinf(object()))
    assert isspecial(3) is False
    assert isspecial(3.5) is False
    assert isspecial(mpf(3.5)) is False
    assert isspecial(0) is True
    assert isspecial(mpf(0)) is True
    assert isspecial(0.0) is True
    assert isspecial(inf) is True
    assert isspecial(-inf) is True
    assert isspecial(nan) is True
    assert isspecial(float(inf)) is True
    assert isspecial(mpc(0, 0)) is True
    assert isspecial(mpc(3, 0)) is False
    assert isspecial(mpc(0, 3)) is False
    assert isspecial(mpc(3, 3)) is False
    assert isspecial(mpc(0, nan)) is True
    assert isspecial(mpc(0, inf)) is True
    assert isspecial(mpc(3, nan)) is True
    assert isspecial(mpc(3, inf)) is True
    assert isspecial(mpc(3, -inf)) is True
    assert isspecial(mpc(nan, 0)) is True
    assert isspecial(mpc(inf, 0)) is True
    assert isspecial(mpc(nan, 3)) is True
    assert isspecial(mpc(inf, 3)) is True
    assert isspecial(mpc(inf, nan)) is True
    assert isspecial(mpc(nan, inf)) is True
    assert isspecial(mpc(nan, nan)) is True
    assert isspecial(mpc(inf, inf)) is True
    assert isspecial(MPQ(3, 2)) is False
    assert isspecial(MPQ(0, 1)) is True
    pytest.raises(TypeError, lambda: isspecial(object()))
    assert isspecial(5e-324) is False  # issue 946
    assert fp.isspecial(5e-324) is False
    assert fp.isspecial(0.0) is True
    assert fp.isspecial(-0.0) is True
    assert isint(3) is True
    assert isint(0) is True
    assert isint(int(3)) is True
    assert isint(int(0)) is True
    assert isint(mpf(3)) is True
    assert isint(mpf(0)) is True
    assert isint(mpf(-3)) is True
    assert isint(mpf(3.2)) is False
    assert isint(3.2) is False
    assert isint(nan) is False
    assert isint(inf) is False
    assert isint(-inf) is False
    assert isint(mpc(0)) is True
    assert isint(mpc(3)) is True
    assert isint(mpc(3.2)) is False
    assert isint(mpc(3, inf)) is False
    assert isint(mpc(inf)) is False
    assert isint(mpc(3, 2)) is False
    assert isint(mpc(0, 2)) is False
    assert isint(mpc(3, 2), gaussian=True) is True
    assert isint(mpc(3, 0), gaussian=True) is True
    assert isint(mpc(0, 3), gaussian=True) is True
    assert isint(3 + 4j) is False
    assert isint(3 + 4j, gaussian=True) is True
    assert isint(3 + 0j) is True
    assert isint(MPQ(3, 2)) is False
    assert isint(MPQ(3, 9)) is False
    assert isint(MPQ(9, 3)) is True
    assert isint(MPQ(0, 4)) is True
    assert isint(MPQ(1, 1)) is True
    assert isint(MPQ(-1, 1)) is True
    pytest.raises(TypeError, lambda: isint(object()))
    assert mp.isnpint(0) is True
    assert mp.isnpint(1) is False
    assert mp.isnpint(-1) is True
    assert mp.isnpint(-1.1) is False
    assert mp.isnpint(-1.0) is True
    assert mp.isnpint(MPQ(1, 2)) is False
    assert mp.isnpint(MPQ(-1, 2)) is False
    assert mp.isnpint(MPQ(-3, 1)) is True
    assert mp.isnpint(MPQ(0, 1)) is True
    assert mp.isnpint(MPQ(1, 1)) is False
    assert mp.isnpint(0 + 0j) is True
    assert mp.isnpint(-1 + 0j) is True
    assert mp.isnpint(-1.1 + 0j) is False
    assert mp.isnpint(-1 + 0.1j) is False
    assert mp.isnpint(0 + 0.1j) is False
    assert mp.isnpint(inf) is False
    with pytest.deprecated_call():
        for ctx in [mp, fp]:
            assert ctx.isnormal(1) is True
            assert ctx.isnormal(0.0) is False
            assert ctx.isnormal(ctx.mpc(0)) is False
            assert ctx.isnormal(ctx.mpc(0, 1)) is True
            assert ctx.isnormal(ctx.mpc(1, inf)) is False


def test_isprime():
    assert isprime(MPZ(2))
    assert not isprime(MPZ(4))


def test_issue_438():
    assert mpf(finf) == mpf('inf')
    assert mpf(fninf) == mpf('-inf')
    assert mpf(fnan)._mpf_ == mpf('nan')._mpf_


def test_ctx_mag():
    assert mp.mag(MPQ(1, 2)) == 0
    assert mp.mag(MPQ(2)) == 2
    assert mp.mag(MPQ(0)) == mpf('-inf')


def test_ctx_mp_mpnumeric():
    with pytest.deprecated_call():
        from mpmath.ctx_mp import mpnumeric

def test_to_man_exp_deprecation():
    with pytest.deprecated_call():
        to_man_exp(fnone)

def test_rational_deprecation():
    with pytest.deprecated_call():
        assert mpmath.rational.mpq(1, 2) == MPQ(1, 2)
    with pytest.deprecated_call():
        pytest.raises(AttributeError, lambda: mpmath.rational.spam)


def test_math2_deprecation():
    with pytest.deprecated_call():
        assert mpmath.math2.log == mpmath.libfp.log


def test_to_from_pickable():
    x = mpf(1.2)._mpf_
    with pytest.deprecated_call():
        assert to_pickable(x) == x
    with pytest.deprecated_call():
        assert from_pickable(x) == x


def test_rand_precision():
    """
    Test precision of rand()
    """
    def get_remainder(x, bits):
        """
        Return ``(x % 2**-bits) * (2**bits)``.
        If this is nonzero, we know for sure that x was generated with a resolution greater than ``bits``.
        """
        x = x * 2 ** bits
        return x - int(x)
    # Python float (to test the tests)
    random.seed(42)
    x = random.random()
    assert x == 0.6394267984578837, "failed to initialize random() reproducibly"
    assert get_remainder(x, 53) == 0
    assert get_remainder(x, 52) != 0 # Note: this is only true for specific random seeds!
    # fp:
    random.seed(42)
    x = fp.rand()
    assert get_remainder(x, 53) == 0
    assert get_remainder(x, 52) != 0 # Note: this is only true for specific random seeds!
    # mp:
    with workprec(123):
        random.seed(43)
        x = mp.rand()
        assert get_remainder(x, 123) == 0
        assert get_remainder(x, 122) != 0  # Note: this is only true for specific random seeds!
    # iv:
    oldprec = iv.prec # REMOVE ME LATER - workaround for the bug that workprec doesn't work for iv
    iv.prec=123 # REMOVE ME LATER - workaround for the bug that workprec doesn't work for iv
    with workprec(123):
        random.seed(43)
        x = iv.rand()
        assert get_remainder(x, 123) == 0
        assert get_remainder(x, 122) != 0  # Note: this is only true for specific random seeds!
    iv.prec = oldprec # REMOVE ME LATER - workaround for the bug that workprec  doesn't work for iv

def test_issue_260():
    assert mpc(str(mpc(1j))) == mpc(1j)


@settings(max_examples=10000)
@given(st.floats(allow_nan=True,
                 allow_infinity=True,
                 allow_subnormal=True),
       st.integers(min_value=0, max_value=15))
@example(0.5, 0)
@example(-0.5, 0)
@example(1.5, 0)
@example(-1.5, 0)
@example(2.675, 2)
@example(math.inf, 3)
@example(-math.inf, 1)
def test_round_bulk(x, n):
    mp.prec = fp.prec
    m = mpf(x)
    mr = round(m, n)
    xr = round(x, n)
    if isnan(x):
        assert isnan(mr)
        assert isnan(xr)
    else:
        assert float(mr) == xr
        # mp context doesn't support negative zero
        if not xr and math.copysign(1., xr) == -1.:
            return
        assert nstr(mr, n=14, base=16, strip_zeros=False,
                    show_zero_exponent=True, binary_exp=True) == xr.hex()
    try:
        xr = round(x)
    except ValueError:
        pytest.raises(ValueError, lambda: round(m))
    except OverflowError:
        pytest.raises(OverflowError, lambda: round(m))
    else:
        mr = round(m)
        assert type(mr) is int
        assert mr == xr


def test_rounding_prop():
    assert mp.rounding == 'n'
    assert mp.sin(1) == mpf('0x1.aed548f090ceep-1')
    mp.rounding = 'u'
    assert mp.rounding == 'u'
    assert mp.sin(1) == mpf('0x1.aed548f090cefp-1')
    with pytest.raises(ValueError):
        mp.rounding = 'x'


def test_from_man_exp():
    with pytest.raises(TypeError):
        mp.mpf(("!", 1))

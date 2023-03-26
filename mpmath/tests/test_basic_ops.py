import mpmath
from mpmath import *
from mpmath.libmp import *
import random
import sys

try:
    long = long
except NameError:
    long = int

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
    mp.dps = 15
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
    # Check that overflow doesn't assign equal hashes to large numbers
    assert hash(mpf('1e1000')) != hash('1e10000')
    assert hash(mpc(100,'1e1000')) != hash(mpc(200,'1e1000'))
    from mpmath.rational import mpq
    assert hash(mp.mpq(1,3))
    assert hash(mp.mpq(0,1)) == 0
    assert hash(mp.mpq(-1,1)) == hash(-1)
    assert hash(mp.mpq(1,1)) == hash(1)
    assert hash(mp.mpq(5,1)) == hash(5)
    assert hash(mp.mpq(1,2)) == hash(0.5)
    if sys.version_info >= (3, 2):
        assert hash(mpf(1)*2**2000) == hash(2**2000)
        assert hash(mpf(1)/2**2000) == hash(mpq(1,2**2000))

# Advanced rounding test
def test_add_rounding():
    mp.dps = 15
    a = from_float(1e-50)
    assert mpf_sub(mpf_add(fone, a, 53, round_up), fone, 53, round_up) == from_float(2.2204460492503131e-16)
    assert mpf_sub(fone, a, 53, round_up) == fone
    assert mpf_sub(fone, mpf_sub(fone, a, 53, round_down), 53, round_down) == from_float(1.1102230246251565e-16)
    assert mpf_add(fone, a, 53, round_down) == fone

def test_almost_equal():
    assert mpf(1.2).ae(mpf(1.20000001), 1e-7)
    assert not mpf(1.2).ae(mpf(1.20000001), 1e-9)
    assert not mpf(-0.7818314824680298).ae(mpf(-0.774695868667929))

def test_arithmetic_functions():
    import operator
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
    mp.dps = 15

def test_exact_integer_arithmetic():
    # XXX: re-fix this so that all operations are tested with all rounding modes
    random.seed(0)
    for prec in [6, 10, 25, 40, 100, 250, 725]:
        for rounding in ['d', 'u', 'f', 'c', 'n']:
            mp.dps = prec
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
    mp.dps = 15

def test_odd_int_bug():
    assert to_int(from_int(3), round_nearest) == 3

def test_str_1000_digits():
    mp.dps = 1001
    # last digit may be wrong
    assert str(mpf(2)**0.5)[-10:-1] == '9518488472'[:9]
    assert str(pi)[-10:-1] == '2164201989'[:9]
    mp.dps = 15

def test_str_10000_digits():
    mp.dps = 10001
    # last digit may be wrong
    assert str(mpf(2)**0.5)[-10:-1] == '5873258351'[:9]
    assert str(pi)[-10:-1] == '5256375678'[:9]
    mp.dps = 15

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

def test_floor_ceil_nint_frac():
    mp.dps = 15
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
    from mpmath.rational import mpq
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
    assert isnan(mpq((3, 2))) is False
    assert isnan(mpq((0, 1))) is False
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
    assert isinf(mpq((3, 2))) is False
    assert isinf(mpq((0, 1))) is False
    assert isnormal(3) is True
    assert isnormal(3.5) is True
    assert isnormal(mpf(3.5)) is True
    assert isnormal(0) is False
    assert isnormal(mpf(0)) is False
    assert isnormal(0.0) is False
    assert isnormal(inf) is False
    assert isnormal(-inf) is False
    assert isnormal(nan) is False
    assert isnormal(float(inf)) is False
    assert isnormal(mpc(0, 0)) is False
    assert isnormal(mpc(3, 0)) is True
    assert isnormal(mpc(0, 3)) is True
    assert isnormal(mpc(3, 3)) is True
    assert isnormal(mpc(0, nan)) is False
    assert isnormal(mpc(0, inf)) is False
    assert isnormal(mpc(3, nan)) is False
    assert isnormal(mpc(3, inf)) is False
    assert isnormal(mpc(3, -inf)) is False
    assert isnormal(mpc(nan, 0)) is False
    assert isnormal(mpc(inf, 0)) is False
    assert isnormal(mpc(nan, 3)) is False
    assert isnormal(mpc(inf, 3)) is False
    assert isnormal(mpc(inf, nan)) is False
    assert isnormal(mpc(nan, inf)) is False
    assert isnormal(mpc(nan, nan)) is False
    assert isnormal(mpc(inf, inf)) is False
    assert isnormal(mpq((3, 2))) is True
    assert isnormal(mpq((0, 1))) is False
    assert isint(3) is True
    assert isint(0) is True
    assert isint(long(3)) is True
    assert isint(long(0)) is True
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
    assert isint(mpq((3, 2))) is False
    assert isint(mpq((3, 9))) is False
    assert isint(mpq((9, 3))) is True
    assert isint(mpq((0, 4))) is True
    assert isint(mpq((1, 1))) is True
    assert isint(mpq((-1, 1))) is True
    assert mp.isnpint(0) is True
    assert mp.isnpint(1) is False
    assert mp.isnpint(-1) is True
    assert mp.isnpint(-1.1) is False
    assert mp.isnpint(-1.0) is True
    assert mp.isnpint(mp.mpq(1, 2)) is False
    assert mp.isnpint(mp.mpq(-1, 2)) is False
    assert mp.isnpint(mp.mpq(-3, 1)) is True
    assert mp.isnpint(mp.mpq(0, 1)) is True
    assert mp.isnpint(mp.mpq(1, 1)) is False
    assert mp.isnpint(0 + 0j) is True
    assert mp.isnpint(-1 + 0j) is True
    assert mp.isnpint(-1.1 + 0j) is False
    assert mp.isnpint(-1 + 0.1j) is False
    assert mp.isnpint(0 + 0.1j) is False


def test_issue_438():
    assert mpf(finf) == mpf('inf')
    assert mpf(fninf) == mpf('-inf')
    assert mpf(fnan)._mpf_ == mpf('nan')._mpf_

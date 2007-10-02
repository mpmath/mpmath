from mpmath.lib import *
from mpmath import *
import random
import time
import math
import cmath

#----------------------------------------------------------------------------
# Low-level tests
#

def test_bitcount():
    assert bitcount(0) == 0
    assert bitcount(1) == 1
    assert bitcount(7) == 3
    assert bitcount(8) == 4
    assert bitcount(2**100) == 101
    assert bitcount(2**100-1) == 100
    assert bitcount(-(2**100)) == 101
    assert bitcount(-(2**100-1)) == 100

def test_trailing_zeros():
    assert trailing_zeros(0) == 0
    assert trailing_zeros(1) == 0
    assert trailing_zeros(2) == 1
    assert trailing_zeros(7) == 0
    assert trailing_zeros(8) == 3
    assert trailing_zeros(2**100) == 100
    assert trailing_zeros(2**100-1) == 0

def test_round_down():
    assert normalize(0, -4, 4, ROUND_DOWN)[:2] == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_DOWN)[:2] == (15, 0)
    assert normalize(0xf1, -4, 4, ROUND_DOWN)[:2] == (15, 0)
    assert normalize(0xff, -4, 4, ROUND_DOWN)[:2] == (15, 0)
    assert normalize(-0xf0, -4, 4, ROUND_DOWN)[:2] == (-15, 0)
    assert normalize(-0xf1, -4, 4, ROUND_DOWN)[:2] == (-15, 0)
    assert normalize(-0xff, -4, 4, ROUND_DOWN)[:2] == (-15, 0)

def test_round_up():
    assert normalize(0, -4, 4, ROUND_UP)[:2] == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_UP)[:2] == (15, 0)
    assert normalize(0xf1, -4, 4, ROUND_UP)[:2] == (1, 4)
    assert normalize(0xff, -4, 4, ROUND_UP)[:2] == (1, 4)
    assert normalize(-0xf0, -4, 4, ROUND_UP)[:2] == (-15, 0)
    assert normalize(-0xf1, -4, 4, ROUND_UP)[:2] == (-1, 4)
    assert normalize(-0xff, -4, 4, ROUND_UP)[:2] == (-1, 4)

def test_round_floor():
    assert normalize(0, -4, 4, ROUND_FLOOR)[:2] == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_FLOOR)[:2] == (15, 0)
    assert normalize(0xf1, -4, 4, ROUND_FLOOR)[:2] == (15, 0)
    assert normalize(0xff, -4, 4, ROUND_FLOOR)[:2] == (15, 0)
    assert normalize(-0xf0, -4, 4, ROUND_FLOOR)[:2] == (-15, 0)
    assert normalize(-0xf1, -4, 4, ROUND_FLOOR)[:2] == (-1, 4)
    assert normalize(-0xff, -4, 4, ROUND_FLOOR)[:2] == (-1, 4)

def test_round_ceiling():
    assert normalize(0, -4, 4, ROUND_CEILING)[:2] == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_CEILING)[:2] == (15, 0)
    assert normalize(0xf1, -4, 4, ROUND_CEILING)[:2] == (1, 4)
    assert normalize(0xff, -4, 4, ROUND_CEILING)[:2] == (1, 4)
    assert normalize(-0xf0, -4, 4, ROUND_CEILING)[:2] == (-15, 0)
    assert normalize(-0xf1, -4, 4, ROUND_CEILING)[:2] == (-15, 0)
    assert normalize(-0xff, -4, 4, ROUND_CEILING)[:2] == (-15, 0)

def test_round_half_up():
    assert normalize(0, -4, 4, ROUND_HALF_UP)[:2] == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_HALF_UP)[:2] == (15, 0)
    assert normalize(0xf7, -4, 4, ROUND_HALF_UP)[:2] == (15, 0)
    assert normalize(0xf8, -4, 4, ROUND_HALF_UP)[:2] == (1, 4)
    assert normalize(0xf9, -4, 4, ROUND_HALF_UP)[:2] == (1, 4)
    assert normalize(0xff, -4, 4, ROUND_HALF_UP)[:2] == (1, 4)
    assert normalize(-0xf0, -4, 4, ROUND_HALF_UP)[:2] == (-15, 0)
    assert normalize(-0xf7, -4, 4, ROUND_HALF_UP)[:2] == (-15, 0)
    assert normalize(-0xf8, -4, 4, ROUND_HALF_UP)[:2] == (-1, 4)
    assert normalize(-0xf9, -4, 4, ROUND_HALF_UP)[:2] == (-1, 4)
    assert normalize(-0xff, -4, 4, ROUND_HALF_UP)[:2] == (-1, 4)

def test_round_half_down():
    assert normalize(0, -4, 4, ROUND_HALF_DOWN)[:2] == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_HALF_DOWN)[:2] == (15, 0)
    assert normalize(0xf7, -4, 4, ROUND_HALF_DOWN)[:2] == (15, 0)
    assert normalize(0xf8, -4, 4, ROUND_HALF_DOWN)[:2] == (15, 0)
    assert normalize(0xf9, -4, 4, ROUND_HALF_DOWN)[:2] == (1, 4)
    assert normalize(0xff, -4, 4, ROUND_HALF_DOWN)[:2] == (1, 4)
    assert normalize(-0xf0, -4, 4, ROUND_HALF_DOWN)[:2] == (-15, 0)
    assert normalize(-0xf7, -4, 4, ROUND_HALF_DOWN)[:2] == (-15, 0)
    assert normalize(-0xf8, -4, 4, ROUND_HALF_DOWN)[:2] == (-15, 0)
    assert normalize(-0xf9, -4, 4, ROUND_HALF_DOWN)[:2] == (-1, 4)
    assert normalize(-0xff, -4, 4, ROUND_HALF_DOWN)[:2] == (-1, 4)

def test_round_half_even():
    assert normalize(0, -4, 4, ROUND_HALF_EVEN)[:2] == (0, 0)
    assert normalize(0xf0, -4, 4, ROUND_HALF_EVEN)[:2] == (15, 0)
    assert normalize(0xf7, -4, 4, ROUND_HALF_EVEN)[:2] == (15, 0)
    assert normalize(0xf8, -4, 4, ROUND_HALF_EVEN)[:2] == (1, 4)    # 1111.1000 -> 10000.0
    assert normalize(0xf9, -4, 4, ROUND_HALF_EVEN)[:2] == (1, 4)    # 1111.1001 -> 10000.0
    assert normalize(0xe8, -4, 4, ROUND_HALF_EVEN)[:2] == (7, 1)    # 1110.1000 -> 1110.0
    assert normalize(0xe9, -4, 4, ROUND_HALF_EVEN)[:2] == (15, 0)     # 1110.1001 -> 1111.0
    assert normalize(-0xf0, -4, 4, ROUND_HALF_EVEN)[:2] == (-15, 0)
    assert normalize(-0xf7, -4, 4, ROUND_HALF_EVEN)[:2] == (-15, 0)
    assert normalize(-0xf8, -4, 4, ROUND_HALF_EVEN)[:2] == (-1, 4)
    assert normalize(-0xf9, -4, 4, ROUND_HALF_EVEN)[:2] == (-1, 4)
    assert normalize(-0xe8, -4, 4, ROUND_HALF_EVEN)[:2] == (-7, 1)
    assert normalize(-0xe9, -4, 4, ROUND_HALF_EVEN)[:2] == (-15, 0)


def test_rounding_bugs():
    # 1 less than power-of-two cases
    assert normalize(72057594037927935, -56, 53, ROUND_UP) == (1, 0, 1)
    assert normalize(31, 0, 4, ROUND_UP) == (1, 5, 1)
    assert normalize(-31, 0, 4, ROUND_FLOOR) == (-1, 5, 1)
    assert normalize(255, 0, 7, ROUND_UP) == (1, 8, 1)
    assert normalize(-255, 0, 7, ROUND_FLOOR) == (-1, 8, 1)


# Advanced rounding test
def test_add_rounding():
    mpf.dps = 15
    mpf.round_up()
    assert (mpf(1) + 1e-50) - 1 == 2.2204460492503131e-16
    assert mpf(1) - 1e-50 == 1.0
    mpf.round_down()
    assert 1 - (mpf(1) - 1e-50) == 1.1102230246251565e-16
    assert mpf(1) + 1e-50 == 1.0
    mpf.round_half_even()

def test_almost_equal():
    assert mpf(1.2).ae(mpf(1.20000001), 1e-7)
    assert not mpf(1.2).ae(mpf(1.20000001), 1e-9)
    assert not mpf(-0.7818314824680298).ae(mpf(-0.774695868667929))


#----------------------------------------------------------------------------
# Test basic arithmetic
#

def test_add():
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

def test_complex():
    # many more tests needed
    assert 1 + mpc(2) == 3
    assert not mpc(2).ae(2 + 1e-13)
    assert mpc(2+1e-15j).ae(2)

# Results should agree with those of small floats at standard precision
def test_standard_float():
    random.seed(0)
    sqrtfail = 0
    for i in range(5000):
        x = (random.random()*1000 - 1000) * 2.0**random.randint(-20, 20)
        y = (random.random()*1000 - 1000) * 2.0**random.randint(-20, 20)
        assert mpf(x) == x
        assert (mpf(x) < mpf(y)) == (x < y)
        assert (mpf(x) > mpf(y)) == (x > y)
        assert (mpf(x) == mpf(y)) == (x == y)
        assert (mpf(x) != mpf(y)) == (x != y)
        assert (mpf(x) <= mpf(y)) == (x <= y)
        assert (mpf(x) >= mpf(y)) == (x >= y)
        assert mpf(x) == mpf(x)
        assert mpf(x) + mpf(y) == x + y
        assert mpf(x) * mpf(y) == x * y
        assert mpf(x) / mpf(y) == x / y
        assert abs(mpf(x)) == abs(x)
        # this fails roughly 1 time out of 1000. not sure whether mpf
        # or float is rounding the wrong way
        sqrtfail += (abs(mpf(x))**0.5 != abs(x)**0.5)
    assert sqrtfail < 5
    # particular bugs
    assert mpf(4.4408920985006262E-16) < mpf(1.7763568394002505E-15)
    assert mpf(-4.4408920985006262E-16) > mpf(-1.7763568394002505E-15)


def test_mixed_types():
    assert 1 + mpf(3) == mpf(3) + 1 == 4
    assert 1 - mpf(3) == -(mpf(3) - 1) == -2
    assert 3 * mpf(2) == mpf(2) * 3 == 6
    assert 6 / mpf(2) == mpf(6) / 2 == 3
    assert 1.0 + mpf(3) == mpf(3) + 1.0 == 4
    assert 1.0 - mpf(3) == -(mpf(3) - 1.0) == -2
    assert 3.0 * mpf(2) == mpf(2) * 3.0 == 6
    assert 6.0 / mpf(2) == mpf(6) / 2.0 == 3


# Test that integer arithmetic is exact
def test_integers():
    random.seed(0)
    for prec in [6, 10, 25, 40, 100, 250, 725]:
        mpf.dps = prec
        M = 10**(prec-2)
        M2 = 10**(prec//2-2)
        for i in range(10):
            a = random.randint(-M, M)
            b = random.randint(-M, M)
            assert mpf(a) == a
            assert int(mpf(a)) == a
            assert int(mpf(str(a))) == a
            assert mpf(a) + mpf(b) == a + b
            assert mpf(a) - mpf(b) == a - b
            assert -mpf(a) == -a
            a = random.randint(-M2, M2)
            b = random.randint(-M2, M2)
            assert mpf(a) * mpf(b) == a*b
    mpf.dps = 15

def test_exact_sqrts():
    for i in range(20000):
        assert sqrt(mpf(i*i)) == i
    random.seed(1)
    for prec in [100, 300, 1000, 10000]:
        mpf.dps = prec
        for i in range(20):
            A = random.randint(10**(prec//2-2), 10**(prec//2-1))
            assert sqrt(mpf(A*A)) == A
    mpf.dps = 15
    for i in range(100):
        for a in [1, 8, 25, 112307]:
            assert sqrt(mpf((a*a, 2*i))) == mpf((a, i))
            assert sqrt(mpf((a*a, -2*i))) == mpf((a, -i))


#----------------------------------------------------------------------------
# Type comparison
#

def test_type_compare():
    assert mpf(2) == mpc(2,0)
    assert mpf(0) == mpc(0)
    assert mpf(2) != mpc(2, 0.00001)
    assert mpf(2) == 2.0
    assert mpf(2) != 3.0
    assert mpf(2) == 2
    assert mpf(2) != '2.0'
    assert mpc(2) != '2.0'


#----------------------------------------------------------------------------
# Constants and functions
#

tpi = "3.1415926535897932384626433832795028841971693993751058209749445923078\
1640628620899862803482534211706798"

tgamma = "0.5772156649015328606065120900824024310421593359399235988057672348\
84867726777664670936947063291746749516"

tlog2 = "0.69314718055994530941723212145817656807550013436025525412068000949\
3393621969694715605863326996418687542"

tlog10 = "2.3025850929940456840179914546843642076011014886287729760333279009\
6757260967735248023599720508959829834"


def test_constants():
    for prec in [3, 7, 10, 15, 20, 37, 80, 100, 29]:
        mpf.dps = prec
        assert pi == mpf(tpi)
        assert cgamma == mpf(tgamma)
        assert clog2 == mpf(tlog2)
        assert clog10 == mpf(tlog10)
    mpf.dps = 15

def test_str_1000_digits():
    mpf.dps = 1001
    assert str(mpf(2)**0.5)[-10:] == '9518488472'
    assert str(pi)[-10:] == '2164201989'
    mpf.dps = 15

def test_str_10000_digits():
    mpf.dps = 10001
    assert str(mpf(2)**0.5)[-10:] == '5873258352'
    assert str(pi)[-10:] == '5256375679'
    mpf.dps = 15

def test_float_sqrt():
    # These should round identically
    for x in [0, 1e-7, 0.1, 0.5, 1, 2, 3, 4, 5, 0.333, 76.19]:
        assert sqrt(mpf(x)) == float(x)**0.5
    assert sqrt(-1) == 1j
    assert sqrt(-2).ae(cmath.sqrt(-2))
    assert sqrt(-3).ae(cmath.sqrt(-3))
    assert sqrt(-100).ae(cmath.sqrt(-100))
    assert sqrt(1j).ae(cmath.sqrt(1j))
    assert sqrt(-1j).ae(cmath.sqrt(-1j))
    assert sqrt(math.pi + math.e*1j).ae(cmath.sqrt(math.pi + math.e*1j))
    assert sqrt(math.pi - math.e*1j).ae(cmath.sqrt(math.pi - math.e*1j))

def test_hypot():
    assert hypot(0, 0) == 0
    assert hypot(0, 0.33) == mpf(0.33)
    assert hypot(0.33, 0) == mpf(0.33)
    assert hypot(-0.33, 0) == mpf(0.33)
    assert hypot(3, 4) == mpf(5)

def test_exp():
    assert exp(0) == 1
    assert exp(10000).ae(mpf('8.8068182256629215873e4342'))
    assert exp(-10000).ae(mpf('1.1354838653147360985e-4343'))
    assert exp(clog2 * 10).ae(1024)
    assert exp(2+2j).ae(cmath.exp(2+2j))

def test_log():
    assert log(1) == 0
    for x in [0.5, 1.5, 2.0, 3.0, 100, 10**50, 1e-50]:
        assert log(x) == math.log(x)
        assert log(x, x) == 1
    assert log(1024, 2) == 10
    assert log(10**1234, 10) == 1234
    assert log(2+2j).ae(cmath.log(2+2j))

def test_trig_basic():
    for x in (range(100) + range(-100,0)):
        t = x / 4.1
        assert cos(mpf(t)).ae(math.cos(t))
        assert sin(mpf(t)).ae(math.sin(t))
        assert tan(mpf(t)).ae(math.tan(t))
    assert sin(1+1j).ae(cmath.sin(1+1j))
    assert sin(-4-3.6j).ae(cmath.sin(-4-3.6j))
    assert cos(1+1j).ae(cmath.cos(1+1j))
    assert cos(-4-3.6j).ae(cmath.cos(-4-3.6j))

def test_complex_functions():
    for x in (range(10) + range(-10,0)):
        for y in (range(10) + range(-10,0)):
            z = complex(x, y)/4.3 + 0.01j
            assert exp(mpc(z)).ae(cmath.exp(z))
            assert log(mpc(z)).ae(cmath.log(z))
            assert cos(mpc(z)).ae(cmath.cos(z))
            assert sin(mpc(z)).ae(cmath.sin(z))
            assert tan(mpc(z)).ae(cmath.tan(z))

def test_complex_powers():
    random.seed(1)
    for i in range(100):
        x1 = random.uniform(-10, 10)
        y1 = random.uniform(-10, 10)
        x2 = random.uniform(-10, 10)
        y2 = random.uniform(-10, 10)
        z1 = complex(x1, y1)
        z2 = complex(x2, y2)
        assert (mpc(z1)**mpc(z2)).ae(z1**z2, 1e-12)
    assert (e**(-pi*1j)).ae(-1)
    mpf.dps = 50
    assert (e**(-pi*1j)).ae(-1)
    mpf.dps = 15

def test_trig_hard():
    mpf.prec = 150
    a = mpf(10**50)
    mpf.prec = 53
    assert sin(a).ae(-0.7896724934293100827)
    assert cos(a).ae(-0.6135286082336635622)
    assert sin(1e-6).ae(9.999999999998333e-007)
    assert cos(1e-6).ae(0.9999999999995)

def test_atan():
    assert atan(-2.3).ae(math.atan(-2.3))
    assert atan2(1,1).ae(math.atan2(1,1))
    assert atan2(1,-1).ae(math.atan2(1,-1))
    assert atan2(-1,-1).ae(math.atan2(-1,-1))
    assert atan2(-1,1).ae(math.atan2(-1,1))
    assert atan2(-1,0).ae(math.atan2(-1,0))
    assert atan2(1,0).ae(math.atan2(1,0))
    assert atan(1e-50).ae(1e-50, abs_eps=0)
    assert atan(1e50).ae(pi/2)
    assert atan(-1e-50).ae(-1e-50, abs_eps=0)
    assert atan(-1e50).ae(-pi/2)


if __name__ == "__main__":
    globals_ = globals().keys()
    t1 = time.time()
    for f in globals_:
        if f.startswith("test_"):
            print "running", f, "..."
            globals()[f]()
    t2 = time.time()
    print "passed all tests in", (t2-t1), "seconds"

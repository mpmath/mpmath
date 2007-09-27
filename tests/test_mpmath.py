from mpmath.lib import *
from mpmath import *
import random
import time

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
    assert normalize(72057594037927935, -56, 53, ROUND_UP) == (1, 0, 1)


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



#----------------------------------------------------------------------------
# Test basic arithmetic
#

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



def test_hypot():
    assert hypot(0, 0) == 0
    assert hypot(0, 0.33) == mpf(0.33)
    assert hypot(0.33, 0) == mpf(0.33)
    assert hypot(-0.33, 0) == mpf(0.33)
    assert hypot(3, 4) == mpf(5)



if __name__ == "__main__":
    globals_ = globals().keys()
    t1 = time.time()
    for f in globals_:
        if f.startswith("test_"):
            print "running", f, "..."
            globals()[f]()
    t2 = time.time()
    print "passed all tests in", (t2-t1), "seconds"

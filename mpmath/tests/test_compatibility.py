from mpmath import *
from random import seed, randint, random
import math

# Test compatibility with Python floats, which are
# IEEE doubles (53-bit)

N = 5000
seed(1)

# Choosing exponents between roughly -140, 140 ensures that
# the Python floats don't overflow or underflow
xs = [(random()-1) * 10**randint(-140, 140) for x in range(N)]
ys = [(random()-1) * 10**randint(-140, 140) for x in range(N)]

# include some equal values
ys[int(N*0.8):] = xs[int(N*0.8):]

def test_double_compatibility():
    mpf.prec = 53
    mpf.round_default()
    for x, y in zip(xs, ys):
        mpx = mpf(x)
        mpy = mpf(y)
        assert mpf(x) == x
        assert (mpx < mpy) == (x < y)
        assert (mpx > mpy) == (x > y)
        assert (mpx == mpy) == (x == y)
        assert (mpx != mpy) == (x != y)
        assert (mpx <= mpy) == (x <= y)
        assert (mpx >= mpy) == (x >= y)
        assert mpx == mpx
        assert mpx + mpy == x + y
        assert mpx * mpy == x * y
        assert mpx / mpy == x / y
        assert mpx % mpy == x % y
        assert abs(mpx) == abs(x)
        assert mpf(repr(x)) == x
        assert ceil(mpx) == math.ceil(x)
        assert floor(mpx) == math.floor(x)

def test_sqrt():
    # this fails quite often. it appers to be float
    # that rounds the wrong way, not mpf
    fail = 0
    mpf.prec = 53
    mpf.round_default()
    for x in xs:
        x = abs(x)
        mpf.prec = 100
        mp_high = mpf(x)**0.5
        mpf.prec = 53
        mp_low = mpf(x)**0.5
        fp = x**0.5
        assert abs(mp_low-mp_high) <= abs(fp-mp_high)
        fail += mp_low != fp
    assert fail < N/10

def test_bugs():
    # particular bugs
    assert mpf(4.4408920985006262E-16) < mpf(1.7763568394002505E-15)
    assert mpf(-4.4408920985006262E-16) > mpf(-1.7763568394002505E-15)

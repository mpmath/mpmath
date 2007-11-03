from mpmath import *
from random import seed, randint, random

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
        assert abs(mpx) == abs(x)
        assert mpf(repr(x)) == x

def test_sqrt():
    # this fails roughly 1 time out of 1000. it appers to be float
    # that rounds the wrong way, not mpf
    fail = 0
    mpf.prec = 53
    mpf.round_default()
    for x in xs:
        fail += (abs(mpf(x))**0.5 != abs(x)**0.5)
    assert fail < 2*(N/1000.0)

def test_bugs():
    # particular bugs
    assert mpf(4.4408920985006262E-16) < mpf(1.7763568394002505E-15)
    assert mpf(-4.4408920985006262E-16) > mpf(-1.7763568394002505E-15)

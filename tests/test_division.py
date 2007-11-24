from mpmath.lib import *
from mpmath import mpf

from random import randint, choice, seed

all_modes = [ROUND_FLOOR, ROUND_CEILING, ROUND_DOWN, ROUND_UP,
  ROUND_HALF_DOWN, ROUND_HALF_UP, ROUND_HALF_EVEN]

fb = from_bstr
fi = from_int
ff = from_float


def test_div_1_3():
    a = fi(1)
    b = fi(3)
    c = fi(-1)

    # Floor rounds down, ceiling rounds up
    assert fdiv(a, b, 7, ROUND_FLOOR)     == fb('0.01010101')
    assert fdiv(a, b, 7, ROUND_CEILING)   == fb('0.01010110')
    assert fdiv(a, b, 7, ROUND_DOWN)      == fb('0.01010101')
    assert fdiv(a, b, 7, ROUND_UP)        == fb('0.01010110')
    assert fdiv(a, b, 7, ROUND_HALF_DOWN) == fb('0.01010101')
    assert fdiv(a, b, 7, ROUND_HALF_UP)   == fb('0.01010101')
    assert fdiv(a, b, 7, ROUND_HALF_EVEN) == fb('0.01010101')

    # Floor rounds up, ceiling rounds down
    assert fdiv(c, b, 7, ROUND_FLOOR)     == fb('-0.01010110')
    assert fdiv(c, b, 7, ROUND_CEILING)   == fb('-0.01010101')
    assert fdiv(c, b, 7, ROUND_DOWN)      == fb('-0.01010101')
    assert fdiv(c, b, 7, ROUND_UP)        == fb('-0.01010110')
    assert fdiv(c, b, 7, ROUND_HALF_DOWN) == fb('-0.01010101')
    assert fdiv(c, b, 7, ROUND_HALF_UP)   == fb('-0.01010101')
    assert fdiv(c, b, 7, ROUND_HALF_EVEN) == fb('-0.01010101')


def test_div_300():

    q = fi(1000000)
    a = fi(300499999)    # a/q is a little less than a half-integer
    b = fi(300500000)    # b/q exactly a half-integer
    c = fi(300500001)    # c/q is a little more than a half-integer

    # Check nearest integer rounding (prec=9 as 2**8 < 300 < 2**9)

    assert fdiv(a, q, 9, ROUND_DOWN) == fi(300)
    assert fdiv(b, q, 9, ROUND_DOWN) == fi(300)
    assert fdiv(c, q, 9, ROUND_DOWN) == fi(300)
    assert fdiv(a, q, 9, ROUND_UP) == fi(301)
    assert fdiv(b, q, 9, ROUND_UP) == fi(301)
    assert fdiv(c, q, 9, ROUND_UP) == fi(301)
    assert fdiv(a, q, 9, ROUND_HALF_UP) == fi(300)
    assert fdiv(b, q, 9, ROUND_HALF_UP) == fi(301)
    assert fdiv(c, q, 9, ROUND_HALF_UP) == fi(301)
    assert fdiv(a, q, 9, ROUND_HALF_DOWN) == fi(300)
    assert fdiv(b, q, 9, ROUND_HALF_DOWN) == fi(300)
    assert fdiv(c, q, 9, ROUND_HALF_DOWN) == fi(301)

    # Nearest even integer is down
    assert fdiv(a, q, 9, ROUND_HALF_EVEN) == fi(300)
    assert fdiv(b, q, 9, ROUND_HALF_EVEN) == fi(300)
    assert fdiv(c, q, 9, ROUND_HALF_EVEN) == fi(301)

    # Nearest even integer is up
    a = fi(301499999)
    b = fi(301500000)
    c = fi(301500001)
    assert fdiv(a, q, 9, ROUND_HALF_EVEN) == fi(301)
    assert fdiv(b, q, 9, ROUND_HALF_EVEN) == fi(302)
    assert fdiv(c, q, 9, ROUND_HALF_EVEN) == fi(302)


def test_tight_integer_division():
    # Test that integer division at tightest possible precision is exact
    N = 100
    seed(1)
    for i in range(N):
        a = choice([1, -1]) * randint(1, 1<<randint(10, 100))
        b = choice([1, -1]) * randint(1, 1<<randint(10, 100))
        p = a * b
        width = bitcount(b) - trailing_zeros(b)
        a = fi(a); b = fi(b); p = fi(p)
        for mode in all_modes:
            assert fdiv(p, a, width, mode) == b


def test_epsilon_rounding():
    # Verify that fdiv uses infinite precision; this result will
    # appear to be exactly 0.101 to a near-sighted algorithm

    a = fb('0.101' + ('0'*200) + '1')
    b = fb('1.10101')
    c = fmul(a, b, 250, ROUND_FLOOR) # Exact
    assert fdiv(c, b, bitcount(a[0]), ROUND_FLOOR) == a # Exact

    # Not exactly half, so both must round up
    assert fdiv(c, b, 2, ROUND_HALF_DOWN) == fb('0.11')
    assert fdiv(c, b, 2, ROUND_HALF_UP) == fb('0.11')

    assert fdiv(c, b, 2, ROUND_DOWN) == fb('0.10')
    assert fdiv(c, b, 3, ROUND_DOWN) == fb('0.101')
    assert fdiv(c, b, 2, ROUND_UP) == fb('0.11')
    assert fdiv(c, b, 3, ROUND_UP) == fb('0.110')
    assert fdiv(c, b, 2, ROUND_FLOOR) == fb('0.10')
    assert fdiv(c, b, 3, ROUND_FLOOR) == fb('0.101')
    assert fdiv(c, b, 2, ROUND_CEILING) == fb('0.11')
    assert fdiv(c, b, 3, ROUND_CEILING) == fb('0.110')

    # The same for negative numbers
    a = fb('-0.101' + ('0'*200) + '1')
    b = fb('1.10101')
    c = fmul(a, b, 250, ROUND_FLOOR)
    assert fdiv(c, b, bitcount(a[0]), ROUND_FLOOR) == a

    # Both must round away from zero
    assert fdiv(c, b, 2, ROUND_HALF_DOWN) == fb('-0.11')
    assert fdiv(c, b, 2, ROUND_HALF_UP) == fb('-0.11')

    assert fdiv(c, b, 2, ROUND_DOWN) == fb('-0.10')
    assert fdiv(c, b, 3, ROUND_UP) == fb('-0.110')

    # Floor goes up, ceiling goes down
    assert fdiv(c, b, 2, ROUND_FLOOR) == fb('-0.11')
    assert fdiv(c, b, 3, ROUND_FLOOR) == fb('-0.110')
    assert fdiv(c, b, 2, ROUND_CEILING) == fb('-0.10')
    assert fdiv(c, b, 3, ROUND_CEILING) == fb('-0.101')


def test_mod():
    assert mpf(234) % 1 == 0
    assert mpf(-3) % 256 == 253
    assert mpf(0.25) % 23490.5 == 0.25
    assert mpf(0.25) % -23490.5 == -23490.25
    assert mpf(-0.25) % 23490.5 == 23490.25
    assert mpf(-0.25) % -23490.5 == -0.25
    # Check that these cases are handled efficiently
    assert mpf('1e10000000000') % 1 == 0
    assert mpf('1.23e-1000000000') % 1 == mpf('1.23e-1000000000')
    # test __rmod__
    assert 3 % mpf('1.75') == 1.25

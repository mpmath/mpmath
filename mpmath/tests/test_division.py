from random import choice, randint, seed

from mpmath import mpf
from mpmath.libmp import (from_int, from_str, mpf_div, mpf_mul, mpf_rdiv_int,
                          round_ceiling, round_down, round_floor,
                          round_nearest, round_up, trailing)


def test_div_1_3():
    a = from_int(1)
    b = from_int(3)
    c = from_int(-1)

    # floor rounds down, ceiling rounds up
    assert mpf_div(a, b, 7, round_floor)   == from_str('0.01010101', base=2)
    assert mpf_div(a, b, 7, round_ceiling) == from_str('0.01010110', base=2)
    assert mpf_div(a, b, 7, round_down)    == from_str('0.01010101', base=2)
    assert mpf_div(a, b, 7, round_up)      == from_str('0.01010110', base=2)
    assert mpf_div(a, b, 7, round_nearest) == from_str('0.01010101', base=2)

    # floor rounds up, ceiling rounds down
    assert mpf_div(c, b, 7, round_floor)   == from_str('-0.01010110', base=2)
    assert mpf_div(c, b, 7, round_ceiling) == from_str('-0.01010101', base=2)
    assert mpf_div(c, b, 7, round_down)    == from_str('-0.01010101', base=2)
    assert mpf_div(c, b, 7, round_up)      == from_str('-0.01010110', base=2)
    assert mpf_div(c, b, 7, round_nearest) == from_str('-0.01010101', base=2)

def test_mpf_divi_1_3():
    a = 1
    b = from_int(3)
    c = -1
    assert mpf_rdiv_int(a, b, 7, round_floor)   == from_str('0.01010101', base=2)
    assert mpf_rdiv_int(a, b, 7, round_ceiling) == from_str('0.01010110', base=2)
    assert mpf_rdiv_int(a, b, 7, round_down)    == from_str('0.01010101', base=2)
    assert mpf_rdiv_int(a, b, 7, round_up)      == from_str('0.01010110', base=2)
    assert mpf_rdiv_int(a, b, 7, round_nearest) == from_str('0.01010101', base=2)
    assert mpf_rdiv_int(c, b, 7, round_floor)   == from_str('-0.01010110', base=2)
    assert mpf_rdiv_int(c, b, 7, round_ceiling) == from_str('-0.01010101', base=2)
    assert mpf_rdiv_int(c, b, 7, round_down)    == from_str('-0.01010101', base=2)
    assert mpf_rdiv_int(c, b, 7, round_up)      == from_str('-0.01010110', base=2)
    assert mpf_rdiv_int(c, b, 7, round_nearest) == from_str('-0.01010101', base=2)


def test_div_300():

    q = from_int(1000000)
    a = from_int(300499999)    # a/q is a little less than a half-integer
    b = from_int(300500000)    # b/q exactly a half-integer
    c = from_int(300500001)    # c/q is a little more than a half-integer

    # Check nearest integer rounding (prec=9 as 2**8 < 300 < 2**9)

    assert mpf_div(a, q, 9, round_down) == from_int(300)
    assert mpf_div(b, q, 9, round_down) == from_int(300)
    assert mpf_div(c, q, 9, round_down) == from_int(300)
    assert mpf_div(a, q, 9, round_up) == from_int(301)
    assert mpf_div(b, q, 9, round_up) == from_int(301)
    assert mpf_div(c, q, 9, round_up) == from_int(301)

    # Nearest even integer is down
    assert mpf_div(a, q, 9, round_nearest) == from_int(300)
    assert mpf_div(b, q, 9, round_nearest) == from_int(300)
    assert mpf_div(c, q, 9, round_nearest) == from_int(301)

    # Nearest even integer is up
    a = from_int(301499999)
    b = from_int(301500000)
    c = from_int(301500001)
    assert mpf_div(a, q, 9, round_nearest) == from_int(301)
    assert mpf_div(b, q, 9, round_nearest) == from_int(302)
    assert mpf_div(c, q, 9, round_nearest) == from_int(302)


def test_tight_integer_division():
    # Test that integer division at tightest possible precision is exact
    N = 100
    seed(1)
    for i in range(N):
        a = choice([1, -1]) * randint(1, 1<<randint(10, 100))
        b = choice([1, -1]) * randint(1, 1<<randint(10, 100))
        p = a * b
        width = b.bit_length() - trailing(b)
        a = from_int(a); b = from_int(b); p = from_int(p)
        for mode in [round_floor, round_ceiling, round_down,
                     round_up, round_nearest]:
            assert mpf_div(p, a, width, mode) == b


def test_epsilon_rounding():
    # Verify that mpf_div uses infinite precision; this result will
    # appear to be exactly 0.101 to a near-sighted algorithm

    a = from_str('0.101' + ('0'*200) + '1', base=2)
    b = from_str('1.10101', base=2)
    c = mpf_mul(a, b, 250, round_floor) # exact
    assert mpf_div(c, b, a[1].bit_length(), round_floor) == a # exact

    assert mpf_div(c, b, 2, round_down) == from_str('0.10', base=2)
    assert mpf_div(c, b, 3, round_down) == from_str('0.101', base=2)
    assert mpf_div(c, b, 2, round_up) == from_str('0.11', base=2)
    assert mpf_div(c, b, 3, round_up) == from_str('0.110', base=2)
    assert mpf_div(c, b, 2, round_floor) == from_str('0.10', base=2)
    assert mpf_div(c, b, 3, round_floor) == from_str('0.101', base=2)
    assert mpf_div(c, b, 2, round_ceiling) == from_str('0.11', base=2)
    assert mpf_div(c, b, 3, round_ceiling) == from_str('0.110', base=2)

    # The same for negative numbers
    a = from_str('-0.101' + ('0'*200) + '1', base=2)
    b = from_str('1.10101', base=2)
    c = mpf_mul(a, b, 250, round_floor)
    assert mpf_div(c, b, a[1].bit_length(), round_floor) == a

    assert mpf_div(c, b, 2, round_down) == from_str('-0.10', base=2)
    assert mpf_div(c, b, 3, round_up) == from_str('-0.110', base=2)

    # Floor goes up, ceiling goes down
    assert mpf_div(c, b, 2, round_floor) == from_str('-0.11', base=2)
    assert mpf_div(c, b, 3, round_floor) == from_str('-0.110', base=2)
    assert mpf_div(c, b, 2, round_ceiling) == from_str('-0.10', base=2)
    assert mpf_div(c, b, 3, round_ceiling) == from_str('-0.101', base=2)


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

def test_div_negative_rnd_bug():
    assert (-3) / mpf('0.1531879017645047') == mpf('-19.583791966887116')
    assert mpf('-2.6342475750861301') / mpf('0.35126216427941814') == mpf('-7.4993775104985909')

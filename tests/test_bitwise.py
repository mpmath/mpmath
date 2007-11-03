"""
Test bit-level integer operations
"""

from mpmath.lib import *

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
    assert normalize(73786976294838205979L, -65, 53, ROUND_HALF_EVEN) == (1, 1, 1)
    assert normalize(31, 0, 4, ROUND_UP) == (1, 5, 1)
    assert normalize(-31, 0, 4, ROUND_FLOOR) == (-1, 5, 1)
    assert normalize(255, 0, 7, ROUND_UP) == (1, 8, 1)
    assert normalize(-255, 0, 7, ROUND_FLOOR) == (-1, 8, 1)

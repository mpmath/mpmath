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

def test_trailing():
    assert trailing(0) == 0
    assert trailing(1) == 0
    assert trailing(2) == 1
    assert trailing(7) == 0
    assert trailing(8) == 3
    assert trailing(2**100) == 100
    assert trailing(2**100-1) == 0

def test_round_down():
    assert from_man_exp(0, -4, 4, round_down)[:2] == (0, 0)
    assert from_man_exp(0xf0, -4, 4, round_down)[:2] == (15, 0)
    assert from_man_exp(0xf1, -4, 4, round_down)[:2] == (15, 0)
    assert from_man_exp(0xff, -4, 4, round_down)[:2] == (15, 0)
    assert from_man_exp(-0xf0, -4, 4, round_down)[:2] == (-15, 0)
    assert from_man_exp(-0xf1, -4, 4, round_down)[:2] == (-15, 0)
    assert from_man_exp(-0xff, -4, 4, round_down)[:2] == (-15, 0)

def test_round_up():
    assert from_man_exp(0, -4, 4, round_up)[:2] == (0, 0)
    assert from_man_exp(0xf0, -4, 4, round_up)[:2] == (15, 0)
    assert from_man_exp(0xf1, -4, 4, round_up)[:2] == (1, 4)
    assert from_man_exp(0xff, -4, 4, round_up)[:2] == (1, 4)
    assert from_man_exp(-0xf0, -4, 4, round_up)[:2] == (-15, 0)
    assert from_man_exp(-0xf1, -4, 4, round_up)[:2] == (-1, 4)
    assert from_man_exp(-0xff, -4, 4, round_up)[:2] == (-1, 4)

def test_round_floor():
    assert from_man_exp(0, -4, 4, round_floor)[:2] == (0, 0)
    assert from_man_exp(0xf0, -4, 4, round_floor)[:2] == (15, 0)
    assert from_man_exp(0xf1, -4, 4, round_floor)[:2] == (15, 0)
    assert from_man_exp(0xff, -4, 4, round_floor)[:2] == (15, 0)
    assert from_man_exp(-0xf0, -4, 4, round_floor)[:2] == (-15, 0)
    assert from_man_exp(-0xf1, -4, 4, round_floor)[:2] == (-1, 4)
    assert from_man_exp(-0xff, -4, 4, round_floor)[:2] == (-1, 4)

def test_round_ceiling():
    assert from_man_exp(0, -4, 4, round_ceiling)[:2] == (0, 0)
    assert from_man_exp(0xf0, -4, 4, round_ceiling)[:2] == (15, 0)
    assert from_man_exp(0xf1, -4, 4, round_ceiling)[:2] == (1, 4)
    assert from_man_exp(0xff, -4, 4, round_ceiling)[:2] == (1, 4)
    assert from_man_exp(-0xf0, -4, 4, round_ceiling)[:2] == (-15, 0)
    assert from_man_exp(-0xf1, -4, 4, round_ceiling)[:2] == (-15, 0)
    assert from_man_exp(-0xff, -4, 4, round_ceiling)[:2] == (-15, 0)

def test_round_half_up():
    assert from_man_exp(0, -4, 4, round_half_up)[:2] == (0, 0)
    assert from_man_exp(0xf0, -4, 4, round_half_up)[:2] == (15, 0)
    assert from_man_exp(0xf7, -4, 4, round_half_up)[:2] == (15, 0)
    assert from_man_exp(0xf8, -4, 4, round_half_up)[:2] == (1, 4)
    assert from_man_exp(0xf9, -4, 4, round_half_up)[:2] == (1, 4)
    assert from_man_exp(0xff, -4, 4, round_half_up)[:2] == (1, 4)
    assert from_man_exp(-0xf0, -4, 4, round_half_up)[:2] == (-15, 0)
    assert from_man_exp(-0xf7, -4, 4, round_half_up)[:2] == (-15, 0)
    assert from_man_exp(-0xf8, -4, 4, round_half_up)[:2] == (-1, 4)
    assert from_man_exp(-0xf9, -4, 4, round_half_up)[:2] == (-1, 4)
    assert from_man_exp(-0xff, -4, 4, round_half_up)[:2] == (-1, 4)

def test_round_half_down():
    assert from_man_exp(0, -4, 4, round_half_down)[:2] == (0, 0)
    assert from_man_exp(0xf0, -4, 4, round_half_down)[:2] == (15, 0)
    assert from_man_exp(0xf7, -4, 4, round_half_down)[:2] == (15, 0)
    assert from_man_exp(0xf8, -4, 4, round_half_down)[:2] == (15, 0)
    assert from_man_exp(0xf9, -4, 4, round_half_down)[:2] == (1, 4)
    assert from_man_exp(0xff, -4, 4, round_half_down)[:2] == (1, 4)
    assert from_man_exp(-0xf0, -4, 4, round_half_down)[:2] == (-15, 0)
    assert from_man_exp(-0xf7, -4, 4, round_half_down)[:2] == (-15, 0)
    assert from_man_exp(-0xf8, -4, 4, round_half_down)[:2] == (-15, 0)
    assert from_man_exp(-0xf9, -4, 4, round_half_down)[:2] == (-1, 4)
    assert from_man_exp(-0xff, -4, 4, round_half_down)[:2] == (-1, 4)

def test_round_half_even():
    assert from_man_exp(0, -4, 4, round_half_even)[:2] == (0, 0)
    assert from_man_exp(0xf0, -4, 4, round_half_even)[:2] == (15, 0)
    assert from_man_exp(0xf7, -4, 4, round_half_even)[:2] == (15, 0)
    assert from_man_exp(0xf8, -4, 4, round_half_even)[:2] == (1, 4)    # 1111.1000 -> 10000.0
    assert from_man_exp(0xf9, -4, 4, round_half_even)[:2] == (1, 4)    # 1111.1001 -> 10000.0
    assert from_man_exp(0xe8, -4, 4, round_half_even)[:2] == (7, 1)    # 1110.1000 -> 1110.0
    assert from_man_exp(0xe9, -4, 4, round_half_even)[:2] == (15, 0)     # 1110.1001 -> 1111.0
    assert from_man_exp(-0xf0, -4, 4, round_half_even)[:2] == (-15, 0)
    assert from_man_exp(-0xf7, -4, 4, round_half_even)[:2] == (-15, 0)
    assert from_man_exp(-0xf8, -4, 4, round_half_even)[:2] == (-1, 4)
    assert from_man_exp(-0xf9, -4, 4, round_half_even)[:2] == (-1, 4)
    assert from_man_exp(-0xe8, -4, 4, round_half_even)[:2] == (-7, 1)
    assert from_man_exp(-0xe9, -4, 4, round_half_even)[:2] == (-15, 0)

def test_rounding_bugs():
    # 1 less than power-of-two cases
    assert from_man_exp(72057594037927935, -56, 53, round_up) == (1, 0, 1)
    assert from_man_exp(73786976294838205979l, -65, 53, round_half_even) == (1, 1, 1)
    assert from_man_exp(31, 0, 4, round_up) == (1, 5, 1)
    assert from_man_exp(-31, 0, 4, round_floor) == (-1, 5, 1)
    assert from_man_exp(255, 0, 7, round_up) == (1, 8, 1)
    assert from_man_exp(-255, 0, 7, round_floor) == (-1, 8, 1)

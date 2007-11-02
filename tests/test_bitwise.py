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

def test_divide_floor():
    assert divide_floor(20, 5) == 4
    assert divide_floor(20, -5) == -4
    assert divide_floor(-20, 5) == -4
    assert divide_floor(-20, -5) == 4
    assert divide_floor(20, 3) == 6
    assert divide_floor(-20, 3) == -7
    assert divide_floor(20, -3) == -7
    assert divide_floor(-20, -3) == 6

def test_divide_ceiling():
    assert divide_ceiling(20, 5) == 4
    assert divide_ceiling(20, -5) == -4
    assert divide_ceiling(-20, 5) == -4
    assert divide_ceiling(-20, -5) == 4
    assert divide_ceiling(20, 3) == 7
    assert divide_ceiling(-20, 3) == -6
    assert divide_ceiling(20, -3) == -6
    assert divide_ceiling(-20, -3) == 7

def test_divide_down():
    assert divide_down(20, 5) == 4
    assert divide_down(20, -5) == -4
    assert divide_down(-20, 5) == -4
    assert divide_down(-20, -5) == 4
    assert divide_down(20, 3) == 6
    assert divide_down(-20, 3) == -6
    assert divide_down(20, -3) == -6
    assert divide_down(-20, -3) == 6

def test_divide_up():
    assert divide_up(20, 5) == 4
    assert divide_up(20, -5) == -4
    assert divide_up(-20, 5) == -4
    assert divide_up(-20, -5) == 4
    assert divide_up(20, 3) == 7
    assert divide_up(-20, 3) == -7
    assert divide_up(20, -3) == -7
    assert divide_up(-20, -3) == 7

def test_divide_half_up():
    assert divide_half_up(5, 1) == 5
    assert divide_half_up(5, 2) == 3
    assert divide_half_up(5, 3) == 2
    assert divide_half_up(5, 4) == 1
    assert divide_half_up(5, 5) == 1
    assert divide_half_up(6, 4) == 2
    assert divide_half_up(100, 2) == 50
    assert divide_half_up(101, 2) == 51
    assert divide_half_up(300499999, 1000000) == 300
    assert divide_half_up(300500000, 1000000) == 301
    assert divide_half_up(300500001, 1000000) == 301
    assert divide_half_up(-300499999, 1000000) == -300
    assert divide_half_up(-300500000, 1000000) == -301
    assert divide_half_up(-300500001, 1000000) == -301

def test_divide_half_down():
    assert divide_half_down(5, 1) == 5
    assert divide_half_down(5, 2) == 2
    assert divide_half_down(5, 3) == 2
    assert divide_half_down(5, 4) == 1
    assert divide_half_down(5, 5) == 1
    assert divide_half_down(6, 4) == 1
    assert divide_half_down(100, 2) == 50
    assert divide_half_down(101, 2) == 50
    assert divide_half_down(300499999, 1000000) == 300
    assert divide_half_down(300500000, 1000000) == 300
    assert divide_half_down(300500001, 1000000) == 301
    assert divide_half_down(-300499999, 1000000) == -300
    assert divide_half_down(-300500000, 1000000) == -300
    assert divide_half_down(-300500001, 1000000) == -301

def test_divide_half_even():
    for j in range(1, 100):
        for k in range(1, 100):
            q = divide_half_even(j, k)
            fq = float(j)/k
            if q*k == j:
                assert q == fq
            else:
                if fq % 1.0 == 0.5: assert not q & 1
                if fq % 1.0 < 0.5:  assert q < fq
                if fq % 1.0 > 0.5:  assert q > fq
    assert divide_half_even(5, 1) == 5
    assert divide_half_even(5, 2) == 2
    assert divide_half_even(5, 3) == 2
    assert divide_half_even(5, 4) == 1
    assert divide_half_even(5, 5) == 1
    assert divide_half_even(6, 4) == 2
    assert divide_half_even(10, 4) == 2
    assert divide_half_even(14, 4) == 4
    assert divide_half_even(18, 4) == 4
    assert divide_half_even(100, 2) == 50
    assert divide_half_even(101, 2) == 50
    assert divide_half_even(102, 2) == 51
    assert divide_half_even(103, 2) == 52
    assert divide_half_even(104, 2) == 52
    assert divide_half_even(105, 2) == 52
    assert divide_half_even(106, 2) == 53
    assert divide_half_even(50479*318, 318) == 50479
    assert divide_half_even(50479*317, 317) == 50479
    assert divide_half_even(300499999, 1000000) == 300
    assert divide_half_even(300500000, 1000000) == 300
    assert divide_half_even(300500001, 1000000) == 301
    assert divide_half_even(300499999, -1000000) == -300
    assert divide_half_even(300500000, -1000000) == -300
    assert divide_half_even(300500001, -1000000) == -301

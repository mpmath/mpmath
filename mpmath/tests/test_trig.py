from mpmath import *

def test_trig_near_zero():
    mp.dps = 15

    mp.rounding = 'nearest'; assert sin(0) == 0 and cos(0) == 1
    mp.rounding = 'down';    assert sin(0) == 0 and cos(0) == 1
    mp.rounding = 'floor';   assert sin(0) == 0 and cos(0) == 1
    mp.rounding = 'up';      assert sin(0) == 0 and cos(0) == 1
    mp.rounding = 'ceiling'; assert sin(0) == 0 and cos(0) == 1
    mp.rounding = 'nearest';

    a = mpf('1e-100')
    b = mpf('-1e-100')

    mp.rounding = 'nearest'; assert sin(a) == a
    mp.rounding = 'down';    assert sin(a) < a
    mp.rounding = 'floor';   assert sin(a) < a
    mp.rounding = 'up';      assert sin(a) == a
    mp.rounding = 'ceiling'; assert sin(a) == a
    mp.rounding = 'nearest'; assert sin(b) == b
    mp.rounding = 'down';    assert sin(b) > b
    mp.rounding = 'floor';   assert sin(b) == b
    mp.rounding = 'up';      assert sin(b) == b
    mp.rounding = 'ceiling'; assert sin(b) > b
    mp.rounding = 'nearest'

    mp.rounding = 'nearest'; assert cos(a) == 1
    mp.rounding = 'down';    assert cos(a) < 1
    mp.rounding = 'floor';   assert cos(a) < 1
    mp.rounding = 'up';      assert cos(a) == 1
    mp.rounding = 'ceiling'; assert cos(a) == 1
    mp.rounding = 'nearest'; assert cos(b) == 1
    mp.rounding = 'down';    assert cos(b) < 1
    mp.rounding = 'floor';   assert cos(b) < 1
    mp.rounding = 'up';      assert cos(b) == 1
    mp.rounding = 'ceiling'; assert cos(b) == 1
    mp.rounding = 'nearest'


def test_trig_near_n_pi():

    mp.dps = 15
    a = [n*pi for n in [1, 2, 6, 11, 100, 1001, 10000, 100001]]
    mp.dps = 135
    a.append(10**100 * pi)
    mp.dps = 15

    assert sin(a[0]) == mpf('1.2246467991473531772e-16')
    assert sin(a[1]) == mpf('-2.4492935982947063545e-16')
    assert sin(a[2]) == mpf('-7.3478807948841190634e-16')
    assert sin(a[3]) == mpf('4.8998251578625894243e-15')
    assert sin(a[4]) == mpf('1.9643867237284719452e-15')
    assert sin(a[5]) == mpf('-8.8632615209684813458e-15')
    assert sin(a[6]) == mpf('-4.8568235395684898392e-13')
    assert sin(a[7]) == mpf('3.9087342299491231029e-11')
    assert sin(a[8]) == mpf('-1.369235466754566993528e-36')

    mp.rounding = 'nearest'
    assert cos(a[0]) == -1
    assert cos(a[1]) == 1
    assert cos(a[2]) == 1
    assert cos(a[3]) == -1
    assert cos(a[4]) == 1
    assert cos(a[5]) == -1
    assert cos(a[6]) == 1
    assert cos(a[7]) == -1
    assert cos(a[8]) == 1

    mp.rounding = 'up'
    assert cos(a[0]) == -1
    assert cos(a[1]) == 1
    assert cos(a[2]) == 1
    assert cos(a[3]) == -1
    assert cos(a[4]) == 1
    assert cos(a[5]) == -1
    assert cos(a[6]) == 1
    assert cos(a[7]) == -1
    assert cos(a[8]) == 1

    mp.rounding = 'down'
    assert cos(a[0]) > -1
    assert cos(a[1]) < 1
    assert cos(a[2]) < 1
    assert cos(a[3]) > -1
    assert cos(a[4]) < 1
    assert cos(a[5]) > -1
    assert cos(a[6]) < 1
    assert cos(a[7]) > -1
    assert cos(a[8]) < 1

    mp.rounding = 'floor'
    assert cos(a[0]) == -1
    assert cos(a[1]) < 1
    assert cos(a[2]) < 1
    assert cos(a[3]) == -1
    assert cos(a[4]) < 1
    assert cos(a[5]) == -1
    assert cos(a[6]) < 1
    assert cos(a[7]) == -1
    assert cos(a[8]) < 1

    mp.rounding = 'ceiling'
    assert cos(a[0]) > -1
    assert cos(a[1]) == 1
    assert cos(a[2]) == 1
    assert cos(a[3]) > -1
    assert cos(a[4]) == 1
    assert cos(a[5]) > -1
    assert cos(a[6]) == 1
    assert cos(a[7]) > -1
    assert cos(a[8]) == 1

    mp.dps = 15
    mp.rounding = 'nearest'

if __name__ == '__main__':
    for f in globals().keys():
        if f.startswith("test_"):
            print f
            globals()[f]()

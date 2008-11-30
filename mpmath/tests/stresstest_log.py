import sys
from gmpy import *
sys.path.insert(0, '../..')
sys.path.insert(0, '..')

from time import time
from mpmath import *

def check_last_20digits(x, sign, p, np, q, nq, res):
    """check that the last 20 digits of x + sign * p**np/q**nq are right
    """
    for dps in res.keys():
        mp.dps = dps
        a = x + sign * p**np / q**nq
        #print nstr(a, 6)
        s = str(log(a))
        assert s[-20:] == res[dps], 'dps=%d a =%s/%s %s != %s' % \
        (mp.dps, p, q, s[-20:], res[dps])

def test_ln2():
    res = {10000: '91340185660135965556', 20000:'20201430377863539655',
           100000: '02205469487696859273', 200000: '80239565939562135199'}
    a = mpf(2)
    check_last_20digits(0, 1, mpf(2), 1, mpf(1), 1, res)

def test_log_large():
    res = {10000: '65959485048291178539', 50000: '07042167126154502209',
           100000: '21539616499670455924'}
    check_last_20digits(0, 1, mpf(10), 4, 3, 1, res)
    res = {10000: '54020631943060007154', 50000: '22874248609981100075',
           100000: '53982029524801368512'}
    check_last_20digits(0, 1, mpf(10), 10, 3, 1, res)
    res = {10000: '36539088351652334665', 50000: '15542472707046584789',
           100000: '88140304764731621231'}
    check_last_20digits(0, 1, mpf(10), 100, 3, 1, res)

def test_log_small():
    res = {10000: '08787418379034837194', 50000: '54342973371838910333',
           100000: '79102579092603240627'}
    check_last_20digits(1, 1, mpf(1), 1, mpf(10), 1, res)
    res = {10000: '86309966505553152412', 50000: '55433962116376422785',
           100000: '50851077929404777824'}
    check_last_20digits(1, -1, mpf(1), 1, mpf(10), 1, res)
    res = {10000: '78315031758836279343', 50000: '29465377538229478283',
           100000: '00573911856777927276'}
    check_last_20digits(1, -1, mpf(1), 1, mpf(10), 5, res)
    res = {10000: '37216724048352523026', 50000: '11655094657459854621',
           100000: '97088358959385804055'}
    check_last_20digits(1, -1, mpf(1), 1, mpf(10), 10, res)

    res = {10000: '27664533076359806938', 50000: '00872849495437295622',
           100000: '89772025034160535926'}
    check_last_20digits(1, 1, mpf(1), 1, mpf(10), 100, res)


if __name__ == '__main__':
    for f in globals().keys():
        if f.startswith("test_"):
            t0 = time()
            globals()[f]()
            print '%s %.2fs' % (globals()[f].__name__, time() - t0)

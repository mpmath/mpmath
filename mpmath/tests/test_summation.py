from mpmath import *

def test_sumem():
    mp.dps = 15
    assert sumem(lambda k: 1/k**2.5, [50, 100]).ae(0.0012524505324784962)
    assert sumem(lambda k: k**4 + 3*k + 1, [10, 100]).ae(2050333103)

def test_sumsh():
    mp.dps = 15
    assert sumsh(lambda k: (-1)**(k+1) / k, [1, inf]).ae(log(2))
    assert sumsh(lambda k: (-1)**(k+1) / k**2, [1, inf]).ae(pi**2 / 12)
    assert sumsh(lambda k: (-1)**k / log(k), [2, inf]).ae(0.9242998972229388)
    assert sumsh(lambda k: 1/factorial(k), [0, inf]).ae(e)

def test_sumrich():
    mp.dps = 15
    assert sumrich(lambda k: 1/k**2, [1, inf]).ae(pi**2 / 6)

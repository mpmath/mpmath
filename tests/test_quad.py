from mpmath import *
from mpmath.quad import *

# This test takes a LOT of time to run without psyco
try:
    import psyco
    psyco.full()
except ImportError:
    pass

def ae(a, b):
    return abs(a-b) < 10**(-mpf.dps+5)

def test_basic_integrals():
    for prec in [15, 30, 100]:
        assert ae(quadts(lambda x: x**3 - 3*x**2, -2, 4), -12)
        assert ae(quadts(sin, 0, pi), 2)
        assert ae(quadts(sin, 0, 2*pi), 0)
        assert ae(quadts(exp, -inf, -1), 1/e)
        assert ae(quadts(lambda x: exp(-x), 0, inf), 1)
        assert ae(quadts(lambda x: exp(-x**2), -inf, inf), sqrt(pi))
        assert ae(quadts(lambda x: 1/(1+x**2), -1, 1), pi/2)
        assert ae(quadts(lambda x: 1/(1+x**2), -inf, inf), pi)
        assert ae(quadts(lambda x: 2*sqrt(1-x**2), -1, 1), pi)

def test_double_integrals():
    for prec in [15, 30]:
        assert ae(quadts(lambda x, y: cos(x+y/2), (-pi/2, pi/2), (0, pi)), 4)
        assert ae(quadts(lambda x, y: (x-1)/((1-x*y)*log(x*y)), (0, 1), (0, 1)), cgamma)
        assert ae(quadts(lambda x, y: 1/sqrt(1+x**2+y**2), (-1, 1), (-1, 1)), 4*log(2+sqrt(3))-2*pi/3)
        assert ae(quadts(lambda x, y: 1/(1-x**2 * y**2), (0, 1), (0, 1)), pi**2 / 8)
        assert ae(quadts(lambda x, y: 1/(1-x*y), (0, 1), (0, 1)), pi**2 / 6)
        assert ae(quadts(lambda x, y: exp(-(x+y)), (0, inf), (0, inf)), 1)
        # very slow
        #assert ae(quadts(lambda x, y: exp(-x**2-y**2), (-inf, inf), (-inf, inf)), pi)

# Test integrals from "Experimentation in Mathematics" by Borwein,
# Bailey & Girgensohn
def test_expmath_integrals():
    for prec in [15, 30, 50]:
        mpf.dps = prec
        assert ae(quadts(lambda x: x/sinh(x), 0, inf),                    pi**2 / 4)
        assert ae(quadts(lambda x: log(x)**2 / (1+x**2), 0, inf),         pi**3 / 8)
        assert ae(quadts(lambda x: (1+x**2)/(1+x**4), 0, inf),            pi/sqrt(2))
        assert ae(quadts(lambda x: log(x)/cosh(x)**2, 0, inf),            log(pi)-2*log(2)-cgamma)
        assert ae(quadts(lambda x: log(1+x**3)/(1-x+x**2), 0, inf),       2*pi*log(3)/sqrt(3))
        assert ae(quadts(lambda x: log(x)**2 / (x**2+x+1), 0, 1),         8*pi**3 / (81*sqrt(3)))
        assert ae(quadts(lambda x: log(cos(x))**2, 0, pi/2),              pi/2 * (log(2)**2+pi**2/12))
        assert ae(quadts(lambda x: x**2 / sin(x)**2, 0, pi/2),            pi*log(2))
        assert ae(quadts(lambda x: x**2/sqrt(exp(x)-1), 0, inf),          4*pi*(log(2)**2 + pi**2/12))
        assert ae(quadts(lambda x: x*exp(-x)*sqrt(1-exp(-2*x)), 0, inf),  pi*(1+2*log(2))/8)
    mpf.dps = 15

# Do not reach full accuracy
def xtest_expmath_fail():
    assert ae(quadts(lambda x: sqrt(tan(x)), 0, pi/2),          pi*sqrt(2)/2)
    assert ae(quadts(lambda x: atan(x)/(x*sqrt(1-x**2)), 0, 1), pi*log(1+sqrt(2))/2)
    assert ae(quadts(lambda x: log(1+x**2)/x**2, 0, 1),         pi/2-log(2))
    assert ae(quadts(lambda x: x**2/((1+x**4)*sqrt(1-x**4)), 0, 1),     pi/8)

import decimal
import random
from decimal import Decimal
from fractions import Fraction

import pytest

from mpmath import inf, isnan, iv, mp, mpc, mpf, mpi, mpmathify, sqrt
from mpmath.libmp import (fhalf, from_float, from_rational, from_str,
                          round_ceiling, round_floor, round_nearest,
                          to_rational, to_str)


def test_basic_string():
    """
    Test basic string conversion
    """
    assert mpf('3') == mpf('3.0') == mpf('0003.') == mpf('0.03e2') == mpf(3.0)
    assert mpf('30') == mpf('30.0') == mpf('00030.') == mpf(30.0)
    for i in range(10):
        for j in range(10):
            assert mpf('%ie%i' % (i,j)) == i * 10**j
    assert str(mpf('25000.0')) == '25000.0'
    assert str(mpf('2500.0')) == '2500.0'
    assert str(mpf('250.0')) == '250.0'
    assert str(mpf('25.0')) == '25.0'
    assert str(mpf('2.5')) == '2.5'
    assert str(mpf('0.25')) == '0.25'
    assert str(mpf('0.025')) == '0.025'
    assert str(mpf('0.0025')) == '0.0025'
    assert str(mpf('0.00025')) == '0.00025'
    assert str(mpf('0.000025')) == '2.5e-5'
    assert str(mpf(0)) == '0.0'
    assert str(mpf('2.5e1000000000000000000000')) == '2.5e+1000000000000000000000'
    assert str(mpf('2.6e-1000000000000000000000')) == '2.6e-1000000000000000000000'
    assert str(mpf(1.23402834e-15)) == '1.23402834e-15'
    assert str(mpf(-1.23402834e-15)) == '-1.23402834e-15'
    assert str(mpf(-1.2344e-15)) == '-1.2344e-15'
    assert repr(mpf(-1.2344e-15)) == "mpf('-1.2343999999999999e-15')"
    assert str(mpf("2163048125L")) == '2163048125.0'
    assert str(mpf("-2163048125l")) == '-2163048125.0'
    assert str(mpf("-2163048125L/1088391168")) == '-1.98738118113799'
    assert str(mpf("2163048125/1088391168l")) == '1.98738118113799'
    assert str(mpf('inf')) == 'inf'

    # issue 613
    assert str(mpf('2_5_0_0.0')) == '2500.0'
    # issue 377
    assert to_str(from_str('1_234.567891', 80), 24) == '1234.567891'
    assert to_str(from_str('1_234.567_891', 80), 24) == '1234.567891'
    assert to_str(from_str('1_234.567_8_9_1', 80), 24) == '1234.567891'
    assert to_str(from_str('1.0_0', 80), 24) == '1.0'
    assert to_str(from_str('.000', 80), 24) == '0.0'

def test_from_str():
    assert mpf(from_str('ABC.ABC', base=16)) == mpf(float.fromhex('ABC.ABC'))
    assert mpf(from_str('0xABC.ABC')) == mpf(float.fromhex('ABC.ABC'))
    assert mpf(from_str('0x3.a7p10')) == mpf(float.fromhex('0x3.a7p10'))
    assert mpf(from_str('0x1.4ace478p+33')) == mpf(float.fromhex('0x1.4ace478p+33'))
    assert mpf(from_str('0x1.4ace478@+33')) == mpf('7.0354608312666732e+39')
    assert mpf(from_str('0b1101.100101')) == mpf('13.578125')
    assert mpf(from_str('0o1101.100101')) == mpf('577.12524795532227')
    assert mpf(from_str('1.99999999', prec=0)) == mpf('1.9999999901046976')

def test_eps_repr():
    mp.dps = 24
    assert repr(mp.eps) == '<epsilon of working precision: 2.06795e-25~>'

def test_to_str():
    assert to_str(from_str('ABC.ABC', base=16), 6, base=16) == '0xabc.abc'
    assert to_str(from_str('0x3.a7p10', base=16), 3, base=16) == '0xe9c.0'
    assert to_str(from_str('0x1.4ace478p+33'), 7, base=16) == '0x2.959c8f@+8'
    assert to_str(from_str('0o1101.100101'), 8, base=8) == '0o1101.1001'
    assert to_str(from_str('0b1101.100101'), 10, base=2) == '0b1101.100101'
    assert to_str(from_str('0x1.4ace478p+33'), 8, base=16, binary_exp=True) == '0x1.4ace478p+33'
    assert to_str(from_str('0x1.4ace478p+33'), 7, base=16, binary_exp=True) == '0x1.4ace48p+33'
    assert to_str(from_str('0x1.4ace478p+33'), 5, base=16, binary_exp=True) == '0x1.4acep+33'
    assert to_str(from_str('1', base=16), 6, base=16, binary_exp=True) == '0x1.0'
    x = mpf('1234.567891')._mpf_
    pytest.raises(ValueError, lambda: to_str(x, 6, binary_exp=True))
    pytest.raises(ValueError, lambda: to_str(x, 6, rnd='Y'))
    pytest.raises(ValueError, lambda: to_str('1e400e2', 6))
    assert to_str(x, 5, rnd='n') == '1234.6'
    assert to_str(x, 5, rnd='d') == '1234.5'
    assert to_str(x, 5, rnd='u') == '1234.6'

def test_pretty():
    mp.pretty = True
    assert repr(mpf(2.5)) == '2.5'
    assert repr(mpc(2.5,3.5)) == '(2.5 + 3.5j)'
    iv.pretty = True
    assert repr(mpi(2.5,3.5)) == '[2.5, 3.5]'

def test_str_whitespace():
    assert mpf('1.26 ') == 1.26

def test_str_format():
    assert to_str(from_float(0.1),15,strip_zeros=False) == '0.100000000000000'
    assert to_str(from_float(0.0),15,show_zero_exponent=True) == '0.0e+0'
    assert to_str(from_float(0.0),0,show_zero_exponent=True) == '.0e+0'
    assert to_str(from_float(0.0),0,show_zero_exponent=False) == '.0'
    assert to_str(from_float(0.0),1,show_zero_exponent=True) == '0.0e+0'
    assert to_str(from_float(0.0),1,show_zero_exponent=False) == '0.0'
    assert to_str(from_float(1.23),3,show_zero_exponent=True) == '1.23e+0'
    assert to_str(from_float(1.23456789000000e-2),15,strip_zeros=False,min_fixed=0,max_fixed=0) == '1.23456789000000e-2'
    assert to_str(from_float(1.23456789000000e+2),15,strip_zeros=False,min_fixed=0,max_fixed=0) == '1.23456789000000e+2'
    assert to_str(from_float(2.1287e14), 15, max_fixed=1000) == '212870000000000.0'
    assert to_str(from_float(2.1287e15), 15, max_fixed=1000) == '2128700000000000.0'
    assert to_str(from_float(2.1287e16), 15, max_fixed=1000) == '21287000000000000.0'
    assert to_str(from_float(2.1287e30), 15, max_fixed=1000) == '2128700000000000000000000000000.0'

def test_tight_string_conversion():
    # In an old version, '0.5' wasn't recognized as representing
    # an exact binary number and was erroneously rounded up or down
    assert from_str('0.5', 10, round_floor) == fhalf
    assert from_str('0.5', 10, round_ceiling) == fhalf

def test_eval_repr_invariant():
    """Test that eval(repr(x)) == x"""
    random.seed(123)
    for dps in [10, 15, 20, 50, 100]:
        mp.dps = dps
        for i in range(1000):
            a = mpf(random.random())**0.5 * 10**random.randint(-100, 100)
            assert eval(repr(a)) == a

def test_str_bugs():
    # Decimal rounding used to give the wrong exponent in some cases
    assert str(mpf('1e600')) == '1.0e+600'
    assert str(mpf('1e10000')) == '1.0e+10000'

def test_str_prec0():
    assert to_str(from_float(1.234), 0) == '.0e+0'
    assert to_str(from_float(1e-15), 0) == '.0e-15'
    assert to_str(from_float(1e+15), 0) == '.0e+15'
    assert to_str(from_float(-1e-15), 0) == '-.0e-15'
    assert to_str(from_float(-1e+15), 0) == '-.0e+15'

def test_convert_rational():
    assert from_rational(30, 5, 53, round_nearest) == (0, 3, 1, 2)
    assert from_rational(-7, 4, 53, round_nearest) == (1, 7, -2, 3)
    assert to_rational((0, 1, -1, 1)) == (1, 2)
    assert to_rational((0, 1, 0, 1)) == (1, 1)
    pytest.raises(ValueError, lambda: to_rational(mpf('nan')._mpf_))
    pytest.raises(OverflowError, lambda: to_rational(mpf('inf')._mpf_))
    pytest.raises(OverflowError, lambda: to_rational(mpf('-inf')._mpf_))

def test_custom_class():
    class mympf:
        @property
        def _mpf_(self):
            return mpf(3.5)._mpf_
    class mympc:
        @property
        def _mpc_(self):
            return mpf(3.5)._mpf_, mpf(2.5)._mpf_
    assert mpf(2) + mympf() == 5.5
    assert mympf() + mpf(2) == 5.5
    assert mpf(mympf()) == 3.5
    assert mympc() + mpc(2) == mpc(5.5, 2.5)
    assert mpc(2) + mympc() == mpc(5.5, 2.5)
    assert mpc(mympc()) == (3.5+2.5j)
    assert mpmathify(mympf()) == mpf(3.5)
    assert mpmathify(mympc()) == mpc(3.5, 2.5)

def test_conversion_methods():
    class SomethingRandom:
        pass
    class SomethingReal:
        def _mpmath_(self, prec, rounding):
            return mp.make_mpf(from_str('1.3', prec, rounding))
    class SomethingComplex:
        def _mpmath_(self, prec, rounding):
            return mp.make_mpc((from_str('1.3', prec, rounding), \
                from_str('1.7', prec, rounding)))
    x = mpf(3)
    z = mpc(3)
    a = SomethingRandom()
    y = SomethingReal()
    w = SomethingComplex()
    for d in [15, 45]:
        mp.dps = d
        assert (x+y).ae(mpf('4.3'))
        assert (y+x).ae(mpf('4.3'))
        assert (x+w).ae(mpc('4.3', '1.7'))
        assert (w+x).ae(mpc('4.3', '1.7'))
        assert (z+y).ae(mpc('4.3'))
        assert (y+z).ae(mpc('4.3'))
        assert (z+w).ae(mpc('4.3', '1.7'))
        assert (w+z).ae(mpc('4.3', '1.7'))
        x-y; y-x; x-w; w-x; z-y; y-z; z-w; w-z
        x*y; y*x; x*w; w*x; z*y; y*z; z*w; w*z
        x/y; y/x; x/w; w/x; z/y; y/z; z/w; w/z
        x**y; y**x; x**w; w**x; z**y; y**z; z**w; w**z
        x==y; y==x; x==w; w==x; z==y; y==z; z==w; w==z
    mp.dps = 15
    assert x.__add__(a) is NotImplemented
    assert x.__radd__(a) is NotImplemented
    assert x.__lt__(a) is NotImplemented
    assert x.__gt__(a) is NotImplemented
    assert x.__le__(a) is NotImplemented
    assert x.__ge__(a) is NotImplemented
    assert x.__eq__(a) is NotImplemented
    assert x.__ne__(a) is NotImplemented
    assert x.__sub__(a) is NotImplemented
    assert x.__rsub__(a) is NotImplemented
    assert x.__mul__(a) is NotImplemented
    assert x.__rmul__(a) is NotImplemented
    assert x.__truediv__(a) is NotImplemented
    assert x.__rtruediv__(a) is NotImplemented
    assert x.__mod__(a) is NotImplemented
    assert x.__rmod__(a) is NotImplemented
    assert x.__pow__(a) is NotImplemented
    assert x.__rpow__(a) is NotImplemented
    assert z.__add__(a) is NotImplemented
    assert z.__radd__(a) is NotImplemented
    assert z.__eq__(a) is NotImplemented
    assert z.__ne__(a) is NotImplemented
    assert z.__sub__(a) is NotImplemented
    assert z.__rsub__(a) is NotImplemented
    assert z.__mul__(a) is NotImplemented
    assert z.__rmul__(a) is NotImplemented
    assert z.__truediv__(a) is NotImplemented
    assert z.__rtruediv__(a) is NotImplemented
    assert z.__pow__(a) is NotImplemented
    assert z.__rpow__(a) is NotImplemented

def test_mpmathify():
    assert mpmathify('1/2') == 0.5
    assert mpmathify('(1.0+1.0j)') == mpc(1, 1)
    assert mpmathify('(1.2e-10 - 3.4e5j)') == mpc('1.2e-10', '-3.4e5')
    assert mpmathify('1j') == mpc(1j)
    assert mpmathify('oo') == mpf('inf')
    assert mpmathify('+oo') == mpf('inf')
    assert mpmathify('-oo') == mpf('-inf')
    assert mpmathify('2+3*I') == mpc(2, 3)
    assert mpmathify('2+3I') == mpc(2, 3)
    assert mpmathify('2/3 + 4/5j') == mpc(2/3, 4/5)

def test_issue548():
    try:
        # This expression is invalid, but may trigger the ReDOS vulnerability
        # in the regular expression for parsing complex numbers.
        mpmathify('(' + '1' * 5000 + '!j')
    except:
        return
    # The expression is invalid and should raise an exception.
    assert False

def test_compatibility():
    from packaging.version import Version, parse
    np = pytest.importorskip("numpy")
    if parse(np.__version__) < Version('2.0.0b1'):
        npcore = np.core
    else:
        npcore = np._core
    # numpy types
    for nptype in npcore.numerictypes.typeDict.values():
        if issubclass(nptype, np.complexfloating):
            x = nptype(complex(0.5, -0.5))
        elif issubclass(nptype, np.floating):
            x = nptype(0.5)
        elif issubclass(nptype, np.integer):
            x = nptype(2)
        # Handle the weird types
        try: diff = np.abs(type(np.sqrt(x))(sqrt(x)) - np.sqrt(x))
        except: continue
        assert diff < np.float64(2.0**-53)
    assert mpf(np.float64('inf')) == inf
    assert isnan(mp.npconvert(np.float64('nan')))
    if hasattr(np, "float128"):
        mp.prec = 64
        assert (mp.npconvert(np.float128('0.841470984807896506652502321630298954')) ==
                mpf('0.841470984807896506653'))
    mp.prec = 53
    # issues 382 and 539
    assert mp.sqrt(np.int64(1)) == mpf('1.0')
    assert mpf(np.int64(1)) == mpf('1.0')
    #Fraction and Decimal
    oldprec = mp.prec
    mp.prec = 1000
    decimal.getcontext().prec = mp.dps
    assert sqrt(Fraction(2, 3)).ae(sqrt(mpf('2/3')))
    assert sqrt(Decimal(2)/Decimal(3)).ae(sqrt(mpf('2/3')))
    mp.prec = oldprec
    assert mpmathify(np.array(123)) == mpf(123)
    assert mpmathify(np.array(1.25)) == mpf(1.25)
    assert mpmathify(np.array(0.5+1j)) == mpc(0.5+1j)
    pytest.raises(TypeError, lambda: mpmathify(np.array([1])))

def test_issue465():
    assert mpf(Fraction(1, 3)) == mpf('0.33333333333333331')

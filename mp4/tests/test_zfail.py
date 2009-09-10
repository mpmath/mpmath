'''
def test_special_nan_cmp():
    assert not (nan < 3)
    assert not (nan > 3)
    assert not (mpf(3) < float_nan)
    assert not (mpf(3) > float_nan)
    assert not (mpf(3) <= float_nan)
    assert not (mpf(3) >= float_nan)
    assert nan != nan
    assert abs(nan) != abs(nan)

@extradps(50)
def test_precision():
    A = randmatrix(10, 10)
    assert mnorm(inverse(inverse(A)) - A, 1) < 1.e-45


def test_arithmetic_functions():
    import operator
    ops = [(operator.add, fadd), (operator.sub, fsub), (operator.mul, fmul),
        (operator.div, fdiv)]
    a = mpf(0.27)
    b = mpf(1.13)
    c = mpc(0.51+2.16j)
    d = mpc(1.08-0.99j)
    for x in [a,b,c,d]:
        for y in [a,b,c,d]:
            for op, fop in ops:
                if fop is not fdiv:
                    mp.prec = 200
                    z0 = op(x,y)
                mp.prec = 60
                z1 = op(x,y)
                mp.prec = 53
                z2 = op(x,y)
                assert fop(x, y, prec=60) == z1
                assert fop(x, y) == z2
                if fop is not fdiv:
                    assert fop(x, y, prec=inf) == z0
                    assert fop(x, y, dps=inf) == z0
                    assert fop(x, y, exact=True) == z0
                assert fneg(fneg(z1, exact=True), prec=inf) == z1
                assert fneg(z1) == -(+z1)
    mp.dps = 15

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

def test_conversion_methods():
    class SomethingRandom:
        pass
    class SomethingReal:
        def _mpmath_(self, prec, rounding):
            return make_mpf(from_str('1.3', prec, rounding))
    class SomethingComplex:
        def _mpmath_(self, prec, rounding):
            return make_mpc((from_str('1.3', prec, rounding), \
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
    assert x.__cmp__(a) is NotImplemented
    assert x.__sub__(a) is NotImplemented
    assert x.__rsub__(a) is NotImplemented
    assert x.__mul__(a) is NotImplemented
    assert x.__rmul__(a) is NotImplemented
    assert x.__div__(a) is NotImplemented
    assert x.__rdiv__(a) is NotImplemented
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
    assert z.__div__(a) is NotImplemented
    assert z.__rdiv__(a) is NotImplemented
    assert z.__pow__(a) is NotImplemented
    assert z.__rpow__(a) is NotImplemented

def test_pretty():
    mp.pretty = True
    assert repr(mpf(2.5)) == '2.5'
    assert repr(mpc(2.5,3.5)) == '(2.5 + 3.5j)'
    assert repr(mpi(2.5,3.5)) == '[2.5, 3.5]'
    mp.pretty = False

def test_atan_inf_sign():
    mp.dps = 15
    pi2 = pi/2
    assert atan(mpc(inf,-1)).ae(pi2)
    assert atan(mpc(inf,0)).ae(pi2)
    assert atan(mpc(inf,1)).ae(pi2)
    assert atan(mpc(1,inf)).ae(pi2)
    assert atan(mpc(0,inf)).ae(pi2)
    assert atan(mpc(-1,inf)).ae(-pi2)
    assert atan(mpc(-inf,1)).ae(-pi2)
    assert atan(mpc(-inf,0)).ae(-pi2)
    assert atan(mpc(-inf,-1)).ae(-pi2)
    assert atan(mpc(-1,-inf)).ae(-pi2)
    assert atan(mpc(0,-inf)).ae(-pi2)
    assert atan(mpc(1,-inf)).ae(pi2)

def test_atanh_inf_sign():
    # Limits at infinity
    jpi2 = j*pi/2
    assert atanh(inf).ae(-jpi2)
    assert atanh(-inf).ae(jpi2)
    assert atanh(mpc(inf,-1)).ae(-jpi2)
    assert atanh(mpc(inf,0)).ae(-jpi2)
    assert atanh(mpc(inf,1)).ae(jpi2)
    assert atanh(mpc(1,inf)).ae(jpi2)
    assert atanh(mpc(0,inf)).ae(jpi2)
    assert atanh(mpc(-1,inf)).ae(jpi2)
    assert atanh(mpc(-inf,1)).ae(jpi2)
    assert atanh(mpc(-inf,0)).ae(jpi2)
    assert atanh(mpc(-inf,-1)).ae(-jpi2)
    assert atanh(mpc(-1,-inf)).ae(-jpi2)
    assert atanh(mpc(0,-inf)).ae(-jpi2)
    assert atanh(mpc(1,-inf)).ae(-jpi2)

def test_cospi_sinpi_rounding():
    mp.dps = 15
    M = 10**15
    assert 0.999 < cospi(x1, rounding='d') < 1
    assert 0.999 < cospi(x2, rounding='d') < 1
    assert 0.999 < sinpi(x3, rounding='d') < 1
    assert 0.999 < sinpi(x4, rounding='d') < 1
    assert -1 < cospi(x5, rounding='d') < -0.999
    assert -1 < cospi(x6, rounding='d') < -0.999
    assert -1 < sinpi(x7, rounding='d') < -0.999
    assert -1 < sinpi(x8, rounding='d') < -0.999
    assert (sinpi(1e-15)*M).ae(pi)
    assert (sinpi(-1e-15)*M).ae(-pi)
    assert cospi(1e-15) == 1
    assert cospi(1e-15, rounding='d') < 1

def test_functions_infs():
    assert ellipk(inf) == 0
    assert isinstance(ellipk(inf), mpc)
    assert ellipe(inf) == mpc(0,inf)
    assert isnan(erfc(nan))
    assert lambertw(inf,1).imag.ae(2*pi)
    assert lambertw(-inf,1).imag.ae(3*pi)
    assert legendre(3.5+1j,-1) == mpc(inf,inf)
    assert legendre(4.5+1j,-1) == mpc(-inf,-inf)
    assert isnan(nthroot(nan, 1))
    assert isnan(nthroot(nan, 0))
    assert isnan(nthroot(nan, -1))
    assert isnan(nthroot(inf, 0))

def test_blah():
    assert (powm1(fadd(1,1e-100,exact=True), 5)*1e100).ae(5)

def test_functions2_erfinv():
    mp.dps = 15
    assert erfinv(0) == 0
    assert erfinv(0.5).ae(0.47693627620446987338)
    assert erfinv(-0.5).ae(-0.47693627620446987338)
    assert erfinv(1) == inf
    assert erfinv(-1) == -inf
    assert erf(erfinv(0.95)).ae(0.95)
    assert erf(erfinv(0.999999999995)).ae(0.999999999995)
    assert erf(erfinv(-0.999999999995)).ae(-0.999999999995)
    mp.dps = 50
    assert erf(erfinv('0.99999999999999999999999999999995')).ae('0.99999999999999999999999999999995')
    assert erf(erfinv('0.999999999999999999999999999999995')).ae('0.999999999999999999999999999999995')
    assert erf(erfinv('-0.999999999999999999999999999999995')).ae('-0.999999999999999999999999999999995')

'''

"""Tests for the general Heun function heung.

The reference values below have 35 decimal places and were computed with
Sage and ore_algebra, for alpha = 1/2, beta = 3/4, gamma = 2/3, delta = 4/5.
The SECOND_* values are the second solution at the origin, which heung gives
through the companion formula of DLMF 31.3.5. Value checks build inputs
from integer ratios, which round the exact reference input at the tested
precision. A decimal float would stay a 53-bit value at every precision.
"""

import pytest

from mpmath import fp, heung, iv, mp, mpc, mpf, sqrt
from mpmath.functions.heun import _series

HL_REAL_A_Z1 = ('1.07141704293306409080371224895332148', '0')
SECOND_REAL_A_Z2 = ('0.72456755998404172628933645189711678', '0.29715591749708703088553718789987941')
HL_REAL_A_Z3 = ('1.38220981766075651657980157047423786', '0.54114764731791657170019805659081215')
SECOND_REAL_A_Z3 = ('1.58450926086134133694131211850893282', '1.39437717856353592467400758045977913')
HL_COMPLEX_A_Z3 = ('-0.60828177744517820856174827377769171', '-2.12603960450399068743598991162540095')
SECOND_COMPLEX_A_Z3 = ('0.51235095261207196978037517983178209', '-2.04215244149463406764329578694999548')

# heung(2, 1/3, 1/2, 3/4, 2/3, 30, 3/2 + i/16): midpoint, to 60 digits, of a
# Sage and ore_algebra enclosure of radius 1.3e-75.
LARGE_DELTA_VALUE = ('0.970717041493602805120284950820025474645044188837025885799404',
      '0.0321915516851356476498338051749911643086247272700423636963543')


def close(got, ref, ulps=64):
    # Reference strings are parsed at the current precision.
    if isinstance(ref, tuple):
        ref = mpc(*ref)
    return abs(got - ref) <= ulps * mpf(2)**(-mp.prec) * abs(ref)


def params():
    return mpf(1)/2, mpf(3)/4, mpf(2)/3, mpf(4)/5


def second_solution(a, q, alpha, beta, gamma, delta, z, **kwargs):
    # DLMF 31.3.5 through heung with exponent 1 - gamma at the origin.
    epsilon = 1 + alpha + beta - gamma - delta
    s = 1 - gamma
    return z**s * heung(a, q + s*(epsilon + a*delta), alpha + s, beta + s,
                        2 - gamma, delta, z, **kwargs)


@pytest.mark.parametrize('dps', [15, 30])
def test_reference_values(dps):
    mp.dps = dps
    third = mpf(1)/3
    a2, ac = mpf(2), mpc(mpf(1)/2, third)
    qc = mpc(third, -mpf(1)/6)
    z1, z2, z3 = mpf(1)/4, mpc(1, 1)/4, mpc(5, 1)/4
    p = params()
    assert close(heung(a2, third, *p, z1), HL_REAL_A_Z1)
    assert close(heung(a2, third, *p, z3), HL_REAL_A_Z3)
    assert close(heung(ac, qc, *p, z3), HL_COMPLEX_A_Z3)
    assert close(second_solution(a2, third, *p, z2), SECOND_REAL_A_Z2)
    assert close(second_solution(a2, third, *p, z3), SECOND_REAL_A_Z3)
    assert close(second_solution(ac, qc, *p, z3), SECOND_COMPLEX_A_Z3)


@pytest.mark.parametrize('dps', [15, 50, 80])
def test_elementary_case(dps):
    # These parameters give 2/((1-z)*sqrt(4-z)).
    mp.dps = dps
    args = (mpf(4), mpf(9)/4, mpf(3)/2, mpf(3)/2, mpf(1)/2, mpf(2))
    for z in [mpf(1)/5, mpc(3, 1)/2, mpc(-3, 2), mpc(6, -1), mpf(-7)]:
        got = heung(*args, z)
        with mp.workdps(dps + 10):
            ref = 2/((1 - z)*sqrt(4 - z))
        assert close(got, ref)


@pytest.mark.parametrize('dps', [15, 50])
def test_hypergeometric_reduction(dps):
    # epsilon = 0 and q = a*alpha*beta reduce Heun to 2F1 (DLMF 31.7.1).
    mp.dps = dps
    A, B, G = mpf(1)/3, mpf(5)/7, mpf(3)/4
    for a in [mpf(3), mpc(mpf(1)/2, mpf(7)/10), mpf(-2)]:
        for z in [mpf(1)/3, mpc(3, 14)/10, mpc(-10, 1)/4]:
            got = heung(a, a*A*B, A, B, G, A + B + 1 - G, z)
            with mp.workdps(dps + 10):
                ref = mp.hyp2f1(A, B, G, z)
            assert close(got, ref)


@pytest.mark.parametrize('dps', [15, 50])
def test_quadratic_reduction(dps):
    # a = -1, q = 0, delta = epsilon: even in z, so a 2F1 in z**2.
    mp.dps = dps
    A, B, G = mpf(1)/3, mpf(5)/7, mpf(3)/4
    for z in [mpf(1)/3, mpc(3, 14)/10, mpc(2, -1)/5]:
        got = heung(-1, 0, A, B, G, (A + B + 1 - G)/2, z)
        with mp.workdps(dps + 10):
            ref = mp.hyp2f1(A/2, B/2, (G + 1)/2, z**2)
        assert close(got, ref)


def test_polynomial_and_trivial_solutions():
    mp.dps = 15
    assert heung(2, 2, -1, 2, 1, -8, -1) == 0
    z = -1 + mpf(2)**-30
    assert close(heung(2, 2, -1, 2, 1, -8, z), 1 + z)
    z = mpc(3, 1)
    assert close(heung(2, 2, -1, 2, 1, -8, z), 1 + z)
    for z in [mpf(1)/4, mpc(5, 1)/4]:
        assert heung(2, 0, 0, *params()[1:], z) == 1


def test_precision_retry():
    # Strong cancellation along the path: 30 guard bits give a wrong
    # value, so the working precision must adapt.
    mp.dps = 15
    args = (2, mpf(1)/3, mpf(1)/2, mpf(3)/4, mpf(2)/3, 30)
    z = mpc(mpf(3)/2, mpf(1)/16)
    got = heung(*args, z)
    assert close(got, LARGE_DELTA_VALUE)
    mp.dps = 40
    # Rebuild the inputs so 1/3 and 2/3 carry 40 digits.
    args = (2, mpf(1)/3, mpf(1)/2, mpf(3)/4, mpf(2)/3, 30)
    hi = heung(*args, z)
    assert close(hi, LARGE_DELTA_VALUE)
    mp.dps = 15
    assert close(got, hi)


@pytest.mark.parametrize('dps', [15, 30])
def test_small_value_near_polynomial_root(dps):
    # 1 + z is the exact solution; just off the root z = -1 the value is
    # small but nonzero, and must not be mistaken for rounding noise.
    mp.dps = dps
    z = mp.fadd(-1, mpf(2)**-(mp.prec + 40), exact=True)
    got = heung(2, 2, -1, 2, 1, -8, z)
    assert got != 0
    assert close(got, 1 + z)


def test_exact_zero_near_singular_point_may_not_converge():
    # H = 1 - z/z0 is an exact solution with H(z0) = 0, but the path passes
    # close to z = 1, where rounding errors grow faster than the noise
    # estimate, so heung raises instead of returning 0.
    mp.dps = 15
    z0 = mpc(1, mpf(1)/8)
    with pytest.raises(mp.NoConvergence):
        heung(2*z0, -2, -1, mpc(42, -160), 1, 20, z0, maxprec=700)


def test_step_below_working_precision_retries():
    # z is closer to the singular point 1 than the first working precision
    # resolves, so a step cannot advance and the precision must increase.
    mp.dps = 15
    with mp.workprec(200):
        z = 1 - mpf(2)**-88
    v = heung(4, mpf(9)/4, mpf(3)/2, mpf(3)/2, mpf(1)/2, 2, z)
    with mp.workprec(300):
        ref = 2/((1 - z)*sqrt(4 - z))
    assert abs(v - ref) <= 64*2**-mp.prec*abs(ref)


def test_symmetries():
    mp.dps = 15
    a, q = mpc(mpf(1)/2, mpf(1)/3), mpc(mpf(1)/3, -mpf(1)/6)
    z = mpc(5, 1)/4
    alpha, beta, gamma, delta = params()
    v = heung(a, q, alpha, beta, gamma, delta, z)
    assert close(heung(a, q, beta, alpha, gamma, delta, z), v)
    c = [mp.conj(x) for x in (a, q, alpha, beta, gamma, delta, z)]
    assert close(heung(*c), v.conjugate())


@pytest.mark.parametrize('args', [
    (0, 1, 1, 1, 1, 1, 0.5),                  # a = 0
    (1, 1, 1, 1, 1, 1, 0.5),                  # a = 1
    (0, 2, 0, 1, 1, 1, 0),                    # a = 0 before the z = 0 shortcut
    (2, 1, 1, 1, 0, 1, 0.5),                  # gamma = 0
    (2, 1, 1, 1, -1, 1, 0.5),
    (2, 1, 1, 1, -2, 1, 0.5),
    (2, 0, 0, 1, -2, 1, 0.5),                 # pole cancels (H = 1), still rejected
    (2, 1, 1, 1, 1, 1, 1),                    # z = 1
    (2, 1, 1, 1, 1, 1, 2),                    # z = a
    (2, 1, 1, 1, 1, 1, 1.5),                  # on the ray [1, inf), before a
    (2, 1, 1, 1, 1, 1, 7),
    (-2, 1, 1, 1, 1, 1, -2),
    (-2, 1, 1, 1, 1, 1, -3),
    (2, 1, 1, 1, 1, 1, float('inf')),
    (2, float('nan'), 1, 1, 1, 1, 0.5),
    (2, 1, 1, 1, 1, float('inf'), 0.5),
])
def test_domain_errors(args):
    mp.dps = 15
    with pytest.raises(ValueError):
        heung(*args)


def test_domain_errors_complex_ray():
    # z = 2*a is exactly on the ray {t*a : t >= 1}.
    mp.dps = 15
    a = mpc(mpf(1)/2, mpf(1)/3)
    with pytest.raises(ValueError):
        heung(a, 1, 1, 1, 1, 1, 2*a)
    with pytest.raises(ValueError):
        heung(a, 1, 1, 1, 1, 1, a)
    # One unit in the last place off the ray is accepted (the unit in the
    # last place of 2*a.imag = 2/3 is 2**-53).
    z = mpc(2*a.real, 2*a.imag + mpf(2)**-53)
    assert abs(heung(a, 1, 1, 1, 1, 1, z)) > 0


def test_near_pole_gamma_accepted():
    mp.dps = 15
    v = heung(2, 1, mpf(1)/2, mpf(1)/3, -2 + mpf(2)**-40, 1, mpf(1)/4)
    assert mp.isfinite(v)


def test_zero_argument():
    mp.dps = 15
    assert heung(2, 1, 1, 1, 1, 1, 0) == 1
    assert heung(2, 1, 1, 1, 1, 1, mpc(0)) == 1


def test_context_and_types():
    mp.dps = 15
    prec = mp.prec
    r = heung(2, mpf(1)/3, *params(), mpf(1)/4)
    assert mp.prec == prec
    assert type(r) is mpf
    c = heung(2, mpf(1)/3, *params(), mpc(1, 1)/4)
    assert type(c) is mpc
    with pytest.raises(ValueError):
        heung(2, 1, 1, 1, 1, 1, 1)
    assert mp.prec == prec
    with pytest.raises(mp.NoConvergence):
        heung(2, mpf(1)/3, *params(), mpf(1)/4, maxprec=60)
    assert mp.prec == prec

    args = (2, mpf(1)/3, *params(), mpf(1)/4)
    ctx = mp.clone()
    ctx.prec = 100
    v = ctx.heung(*args)
    assert mp.prec == prec
    assert ctx.prec == 100
    mp.dps = 40
    ref = heung(*args)
    mp.dps = 15
    assert abs(v - ref) <= 64 * mpf(2)**-100 * abs(ref)
    # v is rounded to 100 bits, so it differs from the 136-bit reference.
    assert abs(v - ref) > 0


def test_maxterms_limit():
    mp.dps = 15
    with pytest.raises(mp.NoConvergence):
        heung(2, mpf(1)/3, *params(), mpf(1)/4, maxterms=3)
    assert mp.prec == 53


def test_series_never_sums_beyond_maxterms():
    # The three tiny terms after 1/2 would end the series at index 4, one
    # above maxterms, so the series must raise before summing index 4.
    mp.dps = 15
    terms = iter([mpf(1)/2] + [mpf(2)**-200]*3)
    with pytest.raises(mp.NoConvergence):
        _series(mp, terms, 0, mpf(1), mpf(0), mpf(2)**-60, 3, 0, 3)


@pytest.mark.parametrize('ctx', [fp, iv])
def test_other_contexts_not_implemented(ctx):
    with pytest.raises(NotImplementedError):
        ctx.heung(2, 1/3, 1/2, 3/4, 2/3, 4/5, 1/4)


def test_mpc_input_with_zero_imaginary_part_returns_mpc():
    mp.dps = 15
    args = (2, mpf(1)/3, *params())
    for z in [mpc(mpf(1)/4, 0), mpc(0)]:
        assert type(heung(*args, z)) is mpc
    assert close(heung(*args, mpc(mpf(1)/4, 0)), heung(*args, mpf(1)/4))

from itertools import product

import pytest

from mpmath import diff, exp, j, jtheta, mp, pi, rtheta, rtheta_jet
from mpmath.functions.riemann_theta import (
    _apply_reduction, _derivative_reduction_threshold, _matrix_tuple,
    _ellipsoid_rows, _multiindices, _partial_inversion, _point_estimate,
    _rtheta_derivatives, _rtheta_sum, _shift_norm_bucket, _transform_vector,
    _upper_gamma_half_integer_bounds,
)


def test_rtheta_genus_one_jtheta_characteristics():
    # Classical-theta correspondence: DLMF 21.2.8-21.2.12.
    # https://dlmf.nist.gov/21.2.iii
    mp.dps = 30
    w = mp.mpc('0.3', '0.1')
    tau = mp.mpc('0.2', '0.9')
    q = exp(pi * j * tau)
    z = w / pi
    cases = [
        (([0.5], [0.5]), -jtheta(1, w, q)),
        (([0.5], [0]), jtheta(2, w, q)),
        (([0], [0]), jtheta(3, w, q)),
        (([0], [0.5]), jtheta(4, w, q)),
    ]
    for characteristic, expected in cases:
        assert mp.almosteq(rtheta([z], [[tau]], characteristic), expected)


def test_rtheta_genus_one_absolute_accuracy_near_theta1_zero():
    # Theta[1/2, 1/2](0 | tau) = 0 by DLMF 21.3.6.  DLMF 21.2.8
    # identifies it with -theta1(pi*z, q), providing an independent check
    # both at the zero and where relative accuracy is lost near it.
    # https://dlmf.nist.gov/21.2.E8
    # https://dlmf.nist.gov/21.3.E6
    mp.dps = 30
    tau = mp.mpc('0.17', '0.003')
    characteristic = ((0.5,), (0.5,))
    arguments = ('0', '1e-5', '1e-20', '1e-40')
    values = [
        rtheta([mp.mpf(w) / pi], [[tau]], characteristic)
        for w in arguments
    ]

    with mp.workdps(80):
        reference_tau = mp.mpc('0.17', '0.003')
        q = exp(pi * j * reference_tau)
        references = [-jtheta(1, mp.mpf(w), q) for w in arguments]

    for value, reference in zip(values, references):
        assert abs(value - reference) < 10 * mp.eps


def test_rtheta_diagonal_factorisation():
    # A diagonal period matrix separates the defining sum in DLMF 21.2.1.
    # https://dlmf.nist.gov/21.2.E1
    mp.dps = 30
    z = [mp.mpc('0.1', '0.03'), mp.mpc('-0.2', '0.04')]
    tau = [[mp.mpc('0.1', '0.8'), 0], [0, mp.mpc('-0.2', '1.1')]]
    value = rtheta(z, tau)
    expected = (rtheta([z[0]], [[tau[0][0]]])
                * rtheta([z[1]], [[tau[1][1]]]))
    assert mp.almosteq(value, expected)


def test_rtheta_block_diagonal_factorisation_with_derivative():
    # Block factorisation follows by separating the sum in DLMF 21.2.1;
    # differentiating the independent factors gives the asserted product.
    # https://dlmf.nist.gov/21.2.E1
    mp.dps = 35
    z = [mp.mpc('0.13', '0.02'), mp.mpc('-0.08', '0.03'),
         mp.mpc('0.06', '-0.01')]
    tau = [
        [mp.mpc('0.15', '0.85'), 0, 0],
        [0, mp.mpc('0.10', '0.95'), mp.mpc('-0.07', '0.04')],
        [0, mp.mpc('-0.07', '0.04'), mp.mpc('-0.12', '1.15')],
    ]
    characteristic = (
        [0.5, 0, 0.5],
        [0, 0.5, 0.5],
    )
    derivative = (1, 0, 1)
    value = rtheta(z, tau, characteristic, derivative)
    expected = rtheta(
        z[:1], [[tau[0][0]]],
        ([characteristic[0][0]], [characteristic[1][0]]), 1,
    ) * rtheta(
        z[1:], [row[1:] for row in tau[1:]],
        (characteristic[0][1:], characteristic[1][1:]), (0, 1),
    )
    assert mp.almosteq(value, expected)


def test_rtheta_parity_and_characteristic_zero():
    # Zero- and half-characteristic parity: DLMF 21.3.1 and 21.3.6.
    # https://dlmf.nist.gov/21.3.E1
    # https://dlmf.nist.gov/21.3.E6
    mp.dps = 30
    tau = [[1j, mp.mpc('0.1', '0.05')],
           [mp.mpc('0.1', '0.05'), mp.mpc('0.2', '1.2')]]
    z = [mp.mpc('0.13', '0.02'), mp.mpc('-0.07', '0.01')]
    assert mp.almosteq(rtheta(z, tau), rtheta([-z[0], -z[1]], tau))
    odd = ([0.5, 0], [0.5, 0])
    assert abs(rtheta([0, 0], tau, odd)) < mp.eps * 10


def test_rtheta_all_genus_two_half_characteristic_parities():
    # Half-characteristic parity: DLMF 21.3.6.
    # https://dlmf.nist.gov/21.3.E6
    mp.dps = 30
    half = 0.5
    tau = [[mp.mpc('0.13', '1.05'), mp.mpc('-0.09', '0.06')],
           [mp.mpc('-0.09', '0.06'), mp.mpc('-0.17', '1.2')]]
    z = [mp.mpc('0.14', '0.03'), mp.mpc('-0.11', '0.02')]
    for bits in product((0, 1), repeat=4):
        a = tuple(half * bits[i] for i in range(2))
        b = tuple(half * bits[i + 2] for i in range(2))
        sign = -1 if sum(bits[i] * bits[i + 2]
                         for i in range(2)) % 2 else 1
        characteristic = (a, b)
        assert mp.almosteq(
            rtheta([-z[0], -z[1]], tau, characteristic),
            sign * rtheta(z, tau, characteristic),
        )
        if sign == -1:
            assert abs(rtheta([0, 0], tau, characteristic)) < 100 * mp.eps


def test_rtheta_absolute_accuracy_near_odd_characteristic_zero():
    # An odd characteristic vanishes at the origin (DLMF 21.3.6).  Close to
    # that zero, test absolute rather than relative accuracy because the
    # defining sum is evaluated by cancellation.
    # https://dlmf.nist.gov/21.3.E6
    mp.dps = 30
    tau = [[mp.mpc('0.13', '1.05'), mp.mpc('-0.09', '0.06')],
           [mp.mpc('-0.09', '0.06'), mp.mpc('-0.17', '1.2')]]
    characteristic = ((0.5, 0), (0.5, 0))
    offsets = ('1e-5', '1e-20', '1e-40')
    values = [
        rtheta([mp.mpf(offset), 2 * mp.mpf(offset)], tau, characteristic)
        for offset in offsets
    ]

    with mp.workdps(80):
        reference_tau = [
            [mp.mpc('0.13', '1.05'), mp.mpc('-0.09', '0.06')],
            [mp.mpc('-0.09', '0.06'), mp.mpc('-0.17', '1.2')],
        ]
        references = [
            rtheta(
                [mp.mpf(offset), 2 * mp.mpf(offset)], reference_tau,
                characteristic,
            )
            for offset in offsets
        ]

    for value, reference in zip(values, references):
        assert abs(value - reference) < 10 * mp.eps


def test_rtheta_quasiperiodicity():
    # Periodicity and quasi-periodicity: DLMF 21.3.2-21.3.3.
    # https://dlmf.nist.gov/21.3.E2
    # https://dlmf.nist.gov/21.3.E3
    mp.dps = 30
    tau = [[mp.mpc('0.1', '0.9'), mp.mpc('0.05', '0.02')],
           [mp.mpc('0.05', '0.02'), mp.mpc('-0.1', '1.1')]]
    z = [mp.mpc('0.12', '0.02'), mp.mpc('-0.08', '0.03')]
    value = rtheta(z, tau)
    assert mp.almosteq(rtheta([z[0] + 1, z[1] - 2], tau), value)

    k = [1, -1]
    shifted = [z[i] + sum(tau[i][j_] * k[j_] for j_ in range(2))
               for i in range(2)]
    quadratic = sum(k[i] * tau[i][j_] * k[j_]
                    for i in range(2) for j_ in range(2))
    linear = sum(k[i] * z[i] for i in range(2))
    multiplier = exp(-pi * j * quadratic - 2 * pi * j * linear)
    assert mp.almosteq(rtheta(shifted, tau), multiplier * value)


def test_rtheta_derivatives():
    mp.dps = 25
    tau = [[mp.mpc('0.1', '0.9'), mp.mpc('0.05', '0.02')],
           [mp.mpc('0.05', '0.02'), mp.mpc('-0.1', '1.1')]]
    z = [mp.mpc('0.12', '0.02'), mp.mpc('-0.08', '0.03')]
    dz0 = rtheta(z, tau, derivative=(1, 0))
    expected = diff(lambda value: rtheta([value, z[1]], tau), z[0])
    assert mp.almosteq(dz0, expected)

    w = mp.mpc('0.2', '0.04')
    tau1 = mp.mpc('0.15', '0.85')
    q = exp(pi * j * tau1)
    # DLMF 21.2.11 with w = pi*z, followed by the chain rule.
    # https://dlmf.nist.gov/21.2.E11
    for order in (1, 2):
        assert mp.almosteq(
            rtheta([w / pi], [[tau1]], derivative=order),
            pi ** order * jtheta(3, w, q, derivative=order))


def test_rtheta_mixed_derivative_order_independence():
    mp.dps = 25
    tau = [[mp.mpc('0.08', '0.92'), mp.mpc('-0.06', '0.04')],
           [mp.mpc('-0.06', '0.04'), mp.mpc('0.12', '1.08')]]
    z = [mp.mpc('0.13', '0.02'), mp.mpc('-0.09', '0.03')]
    actual = rtheta(z, tau, derivative=(1, 1))
    dz0_dz1 = diff(
        lambda x: diff(lambda y: rtheta([x, y], tau), z[1]), z[0])
    dz1_dz0 = diff(
        lambda y: diff(lambda x: rtheta([x, y], tau), z[0]), z[1])
    assert mp.almosteq(actual, dz0_dz1)
    assert mp.almosteq(actual, dz1_dz0)


def test_rtheta_internal_third_order_jet_matches_scalar_calls():
    mp.dps = 25
    tau = ((mp.mpc('0.08', '0.92'), mp.mpc('-0.06', '0.04')),
           (mp.mpc('-0.06', '0.04'), mp.mpc('0.12', '1.08')))
    z = (mp.mpc('0.13', '0.02'), mp.mpc('-0.09', '0.03'))
    a = (0.5, 0)
    b = (0, 0.5)
    derivatives = tuple(_multiindices(2, 3))
    values = _rtheta_derivatives(mp, z, tau, (a, b), derivatives)

    assert len(values) == 10
    for derivative, value in zip(derivatives, values):
        assert mp.almosteq(
            value, rtheta(z, tau, (a, b), derivative=derivative))


def test_rtheta_jet_public_api():
    mp.dps = 25
    tau = ((mp.mpc('0.08', '0.92'), mp.mpc('-0.06', '0.04')),
           (mp.mpc('-0.06', '0.04'), mp.mpc('0.12', '1.08')))
    z = (mp.mpc('0.13', '0.02'), mp.mpc('-0.09', '0.03'))
    characteristic = ((0.5, 0), (0, 0.5))
    derivatives = (
        (0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))

    jet = rtheta_jet(z, tau, 2, characteristic)

    assert tuple(jet) == derivatives
    for derivative in derivatives:
        assert mp.almosteq(
            jet[derivative],
            rtheta(z, tau, characteristic, derivative=derivative))


def test_rtheta_jet_order_zero_and_genus_one():
    mp.dps = 25
    z = [mp.mpc('0.1', '0.02')]
    tau = [[mp.mpc('0.2', '0.9')]]

    zero_jet = rtheta_jet(z, tau, 0)
    second_jet = rtheta_jet(z, tau, 2)

    assert tuple(zero_jet) == ((0,),)
    assert mp.almosteq(zero_jet[(0,)], rtheta(z, tau))
    assert tuple(second_jet) == ((0,), (1,), (2,))
    assert mp.almosteq(
        second_jet[(2,)], rtheta(z, tau, derivative=2))


def test_rtheta_precision_doubling_genus_three():
    mp.dps = 100
    tau = [[1.1j, 0.04j, 0.02j],
           [0.04j, 1.2j, 0.03j],
           [0.02j, 0.03j, 1.3j]]
    z = [mp.mpc('0.1', '0.02'), mp.mpc('-0.15', '0.01'),
         mp.mpc('0.07', '-0.03')]
    value = rtheta(z, tau, derivative=(1, 0, 1))
    with mp.workdps(130):
        reference = rtheta(z, tau, derivative=(1, 0, 1))
    assert mp.almosteq(value, reference, rel_eps=mp.mpf('1e-98'))


def test_rtheta_third_order_jet_with_characteristic_at_100_dps():
    mp.dps = 100
    tau = (
        (mp.mpc('0.10', '1.00'), mp.mpc('-0.12', '0.08'),
         mp.mpc('0.05', '0.04')),
        (mp.mpc('-0.12', '0.08'), mp.mpc('0.20', '1.15'),
         mp.mpc('0.08', '0.06')),
        (mp.mpc('0.05', '0.04'), mp.mpc('0.08', '0.06'),
         mp.mpc('-0.15', '0.90')),
    )
    z = (mp.mpc('0.11', '0.03'), mp.mpc('-0.07', '0.02'),
         mp.mpc('0.05', '-0.01'))
    half = 0.5
    characteristic = ((half, 0, half), (0, half, half))
    derivatives = tuple(_multiindices(3, 3))
    values = _rtheta_derivatives(
        mp, z, tau, characteristic, derivatives)
    with mp.workdps(130):
        references = _rtheta_derivatives(
            mp, z, tau, characteristic, derivatives)
    tolerance = mp.mpf('1e-98')
    for value, reference in zip(values, references):
        assert mp.almosteq(value, reference, rel_eps=tolerance,
                           abs_eps=tolerance)


def test_rtheta_cross_precision_after_reduction():
    mp.dps = 40
    tau = [[mp.mpc('0.2', '0.15'), mp.mpc('0.12', '0.03')],
           [mp.mpc('0.12', '0.03'), mp.mpc('-0.1', '0.2')]]
    z = [mp.mpc('0.13', '0.07'), mp.mpc('-0.21', '0.04')]
    characteristic = ((mp.mpf('0.3'), mp.mpf('-0.2')),
                      (mp.mpf('0.1'), mp.mpf('0.4')))
    value = rtheta(z, tau, characteristic, derivative=(1, 1))
    with mp.workdps(75):
        reference = rtheta(z, tau, characteristic, derivative=(1, 1))
    assert mp.almosteq(value, reference, rel_eps=mp.mpf('1e-38'))


def test_rtheta_truncation_radius_against_enlarged_sum():
    mp.dps = 45
    tau = ((mp.mpc('0.11', '0.9'), mp.mpc('-0.08', '0.05')),
           (mp.mpc('-0.08', '0.05'), mp.mpc('-0.14', '1.1')))
    z = (mp.mpc('0.17', '0.11'), mp.mpc('-0.12', '0.07'))
    zero = (mp.zero, mp.zero)
    tau_data = list(mp._rtheta_tau_data(tau))
    normal = _rtheta_sum(mp, z, tau, zero, zero, ((0, 0),),
                         tuple(tau_data))[0]
    tau_data[5] *= 1.5
    enlarged = _rtheta_sum(mp, z, tau, zero, zero, ((0, 0),),
                           tuple(tau_data))[0]
    assert mp.almosteq(normal, enlarged, rel_eps=mp.mpf('1e-43'))


def test_rtheta_upper_gamma_bounds_are_conservative():
    mp.dps = 50
    for x in (mp.mpf('0.75'), mp.mpf(5), mp.mpf(80)):
        bounds = _upper_gamma_half_integer_bounds(mp, 1, 8, x)
        for twice_s, bound in enumerate(bounds, start=1):
            exact = mp.gammainc(mp.mpf(twice_s) / 2, x, mp.inf)
            assert bound >= exact


def test_rtheta_derivative_radius_cache_is_small_and_reusable():
    mp.dps = 30
    cached_radius = mp._rtheta_derivative_radius
    cached_radius.cache_clear()
    tau_value = mp.mpc('0.2', '1.2')
    tau = [[tau_value]]

    rtheta([mp.mpc('0.1', '0.01')], tau, derivative=1)
    first = cached_radius.cache_info()
    rtheta([mp.mpc('0.2', '0.02')], tau, derivative=1)
    second = cached_radius.cache_info()
    assert first.misses == 1
    assert second.hits == 1
    assert second.maxsize == 32

    tau_data = mp._rtheta_tau_data(((tau_value,),))
    for index in range(40):
        cached_radius(
            1, (1,), tau_data[4], tau_data[3], mp.mpf(index + 1) / 8)
    assert cached_radius.cache_info().currsize == 32

    for shift_norm in (mp.zero, mp.mpf('0.01'), mp.mpf('0.99'),
                       mp.mpf('12.345')):
        assert _shift_norm_bucket(mp, shift_norm) >= shift_norm
    cached_radius.cache_clear()


def test_rtheta_value_fast_path_matches_generic_sum():
    mp.dps = 45
    tau = ((mp.mpc('0.11', '0.9'), mp.mpc('-0.08', '0.05')),
           (mp.mpc('-0.08', '0.05'), mp.mpc('-0.14', '1.1')))
    z = (mp.mpc('0.17', '0.11'), mp.mpc('-0.12', '0.07'))
    a = (mp.mpf('0.3'), mp.mpf('-0.2'))
    b = (mp.mpf('0.1'), mp.mpf('0.4'))
    tau_data = mp._rtheta_tau_data(tau)
    fast = _rtheta_sum(mp, z, tau, a, b, ((0, 0),), tau_data)[0]
    generic = _rtheta_sum(
        mp, z, tau, a, b, ((0, 0), (0, 0)), tau_data)[0]
    assert mp.almosteq(fast, generic, rel_eps=mp.mpf('1e-43'))


def test_ellipsoid_rows_match_brute_force_points():
    mp.dps = 30
    cases = (
        (((mp.mpf('1.2'),),), (mp.mpf('0.3'),), mp.mpf('2.1')),
        (((mp.mpf('1.1'), mp.mpf('0.2'), mp.mpf('-0.1')),
          (mp.zero, mp.mpf('0.9'), mp.mpf('0.15')),
          (mp.zero, mp.zero, mp.mpf('1.3'))),
         (mp.mpf('0.2'), mp.mpf('-0.4'), mp.mpf('0.1')), mp.mpf('2.4')),
    )
    for T, center, radius in cases:
        row_points = [
            (value,) + outer
            for outer, lower, upper in _ellipsoid_rows(
                mp, T, center, radius)
            for value in range(lower, upper + 1)
        ]
        genus = len(center)
        brute_force = []
        for point in product(range(-4, 5), repeat=genus):
            norm_squared = mp.fsum(
                mp.fsum(T[i][j] * (point[j] - center[j])
                        for j in range(genus)) ** 2
                for i in range(genus))
            if norm_squared <= radius ** 2:
                brute_force.append(point)
        brute_force.sort(key=lambda point: tuple(reversed(point)))
        assert row_points == brute_force


def test_rtheta_wolfram_reference_values():
    # Wolfram Engine 14.3 SiegelTheta values generated at 120 digits from
    # exact rational inputs. The tests use 100 digits and retain guard digits
    # in each reference value. The calls were N[SiegelTheta[tau, z], 120]
    # and N[SiegelTheta[{a, b}, tau, z], 120].
    # https://reference.wolfram.com/language/ref/SiegelTheta.html
    mp.dps = 100
    z2 = [mp.mpc('0.11', '0.03'), mp.mpc('-0.07', '0.02')]
    tau2 = [
        [mp.mpc('0.10', '0.90'), mp.mpc('-0.20', '0.12')],
        [mp.mpc('-0.20', '0.12'), mp.mpc('0.05', '1.10')],
    ]
    char2 = ([0.5, 0], [0, 0.5])
    refs2 = [
        mp.mpc(
            '1.15010054896589710068630628822243271125317054007702426004984068418064969671373665362574581188457525423704174753913082529848244697358328',
            '0.029482887366211288410342712159136306372758610990096207370983577503724607487203127164966369240820605787739356900521381741214906114952103'),
        mp.mpc(
            '0.893359834184414572430609669893702469378276216582063538344544618948879436225815023389107207488871647434926394149016401840086595696537675',
            '0.02505879729520504089469665867573950858628706160028878239473482035693999039117935887428157275486164770342340639145800073069780958643882'),
    ]

    z3 = z2 + [mp.mpc('0.05', '-0.01')]
    tau3 = [
        [mp.mpc('0.10', '1.00'), mp.mpc('-0.12', '0.08'),
         mp.mpc('0.05', '0.04')],
        [mp.mpc('-0.12', '0.08'), mp.mpc('0.20', '1.15'),
         mp.mpc('0.08', '0.06')],
        [mp.mpc('0.05', '0.04'), mp.mpc('0.08', '0.06'),
         mp.mpc('-0.15', '0.90')],
    ]
    char3 = ([0.5, 0, 0.5], [0, 0, 0.5])
    refs3 = [
        mp.mpc(
            '1.221609642406339797657102259903691954365704088043749183987260288160622810706259133711255570012034796050287767505312762580074107266344458',
            '-0.007792386568554659454892241427121689475136799802090774325255515968433879603115006454848761576475120048373924432677064251983242053776354'),
        mp.mpc(
            '-0.114796433952467492469811475212021962597569270860961487827600713071632671744755401699850275301415484731581776951960982866615787356410398',
            '0.015318780171864079253785680566109776307328372396748562810674785685663933968972735107125484955849153511806284150913864802098269822109827'),
    ]

    tolerance = mp.mpf('1e-98')
    for z, tau, characteristic, references in [
            (z2, tau2, None, refs2[:1]),
            (z2, tau2, char2, refs2[1:]),
            (z3, tau3, None, refs3[:1]),
            (z3, tau3, char3, refs3[1:])]:
        assert mp.almosteq(rtheta(z, tau, characteristic), references[0],
                           rel_eps=tolerance, abs_eps=tolerance)


def test_rtheta_wolfram_unreduced_reference_values():
    # Fresh Wolfram Engine 14.3 SiegelTheta values generated at 120 digits
    # from exact rational inputs. These exercise reduction, large z and
    # arbitrary real characteristics together. The calls used the same two
    # forms of SiegelTheta recorded in the preceding reference-value test.
    # https://reference.wolfram.com/language/ref/SiegelTheta.html
    mp.dps = 100
    characteristic2 = ((mp.mpf('0.3'), mp.mpf('-0.2')),
                       (mp.mpf('0.1'), mp.mpf('0.4')))
    tau2 = [[mp.mpc('0.2', '0.15'), mp.mpc('0.12', '0.03')],
            [mp.mpc('0.12', '0.03'), mp.mpc('-0.1', '0.2')]]
    z2 = [mp.mpc('2.31', '-1.14'), mp.mpc('-1.72', '0.83')]
    reference2 = mp.mpc(
        '37058424061976446096.7871409113145979695568610946774019090331353765326092054797218949660046676284791419537804054179813145',
        '-8184623208069837737.99363797667810872646677513420698882738727935055524579125490740124342217355390271095127290985384489359')

    characteristic3 = ((mp.mpf('0.3'), mp.mpf('-0.2'), mp.mpf('0.1')),
                       (mp.mpf('0.1'), mp.mpf('0.4'), mp.mpf('-0.15')))
    tau3 = [
        [mp.mpc('0.3', '0.12'), mp.mpc('0.17', '0.03'),
         mp.mpc(0, '0.02')],
        [mp.mpc('0.17', '0.03'), mp.mpc('-0.2', '0.16'),
         mp.mpc(0, '0.02')],
        [mp.mpc(0, '0.02'), mp.mpc(0, '0.02'),
         mp.mpc('0.1', '0.2')],
    ]
    z3 = [mp.mpc('0.11', '0.03'), mp.mpc('-0.07', '0.02'),
          mp.mpc('0.05', '-0.01')]
    reference3 = mp.mpc(
        '3.44042050000036293079504325661815829014001245512124516294720089024113963744406088834494220278232680204967316690887511458',
        '0.182598566348364806829149624629842477324789098204541555796804418499811344983007460195610367554313435739194325595031041258')

    tolerance = mp.mpf('1e-98')
    assert mp.almosteq(rtheta(z2, tau2, characteristic2), reference2,
                       rel_eps=tolerance, abs_eps=tolerance)
    assert mp.almosteq(rtheta(z3, tau3, characteristic3), reference3,
                       rel_eps=tolerance, abs_eps=tolerance)


def test_rtheta_wolfram_first_derivative_reference_values():
    # Wolfram Engine 14.3 NumericalCalculus`ND values, evaluated using
    # Method -> NIntegrate (Cauchy's integral formula) and 65-70 digits of
    # working precision. SiegelTheta has the same normalized z convention.
    # https://reference.wolfram.com/language/NumericalCalculus/ref/ND.html
    # https://reference.wolfram.com/language/ref/SiegelTheta.html
    mp.dps = 40
    z2 = [mp.mpc('0.11', '0.03'), mp.mpc('-0.07', '0.02')]
    tau2 = [
        [mp.mpc('0.10', '0.90'), mp.mpc('-0.20', '0.12')],
        [mp.mpc('-0.20', '0.12'), mp.mpc('0.05', '1.10')],
    ]
    characteristic2 = ([0.5, 0], [0, 0.5])
    references2 = (
        (None, mp.mpc(
            '-0.421288315640741560394787218172215298435192461898651',
            '-0.297478129127372247972205522599835990950319316362461')),
        (characteristic2, mp.mpc(
            '-0.983740135148338116134552482066670136613333239537688',
            '-0.287732601392744703482773208619225173728089061045940')),
    )
    tolerance = mp.mpf('1e-38')
    for characteristic, reference in references2:
        assert mp.almosteq(
            rtheta(z2, tau2, characteristic, derivative=(1, 0)),
            reference, rel_eps=tolerance, abs_eps=tolerance)

    z3 = [mp.mpc('0.11', '0.03'), mp.mpc('-0.07', '0.02'),
          mp.mpc('0.05', '-0.01')]
    tau3 = [
        [mp.mpc('0.10', '1.00'), mp.mpc('-0.12', '0.08'),
         mp.mpc('0.05', '0.04')],
        [mp.mpc('-0.12', '0.08'), mp.mpc('0.20', '1.15'),
         mp.mpc('0.08', '0.06')],
        [mp.mpc('0.05', '0.04'), mp.mpc('0.08', '0.06'),
         mp.mpc('-0.15', '0.90')],
    ]
    reference3 = mp.mpc(
        '-0.349364739121920510015037946041375468350822979142697',
        '-0.218485727335759004831470672865642569157395969625689')
    assert mp.almosteq(
        rtheta(z3, tau3, derivative=(1, 0, 0)), reference3,
        rel_eps=tolerance, abs_eps=tolerance)


def test_rtheta_flint_high_precision_derivative_reference_values():
    # Certified FLINT 3.6.0 acb_theta_jet enclosures computed at 130 dps
    # through Python-FLINT 0.9.0. FLINT's Taylor coefficients were multiplied
    # by the multi-index factorials to obtain ordinary partial derivatives.
    # Every recorded midpoint had enclosure radius below 3.3e-128.
    # https://flintlib.org/doc/acb_theta.html
    mp.dps = 100
    tolerance = mp.mpf('1e-98')
    z2 = [mp.mpc('0.11', '0.03'), mp.mpc('-0.07', '0.02')]
    tau2 = [
        [mp.mpc('0.10', '0.90'), mp.mpc('-0.20', '0.12')],
        [mp.mpc('-0.20', '0.12'), mp.mpc('0.05', '1.10')],
    ]
    half2 = ([0.5, 0], [0, 0.5])
    cases2 = (
        (None, (2, 0), mp.mpc(
            '-3.70302318847962531102746421196741770479519169317250238327399845194544812073725950034281322335010842942379904608583600925',
            '-0.676647983025765584988696701501531718662199744767527809377309914862400209905506350145714074349651391583245694053251869978')),
        (None, (1, 1), mp.mpc(
            '-0.0492830462579949939080327254030849463068181598703018850653065257700124493596665111519633528138121999674893666903240653057',
            '0.188722819532323736052517708362276883159628967942584955480523117746601772368091761871289179463169032107921550117333179792')),
        (None, (2, 1), mp.mpc(
            '0.510109942412506906348039968829021056110472011913214068324156555703573356310679474090641891489414592720300153907825795172',
            '-1.73784821436088196692553878274949762796692250748447504881253235541903673215233093231020937149976587228035056900176806436')),
        (half2, (0, 2), mp.mpc(
            '1.59405366912712602932168341925929329739061553348406695029834721987354070660484011542785642549687401680165762911354469987',
            '0.678197572830669325431970259193745616335781067306349265494919453231503185934986593023788155173334269512663008718046546465')),
        (half2, (1, 1), mp.mpc(
            '0.00506262711244847595711909837342906731732086829953077179875073877425155230227683339892785473221686111182635405587903951584',
            '-0.656735703386161600190210786006010705331432797990974033139682794285353147435261926555954518014351969442099487008921964491')),
        (half2, (1, 2), mp.mpc(
            '-2.56012945425409603208529451407731106831476911040452819698668383790292721185892340100183555947273879000517313127472716906',
            '-3.62621222740679372224268928087397123446120699815806027482289280696049350547713834880575823432416970384201643649467346182')),
    )
    for characteristic, derivative, reference in cases2:
        assert mp.almosteq(
            rtheta(z2, tau2, characteristic, derivative), reference,
            rel_eps=tolerance, abs_eps=tolerance)

    z3 = z2 + [mp.mpc('0.05', '-0.01')]
    tau3 = [
        [mp.mpc('0.10', '1.00'), mp.mpc('-0.12', '0.08'),
         mp.mpc('0.05', '0.04')],
        [mp.mpc('-0.12', '0.08'), mp.mpc('0.20', '1.15'),
         mp.mpc('0.08', '0.06')],
        [mp.mpc('0.05', '0.04'), mp.mpc('0.08', '0.06'),
         mp.mpc('-0.15', '0.90')],
    ]
    half3 = ([0.5, 0, 0.5], [0, 0, 0.5])
    cases3 = (
        (None, (1, 1, 0), mp.mpc(
            '-0.0583203441458341483463356071834041239984343441778997778365789034110166610739832774214620984517732788372057366818649126281',
            '0.0708753457491918488372972088434176429774075409403303751129258604355297050164131636926147821958543195598178704144781753017')),
        (None, (2, 0, 1), mp.mpc(
            '0.114306309858173014349559978028599658075537833007449717601194351437916360428696872814207699015024136431925580013134749270',
            '0.149072636036949478276875893311970273856306502655947226037232741110709312077642678802950638931487633597974738441036317957')),
        (half3, (0, 1, 1), mp.mpc(
            '-0.330573454408454085191250734483989258773110711237355724952771550857516516520533097128909554697460504958496922554130515797',
            '-0.248921244221115960904005325671928134409098763816549339891050443057701465746834442584948227393550587676101381266386567494')),
        (half3, (1, 1, 1), mp.mpc(
            '0.459086351057195620002222892353732574900829787440015815436363687987155908263334103197721895907348597954260934682121834009',
            '-0.924063613793591498369212162932436220909836688736846903748833925742816424588180770490266398471292926543648766005460109829')),
    )
    for characteristic, derivative, reference in cases3:
        assert mp.almosteq(
            rtheta(z3, tau3, characteristic, derivative), reference,
            rel_eps=tolerance, abs_eps=tolerance)


def test_rtheta_wolfram_higher_genus_reference_values():
    # Wolfram Engine 14.3 values from N[SiegelTheta[tau, z], 40], with
    # tau and z constructed from the exact rational counterparts below.
    # https://reference.wolfram.com/language/ref/SiegelTheta.html
    # Strong damping keeps these high-dimensional smoke tests bounded; they
    # do not imply that arbitrary genus-10 or genus-20 inputs are inexpensive.
    mp.dps = 35

    def make_case(genus, scale):
        tau = [[0 for unused in range(genus)] for unused in range(genus)]
        for i in range(genus):
            tau[i][i] = mp.mpc(
                mp.mpf(i % 3) / 20,
                scale + mp.mpf(i % 4) / 10,
            )
            if i + 1 < genus:
                tau[i][i + 1] = tau[i + 1][i] = mp.mpc('0.02', '0.25')
        z = [
            mp.mpc(mp.mpf((i % 5) - 2) / 100,
                   mp.mpf((i % 3) - 1) / 200)
            for i in range(genus)
        ]
        return z, tau

    cases = (
        (10, 8, mp.mpc(
            '1.0000000001671156401975485245765975730135867034602892301525',
            '2.49857014987438782428182159399172658122217981149e-11')),
        (20, 20, mp.mpc(
            '1.0000000000000000000000000133796107582295982965475652944874',
            '2.0101530593695935374230689609606e-27')),
    )
    tolerance = mp.mpf('1e-33')
    for genus, scale, reference in cases:
        z, tau = make_case(genus, scale)
        assert mp.almosteq(rtheta(z, tau), reference,
                           rel_eps=tolerance, abs_eps=tolerance)


def test_rtheta_validation():
    mp.dps = 20
    with pytest.raises(ValueError, match="vector"):
        rtheta(None, [[1j]])
    with pytest.raises(ValueError, match="vector"):
        rtheta(mp.matrix([[0, 0]]), [[1j]])
    with pytest.raises(ValueError, match="finite"):
        rtheta([mp.inf], [[1j]])
    with pytest.raises(ValueError, match="square"):
        rtheta([0], None)
    with pytest.raises(ValueError, match="square"):
        rtheta([0], [[1j, 0]])
    with pytest.raises(ValueError, match="finite"):
        rtheta([0], [[mp.inf + 1j]])
    with pytest.raises(ValueError, match="length"):
        rtheta([0], [[1j, 0], [0, 1j]])
    with pytest.raises(ValueError, match="symmetric"):
        rtheta([0, 0], [[1j, 1], [0, 1j]])
    with pytest.raises(ValueError, match="positive definite"):
        rtheta([0, 0], [[1j, 0], [0, -1j]])
    with pytest.raises(ValueError, match="real"):
        rtheta([0], [[1j]], characteristic=([1j], [0]))
    with pytest.raises(ValueError, match="pair"):
        rtheta([0], [[1j]], characteristic=([0],))
    with pytest.raises(ValueError, match="length"):
        rtheta([0, 0], [[1j, 0], [0, 1j]], derivative=(1,))
    with pytest.raises(ValueError, match="genus 1"):
        rtheta([0, 0], [[1j, 0], [0, 1j]], derivative=1)
    with pytest.raises(ValueError, match="multi-index"):
        rtheta([0], [[1j]], derivative=None)
    with pytest.raises(ValueError, match="nonnegative"):
        rtheta([0], [[1j]], derivative=(1.0,))
    with pytest.raises(ValueError, match="nonnegative"):
        rtheta([0], [[1j]], derivative=(-1,))
    with pytest.raises(ValueError, match="nonempty sequence"):
        _rtheta_derivatives(mp, [0], [[1j]], None, None)
    with pytest.raises(ValueError, match="nonempty sequence"):
        _rtheta_derivatives(mp, [0], [[1j]], None, ())
    with pytest.raises(ValueError, match="nonnegative integer"):
        rtheta_jet([0], [[1j]], -1)
    with pytest.raises(ValueError, match="nonnegative integer"):
        rtheta_jet([0], [[1j]], 1.0)


def test_rtheta_nearly_symmetric_input():
    mp.dps = 20
    delta = mp.eps / 4
    tau = [[1j, mp.mpc('0.1', '0.02')],
           [mp.mpc('0.1', '0.02') + delta, 1.2j]]
    symmetric = [[1j, mp.mpc('0.1', '0.02') + delta / 2],
                 [mp.mpc('0.1', '0.02') + delta / 2, 1.2j]]
    assert mp.almosteq(rtheta([0.1, 0.2], tau),
                       rtheta([0.1, 0.2], symmetric))


def test_rtheta_symplectic_generator_identities():
    # The basis, translation and general modular transformations are DLMF
    # 21.5.5, 21.5.7 and 21.5.4, respectively.
    # https://dlmf.nist.gov/21.5.E5
    # https://dlmf.nist.gov/21.5.E7
    # https://dlmf.nist.gov/21.5.E4
    mp.dps = 35
    tau = mp.matrix([
        [mp.mpc('0.13', '1.1'), mp.mpc('-0.17', '0.08')],
        [mp.mpc('-0.17', '0.08'), mp.mpc('0.21', '1.3')],
    ])
    z = (mp.mpc('0.12', '0.04'), mp.mpc('-0.09', '0.03'))
    value = rtheta(z, tau)

    translation = ((1, 2), (2, -1))
    translated_tau = mp.matrix(tau)
    for row in range(2):
        for column in range(2):
            translated_tau[row, column] += translation[row][column]
    translated_z = (z[0] - 0.5, z[1] + 0.5)
    assert mp.almosteq(rtheta(translated_z, translated_tau), value)

    basis = ((1, 1), (0, 1))
    basis_tau = mp.matrix(2)
    for row in range(2):
        for column in range(2):
            basis_tau[row, column] = mp.fsum(
                basis[k][row] * tau[k, l] * basis[l][column]
                for k in range(2) for l in range(2))
    basis_z = _transform_vector(mp, basis, z)
    assert mp.almosteq(rtheta(basis_z, basis_tau), value)

    combined_z, combined_factor = _apply_reduction(
        mp, z, (("translate", (1, -1)), ("basis", basis)))
    assert combined_z == _transform_vector(
        mp, basis, (z[0] + 0.5, z[1] - 0.5))
    assert combined_factor == 1

    inverted_tau, inversion = _partial_inversion(mp, tau)
    inverted_z, factor = _apply_reduction(mp, z, (("invert", inversion),))
    assert mp.almosteq(factor * rtheta(inverted_z, inverted_tau), value)


def test_rtheta_cost_selected_reduction():
    mp.dps = 35
    tau = (
        (mp.mpc('0.2', '0.15'), mp.mpc('0.12', '0.03')),
        (mp.mpc('0.12', '0.03'), mp.mpc('-0.1', '0.2')),
    )
    z = (mp.mpc('0.13', '0.07'), mp.mpc('-0.21', '0.04'))
    characteristic = ((mp.mpf('0.3'), mp.mpf('-0.2')),
                      (mp.mpf('0.1'), mp.mpf('0.4')))
    unused_reduced_tau, operations, reduced_points = (
        mp._rtheta_reduction_data(_matrix_tuple(mp.matrix(tau))))

    assert _point_estimate(mp, tau) / reduced_points > 4
    assert {kind for kind, unused in operations} == {
        "translate", "basis", "invert"
    }
    for kind, transform in operations:
        if kind == "basis":
            assert abs(transform[0][0] * transform[1][1]
                       - transform[0][1] * transform[1][0]) == 1

    a, b = characteristic
    with mp.extraprec(30):
        direct = _rtheta_sum(
            mp, z, tau, a, b, ((0, 0),)
        )[0]
    assert mp.almosteq(rtheta(z, tau, characteristic), direct)

    derivative = (1, 0)
    assert _point_estimate(mp, tau) / reduced_points > 12
    with mp.extraprec(30):
        direct_derivative = _rtheta_sum(
            mp, z, tau, a, b, (derivative,)
        )[0]
    assert mp.almosteq(
        rtheta(z, tau, characteristic, derivative), direct_derivative)

    # The final check also exercises DLMF 21.3.3 after reduction.
    # https://dlmf.nist.gov/21.3.E3
    shift = (3, -2)
    shifted_z = tuple(
        z[i] + mp.fsum(tau[i][k] * shift[k] for k in range(2))
        for i in range(2))
    quadratic = mp.fsum(shift[i] * tau[i][k] * shift[k]
                        for i in range(2) for k in range(2))
    linear = mp.fsum(shift[i] * z[i] for i in range(2))
    multiplier = exp(-pi * j * quadratic - 2 * pi * j * linear)
    assert mp.almosteq(rtheta(shifted_z, tau),
                       multiplier * rtheta(z, tau))


def test_rtheta_reduction_is_stable():
    mp.dps = 35
    tau = (
        (mp.mpc('0.2', '0.15'), mp.mpc('0.12', '0.03')),
        (mp.mpc('0.12', '0.03'), mp.mpc('-0.1', '0.2')),
    )
    reduced, operations, _ = mp._rtheta_reduction_data(tau)
    reduced_again, further_operations, _ = (
        mp._rtheta_reduction_data(reduced))
    assert operations
    assert not further_operations
    for row in range(2):
        for column in range(2):
            assert mp.almosteq(reduced[row][column],
                               reduced_again[row][column])


def test_rtheta_cost_selected_reduction_genus_three():
    mp.dps = 25
    tau = (
        (mp.mpc('0.3', '0.12'), mp.mpc('0.17', '0.03'),
         mp.mpc(0, '0.02')),
        (mp.mpc('0.17', '0.03'), mp.mpc('-0.2', '0.16'),
         mp.mpc(0, '0.02')),
        (mp.mpc(0, '0.02'), mp.mpc(0, '0.02'),
         mp.mpc('0.1', '0.2')),
    )
    z = (mp.mpc('0.13', '0.07'), mp.mpc('-0.21', '0.04'),
         mp.mpc('0.06', '-0.03'))
    zero = (mp.zero,) * 3
    reduced_tau, operations, reduced_points = mp._rtheta_reduction_data(tau)

    assert _point_estimate(mp, tau) / reduced_points > 4
    assert len(reduced_tau) == 3
    assert any(kind == "basis" for kind, unused in operations)
    for kind, transform in operations:
        if kind == "basis":
            assert abs(mp.det(mp.matrix(transform))) == 1

    with mp.extraprec(30):
        direct = _rtheta_sum(
            mp, z, tau, zero, zero, ((0, 0, 0),)
        )[0]
    assert mp.almosteq(rtheta(z, tau), direct)


def test_rtheta_cost_selected_genus_three_second_derivative():
    mp.dps = 20
    tau = (
        (mp.mpc('0.3', '0.12'), mp.mpc('0.17', '0.03'),
         mp.mpc(0, '0.02')),
        (mp.mpc('0.17', '0.03'), mp.mpc('-0.2', '0.16'),
         mp.mpc(0, '0.02')),
        (mp.mpc(0, '0.02'), mp.mpc(0, '0.02'),
         mp.mpc('0.1', '0.2')),
    )
    z = (mp.mpc('0.11', '0.03'), mp.mpc('-0.07', '0.02'),
         mp.mpc('0.05', '-0.01'))
    derivative = (1, 1, 0)
    zero = (mp.zero,) * 3
    unused_tau, unused_operations, reduced_points = (
        mp._rtheta_reduction_data(tau))
    point_ratio = _point_estimate(mp, tau) / reduced_points
    threshold = _derivative_reduction_threshold(
        mp, 3, 2, (derivative,))

    # The old threshold charged ten complete traversals and was 40. The
    # shared-work estimate is 13, so this reduction is no longer missed.
    assert threshold < point_ratio < 40
    with mp.extraprec(25):
        direct = _rtheta_sum(
            mp, z, tau, zero, zero, (derivative,)
        )[0]
    assert mp.almosteq(rtheta(z, tau, derivative=derivative), direct)


def test_rtheta_reduced_mixed_derivative():
    mp.dps = 25
    tau = (
        (mp.mpc('0.4', '0.03'), mp.mpc('0.1', '0.008')),
        (mp.mpc('0.1', '0.008'), mp.mpc('-0.3', '0.05')),
    )
    z = (mp.mpc('0.13', '0.07'), mp.mpc('-0.21', '0.04'))
    a = (0.5, 0)
    b = (0, 0.5)
    derivative = (1, 1)
    unused_tau, unused_operations, reduced_points = (
        mp._rtheta_reduction_data(tau))

    # Six derivatives of total degree at most two are combined, so the
    # selection threshold is four times six.
    assert _point_estimate(mp, tau) / reduced_points > 24
    with mp.extraprec(30):
        direct = _rtheta_sum(
            mp, z, tau, a, b, (derivative,)
        )[0]
    assert mp.almosteq(rtheta(z, tau, (a, b), derivative), direct)


def test_rtheta_reduced_derivative_jet_matches_direct_sum():
    mp.dps = 25
    tau = (
        (mp.mpc('0.4', '0.03'), mp.mpc('0.1', '0.008')),
        (mp.mpc('0.1', '0.008'), mp.mpc('-0.3', '0.05')),
    )
    z = (mp.mpc('0.13', '0.07'), mp.mpc('-0.21', '0.04'))
    a = (0.5, 0)
    b = (0, 0.5)
    derivatives = tuple(_multiindices(2, 2))

    values = _rtheta_derivatives(mp, z, tau, (a, b), derivatives)
    with mp.extraprec(30):
        direct = _rtheta_sum(mp, z, tau, a, b, derivatives)
    for value, reference in zip(values, direct):
        assert mp.almosteq(value, reference)


def test_rtheta_unreduced_genus_one_against_jtheta():
    # Classical-theta correspondence: DLMF 21.2.8-21.2.12.
    # https://dlmf.nist.gov/21.2.iii
    mp.dps = 35
    tau = mp.mpc('0.1', '0.1')
    w = mp.mpc('0.3', '0.17')
    q = exp(pi * j * tau)
    cases = [
        (([0.5], [0.5]), -jtheta(1, w, q)),
        (([0.5], [0]), jtheta(2, w, q)),
        (([0], [0]), jtheta(3, w, q)),
        (([0], [0.5]), jtheta(4, w, q)),
    ]
    for characteristic, expected in cases:
        assert mp.almosteq(rtheta([w / pi], [[tau]], characteristic),
                           expected)

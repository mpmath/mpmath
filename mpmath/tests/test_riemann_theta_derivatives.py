"""Tests for Riemann theta derivatives and jets."""

import pytest

from mpmath import diff, exp, j, jtheta, mp, pi, rtheta, rtheta_jet
from mpmath.functions.riemann_theta import (
    _derivative_reduction_threshold, _matrix_tuple, _multiindices,
    _point_estimate, _reduced_jet_data, _rtheta_derivatives, _rtheta_sum,
    _shift_norm_bucket,
)


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


def test_rtheta_jet_reduction_cache_reuse_and_precision():
    with mp.workdps(30):
        cached = mp._rtheta_jet_reduction_data
        cached.cache_clear()
        tau = ((mp.mpc('0.2', '0.15'), mp.mpc('0.12', '0.03')),
               (mp.mpc('0.12', '0.03'), mp.mpc('-0.1', '0.2')))
        reduced, operations, _ = mp._rtheta_reduction_data(tau)
        steps, jacobian, quadratic = cached(2, operations)
        bases = tuple(_matrix_tuple(data) for kind, data in steps
                      if kind == "basis")
        assert bases
        assert quadratic
        z = (mp.mpc('0.13', '0.07'), mp.mpc('-0.21', '0.04'))
        zero = (mp.zero, mp.zero)
        first = _reduced_jet_data(
            mp, z, tau, zero, zero, operations, reduced)
        second = _reduced_jet_data(
            mp, (z[0] + 1, z[1] - 2), tau, (mp.mpf('0.25'), mp.zero),
            zero, operations, reduced)
        assert first[2] is second[2] is jacobian
        assert cached.cache_info().misses == 1
        assert cached.cache_info().hits == 2
        assert bases == tuple(_matrix_tuple(data) for kind, data in steps
                              if kind == "basis")
        with mp.extraprec(20):
            assert cached(2, operations)[1] is not jacobian
        assert cached.cache_info().misses == 2
        assert cached(2, operations)[1] is jacobian
        for i in range(20):
            cached(1, (("translate", (i,)),))
        assert cached.cache_info().maxsize == 16
        assert cached.cache_info().currsize == 16
        cached.cache_clear()


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


def test_rtheta_derivative_validation():
    mp.dps = 20
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


def test_rtheta_cost_selected_derivative_reduction():
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

    a, b = characteristic
    derivative = (1, 0)
    assert _point_estimate(mp, tau) / reduced_points > 12
    with mp.extraprec(30):
        direct_derivative = _rtheta_sum(
            mp, z, tau, a, b, (derivative,)
        )[0]
    assert mp.almosteq(
        rtheta(z, tau, characteristic, derivative), direct_derivative)


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

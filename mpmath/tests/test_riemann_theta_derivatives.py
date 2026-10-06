"""Tests for Riemann theta derivatives and jets."""

import pytest

from mpmath import diff, exp, j, jtheta, mp, mpc, mpf, pi, rtheta, rtheta_jet


def test_rtheta_block_diagonal_factorisation_with_derivative():
    # Block factorisation follows by separating the sum in DLMF 21.2.1;
    # differentiating the independent factors gives the asserted product.
    # https://dlmf.nist.gov/21.2.E1
    mp.dps = 35
    z = [mpc('0.13', '0.02'), mpc('-0.08', '0.03'),
         mpc('0.06', '-0.01')]
    tau = [
        [mpc('0.15', '0.85'), 0, 0],
        [0, mpc('0.10', '0.95'), mpc('-0.07', '0.04')],
        [0, mpc('-0.07', '0.04'), mpc('-0.12', '1.15')],
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
    tau = [[mpc('0.1', '0.9'), mpc('0.05', '0.02')],
           [mpc('0.05', '0.02'), mpc('-0.1', '1.1')]]
    z = [mpc('0.12', '0.02'), mpc('-0.08', '0.03')]
    dz0 = rtheta(z, tau, derivative=(1, 0))
    expected = diff(lambda value: rtheta([value, z[1]], tau), z[0])
    assert mp.almosteq(dz0, expected)

    w = mpc('0.2', '0.04')
    tau1 = mpc('0.15', '0.85')
    q = exp(pi * j * tau1)
    # DLMF 21.2.11 with w = pi*z, followed by the chain rule.
    # https://dlmf.nist.gov/21.2.E11
    for order in (1, 2):
        assert mp.almosteq(
            rtheta([w / pi], [[tau1]], derivative=order),
            pi ** order * jtheta(3, w, q, derivative=order))


def test_rtheta_mixed_derivative_order_independence():
    mp.dps = 25
    tau = [[mpc('0.08', '0.92'), mpc('-0.06', '0.04')],
           [mpc('-0.06', '0.04'), mpc('0.12', '1.08')]]
    z = [mpc('0.13', '0.02'), mpc('-0.09', '0.03')]
    actual = rtheta(z, tau, derivative=(1, 1))
    dz0_dz1 = diff(
        lambda x: diff(lambda y: rtheta([x, y], tau), z[1]), z[0])
    dz1_dz0 = diff(
        lambda y: diff(lambda x: rtheta([x, y], tau), z[0]), z[1])
    assert mp.almosteq(actual, dz0_dz1)
    assert mp.almosteq(actual, dz1_dz0)


def test_rtheta_third_order_jet_matches_scalar_calls():
    mp.dps = 25
    tau = ((mpc('0.08', '0.92'), mpc('-0.06', '0.04')),
           (mpc('-0.06', '0.04'), mpc('0.12', '1.08')))
    z = (mpc('0.13', '0.02'), mpc('-0.09', '0.03'))
    a = (0.5, 0)
    b = (0, 0.5)
    jet = rtheta_jet(z, tau, 3, (a, b))

    assert tuple(jet) == (
        (0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2),
        (3, 0), (2, 1), (1, 2), (0, 3))
    for derivative, value in jet.items():
        assert mp.almosteq(
            value, rtheta(z, tau, (a, b), derivative=derivative))


def test_rtheta_jet_public_api():
    mp.dps = 25
    tau = ((mpc('0.08', '0.92'), mpc('-0.06', '0.04')),
           (mpc('-0.06', '0.04'), mpc('0.12', '1.08')))
    z = (mpc('0.13', '0.02'), mpc('-0.09', '0.03'))
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
    z = [mpc('0.1', '0.02')]
    tau = [[mpc('0.2', '0.9')]]

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
    z = [mpc('0.1', '0.02'), mpc('-0.15', '0.01'),
         mpc('0.07', '-0.03')]
    value = rtheta(z, tau, derivative=(1, 0, 1))
    with mp.workdps(130):
        reference = rtheta(z, tau, derivative=(1, 0, 1))
    assert mp.almosteq(value, reference, rel_eps=mpf('1e-98'))


def test_rtheta_third_order_jet_with_characteristic_at_100_dps():
    mp.dps = 100
    tau = (
        (mpc('0.10', '1.00'), mpc('-0.12', '0.08'),
         mpc('0.05', '0.04')),
        (mpc('-0.12', '0.08'), mpc('0.20', '1.15'),
         mpc('0.08', '0.06')),
        (mpc('0.05', '0.04'), mpc('0.08', '0.06'),
         mpc('-0.15', '0.90')),
    )
    z = (mpc('0.11', '0.03'), mpc('-0.07', '0.02'),
         mpc('0.05', '-0.01'))
    half = 0.5
    characteristic = ((half, 0, half), (0, half, half))
    jet = rtheta_jet(z, tau, 3, characteristic)
    with mp.workdps(130):
        reference_jet = rtheta_jet(z, tau, 3, characteristic)
    tolerance = mpf('1e-98')
    assert jet.keys() == reference_jet.keys()
    for derivative, value in jet.items():
        assert mp.almosteq(value, reference_jet[derivative], rel_eps=tolerance,
                           abs_eps=tolerance)


def test_rtheta_cross_precision_after_reduction():
    mp.dps = 40
    tau = [[mpc('0.2', '0.15'), mpc('0.12', '0.03')],
           [mpc('0.12', '0.03'), mpc('-0.1', '0.2')]]
    z = [mpc('0.13', '0.07'), mpc('-0.21', '0.04')]
    characteristic = ((mpf('0.3'), mpf('-0.2')),
                      (mpf('0.1'), mpf('0.4')))
    value = rtheta(z, tau, characteristic, derivative=(1, 1))
    with mp.workdps(75):
        reference = rtheta(z, tau, characteristic, derivative=(1, 1))
    assert mp.almosteq(value, reference, rel_eps=mpf('1e-38'))


def test_rtheta_derivatives_after_varying_arguments():
    mp.dps = 30
    tau = [[mpc('0.2', '1.2')]]
    z = [mpc('0.1', '0.01')]
    initial = rtheta(z, tau, derivative=1)
    for i in range(40):
        argument = mpc(mpf(i + 1) / 8, mpf(i + 1) / 100)
        value = rtheta([argument], tau, derivative=1)
        with mp.workdps(60):
            q = exp(pi * j * tau[0][0])
            reference = pi * jtheta(3, pi * argument, q, derivative=1)
        assert mp.almosteq(value, reference,
                           rel_eps=mpf('1e-28'), abs_eps=mpf('1e-28'))
    assert mp.almosteq(rtheta(z, tau, derivative=1), initial)


def test_rtheta_jets_preserve_inputs_and_precision():
    mp.dps = 30
    tau = mp.matrix([[mpc('0.2', '0.15'), mpc('0.12', '0.03')],
                     [mpc('0.12', '0.03'), mpc('-0.1', '0.2')]])
    z = mp.matrix([mpc('0.13', '0.07'), mpc('-0.21', '0.04')])
    characteristic = (mp.matrix([0.5, 0]), mp.matrix([0, 0.5]))
    original_tau, original_z = tau.copy(), z.copy()
    original_characteristic = tuple(v.copy() for v in characteristic)
    initial = rtheta_jet(z, tau, 2, characteristic)
    for argument, char in ((z + mp.matrix([1, -2]), characteristic),
                           (z, ([0.25, 0], [0, 0.5]))):
        jet = rtheta_jet(argument, tau, 2, char)
        with mp.workdps(65):
            reference = rtheta_jet(argument, tau, 2, char)
        for derivative, value in jet.items():
            assert mp.almosteq(value, reference[derivative],
                               rel_eps=mpf('1e-28'), abs_eps=mpf('1e-28'))
    repeated = rtheta_jet(z, tau, 2, characteristic)
    for derivative, value in repeated.items():
        assert mp.almosteq(value, initial[derivative])
    assert mp.dps == 30
    assert tau == original_tau and z == original_z
    assert characteristic == original_characteristic


def test_rtheta_wolfram_first_derivative_reference_values():
    # Wolfram Engine 14.3 NumericalCalculus`ND values, evaluated using
    # Method -> NIntegrate (Cauchy's integral formula) and 65-70 digits of
    # working precision. SiegelTheta has the same normalized z convention.
    # https://reference.wolfram.com/language/NumericalCalculus/ref/ND.html
    # https://reference.wolfram.com/language/ref/SiegelTheta.html
    mp.dps = 40
    z2 = [mpc('0.11', '0.03'), mpc('-0.07', '0.02')]
    tau2 = [
        [mpc('0.10', '0.90'), mpc('-0.20', '0.12')],
        [mpc('-0.20', '0.12'), mpc('0.05', '1.10')],
    ]
    characteristic2 = ([0.5, 0], [0, 0.5])
    references2 = (
        (None, mpc(
            '-0.421288315640741560394787218172215298435192461898651',
            '-0.297478129127372247972205522599835990950319316362461')),
        (characteristic2, mpc(
            '-0.983740135148338116134552482066670136613333239537688',
            '-0.287732601392744703482773208619225173728089061045940')),
    )
    tolerance = mpf('1e-38')
    for characteristic, reference in references2:
        assert mp.almosteq(
            rtheta(z2, tau2, characteristic, derivative=(1, 0)),
            reference, rel_eps=tolerance, abs_eps=tolerance)

    z3 = [mpc('0.11', '0.03'), mpc('-0.07', '0.02'),
          mpc('0.05', '-0.01')]
    tau3 = [
        [mpc('0.10', '1.00'), mpc('-0.12', '0.08'),
         mpc('0.05', '0.04')],
        [mpc('-0.12', '0.08'), mpc('0.20', '1.15'),
         mpc('0.08', '0.06')],
        [mpc('0.05', '0.04'), mpc('0.08', '0.06'),
         mpc('-0.15', '0.90')],
    ]
    reference3 = mpc(
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
    tolerance = mpf('1e-98')
    z2 = [mpc('0.11', '0.03'), mpc('-0.07', '0.02')]
    tau2 = [
        [mpc('0.10', '0.90'), mpc('-0.20', '0.12')],
        [mpc('-0.20', '0.12'), mpc('0.05', '1.10')],
    ]
    half2 = ([0.5, 0], [0, 0.5])
    cases2 = (
        (None, (2, 0), mpc(
            '-3.70302318847962531102746421196741770479519169317250238327399845194544812073725950034281322335010842942379904608583600925',
            '-0.676647983025765584988696701501531718662199744767527809377309914862400209905506350145714074349651391583245694053251869978')),
        (None, (1, 1), mpc(
            '-0.0492830462579949939080327254030849463068181598703018850653065257700124493596665111519633528138121999674893666903240653057',
            '0.188722819532323736052517708362276883159628967942584955480523117746601772368091761871289179463169032107921550117333179792')),
        (None, (2, 1), mpc(
            '0.510109942412506906348039968829021056110472011913214068324156555703573356310679474090641891489414592720300153907825795172',
            '-1.73784821436088196692553878274949762796692250748447504881253235541903673215233093231020937149976587228035056900176806436')),
        (half2, (0, 2), mpc(
            '1.59405366912712602932168341925929329739061553348406695029834721987354070660484011542785642549687401680165762911354469987',
            '0.678197572830669325431970259193745616335781067306349265494919453231503185934986593023788155173334269512663008718046546465')),
        (half2, (1, 1), mpc(
            '0.00506262711244847595711909837342906731732086829953077179875073877425155230227683339892785473221686111182635405587903951584',
            '-0.656735703386161600190210786006010705331432797990974033139682794285353147435261926555954518014351969442099487008921964491')),
        (half2, (1, 2), mpc(
            '-2.56012945425409603208529451407731106831476911040452819698668383790292721185892340100183555947273879000517313127472716906',
            '-3.62621222740679372224268928087397123446120699815806027482289280696049350547713834880575823432416970384201643649467346182')),
    )
    for characteristic, derivative, reference in cases2:
        assert mp.almosteq(
            rtheta(z2, tau2, characteristic, derivative), reference,
            rel_eps=tolerance, abs_eps=tolerance)

    z3 = z2 + [mpc('0.05', '-0.01')]
    tau3 = [
        [mpc('0.10', '1.00'), mpc('-0.12', '0.08'),
         mpc('0.05', '0.04')],
        [mpc('-0.12', '0.08'), mpc('0.20', '1.15'),
         mpc('0.08', '0.06')],
        [mpc('0.05', '0.04'), mpc('0.08', '0.06'),
         mpc('-0.15', '0.90')],
    ]
    half3 = ([0.5, 0, 0.5], [0, 0, 0.5])
    cases3 = (
        (None, (1, 1, 0), mpc(
            '-0.0583203441458341483463356071834041239984343441778997778365789034110166610739832774214620984517732788372057366818649126281',
            '0.0708753457491918488372972088434176429774075409403303751129258604355297050164131636926147821958543195598178704144781753017')),
        (None, (2, 0, 1), mpc(
            '0.114306309858173014349559978028599658075537833007449717601194351437916360428696872814207699015024136431925580013134749270',
            '0.149072636036949478276875893311970273856306502655947226037232741110709312077642678802950638931487633597974738441036317957')),
        (half3, (0, 1, 1), mpc(
            '-0.330573454408454085191250734483989258773110711237355724952771550857516516520533097128909554697460504958496922554130515797',
            '-0.248921244221115960904005325671928134409098763816549339891050443057701465746834442584948227393550587676101381266386567494')),
        (half3, (1, 1, 1), mpc(
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
    with pytest.raises(TypeError):
        rtheta([0], [[1j]], derivative=None)
    with pytest.raises(ValueError, match="nonnegative"):
        rtheta([0], [[1j]], derivative=(1.0,))
    with pytest.raises(ValueError, match="nonnegative"):
        rtheta([0], [[1j]], derivative=(-1,))
    with pytest.raises(ValueError, match="nonnegative integer"):
        rtheta_jet([0], [[1j]], -1)
    with pytest.raises(ValueError, match="nonnegative integer"):
        rtheta_jet([0], [[1j]], 1.0)


def test_rtheta_derivative_with_small_imaginary_periods():
    mp.dps = 35
    tau = ((mpc('0.2', '0.15'), mpc('0.12', '0.03')),
           (mpc('0.12', '0.03'), mpc('-0.1', '0.2')))
    z = (mpc('0.13', '0.07'), mpc('-0.21', '0.04'))
    characteristic = ((mpf('0.3'), mpf('-0.2')), (mpf('0.1'), mpf('0.4')))
    value = rtheta(z, tau, characteristic, derivative=(1, 0))
    with mp.workdps(60):
        reference = diff(lambda x: rtheta([x, z[1]], tau, characteristic), z[0])
    assert mp.almosteq(value, reference,
                       rel_eps=mpf('1e-33'), abs_eps=mpf('1e-33'))


def test_rtheta_genus_three_mixed_derivative_with_small_imaginary_periods():
    mp.dps = 20
    tau = (
        (mpc('0.3', '0.12'), mpc('0.17', '0.03'), mpc(0, '0.02')),
        (mpc('0.17', '0.03'), mpc('-0.2', '0.16'), mpc(0, '0.02')),
        (mpc(0, '0.02'), mpc(0, '0.02'), mpc('0.1', '0.2')),
    )
    z = (mpc('0.11', '0.03'), mpc('-0.07', '0.02'), mpc('0.05', '-0.01'))
    value = rtheta(z, tau, derivative=(1, 1, 0))
    with mp.workdps(45):
        reference = diff(lambda x: diff(
            lambda y: rtheta([x, y, z[2]], tau), z[1]), z[0])
    assert mp.almosteq(value, reference,
                       rel_eps=mpf('1e-18'), abs_eps=mpf('1e-18'))


def test_rtheta_mixed_derivative_with_small_imaginary_periods():
    mp.dps = 25
    tau = ((mpc('0.4', '0.03'), mpc('0.1', '0.008')),
           (mpc('0.1', '0.008'), mpc('-0.3', '0.05')))
    z = (mpc('0.13', '0.07'), mpc('-0.21', '0.04'))
    characteristic = ((0.5, 0), (0, 0.5))
    value = rtheta(z, tau, characteristic, derivative=(1, 1))
    with mp.workdps(50):
        reference = diff(lambda x: diff(
            lambda y: rtheta([x, y], tau, characteristic), z[1]), z[0])
    assert mp.almosteq(value, reference,
                       rel_eps=mpf('1e-23'), abs_eps=mpf('1e-23'))


def test_rtheta_jet_with_small_imaginary_periods_matches_scalar_calls():
    mp.dps = 25
    tau = ((mpc('0.4', '0.03'), mpc('0.1', '0.008')),
           (mpc('0.1', '0.008'), mpc('-0.3', '0.05')))
    z = (mpc('0.13', '0.07'), mpc('-0.21', '0.04'))
    characteristic = ((0.5, 0), (0, 0.5))
    jet = rtheta_jet(z, tau, 2, characteristic)
    for derivative, value in jet.items():
        with mp.workdps(55):
            reference = rtheta(z, tau, characteristic, derivative=derivative)
        assert mp.almosteq(value, reference,
                           rel_eps=mpf('1e-23'), abs_eps=mpf('1e-23'))

"""Tests for Riemann theta values and period geometry."""

from itertools import product

import pytest

from mpmath import exp, j, jtheta, mp, mpc, mpf, pi, rtheta


def test_rtheta_genus_one_jtheta_characteristics():
    # Classical-theta correspondence: DLMF 21.2.8-21.2.12.
    # https://dlmf.nist.gov/21.2.iii
    mp.dps = 30
    w = mpc('0.3', '0.1')
    tau = mpc('0.2', '0.9')
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
    tau = mpc('0.17', '0.003')
    characteristic = ((0.5,), (0.5,))
    arguments = ('0', '1e-5', '1e-20', '1e-40')
    values = [
        rtheta([mpf(w) / pi], [[tau]], characteristic)
        for w in arguments
    ]

    with mp.workdps(80):
        reference_tau = mpc('0.17', '0.003')
        q = exp(pi * j * reference_tau)
        references = [-jtheta(1, mpf(w), q) for w in arguments]

    for value, reference in zip(values, references):
        assert abs(value - reference) < 10 * mp.eps


def test_rtheta_diagonal_factorisation():
    # A diagonal period matrix separates the defining sum in DLMF 21.2.1.
    # https://dlmf.nist.gov/21.2.E1
    mp.dps = 30
    z = [mpc('0.1', '0.03'), mpc('-0.2', '0.04')]
    tau = [[mpc('0.1', '0.8'), 0], [0, mpc('-0.2', '1.1')]]
    value = rtheta(z, tau)
    expected = (rtheta([z[0]], [[tau[0][0]]])
                * rtheta([z[1]], [[tau[1][1]]]))
    assert mp.almosteq(value, expected)


def test_rtheta_parity_and_characteristic_zero():
    # Zero- and half-characteristic parity: DLMF 21.3.1 and 21.3.6.
    # https://dlmf.nist.gov/21.3.E1
    # https://dlmf.nist.gov/21.3.E6
    mp.dps = 30
    tau = [[1j, mpc('0.1', '0.05')],
           [mpc('0.1', '0.05'), mpc('0.2', '1.2')]]
    z = [mpc('0.13', '0.02'), mpc('-0.07', '0.01')]
    assert mp.almosteq(rtheta(z, tau), rtheta([-z[0], -z[1]], tau))
    odd = ([0.5, 0], [0.5, 0])
    assert abs(rtheta([0, 0], tau, odd)) < mp.eps * 10


def test_rtheta_all_genus_two_half_characteristic_parities():
    # Half-characteristic parity: DLMF 21.3.6.
    # https://dlmf.nist.gov/21.3.E6
    mp.dps = 30
    half = 0.5
    tau = [[mpc('0.13', '1.05'), mpc('-0.09', '0.06')],
           [mpc('-0.09', '0.06'), mpc('-0.17', '1.2')]]
    z = [mpc('0.14', '0.03'), mpc('-0.11', '0.02')]
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
    tau = [[mpc('0.13', '1.05'), mpc('-0.09', '0.06')],
           [mpc('-0.09', '0.06'), mpc('-0.17', '1.2')]]
    characteristic = ((0.5, 0), (0.5, 0))
    offsets = ('1e-5', '1e-20', '1e-40')
    values = [
        rtheta([mpf(offset), 2 * mpf(offset)], tau, characteristic)
        for offset in offsets
    ]

    with mp.workdps(80):
        reference_tau = [
            [mpc('0.13', '1.05'), mpc('-0.09', '0.06')],
            [mpc('-0.09', '0.06'), mpc('-0.17', '1.2')],
        ]
        references = [
            rtheta(
                [mpf(offset), 2 * mpf(offset)], reference_tau,
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
    tau = [[mpc('0.1', '0.9'), mpc('0.05', '0.02')],
           [mpc('0.05', '0.02'), mpc('-0.1', '1.1')]]
    z = [mpc('0.12', '0.02'), mpc('-0.08', '0.03')]
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


@pytest.mark.parametrize("genus,radius", [(1, 8), (2, 8), (3, 6)])
def test_rtheta_matches_defining_sum(genus, radius):
    # Independently sum a rectangular box in the defining series. These
    # strongly damped periods make omitted terms negligible at 40 digits.
    mp.dps = 40
    tau = mp.matrix([
        [mpc('0.11', '1.1'), mpc('-0.08', '0.05'), mpc('0.04', '0.03')],
        [mpc('-0.08', '0.05'), mpc('-0.14', '1.2'), mpc('0.02', '0.04')],
        [mpc('0.04', '0.03'), mpc('0.02', '0.04'), mpc('0.09', '1.3')],
    ])[:genus, :genus]
    z = [mpc('0.17', '0.11'), mpc('-0.12', '0.07'),
         mpc('0.05', '-0.02')][:genus]
    a = [mpf('0.3'), mpf('-0.2'), mpf('0.1')][:genus]
    b = [mpf('0.1'), mpf('0.4'), mpf('-0.3')][:genus]
    with mp.workdps(65):
        terms = []
        for point in product(range(-radius, radius + 1), repeat=genus):
            n = mp.matrix([point[i] + a[i] for i in range(genus)])
            terms.append(exp(pi * j * (n.T * tau * n)[0]
                             + 2 * pi * j * (n.T * (mp.matrix(z) + mp.matrix(b)))[0]))
        reference = mp.fsum(terms)
    assert mp.almosteq(rtheta(z, tau, (a, b)), reference,
                       rel_eps=mpf('1e-38'), abs_eps=mpf('1e-38'))


def test_rtheta_wolfram_reference_values():
    # Wolfram Engine 14.3 SiegelTheta values generated at 120 digits from
    # exact rational inputs. The tests use 100 digits and retain guard digits
    # in each reference value. The calls were N[SiegelTheta[tau, z], 120]
    # and N[SiegelTheta[{a, b}, tau, z], 120].
    # https://reference.wolfram.com/language/ref/SiegelTheta.html
    mp.dps = 100
    z2 = [mpc('0.11', '0.03'), mpc('-0.07', '0.02')]
    tau2 = [
        [mpc('0.10', '0.90'), mpc('-0.20', '0.12')],
        [mpc('-0.20', '0.12'), mpc('0.05', '1.10')],
    ]
    char2 = ([0.5, 0], [0, 0.5])
    refs2 = [
        mpc(
            '1.15010054896589710068630628822243271125317054007702426004984068418064969671373665362574581188457525423704174753913082529848244697358328',
            '0.029482887366211288410342712159136306372758610990096207370983577503724607487203127164966369240820605787739356900521381741214906114952103'),
        mpc(
            '0.893359834184414572430609669893702469378276216582063538344544618948879436225815023389107207488871647434926394149016401840086595696537675',
            '0.02505879729520504089469665867573950858628706160028878239473482035693999039117935887428157275486164770342340639145800073069780958643882'),
    ]

    z3 = z2 + [mpc('0.05', '-0.01')]
    tau3 = [
        [mpc('0.10', '1.00'), mpc('-0.12', '0.08'),
         mpc('0.05', '0.04')],
        [mpc('-0.12', '0.08'), mpc('0.20', '1.15'),
         mpc('0.08', '0.06')],
        [mpc('0.05', '0.04'), mpc('0.08', '0.06'),
         mpc('-0.15', '0.90')],
    ]
    char3 = ([0.5, 0, 0.5], [0, 0, 0.5])
    refs3 = [
        mpc(
            '1.221609642406339797657102259903691954365704088043749183987260288160622810706259133711255570012034796050287767505312762580074107266344458',
            '-0.007792386568554659454892241427121689475136799802090774325255515968433879603115006454848761576475120048373924432677064251983242053776354'),
        mpc(
            '-0.114796433952467492469811475212021962597569270860961487827600713071632671744755401699850275301415484731581776951960982866615787356410398',
            '0.015318780171864079253785680566109776307328372396748562810674785685663933968972735107125484955849153511806284150913864802098269822109827'),
    ]

    tolerance = mpf('1e-98')
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
    characteristic2 = ((mpf('0.3'), mpf('-0.2')),
                       (mpf('0.1'), mpf('0.4')))
    tau2 = [[mpc('0.2', '0.15'), mpc('0.12', '0.03')],
            [mpc('0.12', '0.03'), mpc('-0.1', '0.2')]]
    z2 = [mpc('2.31', '-1.14'), mpc('-1.72', '0.83')]
    reference2 = mpc(
        '37058424061976446096.7871409113145979695568610946774019090331353765326092054797218949660046676284791419537804054179813145',
        '-8184623208069837737.99363797667810872646677513420698882738727935055524579125490740124342217355390271095127290985384489359')

    characteristic3 = ((mpf('0.3'), mpf('-0.2'), mpf('0.1')),
                       (mpf('0.1'), mpf('0.4'), mpf('-0.15')))
    tau3 = [
        [mpc('0.3', '0.12'), mpc('0.17', '0.03'),
         mpc(0, '0.02')],
        [mpc('0.17', '0.03'), mpc('-0.2', '0.16'),
         mpc(0, '0.02')],
        [mpc(0, '0.02'), mpc(0, '0.02'),
         mpc('0.1', '0.2')],
    ]
    z3 = [mpc('0.11', '0.03'), mpc('-0.07', '0.02'),
          mpc('0.05', '-0.01')]
    reference3 = mpc(
        '3.44042050000036293079504325661815829014001245512124516294720089024113963744406088834494220278232680204967316690887511458',
        '0.182598566348364806829149624629842477324789098204541555796804418499811344983007460195610367554313435739194325595031041258')

    tolerance = mpf('1e-98')
    assert mp.almosteq(rtheta(z2, tau2, characteristic2), reference2,
                       rel_eps=tolerance, abs_eps=tolerance)
    assert mp.almosteq(rtheta(z3, tau3, characteristic3), reference3,
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
            tau[i][i] = mpc(
                mpf(i % 3) / 20,
                scale + mpf(i % 4) / 10,
            )
            if i + 1 < genus:
                tau[i][i + 1] = tau[i + 1][i] = mpc('0.02', '0.25')
        z = [
            mpc(mpf((i % 5) - 2) / 100,
                   mpf((i % 3) - 1) / 200)
            for i in range(genus)
        ]
        return z, tau

    cases = (
        (10, 8, mpc(
            '1.0000000001671156401975485245765975730135867034602892301525',
            '2.49857014987438782428182159399172658122217981149e-11')),
        (20, 20, mpc(
            '1.0000000000000000000000000133796107582295982965475652944874',
            '2.0101530593695935374230689609606e-27')),
    )
    tolerance = mpf('1e-33')
    for genus, scale, reference in cases:
        z, tau = make_case(genus, scale)
        assert mp.almosteq(rtheta(z, tau), reference,
                           rel_eps=tolerance, abs_eps=tolerance)


def test_rtheta_validation():
    mp.dps = 20
    with pytest.raises(TypeError):
        rtheta(None, [[1j]])
    with pytest.raises(ValueError, match="vector"):
        rtheta(mp.matrix([[0, 0]]), [[1j]])
    with pytest.raises(ValueError, match="finite"):
        rtheta([mp.inf], [[1j]])
    with pytest.raises(TypeError):
        rtheta([0], None)
    with pytest.raises(ValueError, match="square"):
        rtheta([0], [[1j, 0]])
    with pytest.raises(ValueError, match="finite"):
        rtheta([0], [[mp.inf + 1j]])
    with pytest.raises(ValueError, match="length"):
        rtheta([0], [[1j, 0], [0, 1j]])
    with pytest.raises(ValueError, match="symmetric"):
        rtheta([0, 0], [[1j, 1], [0, 1j]])
    for tau in ([[0]], [[1j, 0], [0, -1j]], [[1j, 1j], [1j, 1j]]):
        with pytest.raises(ValueError, match="positive-definite"):
            rtheta([0] * len(tau), tau)
    with pytest.raises(ValueError, match="real"):
        rtheta([0], [[1j]], characteristic=([1j], [0]))
    with pytest.raises(ValueError, match="pair"):
        rtheta([0], [[1j]], characteristic=([0],))


def test_rtheta_nearly_symmetric_input():
    mp.dps = 20
    delta = mp.eps / 4
    tau = [[1j, mpc('0.1', '0.02')],
           [mpc('0.1', '0.02') + delta, 1.2j]]
    symmetric = [[1j, mpc('0.1', '0.02') + delta / 2],
                 [mpc('0.1', '0.02') + delta / 2, 1.2j]]
    assert mp.almosteq(rtheta([0.1, 0.2], tau),
                       rtheta([0.1, 0.2], symmetric))


def test_rtheta_symplectic_generator_identities():
    # DLMF 21.5.5, 21.5.7 and 21.5.4, checked through public values.
    # https://dlmf.nist.gov/21.5.E5
    # https://dlmf.nist.gov/21.5.E7
    # https://dlmf.nist.gov/21.5.E4
    mp.dps = 35
    tau = mp.matrix([
        [mpc('0.13', '1.1'), mpc('-0.17', '0.08')],
        [mpc('-0.17', '0.08'), mpc('0.21', '1.3')],
    ])
    z = mp.matrix([mpc('0.12', '0.04'), mpc('-0.09', '0.03')])
    value = rtheta(z, tau)
    translation = mp.matrix([[1, 2], [2, -1]])
    translated_tau = tau + translation
    translated_z = z - mp.matrix([0.5, -0.5])
    assert mp.almosteq(rtheta(translated_z, translated_tau), value)

    U = mp.matrix([[1, 1], [0, 1]])
    assert mp.almosteq(rtheta(U.T * z, U.T * tau * U), value)
    assert mp.almosteq(
        rtheta(U.T * translated_z, U.T * translated_tau * U), value)

    t, c = tau[0, 0], tau[0, 1]
    inverted_tau = [[-1 / t, c / t], [c / t, tau[1, 1] - c**2 / t]]
    inverted_z = [z[0] / t, z[1] - c * z[0] / t]
    factor = exp(-pi * j * z[0]**2 / t) / mp.sqrt(-j * t)
    assert mp.almosteq(factor * rtheta(inverted_z, inverted_tau), value)


def test_rtheta_skewed_basis_invariance():
    # A unimodular change of basis reindexes the defining lattice sum.
    # Dyadic entries avoid rounding the input while constructing skewed
    # periods; the function must still agree at ordinary working precision.
    mp.dps = 15
    tau = mp.diag([1j, 1.5j, 2j])
    z = mp.matrix([mpc(0.125, 0.0625), mpc(-0.0625, 0.03125),
                   mpc(0.03125, -0.015625)])
    U = mp.matrix([[1, 20, 100], [0, 1, 20], [0, 0, 1]])
    value = rtheta(U.T * z, U.T * tau * U)
    with mp.workdps(50):
        reference = mp.fprod(jtheta(3, pi * z[i], exp(pi * j * tau[i, i]))
                             for i in range(3))
    assert mp.almosteq(value, reference,
                       rel_eps=mpf('1e-12'), abs_eps=mpf('1e-12'))


def test_rtheta_quasiperiodicity_with_small_imaginary_periods():
    mp.dps = 35
    tau = mp.matrix([[mpc('0.2', '0.15'), mpc('0.12', '0.03')],
                     [mpc('0.12', '0.03'), mpc('-0.1', '0.2')]])
    z = mp.matrix([mpc('0.13', '0.07'), mpc('-0.21', '0.04')])
    characteristic = ([mpf('0.3'), mpf('-0.2')], [mpf('0.1'), mpf('0.4')])
    value = rtheta(z, tau, characteristic)
    with mp.workdps(70):
        reference = rtheta(z, tau, characteristic)
    assert mp.almosteq(value, reference, rel_eps=mpf('1e-33'))

    shift = mp.matrix([3, -2])
    quadratic = (shift.T * tau * shift)[0]
    linear = (shift.T * z)[0]
    multiplier = exp(-pi * j * quadratic - 2 * pi * j * linear)
    assert mp.almosteq(rtheta(z + tau * shift, tau),
                       multiplier * rtheta(z, tau))


def test_rtheta_small_imaginary_periods_genus_three():
    mp.dps = 25
    tau = (
        (mpc('0.3', '0.12'), mpc('0.17', '0.03'), mpc(0, '0.02')),
        (mpc('0.17', '0.03'), mpc('-0.2', '0.16'), mpc(0, '0.02')),
        (mpc(0, '0.02'), mpc(0, '0.02'), mpc('0.1', '0.2')),
    )
    z = (mpc('0.13', '0.07'), mpc('-0.21', '0.04'), mpc('0.06', '-0.03'))
    value = rtheta(z, tau)
    with mp.workdps(60):
        reference = rtheta(z, tau)
    assert mp.almosteq(value, reference,
                       rel_eps=mpf('1e-23'), abs_eps=mpf('1e-23'))


def test_rtheta_unreduced_genus_one_against_jtheta():
    # Classical-theta correspondence: DLMF 21.2.8-21.2.12.
    # https://dlmf.nist.gov/21.2.iii
    mp.dps = 35
    tau = mpc('0.1', '0.1')
    w = mpc('0.3', '0.17')
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

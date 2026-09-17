"""Public LLL interface checks using exact rational arithmetic."""

from fractions import Fraction
from itertools import permutations
import random
import math

import pytest

from mpmath import fp, lll_gram, mp
from mpmath.matrices import lattice


@pytest.fixture(params=['automatic', 'numerical'], autouse=True)
def reduction_path(request, monkeypatch):
    """Keep the complete public-contract corpus exercising both implementations."""
    if request.param == 'numerical':
        monkeypatch.setattr(lattice, '_use_exact', lambda *args, **kwargs: False)


def determinant(A):
    n = len(A)
    total = 0
    for p in permutations(range(n)):
        sign = (-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = sign
        for i in range(n):
            term *= A[i][p[i]]
        total += term
    return total


def assert_reduced(Y, U, delta=Fraction(3, 4)):
    n = len(Y)
    assert isinstance(U, tuple) and all(isinstance(row, tuple) for row in U)
    assert all(len(row) == n for row in U)
    assert all(type(x) is int for row in U for x in row)
    assert determinant(U) == 1
    G = [[sum(U[k][i] * Fraction(Y[k][l]) * U[l][j]
              for k in range(n) for l in range(n))
          for j in range(n)] for i in range(n)]
    # Orthogonalize coordinate vectors in the inner product defined by G.
    def inner(v, w):
        return sum(v[i] * G[i][j] * w[j] for i in range(n) for j in range(n))
    orthogonal = []
    lengths = []
    for i in range(n):
        e = [Fraction(int(i == j)) for j in range(n)]
        v = e[:]
        coefficients = []
        for j, b in enumerate(orthogonal):
            mu = inner(e, b) / lengths[j]
            assert abs(mu) <= Fraction(1, 2)
            coefficients.append(mu)
            v = [x - mu * y for x, y in zip(v, b)]
        length = inner(v, v)
        assert length > 0
        if i:
            assert length >= (delta - coefficients[-1] ** 2) * lengths[-1]
        orthogonal.append(v)
        lengths.append(length)


def gram(B):
    return [[sum(row[i] * row[j] for row in B)
             for j in range(len(B[0]))] for i in range(len(B[0]))]


@pytest.mark.parametrize('ctx', [mp, fp])
def test_lll_gram_examples_and_orientation(ctx):
    lll_gram = ctx.lll_gram
    assert lll_gram([[7]]) == ((1,),)
    assert lll_gram([[1, 10], [10, 101]]) == ((1, -10), (0, 1))
    assert lll_gram([[4, 0], [0, 1]]) == ((0, -1), (1, 0))
    assert_reduced([[5, 3], [3, 2]], lll_gram([[5, 3], [3, 2]]))


@pytest.mark.parametrize('ctx', [mp, fp])
@pytest.mark.parametrize('delta', ['0.26', '0.5', '0.75', '0.99'])
def test_lll_gram_random_bases(ctx, delta):
    rng = random.Random(1182)
    d = ctx.mpf(delta)
    exact_delta = Fraction(*d.as_integer_ratio())
    for n in range(2, 6):
        for unused in range(3):
            B = [[rng.randint(-8, 8) for j in range(n)] for i in range(n)]
            # Appending I ensures independent columns without a rank filter.
            B += [[int(i == j) for j in range(n)] for i in range(n)]
            Y = gram(B)
            assert_reduced(Y, ctx.lll_gram(Y, delta=d), exact_delta)


@pytest.mark.parametrize('dps', [15, 30, 60])
def test_lll_gram_noninteger_inputs(dps):
    mp.dps = dps
    B = [[Fraction(1, 8), Fraction(3, 2), Fraction(5, 4)],
         [0, Fraction(1, 4), Fraction(7, 8)], [0, 0, Fraction(1, 2)]]
    Y = gram(B)
    assert_reduced(Y, lll_gram([[mp.mpf(x.numerator) / x.denominator for x in row] for row in Y]))


@pytest.mark.parametrize('power', [60, 200])
def test_lll_gram_inequality_boundaries(power):
    a = 2 ** power
    for offset in (-1, 0, 1):
        for Y in ([[a, a // 2 + offset], [a // 2 + offset, a]],
                  [[a, 0], [0, 3 * a // 4 + offset]]):
            assert_reduced(Y, lll_gram(Y))


@pytest.mark.parametrize('rounding', ['n', 'f', 'c', 'd', 'u'])
def test_lll_gram_signed_half_integer_coefficients(rounding):
    ctx = mp.clone()
    for numerator in (-5, -3, -1, 1, 3, 5):
        for offset in (-1, 0, 1):
            with ctx.workprec(180):
                q = ctx.mpf(numerator) / 2 + ctx.mpf(offset) * ctx.mpf(2) ** -70
                Y = ctx.matrix([[1, q], [q, q*q + 1]])
            exact = [[Fraction(*x.as_integer_ratio()) for x in row] for row in Y.tolist()]
            ctx.rounding = rounding
            assert_reduced(exact, ctx.lll_gram(Y))


@pytest.mark.parametrize('ctx', [mp, fp])
@pytest.mark.parametrize('power', [-100, 0, 100])
def test_lll_gram_mixed_binary_denominators(ctx, power):
    B = [[Fraction(1, 32), Fraction(-3, 16), Fraction(5, 8)],
         [0, Fraction(1, 4), Fraction(-7, 2)], [0, 0, 16]]
    scale = Fraction(2) ** power
    Y = [[x * scale for x in row] for row in gram(B)]
    converted = [[ctx.mpf(x.numerator) / x.denominator for x in row] for row in Y]
    assert_reduced(Y, ctx.lll_gram(converted))


@pytest.mark.parametrize('n', [3, 4, 5])
@pytest.mark.parametrize('s', [10 ** 8, 10 ** 12])
def test_lll_gram_skew_and_large_integer_transform(n, s):
    B = [[int(i == j) + (s if j == i + 1 else 0) for j in range(n)] for i in range(n)]
    Y = gram(B)
    U = lll_gram(Y)
    assert_reduced(Y, U)
    if n >= 4:
        assert max(abs(x) for row in U for x in row) > 2 ** mp.prec


def test_lll_gram_preserves_context_and_input():
    ctx = mp.clone()
    ctx.dps = 80
    Y = ctx.matrix([[1, ctx.mpf('0.1')], [ctx.mpf('0.1'), 2]])
    before = Y.tolist()
    ctx.dps = 15
    for rounding in ('n', 'f', 'c', 'd', 'u'):
        ctx.rounding = rounding
        U = ctx.lll_gram(Y)
        assert_reduced([[Fraction(*x.as_integer_ratio()) for x in row] for row in before], U)
        assert Y.tolist() == before
        assert ctx.dps == 15 and ctx.rounding == rounding
    assert mp.dps == 15 and mp.rounding == 'n'


@pytest.mark.parametrize('power', [-2000, 2000])
def test_lll_gram_outside_float_range(power):
    ctx = mp.clone()
    scale = ctx.mpf(2) ** power
    Y = ctx.matrix([[2*scale, scale], [scale, 2*scale]])
    exact = [[Fraction(*x.as_integer_ratio()) for x in row] for row in Y.tolist()]
    before = ctx.prec
    assert_reduced(exact, ctx.lll_gram(Y))
    assert ctx.prec == before


def test_lll_gram_high_precision_input_at_low_working_precision():
    ctx = mp.clone()
    s = 2 ** 2000
    with ctx.workprec(4030):
        Y = ctx.matrix([[1, s], [s, s*s + 1]])
    before = ctx.prec
    assert_reduced([[1, s], [s, s*s + 1]], ctx.lll_gram(Y))
    assert ctx.prec == before


@pytest.mark.parametrize('Y', [[], [1, 2], [[1, 2, 3], [2, 4, 5]],
                               [[1, 0], [0]], None])
@pytest.mark.parametrize('ctx', [mp, fp])
def test_lll_gram_invalid_shape(ctx, Y):
    with pytest.raises(ValueError, match='square'):
        ctx.lll_gram(Y)


@pytest.mark.parametrize('Y, message', [
    ([[1, 2], [0, 1]], 'symmetric'),
    ([[1, 0], [0, 0]], 'positive definite'),
    ([[1, 1], [1, 1]], 'positive definite'),
    ([[1, 2], [2, 1]], 'positive definite'),
    ([[-1]], 'positive definite'),
    ([[mp.inf]], 'finite and real'),
    ([[mp.nan]], 'finite and real'),
    ([[1j]], 'finite and real'),
])
@pytest.mark.parametrize('ctx', [mp, fp])
def test_lll_gram_invalid_entries(ctx, Y, message):
    with pytest.raises(ValueError, match=message):
        ctx.lll_gram(Y)


@pytest.mark.parametrize('delta', [0.25, 1, 2, -1, mp.nan, mp.inf, 1j, 'bad', None])
@pytest.mark.parametrize('ctx', [mp, fp])
def test_lll_gram_invalid_delta(ctx, delta):
    with pytest.raises(ValueError, match='delta'):
        ctx.lll_gram([[1]], delta=delta)


@pytest.mark.parametrize('limit', [0, -1, 1.5, None])
@pytest.mark.parametrize('ctx', [mp, fp])
def test_lll_gram_invalid_step_limit(ctx, limit):
    with pytest.raises(ValueError, match='maxsteps'):
        ctx.lll_gram([[1]], maxsteps=limit)


def test_lll_gram_step_limit_restores_precision():
    before = mp.prec
    with pytest.raises(mp.NoConvergence, match='maxsteps'):
        lll_gram([[4, 0], [0, 1]], maxsteps=1)
    assert mp.prec == before


def test_lll_gram_persistent_numerical_failure(monkeypatch):
    monkeypatch.setattr(lattice, '_use_exact', lambda *args, **kwargs: False)
    ctx = mp.clone()
    before = ctx.prec
    precisions = []
    # Float conversion loses positive definiteness, forcing the mp fallback.
    s = 2 ** 100
    with ctx.workprec(240):
        Y = ctx.matrix([[1, s], [s, s*s + 1]])

    def failing_sum(*args, **kwargs):
        precisions.append(ctx.prec)
        raise ZeroDivisionError("simulated numerical failure")

    # Exercise the public failure contract without a pathological large matrix.
    monkeypatch.setattr(ctx, 'fsum', failing_sum)
    with pytest.raises(ctx.NoConvergence, match='five precision attempts'):
        ctx.lll_gram(Y)
    assert len(precisions) == 5
    assert all(b > a for a, b in zip(precisions, precisions[1:]))
    assert ctx.prec == before


@pytest.mark.parametrize('scale', [1e-300, 1e-200, 1.0, 1e200, 1e300])
def test_lll_gram_fp_scaled_inputs(scale):
    Y = [[2 * scale, scale], [scale, 2 * scale]]
    assert_reduced(Y, fp.lll_gram(Y))


def test_lll_gram_fp_boundaries():
    for value in (math.nextafter(0.5, 0), 0.5, math.nextafter(0.5, 1)):
        Y = [[1.0, value], [value, 1.0]]
        assert_reduced(Y, fp.lll_gram(Y))
    for value in (math.nextafter(0.75, 0), 0.75, math.nextafter(0.75, 1)):
        Y = [[1.0, 0.0], [0.0, value]]
        assert_reduced(Y, fp.lll_gram(Y))


@pytest.mark.parametrize('n, s', [(3, 10 ** 6), (5, 10 ** 5)])
def test_lll_gram_fp_skewed_inputs(n, s):
    B = [[int(i == j) + (s if j == i + 1 else 0) for j in range(n)] for i in range(n)]
    Y = gram(B)
    before = fp.prec
    assert_reduced(Y, fp.lll_gram(Y))
    assert fp.prec == before
    # The input survived conversion; mp can resolve the reduction itself.
    assert_reduced(Y, mp.lll_gram(fp.matrix(Y)))


def test_lll_gram_fp_input_rounding():
    s = 10 ** 8
    Y = [[1, s], [s, s * s + 1]]
    with pytest.raises(ValueError, match='positive definite'):
        fp.lll_gram(Y)
    assert_reduced(Y, mp.lll_gram(Y))


def test_lll_gram_fp_rejects_unverified_reduction(monkeypatch):
    s = 2 ** 27
    Y = [[s, 1, 0], [1, s, s // 2], [0, s // 2, s]]
    # mu[2,1] = s**2 / (2*(s**2-1)) is just above 1/2, but fp
    # rounds it to 1/2 and leaves the basis unchanged. Verification must
    # reject it even though no arithmetic operation raised an exception.
    monkeypatch.setattr(lattice, '_use_exact', lambda *args, **kwargs: False)
    with pytest.raises(fp.NoConvergence, match='could not be verified'):
        fp.lll_gram(Y)
    assert_reduced(Y, mp.lll_gram(Y))


@pytest.mark.parametrize('ctx', [mp, fp])
def test_lll_gram_large_dimension(ctx):
    assert ctx.lll_gram(ctx.eye(17)) == tuple(tuple(int(i == j) for j in range(17))
                                           for i in range(17))


def test_lll_gram_exact_normalization_sign():
    for x in (0, -1, -2, -mp.mpf(2)**-2000):
        with pytest.raises(ValueError, match='positive definite'):
            mp.lll_gram([[x]])


def test_lll_gram_exact_fallback():
    ctx = mp.clone()
    s = 2 ** 1000
    with ctx.workprec(2030):
        Y = ctx.matrix([[1, s], [s, s*s + 1]])
    assert_reduced([[1, s], [s, s*s + 1]], ctx.lll_gram(Y))


def test_lll_gram_fp_overflow_and_step_limit():
    Y = [[1e-300, 1.0], [1.0, 2e300]]
    with pytest.raises(fp.NoConvergence, match='fixed precision'):
        fp.lll_gram(Y)
    assert_reduced(Y, mp.lll_gram(Y))
    with pytest.raises(fp.NoConvergence, match='maxsteps'):
        fp.lll_gram([[4, 0], [0, 1]], maxsteps=1)

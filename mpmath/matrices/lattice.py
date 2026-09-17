"""Lattice reduction from positive-definite Gram matrices."""

from functools import reduce
from math import gcd
from operator import index

from ..libmp.backend import MPQ, MPZ


# Implementation overview
# -----------------------
# lll_gram validates and converts the input once, then _reduce_gram owns method
# selection. Both algorithms return a unimodular integer change of basis U for
# the column convention U.T * Y * U, with orientation adjusted to det(U) = +1.
#
# Small inputs use fraction-free exact reduction. Clearing binary denominators
# and removing a positive common factor preserves the LLL conditions. The
# algorithm maintains determinant-scaled Gram-Schmidt coefficients, so its
# positivity and reduction tests use integer arithmetic throughout.
#
# Other inputs first get a numerical candidate. An mp call tries machine
# precision before exact fallback or increasing-precision retries; fp has one
# numerical attempt. Every numerical candidate is checked against the original
# converted input using exact rational arithmetic. Dispatch thresholds affect
# speed only; they do not relax the reduction conditions.
#
# Lists retain Python/backend integers and rationals without context rounding.
# ctx.matrix, Cholesky and QR are numerical operations, so they cannot replace
# the exact transforms and Gram-Schmidt checks below. We reuse ctx.fsum and
# ctx.nint for numerical arithmetic and MPZ/MPQ for backend-aware exact arithmetic.
# Algorithm-specific helpers stay with their owning implementation below.


class LatticeMethods:
    def lll_gram(ctx, Y, delta=0.75, *, maxsteps=10000):
        r"""
        Return an LLL-reducing change of basis for the Gram matrix *Y*.

        *Y* must be a nonempty, real, symmetric positive-definite matrix or
        nested sequence. Symmetry is required exactly, rather than to a
        tolerance. Invalid inputs raise ``ValueError``.

        The result *U* is a tuple of row tuples of Python integers with
        determinant +1, such that :math:`U^T Y U` is reduced. For a basis *B*
        stored in columns, the reduced basis is :math:`B U`. Reduced bases
        are not unique. Converting *U* to an mpmath matrix makes subsequent
        arithmetic subject to the working precision.

        The parameter ``delta`` specifies :math:`\delta`, with
        :math:`1/4 < \delta < 1` and default :math:`\delta = 3/4`.
        Larger values demand stronger reduction and may take longer.
        The Gram-Schmidt coefficients and vectors satisfy

        .. math ::

            |\mu_{i,j}| \leq \tfrac12 \quad (j < i),
            \qquad
            \|b_i^*\|^2 \geq
            (\delta - \mu_{i,i-1}^2)\|b_{i-1}^*\|^2 \quad (i > 0).

        These conditions hold for the numerical inputs after conversion to
        the context. ``maxsteps`` must be a positive integer and limits
        iterations per reduction attempt. Failure to obtain a verified
        reduction within the work limits raises ``mp.NoConvergence``.
        The ``mp`` and ``fp`` contexts are supported. If fixed precision is
        insufficient, ``fp.lll_gram`` raises ``fp.NoConvergence``; use
        ``mp.lll_gram`` with the original inputs to compute at higher precision.

        **Examples**

        Reduce a skewed basis::

            >>> from mpmath import mp, lll_gram
            >>> Y = mp.matrix([[1, 10], [10, 101]])
            >>> U = lll_gram(Y)
            >>> U
            ((1, -10), (0, 1))
            >>> V = mp.matrix(U)
            >>> V.T * Y * V == mp.eye(2)
            True

        **References**

        A. K. Lenstra, H. W. Lenstra Jr. and L. Lovasz,
        "Factoring polynomials with rational coefficients", Mathematische
        Annalen 261 (1982), 515-534.
        """
        try:
            maxsteps = index(maxsteps)
        except TypeError:
            raise ValueError("maxsteps must be a positive integer") from None
        if maxsteps <= 0:
            raise ValueError("maxsteps must be a positive integer")
        try:
            delta = ctx.convert(delta)
            if not ctx.isfinite(delta) or ctx.im(delta):
                raise ValueError
            delta = delta.real
            if not ctx.mpf('0.25') < delta < 1:
                raise ValueError
        except (TypeError, ValueError):
            raise ValueError("delta must be real and satisfy 0.25 < delta < 1") from None
        if isinstance(Y, (list, tuple)) and any(
                not isinstance(row, (list, tuple)) or len(row) != len(Y) for row in Y):
            raise ValueError("Y must be a nonempty square matrix")
        try:
            Y = ctx.matrix(Y)
        except (TypeError, ValueError, IndexError):
            raise ValueError("Y must be a nonempty square matrix") from None
        n = Y.rows
        if not n or n != Y.cols:
            raise ValueError("Y must be a nonempty square matrix")
        if any(not ctx.isfinite(x) or ctx.im(x) for x in Y):
            raise ValueError("Y entries must be finite and real")
        # Avoid ctx.re here: unary operations can round an existing mpf.
        gram = [[Y[i, j].real for j in range(n)] for i in range(n)]
        if any(gram[i][j] != gram[j][i] for i in range(n) for j in range(i)):
            raise ValueError("Y must be symmetric")
        return _reduce_gram(ctx, gram, delta, maxsteps)


# Method selection and numerical recovery.


def _reduce_gram(ctx, gram, delta, maxsteps):
    """Choose a reduction method; numerical candidates need exact verification."""
    n = len(gram)
    # Only small matrices are eligible for exact reduction. Normalize once,
    # so the initial choice and the fallback use the same size estimate.
    integral = _integer_gram(gram) if n <= 16 else None
    bits = (max(abs(x).bit_length() for row in integral for x in row)
            if integral is not None else 0)
    if _use_exact(n, bits):
        return _lll_exact(ctx, integral, delta, maxsteps)

    exact = [[MPQ(*x.as_integer_ratio()) for x in row] for row in gram]
    exact_delta = MPQ(*delta.as_integer_ratio())
    if ctx._fixed_precision:
        return _reduce_with_precision(ctx, gram, delta, maxsteps, exact, exact_delta)

    # A cheap candidate may suffice even when the input carries many digits.
    transform = _try_float_reduction(gram, delta, maxsteps, exact, exact_delta)
    if transform is not None:
        return transform
    if _use_exact(n, bits, fallback=True):
        return _lll_exact(ctx, integral, delta, maxsteps)
    return _reduce_with_precision(ctx, gram, delta, maxsteps, exact, exact_delta)



def _use_exact(n, bits, fallback=False):
    # Conservative timing crossover for the determinant-sized intermediates.
    # These bounds affect performance only, not the reduction conditions.
    if fallback:
        return n <= 8 and n * bits <= 32768
    return n <= 16 and n * bits <= 1536


def _try_float_reduction(gram, delta, maxsteps, exact, exact_delta):
    """Return a verified machine-precision candidate, or None to try another method."""
    from .. import fp
    try:
        transform = _lll_numerical(fp, [[float(x) for x in row] for row in gram],
                                   float(delta), maxsteps)
        # Check the original input, not its float approximation. Since U is
        # unimodular, successful verification also proves positive definiteness.
        if _is_reduced(exact, transform, exact_delta):
            return transform
    except (ValueError, ZeroDivisionError, OverflowError, fp.NoConvergence):
        pass
    return None


def _reduce_with_precision(ctx, gram, delta, maxsteps, exact, exact_delta):
    """Use the requested context, retrying with more precision when available."""
    # Distinguish invalid input from a valid matrix that needs more precision.
    _gram_schmidt(exact)
    if ctx._fixed_precision:
        try:
            transform = _lll_numerical(ctx, gram, delta, maxsteps)
        except (ValueError, ZeroDivisionError, OverflowError):
            raise ctx.NoConvergence(
                "LLL reduction failed at fixed precision; use mp.lll_gram") from None
        if _is_reduced(exact, transform, exact_delta):
            return transform
        raise ctx.NoConvergence(
            "LLL reduction could not be verified at fixed precision; use mp.lll_gram")

    # Retain precision carried by the input, even if ctx.prec was lowered.
    # Guard bits reduce retries; exact verification decides acceptance.
    precision = max([ctx.prec, delta._mpf_[3]]
                    + [x._mpf_[3] for row in gram for x in row]) + 20
    for attempt in range(5):
        try:
            with ctx.workprec(precision << attempt):
                transform = _lll_numerical(ctx, gram, delta, maxsteps)
        except (ValueError, ZeroDivisionError):
            continue
        if _is_reduced(exact, transform, exact_delta):
            return transform
    raise ctx.NoConvergence("LLL reduction could not be verified after five precision attempts")


# Fraction-free exact reduction.


def _integer_gram(gram):
    """Clear binary denominators and remove a positive common integer factor."""
    ratios = [[x.as_integer_ratio() for x in row] for row in gram]
    scale = max(d for row in ratios for _, d in row)
    A = [[MPZ(p) * (scale // d) for p, d in row] for row in ratios]
    common = reduce(gcd, (x for row in A for x in row), 0)
    return [[x // common for x in row] for row in A] if common else A


def _lll_exact(ctx, A, delta, maxsteps):
    """Reduce an integer Gram matrix using determinant-scaled coefficients."""
    def _exact_div(a, b):
        q, r = divmod(a, b)
        assert not r, 'fraction-free invariant failed'
        return q

    n = len(A)
    # D[j] is the determinant of the leading j by j Gram submatrix;
    # L[i][j] = mu[i][j] * D[j+1]. All these quantities are integers.
    D = [MPZ(1)]
    L = [[MPZ(0)] * n for _ in range(n)]
    for j in range(n):
        for i in range(j, n):
            value = A[i][j]
            for h in range(j):
                value = _exact_div(D[h+1]*value - L[i][h]*L[j][h], D[h])
            L[i][j] = value
        if L[j][j] <= 0:
            raise ValueError('Y must be positive definite')
        D.append(L[j][j])
    U = [[int(i == j) for j in range(n)] for i in range(n)]
    numerator, denominator = delta.as_integer_ratio()
    k = 1
    steps = 0
    orientation = 1
    while k < n:
        if steps == maxsteps:
            raise ctx.NoConvergence('LLL reduction exceeded maxsteps')
        steps += 1
        # Size-reduce column k, rounding half-integers to the nearest even integer.
        for j in range(k-1, -1, -1):
            q, remainder = divmod(L[k][j], D[j+1])
            if 2*remainder > D[j+1] or (2*remainder == D[j+1] and q % 2):
                q += 1
            if q:
                q = int(q)
                for i in range(n):
                    U[i][k] -= q*U[i][j]
                for i in range(j):
                    L[k][i] -= q*L[j][i]
                L[k][j] -= q*D[j+1]
        # Test the Lovasz condition after clearing its positive denominators.
        a, b, c, ell = D[k-1], D[k], D[k+1], L[k][k-1]
        combined = a*c + ell*ell
        if denominator*combined >= numerator*b*b:
            k += 1
        else:
            # Swap adjacent columns and update the determinant-scaled projections.
            D[k] = _exact_div(combined, b)
            for i in range(k-1):
                L[k][i], L[k-1][i] = L[k-1][i], L[k][i]
            for i in range(k+1, n):
                x, y = L[i][k-1], L[i][k]
                L[i][k-1] = _exact_div(a*y + ell*x, b)
                L[i][k] = _exact_div(c*x - ell*y, b)
            for row in U:
                row[k], row[k-1] = row[k-1], row[k]
            orientation = -orientation
            k = max(1, k-1)
    if orientation < 0:
        for row in U:
            row[-1] = -row[-1]
    return tuple(tuple(row) for row in U)


# Numerical reduction and exact verification.


def _lll_numerical(ctx, gram, delta, maxsteps):
    """Generate a candidate using incremental numerical Gram-Schmidt updates."""
    n = len(gram)
    transform = [[int(i == j) for j in range(n)] for i in range(n)]
    orientation = 1
    mu, lengths = _gram_schmidt(gram, ctx.fsum)
    k = 1
    steps = 0
    while k < n:
        if steps == maxsteps:
            raise ctx.NoConvergence("LLL reduction exceeded maxsteps")
        steps += 1
        for j in range(k - 1, -1, -1):
            # prec=0 preserves the entire integer, independently of ctx.prec.
            q = (round(mu[k][j]) if ctx._fixed_precision
                 else int(ctx.nint(mu[k][j], prec=0)))
            if q:
                for i in range(n):
                    transform[i][k] -= q * transform[i][j]
                # Subtracting an earlier basis vector leaves the orthogonal
                # lengths unchanged; only these projection coefficients move.
                for i in range(j):
                    mu[k][i] -= q * mu[j][i]
                mu[k][j] -= q
        coefficient = mu[k][k - 1]
        if lengths[k] >= (delta - coefficient ** 2) * lengths[k - 1]:
            k += 1
        else:
            # Update the two swapped orthogonal vectors and the projections
            # of later vectors, rather than rebuilding the entire Gram matrix.
            previous, current = lengths[k - 1], lengths[k]
            combined = ctx.fsum((current, coefficient ** 2 * previous))
            new_coefficient = coefficient * previous / combined
            lengths[k] = previous * (current / combined)
            lengths[k - 1] = combined
            for i in range(k - 1):
                mu[k][i], mu[k - 1][i] = mu[k - 1][i], mu[k][i]
            for i in range(k + 1, n):
                old = mu[i][k]
                mu[i][k] = mu[i][k - 1] - coefficient * old
                mu[i][k - 1] = old + new_coefficient * mu[i][k]
            mu[k][k - 1] = new_coefficient
            for row in transform:
                row[k], row[k - 1] = row[k - 1], row[k]
            orientation = -orientation
            k = max(1, k - 1)
    if orientation < 0:
        for row in transform:
            row[-1] = -row[-1]
    return tuple(tuple(row) for row in transform)


def _gram_schmidt(gram, summation=sum):
    """Compute Gram-Schmidt coefficients and squared lengths from a Gram matrix."""
    n = len(gram)
    mu = [[0] * n for _ in range(n)]
    lengths = []
    for i in range(n):
        length = gram[i][i] - summation(
            mu[i][j] ** 2 * lengths[j] for j in range(i))
        if length <= 0:
            raise ValueError("Y must be positive definite")
        lengths.append(length)
        for k in range(i + 1, n):
            mu[k][i] = (gram[k][i] - summation(
                mu[k][j] * mu[i][j] * lengths[j]
                for j in range(i))) / length
    return mu, lengths


def _is_reduced(gram, transform, delta):
    # Binary inputs have power-of-two denominators. A common positive scale
    # leaves the reduction conditions unchanged and makes the transform integral.
    denominator = max(x.denominator for row in gram for x in row)
    integral = [[int(x.numerator) * (int(denominator) // int(x.denominator))
                 for x in row] for row in gram]
    reduced = _congruence(integral, transform)
    mu, lengths = _gram_schmidt([[MPQ(x) for x in row] for row in reduced])
    n = len(gram)
    return (all(abs(mu[i][j]) <= MPQ(1, 2)
                for i in range(n) for j in range(i))
            and all(lengths[i] >= (delta - mu[i][i - 1] ** 2) * lengths[i - 1]
                    for i in range(1, n)))


def _congruence(gram, transform, summation=sum):
    """Apply an integer change of basis without converting its coefficients."""
    n = len(gram)
    right = [[summation(gram[i][k] * transform[k][j] for k in range(n))
              for j in range(n)] for i in range(n)]
    return [[summation(transform[k][i] * right[k][j] for k in range(n))
             for j in range(n)] for i in range(n)]

"""
Filename: eigen_markov.py

Author: Daisuke Oyama

This file contains some routines specialized for stochastic (Markov) matrices.

stoch_eig : returns the stochastic eigenvector (stationary distribution)
    of an irreducible stochastic matrix P, i.e., the stochastic vector x
    with x P = x
gth_solve : returns the (normalized) nontrivial solution to x A = 0 for
    an irreducible transition rate matrix A


Stochastic matrices
...................

The routine ``stoch_eig`` returns the stochastic eigenvector of an
irreducible stochastic matrix *P*, i.e., the stochastic vector `x` with
`x P = x` (a stochastic matrix, or Markov transition matrix, is a real
square matrix whose entries are nonnegative and rows sum to one).
Internally, the routine passes the input to the ``gth_solve`` routine.

The routine ``gth_solve`` solves for a (normalized) nontrivial solution
to `x A = 0` for an irreducible transition rate matrix *A* (a transition
rate matrix is a real square matrix whose off-diagonal entries are
nonnegative and rows sum to zero), by using the Grassmann-Taksar-Heyman
(GTH) algorithm [1]_, a numerically stable variant of Gaussian
elimination.


Notes
-----
For stochastic matrices, see en.wikipedia.org/wiki/Stochastic_matrix
For transition rate matrices, see en.wikipedia.org/wiki/Transition_rate_matrix
For irreducible matrices, see
en.wikipedia.org/wiki/Perron-Frobenius_theorem#Classification_of_matrices

If the input matrix is not irreducible, ``stoch_eig`` (``gth_solve``)
returns the solution corresponding to the irreducible class of indices
that contains the first recurrent index.

For the GTH algorithm, see the excellent lecture notes [2]_.

References
----------
.. [1] W. K. Grassmann, M. I. Taksar and D. P. Heyman, "Regenerative
   Analysis and Steady State Distributions for Markov Chains,"
   Operations Research (1985), 1107-1116.
.. [2] W. J. Stewart, "Performance Modelling and Markov Chains,"
   www.sti.uniurb.it/events/sfm07pe/slides/Stewart_1.pdf

"""
from ..libmp.backend import xrange
from .eigen import defun


@defun
def stoch_eig(ctx, P, overwrite=False):
    r"""
    This routine returns the stochastic eigenvector (stationary
    probability distribution vector) of an irreducible stochastic matrix
    *P*, i.e., the solution to `x P = x`, normalized so that its 1-norm
    equals one. Internally, the routine passes the input to the
    ``gth_solve`` routine.

    Parameters
    ----------
    P : matrix of shape (n, n)
        Stochastic matrix.
    overwrite : bool, optional(default=False)
        If True, allows modification of P which may improve performance;
        if False, P is not modified.

    Returns
    -------
    matrix of shape (1, n)
        Stochastic eigenvalue (stationary distribution) of P, i.e., the
        solution to x P = x, normalized so that its 1-norm equals one.

    Examples
    --------
    >>> import mpmath as mp
    >>> from eigen_markov import stoch_eig
    >>> P = mp.matrix([[0.9, 0.075, 0.025], [0.15, 0.8, 0.05], [0.25, 0.25, 0.5]])
    >>> x = mp.mp.stoch_eig(P)
    >>> print x
    [0.625  0.3125  0.0625]
    >>> print x * P
    [0.625  0.3125  0.0625]

    """
    # In fact, stoch_eig, which for the user is a routine to solve
    # x P = x, or x (P - I) = 0, is just another name of the function
    # gth_solve, which solves x A = 0, where the GTH algorithm,
    # the algorithm used there, does not use the actual values of
    # the diagonals of A, under the assumption that
    # A_{ii} = \sum_{j \neq i} A_{ij}, and therefore,
    # gth_solve(P-I) = gth_solve(P), so that it is irrelevant whether to
    # pass P or P - I to gth_solve.
    return ctx.gth_solve(P, overwrite=overwrite)


@defun
def gth_solve(ctx, A, overwrite=False):
    r"""
    This routine computes a nontrivial solution of a linear equation
    system of the form `x A = 0`, where *A* is an irreducible transition
    rate matrix, by using the Grassmann-Taksar-Heyman (GTH) algorithm, a
    variant of Gaussian elimination. The solution is normalized so that
    its 1-norm equals one.

    Parameters
    ----------
    A : matrix of shape (n, n)
        Transition rate matrix.
    overwrite : bool, optional(default=False)
        If True, allows modification of A which may improve performance;
        if False, A is not modified.

    Returns
    -------
    matrix of shape (1, n)
        Non-zero solution to x A = 0, normalized so that its 1-norm
        equals one.

    Examples
    --------
    >>> import mpmath as mp
    >>> from eigen_markov import gth_solve
    >>> A = mp.matrix([[-0.1, 0.075, 0.025], [0.15, -0.2, 0.05], [0.25, 0.25, -0.5]])
    >>> x = mp.mp.gth_solve(A)
    >>> print x
    [0.625  0.3125  0.0625]
    >>> print mp.chop(x*A)
    [0.0  0.0  0.0]

    """
    if not isinstance(A, ctx.matrix):
        A = ctx.matrix(A)
    elif not overwrite:
        A = A.copy()

    n, m = A.rows, A.cols

    if n != m:
        raise ValueError('matrix must be square')

    x = ctx.zeros(1, n)

    # === Reduction === #
    for j in xrange(n-1):
        scale = ctx.fsum(A[j, j+1:n])
        if scale <= 0:
            # Only consider the leading principal minor of size j+1,
            # which is irreducible
            n = j+1
            break
        A[j+1:n, j] /= scale

        for i in xrange(j+1, n):
            for k in xrange(j+1, n):
                A[k, i] += A[j, i] * A[k, j]

    # === Backward substitution === #
    x[n-1] = 1
    for i in xrange(n-2, -1, -1):
        x[i] = ctx.fsum((x[j] * A[j, i] for j in xrange(i+1, n)))

    # === Normalization === #
    x /= ctx.fsum(x)

    return x

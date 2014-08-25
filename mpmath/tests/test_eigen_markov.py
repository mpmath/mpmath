from __future__ import division, print_function

import mpmath as mp

VERBOSE = 0


# Generate test cases
def KMR_Markov_matrix_sequential(N, p, epsilon):
    """
    Generate a Markov matrix arising from a certain game-theoretic model

    Parameters
    ----------
    N : int

    p : float
        Between 0 and 1

    epsilon : float
        Between 0 and 1

    Returns
    -------
    P : matrix of shape (N+1, N+1)

    """
    P = mp.zeros(N+1, N+1)
    P[0, 0], P[0, 1] = 1 - epsilon * (1/2), epsilon * (1/2)
    for n in range(1, N):
        P[n, n-1] = \
            (n/N) * (epsilon * (1/2) +
                     (1 - epsilon) * (((n-1)/(N-1) < p) + ((n-1)/(N-1) == p) * (1/2))
                     )
        P[n, n+1] = \
            ((N-n)/N) * (epsilon * (1/2) +
                         (1 - epsilon) * ((n/(N-1) > p) + (n/(N-1) == p) * (1/2))
                         )
        P[n, n] = 1 - P[n, n-1] - P[n, n+1]
    P[N, N-1], P[N, N] = epsilon * (1/2), 1 - epsilon * (1/2)
    return P


# Stochastic matrix instances
Ps = []

Ps.append(
    mp.matrix([[0.9 , 0.075, 0.025],
               [0.15, 0.8  , 0.05 ],
               [0.25, 0.25 , 0.5  ]])
)
Ps.append(KMR_Markov_matrix_sequential(N=3, p=1./3, epsilon=1e-14))
Ps.append(KMR_Markov_matrix_sequential(N=27, p=1./3, epsilon=1e-2))

# Transition rate matrix instances
As = [P.copy() for P in Ps]
for A in As:
    for i in range(A.rows):
        A[i, i] = -mp.fsum((A[i, j] for j in range(A.cols) if j != i))


def run_stoch_eig(P, verbose=0):
    """
    stoch_eig returns a stochastic vector x such that x P = x
    for an irreducible stochstic matrix P.
    """
    if verbose > 1:
        print("original matrix (stoch_eig):\n", P)

    x = mp.stoch_eig(P)

    if verbose > 1:
        print("x\n", x)

    eps = mp.exp(0.8 * mp.log(mp.eps))  # From test_eigen.py

    # x is a left eigenvector of P with eigenvalue unity
    err0 = mp.norm(x*P-x, p=1)
    if verbose > 0:
        print("|xP - x| (stoch_eig):", err0)
    assert err0 < eps

    # x is a nonnegative vector
    if verbose > 0:
        print("min(x) (stoch_eig):", min(x))
    assert min(x) >= 0 - eps

    # 1-norm of x is one
    err1 = mp.fabs(mp.norm(x, p=1) - 1)
    if verbose > 0:
        print("||x| - 1| (stoch_eig):", err1)
    assert err1 < eps


def run_gth_solve(A, verbose=0):
    """
    gth_solve returns a stochastic vector x such that x A = 0
    for an irreducible transition rate matrix A.
    """
    if verbose > 1:
        print("original matrix (gth_solve):\n", A)

    x = mp.gth_solve(A)

    if verbose > 1:
        print("x\n", x)

    eps = mp.exp(0.8 * mp.log(mp.eps))  # test_eigen.py

    # x is a solution to x A = 0
    err0 = mp.norm(x*A, p=1)
    if verbose > 0:
        print("|xA| (gth_solve):", err0)
    assert err0 < eps

    # x is a nonnegative vector
    if verbose > 0:
        print("min(x) (gth_solve):", min(x))
    assert min(x) >= 0 - eps

    # 1-norm of x is one
    err1 = mp.fabs(mp.norm(x, p=1) - 1)
    if verbose > 0:
        print("||x| - 1| (gth_solve):", err1)
    assert err1 < eps


#######################


def test_stoch_eig_fixed_matrix():
    for P in Ps:
        run_stoch_eig(P, verbose=VERBOSE)


def test_gth_solve_fixed_matrix():
    for A in As:
        run_gth_solve(A, verbose=VERBOSE)


def test_stoch_eig_randmatrix():
    N = 5

    for j in range(10):
        P = mp.randmatrix(N, N)

        for i in range(N):
            P[i, :] /= mp.fsum(P[i, :])

        run_stoch_eig(P, verbose=VERBOSE)


def test_gth_solve_randmatrix():
    N = 5

    for j in range(10):
        A = mp.randmatrix(N, N)

        for i in range(N):
            A[i, :] /= mp.fsum(A[i, :])
            A[i, i] = -mp.fsum((A[i, j] for j in range(N) if j != i))

        run_gth_solve(A, verbose=VERBOSE)


def test_stoch_eig_high_prec():
    n = 1e-100
    with mp.workdps(100):
        P = mp.matrix([[1-3*(mp.exp(n)-1), 3*(mp.exp(n)-1)],
                       [mp.exp(n)-1      , 1-(mp.exp(n)-1)]])

    run_stoch_eig(P, verbose=VERBOSE)


def test_gth_solve_high_prec():
    n = 1e-100
    with mp.workdps(100):
        P = mp.matrix([[-3*(mp.exp(n)-1), 3*(mp.exp(n)-1)],
                       [mp.exp(n)-1     , -(mp.exp(n)-1) ]])

    run_gth_solve(P, verbose=VERBOSE)


def test_stoch_eig_fp():
    P = mp.fp.matrix([[0.9 , 0.075, 0.025],
                      [0.15, 0.8  , 0.05 ],
                      [0.25, 0.25 , 0.5  ]])
    x_expected = mp.fp.matrix([[0.625, 0.3125, 0.0625]])
    x = mp.fp.stoch_eig(P)
    eps = mp.exp(0.8 * mp.log(mp.eps))  # test_eigen.py
    err0 = mp.norm(x-x_expected, p=1)
    assert err0 < eps


def test_gth_solve_fp():
    P = mp.fp.matrix([[-0.1, 0.075, 0.025],
                      [0.15, -0.2 , 0.05 ],
                      [0.25, 0.25 , -0.5 ]])
    x_expected = mp.fp.matrix([[0.625, 0.3125, 0.0625]])
    x = mp.fp.gth_solve(P)
    eps = mp.exp(0.8 * mp.log(mp.eps))  # test_eigen.py
    err0 = mp.norm(x-x_expected, p=1)
    assert err0 < eps


def test_stoch_eig_iv():
    P = mp.iv.matrix([[0.9 , 0.075, 0.025],
                      [0.15, 0.8  , 0.05 ],
                      [0.25, 0.25 , 0.5  ]])
    x_expected = mp.matrix([[0.625, 0.3125, 0.0625]])
    x_iv = mp.iv.stoch_eig(P)
    for value, interval in zip(x_expected, x_iv):
        assert value in interval


def test_gth_solve_iv():
    P = mp.iv.matrix([[-0.1, 0.075, 0.025],
                      [0.15, -0.2 , 0.05 ],
                      [0.25, 0.25 , -0.5 ]])
    x_expected = mp.matrix([[0.625, 0.3125, 0.0625]])
    x_iv = mp.iv.gth_solve(P)
    for value, interval in zip(x_expected, x_iv):
        assert value in interval

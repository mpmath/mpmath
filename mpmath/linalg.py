# TODO:
# *implement high-level qr() and lu()
# *cache LU decompositon
# *test unitvector
# *iterative solving
# *condition numbers
# *iterative improving of solution

from __future__ import division

from mptypes import *
from matrices import matrix, swap_row, extend, mnorm_1, norm_p
from copy import copy

def LU_decomp(A):
    """
    LU-factorization of a n*n matrix using the Gauss algorithm.
    Returns L and U in one matrix and the pivot indices.
    """
    if not A.rows == A.cols:
        raise ValueError('need n*n matrix')
    tol = mnorm_1(A) * eps # each pivot element has to be bigger than this
    n = A.rows
    p = [None]*(n - 1)
    for j in xrange(n - 1):
        # pivoting, choose max(abs(reciprocal row sum)*abs(pivot element))
        biggest = 0
        for k in xrange(j, n):
            current = 1/sum([abs(A[k,l]) for l in xrange(j, n)]) * abs(A[k,j])
            if current > biggest: # TODO: what if equal?
                biggest = current
                p[j] = k
        # swap rows according to p
        swap_row(A, j, p[j])
        if abs(A[j,j]) < tol:
            raise ZeroDivisionError('matrix is numerically singular')
        # calculate elimination factors and add rows
        for i in xrange(j + 1, n):
            A[i,j] /= A[j,j]
            for k in xrange(j + 1, n):
                A[i,k] -= A[i,j]*A[j,k]
    return A, p

def L_solve(L, b, p=None):
    """
    Solve the lower part of a LU factorized matrix for y.
    """
    L.rows == L.cols, 'need n*n matrix'
    n = L.rows
    assert len(b) == n
    if p: # swap b according to p
        for k in xrange(0, len(p)):
            swap_row(b, k, p[k])
    # solve
    for i in xrange(1, n):
        for j in xrange(i):
            b[i] -= L[i,j] * b[j]
    return b

def U_solve(U, y):
    """
    Solve the upper part of a LU factorized matrix for x.
    """
    assert U.rows == U.cols, 'need n*n matrix'
    n = U.rows
    assert len(y) == n
    x = copy(y)
    for i in xrange(n - 1, -1, -1):
        for j in xrange(i + 1, n):
            x[i] -= U[i,j] * x[j]
        x[i] /= U[i,i]
    return x

@extraprec(10)
def lu_solve(A, b):
    """
    Ax = b => x

    Solve a determined or overdetermined linear equations system.
    Fast LU decomposition is used, which is less accurate than QR decomposition
    (especially for overdetermined systems), but it's twice as efficient.
    Use qr_solve if you want more precison or have to solve a very ill-
    conditioned system.
    """
    # do not overwrite A nor b
    A, b = matrix(A).copy(), matrix(b).copy()
    if A.rows < A.cols:
        raise ValueError('cannot solve underdetermined system')
    if A.rows > A.cols:
        # use least-squares method if overdetermined
        # (this increases errors)
        AT = A.T
        A = AT * A
        b = AT * b
        return cholesky_solve(A, b)
    else:
        # LU factorization
        A, p = LU_decomp(A)
        b = L_solve(A, b, p)
        x = U_solve(A, b)
        return x

def unitvector(n, i):
    """
    Return the i-th n-dimensional unit vector.
    """
    assert 0 < i <= n, 'this unit vector does not exist'
    return [0]*(i-1) + [1] + [0]*(n-i)

@extraprec(10)
def inverse(A):
    """
    Calculate the inverse of a matrix.

    If you want to solve an equation system Ax = b, it's recommended to use
    solve(A, b) instead, it's about 3 times more efficient.
    """
    # do not overwrite A
    A = matrix(A).copy()
    n = A.rows
    # get LU factorisation
    A, p = LU_decomp(A)
    cols = []
    # calculate unit vectors and solve corresponding system to get columns
    for i in xrange(1, n + 1):
        e = unitvector(n, i)
        y = L_solve(A, e, p)
        cols.append(U_solve(A, y))
    # convert columns to matrix
    inv = []
    for i in xrange(n):
        row = []
        for j in xrange(n):
            row.append(cols[j][i])
        inv.append(row)
    return matrix(inv)

def householder(A):
    """
    (A|b) -> H, p, x, res
    """
    assert isinstance(A, matrix)
    m = A.rows
    n = A.cols
    assert m >= n - 1
    # calculate Householder matrix
    p = []
    for j in xrange(0, n - 1):
        s = 0.
        for i in xrange(j, m):
            s += (A[i,j])**2
        if not abs(s) > eps:
            raise ValueError('matrix is numerically singular')
        if A[j,j] < 0:
            p.append(sqrt(s))
        else:
            p.append(-sqrt(s))
        kappa = s - p[j] * A[j,j]
        A[j,j] -= p[j]
        for k in xrange(j+1, n):
            y = 0.
            for i in xrange(j, m):
                y += A[i,j] * A[i,k]
            y /= kappa
            for i in xrange(j, m):
                A[i,k] -= A[i,j] * y
    # solve Rx = c1
    x = []
    for i in xrange(n - 1):
        x.append(A[i,n - 1])
    for i in xrange(n - 2, -1, -1):
        for j in xrange(i + 1, n - 1):
            x[i] -= A[i,j] * x[j]
        x[i] /= p[i]
    # calculate residual
    if not m == n - 1:
        r = []
        for i in xrange(m - n + 1):
            r.append(A[m-1-i, n-1])
    else:
        # determined system, residual should be 0
        r = [0]*m
    return A, p, x, r

def residual(A, x, b):
    """
    Calculate the residual of a solution to a linear equation system.

    r = A*x - b for A*x = b
    """
    oldprec = mp.prec
    try:
        mp.prec *= 2
        A, x, b = matrix(A), matrix(x), matrix(b)
        return A*x - b
    finally:
        mp.prec = oldprec

@extraprec(10)
def qr_solve(A, b, norm=lambda x: norm_p(x, 2)):
    """
    Ax = b => x, ||Ax - b||

    Solve a determined or overdetermined linear equations system and
    calculate the norm of the residual (error).
    QR decompostion using Householder factorization is applied, which gives very
    accurate results even for ill-conditioned matrices. qr_solve is twice as
    efficient.
    """
    # do not overwrite A nor b
    A, b = matrix(A).copy(), matrix(b).copy()
    if A.rows < A.cols:
        raise ValueError('cannot solve underdetermined system')
    H, p, x, r = householder(extend(A, b))
    res = norm(r)
    # calculate residual "manually" for determined systems
    if res == 0:
        res = norm(residual(A, x, b))
    return matrix(x), res

def cholesky(A):
    """
    Cholesky decompositon of a symmetric positive-definite matrix.

    Can be used to solve linear equation systems twice as efficient compared
    to LU decomposition or to test whether A is positive-definite.

    A = L * L.T
    Only L (the lower part) is returned.
    """
    assert isinstance(A, matrix)
    if not A.rows == A.cols:
        raise ValueError('need n*n matrix')
    n = A.rows
    L = matrix(n)
    for j in xrange(n):
        s = A[j,j] - sum((L[j,k]**2 for k in xrange(j)))
        if s < eps:
            raise ValueError('matrix not positive-definite')
        L[j,j] = sqrt(s)
        for i in xrange(j, n):
            L[i,j] = (A[i,j] - sum((L[i,k] * L[j,k] for k in xrange(j)))) \
                     / L[j,j]
    return L

@extraprec(10)
def cholesky_solve(A, b):
    """
    Ax = b => x

    Solve a symmetric positive-definite linear equation system.
    This is twice as efficient as lu_solve.

    Typical use cases:
    * A.T*A
    * Hessian matrix
    * differential equations
    """
    # do not overwrite A nor b
    A, b = matrix(A).copy(), matrix(b).copy()
    if A.rows !=  A.cols:
        raise ValueError('can only solve determined system')
    # Cholesky factorization
    L = cholesky(A)
    # solve
    n = L.rows
    assert len(b) == n
    for i in xrange(n):
        b[i] -= sum((L[i,j] * b[j] for j in xrange(i)))
        b[i] /= L[i,i]
    x = U_solve(L.T, b)
    return x

@extraprec(10)
def det(A):
    """
    Calculate the determinant of a matrix.
    """
    # do not overwrite A
    A = matrix(A).copy()
    # use LU factorization to calculate determinant
    try:
        R, p = LU_decomp(A)
    except ZeroDivisionError:
        return 0
    z = 1
    for i, e in enumerate(p):
        if i != e:
            z *= -1
    for i in xrange(A.rows):
        z *= R[i,i]
    return z

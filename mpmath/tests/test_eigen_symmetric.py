#!/usr/bin/python
# -*- coding: utf-8 -*-

from mpmath import mp
from mpmath import libmp

xrange = libmp.backend.xrange

def run_eigsy(A, verbose = False):
    if verbose:
        print("original matrix:\n", str(A))

    D, Q = mp.eigsy(A)
    B = Q * mp.diag(D) * Q.transpose()
    C = A - B
    E = Q * Q.transpose() - mp.eye(A.rows)

    if verbose:
        print("eigenvalues:\n", D)
        print("eigenvectors:\n", Q)

    NC = mp.mnorm(C)
    NE = mp.mnorm(E)

    if verbose:
        print("difference:", NC, "\n", C, "\n")
        print("difference:", NE, "\n", E, "\n")

    assert NC < mp.exp( 0.8 * mp.log(mp.eps))
    assert NE < mp.exp( 0.8 * mp.log(mp.eps))

    return NC

def run_eighe(A, verbose = False):
    if verbose:
        print("original matrix:\n", str(A))

    D, Q = mp.eighe(A)
    B = Q * mp.diag(D) * Q.transpose_conj()
    C = A - B
    E = Q * Q.transpose_conj() - mp.eye(A.rows)

    if verbose:
        print("eigenvalues:\n", D)
        print("eigenvectors:\n", Q)

    NC = mp.mnorm(C)
    NE = mp.mnorm(E)

    if verbose:
        print("difference:", NC, "\n", C, "\n")
        print("difference:", NE, "\n", E, "\n")

    assert NC < mp.exp( 0.8 * mp.log(mp.eps))
    assert NE < mp.exp( 0.8 * mp.log(mp.eps))

    return NC

def irandmatrix(n, range = 10):
    """
    random matrix with integer entries
    """
    A = mp.matrix(n, n)
    for i in xrange(n):
        for j in xrange(n):
            A[i,j]=int( (2 * mp.rand() - 1) * range)
    return A

#######################

def test_eigsy_fixed_matrix():
    A = mp.matrix([[2, 3], [3, 5]])
    run_eigsy(A)

    A = mp.matrix([[7, -11], [-11, 13]])
    run_eigsy(A)

    A = mp.matrix([[2, 11, 7], [11, 3, 13], [7, 13, 5]])
    run_eigsy(A)

    A = mp.matrix([[2, 0, 7], [0, 3, 1], [7, 1, 5]])
    run_eigsy(A)


def test_eighe_fixed_matrix():
    A = mp.matrix([[2, 3], [3, 5]])
    run_eighe(A)

    A = mp.matrix([[7, -11], [-11, 13]])
    run_eighe(A)

    A = mp.matrix([[2, 11, 7], [11, 3, 13], [7, 13, 5]])
    run_eighe(A)

    A = mp.matrix([[2, 0, 7], [0, 3, 1], [7, 1, 5]])
    run_eighe(A)

    #

    A = mp.matrix([[2, 3+7j], [3-7j, 5]])
    run_eighe(A)

    A = mp.matrix([[2, -11j, 0], [+11j, 3, 29j], [0, -29j, 5]])
    run_eighe(A)

    A = mp.matrix([[2, 11 + 17j, 7 + 19j], [11 - 17j, 3, -13 + 23j], [7 - 19j, -13 - 23j, 5]])
    run_eighe(A)

def test_eigsy_randmatrix():
    N = 5

    for a in xrange(10):
        A=mp.randmatrix(N, N)

        for i in xrange(0, N):
            for j in xrange(i + 1, N):
                A[j,i] = A[i,j]

        run_eigsy(A)

def test_eighe_randmatrix():
    N = 5

    for a in xrange(10):
        A = mp.randmatrix(N, N) + 1j * mp.randmatrix(N, N)

        for i in xrange(0, N):
            A[i,i] = mp.re(A[i,i])
            for j in xrange(i + 1, N):
                A[j,i] = mp.conj(A[i,j])

        run_eighe(A)

def test_eigsy_irandmatrix():
    N = 4
    R = 4

    for a in xrange(10):
        A=irandmatrix(N, R)

        for i in xrange(0, N):
            for j in xrange(i + 1, N):
                A[j,i] = A[i,j]

        run_eigsy(A)

def test_eighe_irandmatrix():
    N = 4
    R = 4

    for a in xrange(10):
        A=irandmatrix(N, R) + 1j * irandmatrix(N, R)

        for i in xrange(0, N):
            A[i,i] = mp.re(A[i,i])
            for j in xrange(i + 1, N):
                A[j,i] = mp.conj(A[i,j])

        run_eighe(A)


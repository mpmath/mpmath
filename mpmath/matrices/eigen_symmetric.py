#!/usr/bin/python
# -*- coding: utf-8 -*-

##################################################################################################
#     module for the symmetric eigenvalue problem
#       Copyright 2013 Timo Hartmann (thartmann15 at googlemail.com)
#
# todo:
#  - implement balancing
#  - svd
#
##################################################################################################

"""
  The symmetric eigenvalue problem.
  ---------------------------------

  This file contains routines for the symmetric eigenvalue problem.

  high level routines:

    eigsy : real symmetric (ordinary) eigenvalue problem
    eighe : complex hermitian (ordinary) eigenvalue problem
    eigh  : unified interface for eigsy and eighe

  low level routines:

    r_sy_tridiag : reduction of real symmetric matrix to real symmetric tridiagonal matrix
    c_he_tridiag_0 : reduction of complex hermitian matrix to real symmetric tridiagonal matrix
    c_he_tridiag_1 : auxiliary routine to c_he_tridiag_0
    c_he_tridiag_2 : auxiliary routine to c_he_tridiag_0
    tridiag_eigen : solves the real symmetric tridiagonal matrix eigenvalue problem
"""

from ..libmp.backend import xrange

class EigenSymmetric(object):
  pass

def defun(f):
  setattr(EigenSymmetric,f.__name__,f)


def r_sy_tridiag(ctx,A,D,E,calc_ev=True):
  """
  This routine transforms a real symmetric matrix A to a real symmetric
  tridiagonal matrix T using an orthogonal similarity transformation:
        Q' * A * Q = T     (here ' denotes the matrix transpose).
  The orthogonal matrix Q is build up from Householder reflectors.

  parameters:
    A         (input/output) On input, A contains the real symmetric matrix of
              dimension (n,n). On output, if calc_ev is true, A contains the
              orthogonal matrix Q, otherwise A is destroyed.

    D         (output) real array of length n, contains the diagonal elements
              of the tridiagonal matrix

    E         (output) real array of length n, contains the offdiagonal elements
              of the tridiagonal matrix in E[0:(n-1)] where is the dimension of
              the matrix A. E[n-1] is undefined.

    calc_ev   (input) If calc_ev is true, this routine explicitly calculates the
              orthogonal matrix Q which is then returned in A. If calc_ev is
              false, Q is not explicitly calculated resulting in a shorter run time.

  This routine is a python translation of the fortran routine tred2.f in the
  software library EISPACK (see netlib.org) which itself is based on the algol
  procedure tred2 described in:
    - Num. Math. 11, p.181-195 (1968) by Martin, Reinsch and Wilkonson
    - Handbook for auto. comp., Vol II, Linear Algebra, p.212-226 (1971)

  For a good introduction to Householder reflections, see also
    Stoer, Bulirsch - Introduction to Numerical Analysis.
  """

  # note : the vector v of the i-th houshoulder reflector is stored in a[(i+1):,i]
  #        whereas v/<v,v> is stored in a[i,(i+1):]

  n=A.rows
  for i in xrange(n-1,0,-1):

                                    # scale the vector

    scale=0
    for k in xrange(0,i):
      scale+=abs(A[k,i])

    scale_inv=0
    if scale!=0:
      scale_inv=1/scale

                                    # sadly there are floating point numbers not equal to zero whose reciprocal is infinity

    if i==1 or scale==0 or ctx.isinf(scale_inv):
      E[i]=A[i-1,i]                 # nothing to do
      D[i]=0
      continue

                                    # calculate parameters for housholder transformation

    H=0
    for k in xrange(0,i):
      A[k,i]*=scale_inv
      H+=A[k,i]*A[k,i]

    F=A[i-1,i]
    G=ctx.sqrt(H)
    if F>0:
      G=-G
    E[i]=scale*G
    H-=F*G
    A[i-1,i]=F-G
    F=0

                                   # apply housholder transformation

    for j in xrange(0,i):
      if calc_ev:
        A[i,j]=A[j,i]/H

      G=0                          # calculate A*U
      for k in xrange(0,j+1):
        G+=A[k,j]*A[k,i]
      for k in xrange(j+1,i):
        G+=A[j,k]*A[k,i]

      E[j]=G/H                     # calculate P
      F+=E[j]*A[j,i]

    HH=F/(2*H)

    for j in xrange(0,i):          # calculate reduced A
      F=A[j,i]
      G=E[j]-HH*F                  # calculate Q
      E[j]=G

      for k in xrange(0,j+1):
        A[k,j]-=F*E[k]+G*A[k,i]

    D[i]=H

  for i in xrange(1,n):            # better for compatibility
    E[i-1]=E[i]
  E[n-1]=0

  if calc_ev:
    D[0]=0
    for i in xrange(0,n):
      if D[i]!=0:
        for j in xrange(0,i):      # accumulate transformation matrices
          G=0
          for k in xrange(0,i):
            G+=A[i,k]*A[k,j]
          for k in xrange(0,i):
            A[k,j]-=G*A[k,i]

      D[i]=A[i,i]
      A[i,i]=1

      for j in xrange(0,i):
        A[j,i]=A[i,j]=0
  else:
    for i in xrange(0,n):
      D[i]=A[i,i]










def c_he_tridiag_0(ctx,A,D,E,T):
  """
  This routine transforms a complex hermitian matrix A to a real symmetric
  tridiagonal matrix T using an unitary similarity transformation:
        Q' * A * Q = T     (here ' denotes the hermitian matrix transpose,
                            i.e. transposition und conjugation).
  The unitary matrix Q is build up from Householder reflectors and
  an unitary diagonal matrix.

  parameters:
    A         (input/output) On input, A contains the complex hermitian matrix
              of dimension (n,n). On output, A contains the unitary matrix Q
              in compressed form.

    D         (output) real array of length n, contains the diagonal elements
              of the tridiagonal matrix.

    E         (output) real array of length n, contains the offdiagonal elements
              of the tridiagonal matrix in E[0:(n-1)] where is the dimension of
              the matrix A. E[n-1] is undefined.

    T         (output) complex array of length n, contains a unitary diagonal
              matrix.

  This routine is a python translation (in slightly modified form) of the fortran
  routine htridi.f in the software library EISPACK (see netlib.org) which itself
  is a complex version of the algol procedure tred1 described in:
    - Num. Math. 11, p.181-195 (1968) by Martin, Reinsch and Wilkonson
    - Handbook for auto. comp., Vol II, Linear Algebra, p.212-226 (1971)

  For a good introduction to Householder reflections, see also
    Stoer, Bulirsch - Introduction to Numerical Analysis.
  """

  n=A.rows
  T[n-1]=1
  for i in xrange(n-1,0,-1):

                                      # scale the vector

    scale=0
    for k in xrange(0,i):
      scale+=abs(ctx.re(A[k,i]))+abs(ctx.im(A[k,i]))

    scale_inv=0
    if scale!=0:
      scale_inv=1/scale

                                      # sadly there are floating point numbers not equal to zero whose reciprocal is infinity

    if scale==0 or ctx.isinf(scale_inv):
      E[i]=0
      D[i]=0
      T[i-1]=1
      continue

    if i==1:
      F=A[i-1,i]
      f=abs(F)
      E[i]=f
      D[i]=0
      if f!=0:
        T[i-1]=T[i]*F/f
      else:
        T[i-1]=T[i]
      continue

                                      # calculate parameters for housholder transformation

    H=0
    for k in xrange(0,i):
      A[k,i]*=scale_inv
      rr=ctx.re(A[k,i])
      ii=ctx.im(A[k,i])
      H+=rr*rr+ii*ii

    F=A[i-1,i]
    f=abs(F)
    G=ctx.sqrt(H)
    H+=G*f
    E[i]=scale*G
    if f!=0:
      F=F/f
      TZ=-T[i]*F      # T[i-1]=-T[i]*F, but we need T[i-1] as temporary storage
      G*=F
    else:
      TZ=-T[i]        # T[i-1]=-T[i]
    A[i-1,i]+=G
    F=0

                                     # apply housholder transformation

    for j in xrange(0,i):
      A[i,j]=A[j,i]/H

      G=0                            # calculate A*U
      for k in xrange(0,j+1):
        G+=ctx.conj(A[k,j])*A[k,i]
      for k in xrange(j+1,i):
        G+=A[j,k]*A[k,i]

      T[j]=G/H                       # calculate P
      F+=ctx.conj(T[j])*A[j,i]

    HH=F/(2*H)

    for j in xrange(0,i):            # calculate reduced A
      F=A[j,i]
      G=T[j]-HH*F                    # calculate Q
      T[j]=G

      for k in xrange(0,j+1):
        A[k,j]-=ctx.conj(F)*T[k]+ctx.conj(G)*A[k,i]
                                     # as we use the lower left part for storage
                                     # we have to use the transpose of the normal formula

    T[i-1]=TZ
    D[i]=H

  for i in xrange(1,n):              # better for compatibility
    E[i-1]=E[i]
  E[n-1]=0

  D[0]=0
  for i in xrange(0,n):
    zw=D[i]
    D[i]=ctx.re(A[i,i])
    A[i,i]=zw







def c_he_tridiag_1(ctx,A,T):
  """
  This routine forms the unitary matrix Q described in c_he_tridiag_0.

  parameters:
    A    (input/output) On input, A is the same matrix as delivered by
         c_he_tridiag_0. On output, A is set to Q.

    T    (input) On input, T is the same array as delivered by c_he_tridiag_0.

  """

  n=A.rows

  for i in xrange(0,n):
    if A[i,i]!=0:
      for j in xrange(0,i):
        G=0
        for k in xrange(0,i):
          G+=ctx.conj(A[i,k])*A[k,j]
        for k in xrange(0,i):
          A[k,j]-=G*A[k,i]

    A[i,i]=1

    for j in xrange(0,i):
      A[j,i]=A[i,j]=0

  for i in xrange(0,n):
    for k in xrange(0,n):
      A[i,k]*=T[k]




def c_he_tridiag_2(ctx,A,T,B):
  """
  This routine applied the unitary matrix Q described in c_he_tridiag_0
  onto the the matrix B, i.e. it forms Q*B.

  parameters:
    A    (input) On input, A is the same matrix as delivered by c_he_tridiag_0.

    T    (input) On input, T is the same array as delivered by c_he_tridiag_0.

    B    (input/output) On input, B is a complex matrix. On output B is replaced
         by Q*B.

  This routine is a python translation of the fortran routine htribk.f in the
  software library EISPACK (see netlib.org). See c_he_tridiag_0 for more
  references.
  """
  n=A.rows

  for i in xrange(0,n):
    for k in xrange(0,n):
      B[k,i]*=T[k]

  for i in xrange(0,n):
    if A[i,i]!=0:
      for j in xrange(0,n):
        G=0
        for k in xrange(0,i):
          G+=ctx.conj(A[i,k])*B[k,j]
        for k in xrange(0,i):
          B[k,j]-=G*A[k,i]





def tridiag_eigen(ctx,d,e,z=False):
  """
  This subroutine find the eigenvalues and the first components of the
  eigenvectors of a real symmetric tridiagonal matrix using the implicit
  QL method.

  parameters:

    d (input/output) real array of length n. on input, d contains the diagonal
      elements of the input matrix. on output, d contains the eigenvalues in
      ascending order.

    e (input) real array of length n. on input, e contains the offdiagonal
      elements of the input matrix in e[0:(n-1)]. On output, e has been
      destroyed.

    z (input/output) If z is equal to False, no eigenvectors will be computed.
      Otherwise on input z should have the format z[0:m,0:n] (i.e. a real or
      complex matrix of dimension (m,n) ). On output this matrix will be
      multiplied by the matrix of the eigenvectors (i.e. the columns of this
      matrix are the eigenvectors): z --> z*EV
      That means if z[i,j]={1 if j==j; 0 otherwise} on input, then on output
      z will contain the first m components of the eigenvectors. That means
      if m is equal to n, the i-th eigenvector will be z[:,i].

  This routine is a python translation (in slightly modified form) of the
  fortran routine imtql2.f in the software library EISPACK (see netlib.org)
  which itself is based on the algol procudure imtql2 desribed in:
   - num. math. 12, p. 377-383(1968) by matrin and wilkinson
   - modified in num. math. 15, p. 450(1970) by dubrulle
   - handbook for auto. comp., vol. II-linear algebra, p. 241-248 (1971)
  See also the routine gaussq.f in netlog.org or acm algorithm 726.
  """
  n=len(d)
  e[n-1]=0
  iterlim=2*ctx.dps

  for l in xrange(n):
    j=0
    while 1:
      m=l
      while 1:
        # look for a small subdiagonal element
        if m+1==n:
          break
        if abs(e[m])<=ctx.eps*(abs(d[m])+abs(d[m+1])):
          break
        m=m+1
      if m==l:
        break

      if j>=iterlim:
        raise RuntimeError("tridiag_eigen: no convergence to an eigenvalue after %d iterations" % iterlim)

      j+=1

      # form shift

      p=d[l]
      g=(d[l+1]-p)/(2*e[l])
      r=ctx.hypot(g,1)

      if(g<0):
        s=g-r
      else:
        s=g+r

      g=d[m]-p+e[l]/s

      s,c,p=1,1,0

      for i in xrange(m-1,l-1,-1):
        f=s*e[i]
        b=c*e[i]
        if abs(f)>abs(g):       # this here is a slight improvement also used in gaussq.f or acm algorithm 726.
          c=g/f
          r=ctx.hypot(c,1)
          e[i+1]=f*r
          s=1/r
          c=c*s
        else:
          s=f/g
          r=ctx.hypot(s,1)
          e[i+1]=g*r
          c=1/r
          s=s*c
        g=d[i+1]-p
        r=(d[i]-g)*s+2*c*b
        p=s*r
        d[i+1]=g+p
        g=c*r-b

        if not isinstance(z,bool):
          # calculate eigenvectors
          for w in xrange(z.rows):
            f=         z[w,i+1]
            z[w,i+1]=s*z[w,i  ]+c*f
            z[w,i  ]=c*z[w,i  ]-s*f

      d[l]=d[l]-p
      e[l]=g
      e[m]=0

  for ii in xrange(1,n):
    # order eigenvalues and eigenvectors
    i=ii-1
    k=i
    p=d[i]
    for j in xrange(ii,n):
      if d[j]>=p:
        continue
      k=j
      p=d[k]
    if k==i:
      continue
    d[k]=d[i]
    d[i]=p

    if not isinstance(z,bool):
      for w in xrange(z.rows):
        p=     z[w,i];
        z[w,i]=z[w,k];
        z[w,k]=p;

########################################################################################

@defun
def eigsy(ctx,A,eigvals_only=False,overwrite_a=False):
  """
    This routine solves the (ordinary) eigenvalue problem for a real symmetric
    square matrix A. Given A, an orthogonal matrix Q is calculated which
    diagonalizes A:

          Q' A Q = diag(E)

    Here diag(E) is a diagonal matrix whose diagonal is E.
    ' denotes the transpose.

    The columns of Q are the eigenvectors of A and E contains the eigenvalues:

          A Q[:,i]=E[i] Q[:,i]


    input:

      A: real matrix of format (n,n) which is symmetric
         (i.e. A=A' or A[i,j]=A[j,i])

      eigvals_only: if true, calculates only the eigenvalues E.
                    if false, calculates both eigenvectors and eigenvalues.

      overwrite_a: if true, allows modification of A which may improve
                   performance. if false, A is not modified.

    output:

      E: vector of format (n). contains the eigenvalues of A.

      Q: orthogonal matrix of format (n,n). contains the eigenvectors
         of A as columns.

    return value:

           E        if eigvals_only is true
         (E,Q)      if eigvals_only is false

    example:
      >>> from mpmath import mp
      >>> A = mp.matrix([[1,2],[2,3]])
      >>> E,Q = mp.eigsy(A)
      >>> mp.chop(A*Q[:,0]-E[0]*Q[:,0])
      matrix([['0.0'],['0.0']])


    see also: eighe, eigh
  """
  if not overwrite_a: A=A.copy()

  d=ctx.zeros(A.rows,1)
  e=ctx.zeros(A.rows,1)

  if eigvals_only:
    r_sy_tridiag(ctx,A,d,e,calc_ev=False)
    tridiag_eigen(ctx,d,e,False)
    return d
  else:
    r_sy_tridiag(ctx,A,d,e,calc_ev=True)
    tridiag_eigen(ctx,d,e,A)
    return (d,A)


@defun
def eighe(ctx,A,eigvals_only=False,overwrite_a=False):
  """
    This routine solves the (ordinary) eigenvalue problem for a complex
    hermitian square matrix A. Given A, an unitary matrix Q is calculated which
    diagonalizes A:

        Q' A Q = diag(E)

    Here diag(E) a is diagonal matrix whose diagonal is E.
    ' denotes the hermitian transpose (i.e. ordinary transposition and
    complex conjugation).

    The columns of Q are the eigenvectors of A and E contains the eigenvalues:

        A Q[:,i]=E[i] Q[:,i]


    input:

      A: complex matrix of format (n,n) which is hermitian
         (i.e. A=A' or A[i,j]=conj(A[j,i]))

      eigvals_only: if true, calculates only the eigenvalues E.
                    if false, calculates both eigenvectors and eigenvalues.

      overwrite_a: if true, allows modification of A which may improve
                   performance. if false, A is not modified.

    output:

      E: vector of format (n). contains the eigenvalues of A.

      Q: unitary matrix of format (n,n). contains the eigenvectors
         of A as columns.

    return value:

           E        if eigvals_only is true
         (E,Q)      if eigvals_only is false

    example:
      >>> from mpmath import mp
      >>> A = mp.matrix([[1,2+5j],[2-5j,3]])
      >>> E,Q = mp.eighe(A)
      >>> mp.chop(A*Q[:,0]-E[0]*Q[:,0])
      matrix([['0.0'],['0.0']])

    see also: eigsy, eigh
  """
  if not overwrite_a: A=A.copy()

  d=ctx.zeros(A.rows,1)
  e=ctx.zeros(A.rows,1)
  t=ctx.zeros(A.rows,1)

  if eigvals_only:
    c_he_tridiag_0(ctx,A,d,e,t)
    tridiag_eigen(ctx,d,e,False)
    return d
  else:
    c_he_tridiag_0(ctx,A,d,e,t)
    B=ctx.eye(A.rows)
    tridiag_eigen(ctx,d,e,B)
    c_he_tridiag_2(ctx,A,t,B)
    return (d,B)

@defun
def eigh(ctx,A,eigvals_only=False,overwrite_a=False):
  """
    "eigh" is a unified interface for "eigsy" and "eighe". Depending on
    whether A is real or complex the appropriate function is called.

    This routine solves the (ordinary) eigenvalue problem for a real symmetric
    or complex hermitian square matrix A. Given A, an orthogonal (A real) or
    unitary (A complex) matrix Q is calculated which diagonalizes A:

        Q' A Q = diag(E)

    Here diag(E) a is diagonal matrix whose diagonal is E.
    ' denotes the hermitian transpose (i.e. ordinary transposition and
    complex conjugation).

    The columns of Q are the eigenvectors of A and E contains the eigenvalues:

        A Q[:,i]=E[i] Q[:,i]

    input:

      A: a real or complex square matrix of format (n,n) which is symmetric
         (i.e. A[i,j]=A[j,i]) or hermitian (i.e. A[i,j]=conj(A[j,i])).

      eigvals_only: if true, calculates only the eigenvalues E.
                    if false, calculates both eigenvectors and eigenvalues.

      overwrite_a: if true, allows modification of A which may improve
                   performance. if false, A is not modified.

    output:

      E: vector of format (n). contains the eigenvalues of A.

      Q: an orthogonal or unitary matrix of format (n,n). contains the
         eigenvectors of A as columns.

    return value:

           E        if eigvals_only is true
         (E,Q)      if eigvals_only is false

    example:
      >>> from mpmath import mp
      >>> A = mp.matrix([[1,2],[2,3]])
      >>> E,Q = mp.eigh(A)
      >>> mp.chop(A*Q[:,0]-E[0]*Q[:,0])
      matrix([['0.0'],['0.0']])

      >>> A = mp.matrix([[1,2+5j],[2-5j,3]])
      >>> E,Q = mp.eigh(A)
      >>> mp.chop(A*Q[:,0]-E[0]*Q[:,0])
      matrix([['0.0'],['0.0']])

    see also: eigsy, eighe
  """

  iscomplex = any(type(x) is ctx.mpc for x in A)

  if iscomplex:
    return ctx.eighe(A,eigvals_only=eigvals_only,overwrite_a=overwrite_a)
  else:
    return ctx.eigsy(A,eigvals_only=eigvals_only,overwrite_a=overwrite_a)

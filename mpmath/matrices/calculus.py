# TODO: should use diagonalization-based algorithms

class MatrixCalculusMethods:

    def _exp_pade(ctx, a):
        """
        Exponential of a matrix using Pade approximants.

        See G. H. Golub, C. F. van Loan 'Matrix Computations',
        third Ed., page 572

        TODO:
         - find a good estimate for q
         - reduce the number of matrix multiplications to improve
           performance
        """
        def eps_pade(p):
            return ctx.mpf(2)**(3-2*p) * \
                ctx.factorial(p)**2/(ctx.factorial(2*p)**2 * (2*p + 1))
        q = 4
        extraq = 8
        while 1:
            if eps_pade(q) < ctx.eps:
                break
            q += 1
        q += extraq
        j = int(max(1, ctx.mag(ctx.mnorm(a,'inf'))))
        extra = q
        prec = ctx.prec
        ctx.dps += extra + 3
        try:
            a = a/2**j
            na = a.rows
            den = ctx.eye(na)
            num = ctx.eye(na)
            x = ctx.eye(na)
            c = ctx.mpf(1)
            for k in range(1, q+1):
                c *= ctx.mpf(q - k + 1)/((2*q - k + 1) * k)
                x = a*x
                cx = c*x
                num += cx
                den += (-1)**k * cx
            f = ctx.lu_solve_mat(den, num)
            for k in range(j):
                f = f*f
        finally:
            ctx.prec = prec
        return f*1

    def expm(ctx, A, method='taylor'):
        r"""
        Computes the matrix exponential of a square matrix `A`, which is defined
        by the power series

        .. math ::

            \exp(A) = I + A + \frac{A^2}{2!} + \frac{A^3}{3!} + \ldots

        With method='taylor', the matrix exponential is computed
        using the Taylor series. With method='pade', Pade approximants
        are used instead.

        **Examples**

        Basic examples::

            >>> from mpmath import (mp, expm, zeros, eye, j, hilbert, chop,
            ...                     mnorm, ones, matrix)
            >>> mp.pretty = True
            >>> expm(zeros(3))
            [1.0  0.0  0.0]
            [0.0  1.0  0.0]
            [0.0  0.0  1.0]
            >>> expm(eye(3))
            [2.718281828459045                0.0                0.0]
            [              0.0  2.718281828459045                0.0]
            [              0.0                0.0  2.718281828459045]
            >>> expm([[1,1,0],[1,0,1],[0,1,0]])
            [ 3.868145006154145  2.268128708521446  0.8411308412301961]
            [ 2.268128708521446  2.441147138862895    1.42699786729125]
            [0.8411308412301961   1.42699786729125  1.6000162976326988]
            >>> expm([[1,1,0],[1,0,1],[0,1,0]], method='pade')
            [ 3.868145006154145  2.268128708521446  0.8411308412301961]
            [ 2.268128708521446  2.441147138862895    1.42699786729125]
            [0.8411308412301961   1.42699786729125  1.6000162976326988]
            >>> expm([[1+j, 0], [1+j,1]])
            [(1.4686939399158851+2.2873552871788423j)                0.0]
            [(1.0377673986356823+3.5369431757220027j)  2.718281828459045]

        Matrices with large entries are allowed::

            >>> expm(matrix([[1,2],[2,3]])**25)
            [5.650240640484147e+2050488462815550  9.142281400919324e+2050488462815550]
            [9.142281400919324e+2050488462815550  1.479252204140347e+2050488462815551]

        The identity `\exp(A+B) = \exp(A) \exp(B)` does not hold for
        noncommuting matrices::

            >>> A = hilbert(3)
            >>> B = A + eye(3)
            >>> chop(mnorm(A*B - B*A))
            0.0
            >>> chop(mnorm(expm(A+B) - expm(A)*expm(B)))
            0.0
            >>> B = A + ones(3)
            >>> mnorm(A*B - B*A)
            1.8000000000000003
            >>> mnorm(expm(A+B) - expm(A)*expm(B))
            42.0927851137247

        """
        A = ctx.matrix(A)

        if method == 'pade':
            prec = ctx.prec
            try:
                ctx.prec += 2*A.rows
                res = ctx._exp_pade(A)
            finally:
                ctx.prec = prec
            return res

        prec = ctx.prec
        j = int(max(1, ctx.mag(ctx.mnorm(A,'inf'))))
        j += int(0.5*prec**0.5)
        try:
            ctx.prec += 10 + 2*j
            tol = +ctx.eps
            A = A/2**j
            T = A
            Y = A**0 + A
            k = 2
            while 1:
                T *= A * (1/ctx.mpf(k))
                if ctx.mnorm(T, 'inf') < tol:
                    break
                Y += T
                k += 1
            for k in range(j):
                Y = Y*Y
        finally:
            ctx.prec = prec
        Y *= 1
        return Y

    def cosm(ctx, A):
        r"""
        Gives the cosine of a square matrix `A`, defined in analogy
        with the matrix exponential.

        Examples::

            >>> from mpmath import mp, eye, cosm, hilbert, j, matrix
            >>> mp.pretty = True
            >>> X = eye(3)
            >>> cosm(X)
            [0.5403023058681398                 0.0                 0.0]
            [               0.0  0.5403023058681398                 0.0]
            [               0.0                 0.0  0.5403023058681398]
            >>> X = hilbert(3)
            >>> cosm(X)
            [ 0.42440383456955494  -0.3166434130471669  -0.22147494594929265]
            [ -0.3166434130471669    0.820646708837824   -0.1271836947700393]
            [-0.22147494594929265  -0.1271836947700393     0.909236687217541]
            >>> X = matrix([[1+j,-2],[0,-j]])
            >>> cosm(X)
            [(0.833730025131149-0.9888977057628651j)  (1.07485840848393-0.17192140544212975j)]
            [                                    0.0                  (1.5430806348152437+0j)]
        """
        A = ctx.matrix(A)
        B = 0.5 * (ctx.expm(A*ctx.j) + ctx.expm(A*(-ctx.j)))
        if not sum(A.apply(ctx.im).apply(abs)):
            B = B.apply(ctx.re)
        return B

    def sinm(ctx, A):
        r"""
        Gives the sine of a square matrix `A`, defined in analogy
        with the matrix exponential.

        Examples::

            >>> from mpmath import mp, eye, sinm, hilbert, matrix, j
            >>> mp.pretty = True
            >>> X = eye(3)
            >>> sinm(X)
            [0.8414709848078965                 0.0                 0.0]
            [               0.0  0.8414709848078965                 0.0]
            [               0.0                 0.0  0.8414709848078965]
            >>> X = hilbert(3)
            >>> sinm(X)
            [ 0.7116085121509942   0.3397839132474387  0.22074283731474142]
            [ 0.3397839132474387  0.24411386569553178  0.18723127117437235]
            [0.22074283731474142  0.18723127117437235  0.15581673076963523]
            >>> X = matrix([[1+j,-2],[0,-j]])
            >>> sinm(X)
            [(1.2984575814159773+0.6349639147847361j)  (-1.9675151193092209+0.3147000217613668j)]
            [                                     0.0                       -1.1752011936438014j]
        """
        A = ctx.matrix(A)
        B = (-0.5j) * (ctx.expm(A*ctx.j) - ctx.expm(A*(-ctx.j)))
        if not sum(A.apply(ctx.im).apply(abs)):
            B = B.apply(ctx.re)
        return B

    def _sqrtm_rot(ctx, A, _may_rotate):
        # If the iteration fails to converge, cheat by performing
        # a rotation by a complex number
        u = ctx.j**0.3
        return ctx.sqrtm(u*A, _may_rotate) / ctx.sqrt(u)

    def sqrtm(ctx, A, _may_rotate=2):
        r"""
        Computes a square root of the square matrix `A`, i.e. returns
        a matrix `B = A^{1/2}` such that `B^2 = A`. The square root
        of a matrix, if it exists, is not unique.

        **Examples**

        Square roots of some simple matrices::

            >>> from mpmath import mp, sqrtm, j, matrix, cos, sin, chop, mnorm
            >>> mp.pretty = True
            >>> sqrtm([[1,0], [0,1]])
            [1.0  0.0]
            [0.0  1.0]
            >>> sqrtm([[0,0], [0,0]])
            [0.0  0.0]
            [0.0  0.0]
            >>> sqrtm([[2,0],[0,1]])
            [1.4142135623730951  0.0]
            [               0.0  1.0]
            >>> sqrtm([[1,1],[1,0]])
            [(0.920442065259926-0.21728689675164015j)  (0.5688644810057831+0.3515775842541429j)]
            [(0.5688644810057831+0.3515775842541429j)  (0.3515775842541429-0.5688644810057831j)]
            >>> sqrtm([[1,0],[0,1]])
            [1.0  0.0]
            [0.0  1.0]
            >>> sqrtm([[-1,0],[0,1]])
            [-1j     0.0]
            [0.0  (1+0j)]
            >>> sqrtm([[j,0],[0,j]])
            [(0.7071067811865475+0.7071067811865475j)                                       0.0]
            [                                     0.0  (0.7071067811865475+0.7071067811865475j)]

        A square root of a rotation matrix, giving the corresponding
        half-angle rotation matrix::

            >>> t1 = 0.75
            >>> t2 = t1 * 0.5
            >>> A1 = matrix([[cos(t1), -sin(t1)], [sin(t1), cos(t1)]])
            >>> A2 = matrix([[cos(t2), -sin(t2)], [sin(t2), cos(t2)]])
            >>> sqrtm(A1)
            [0.9305076219123143  -0.3662725290860475]
            [0.3662725290860475   0.9305076219123143]
            >>> A2
            [ 0.9305076219123143  -0.36627252908604757]
            [0.36627252908604757    0.9305076219123143]

        The identity `(A^2)^{1/2} = A` does not necessarily hold::

            >>> A = matrix([[4,1,4],[7,8,9],[10,2,11]])
            >>> sqrtm(A**2)
            [ 4.0  1.0   4.0]
            [ 7.0  8.0   9.0]
            [10.0  2.0  11.0]
            >>> sqrtm(A)**2
            [3.9999999999999996  1.0  3.9999999999999996]
            [               7.0  8.0                 9.0]
            [ 9.999999999999998  2.0  10.999999999999998]
            >>> A = matrix([[-4,1,4],[7,-8,9],[10,2,11]])
            >>> sqrtm(A**2)
            [  7.4371511219499515  -0.32412756998547376  1.8481718827526006]
            [-0.25154971571694174     9.326997659004016  2.4822118098514743]
            [   4.116093888336158    0.7757518770982581   13.01795569734202]
            >>> chop(sqrtm(A)**2)
            [-4.0   1.0   4.0]
            [ 7.0  -8.0   9.0]
            [10.0   2.0  11.0]

        For some matrices, a square root does not exist::

            >>> sqrtm([[0,1], [0,0]])
            Traceback (most recent call last):
              ...
            ZeroDivisionError: matrix is numerically singular

        Two examples from the documentation for Matlab's ``sqrtm``::

            >>> mp.pretty = True
            >>> sqrtm([[7,10],[15,22]])
            [1.5666989036012806  1.7407765595569784]
            [2.6111648393354674   4.177863742936748]
            >>>
            >>> X = matrix(\
            ...   [[5,-4,1,0,0],
            ...   [-4,6,-4,1,0],
            ...   [1,-4,6,-4,1],
            ...   [0,1,-4,6,-4],
            ...   [0,0,1,-4,5]])
            >>> Y = matrix(\
            ...   [[2,-1,-0,-0,-0],
            ...   [-1,2,-1,0,-0],
            ...   [0,-1,2,-1,0],
            ...   [-0,0,-1,2,-1],
            ...   [-0,-0,-0,-1,2]])
            >>> mnorm(sqrtm(X) - Y)
            4.531553283261138e-19

        """
        A = ctx.matrix(A)
        # Trivial
        if A*0 == A:
            return A
        prec = ctx.prec
        if _may_rotate:
            d = ctx.det(A)
            if abs(ctx.im(d)) < 16*ctx.eps and ctx.re(d) < 0:
                return ctx._sqrtm_rot(A, _may_rotate-1)
        try:
            ctx.prec += 10
            tol = ctx.eps * 128
            Y = A
            Z = I = A**0
            k = 0
            # Denman-Beavers iteration
            while 1:
                Yprev = Y
                try:
                    Y, Z = 0.5*(Y+ctx.inverse(Z)), 0.5*(Z+ctx.inverse(Y))
                except ZeroDivisionError:
                    if _may_rotate:
                        Y = ctx._sqrtm_rot(A, _may_rotate-1)
                        break
                    else:
                        raise
                mag1 = ctx.mnorm(Y-Yprev, 'inf')
                mag2 = ctx.mnorm(Y, 'inf')
                if mag1 <= mag2*tol:
                    break
                if _may_rotate and k > 6 and not mag1 < mag2 * 0.001:
                    return ctx._sqrtm_rot(A, _may_rotate-1)
                k += 1
                if k > ctx.prec:
                    raise ctx.NoConvergence
        finally:
            ctx.prec = prec
        Y *= 1
        return Y

    def logm(ctx, A):
        r"""
        Computes a logarithm of the square matrix `A`, i.e. returns
        a matrix `B = \log(A)` such that `\exp(B) = A`. The logarithm
        of a matrix, if it exists, is not unique.

        **Examples**

        Logarithms of some simple matrices::

            >>> from mpmath import (mp, eye, logm, expm, matrix, j, nprint,
            ...                     chop, hilbert, cos, sin, pi, re)
            >>> mp.pretty = True
            >>> X = eye(3)
            >>> logm(X)
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            [0.0  0.0  0.0]
            >>> logm(2*X)
            [0.6931471805599453                 0.0                 0.0]
            [               0.0  0.6931471805599453                 0.0]
            [               0.0                 0.0  0.6931471805599453]
            >>> logm(expm(X))
            [1.0  0.0  0.0]
            [0.0  1.0  0.0]
            [0.0  0.0  1.0]

        A logarithm of a complex matrix::

            >>> X = matrix([[2+j, 1, 3], [1-j, 1-2*j, 1], [-4, -5, j]])
            >>> B = logm(X)
            >>> nprint(B)
            [ (0.808757 + 0.107759j)    (2.20752 + 0.202762j)   (1.07376 - 0.773874j)]
            [ (0.905709 - 0.107795j)  (0.0287395 - 0.824993j)  (0.111619 + 0.514272j)]
            [(-0.930151 + 0.399512j)   (-2.06266 - 0.674397j)  (0.791552 + 0.519839j)]
            >>> chop(expm(B))
            [(2+0.9999999999999999j)                      1.0  2.9999999999999996]
            [(0.9999999999999999-1j)  (0.9999999999999999-2j)  0.9999999999999999]
            [                   -4.0                     -5.0                  1j]

        A matrix `X` close to the identity matrix, for which
        `\log(\exp(X)) = \exp(\log(X)) = X` holds::

            >>> X = eye(3) + hilbert(3)/4
            >>> X
            [               1.25               0.125  0.08333333333333333]
            [              0.125  1.0833333333333333               0.0625]
            [0.08333333333333333              0.0625                 1.05]
            >>> logm(expm(X))
            [               1.25               0.125  0.08333333333333333]
            [              0.125  1.0833333333333333               0.0625]
            [0.08333333333333333              0.0625                 1.05]
            >>> expm(logm(X))
            [               1.25               0.125  0.08333333333333334]
            [              0.125  1.0833333333333333               0.0625]
            [0.08333333333333334              0.0625                 1.05]

        A logarithm of a rotation matrix, giving back the angle of
        the rotation::

            >>> t = 3.7
            >>> A = matrix([[cos(t),sin(t)],[-sin(t),cos(t)]])
            >>> chop(logm(A))
            [              0.0  -2.583185307179586]
            [2.583185307179586                 0.0]
            >>> (2*pi-t)
            2.583185307179586

        For some matrices, a logarithm does not exist::

            >>> logm([[1,0], [0,0]])
            Traceback (most recent call last):
              ...
            ZeroDivisionError: matrix is numerically singular

        Logarithm of a matrix with large entries::

            >>> logm(hilbert(3) * 10**20).apply(re)
            [  45.55975135934325  1.277210060427989  0.31766268771797823]
            [  1.277210060427989  42.52227789735423    2.240037087916037]
            [0.31766268771797823  2.240037087916037    42.39521282226704]

        """
        A = ctx.matrix(A)
        if ctx.mnorm(A, 'inf') == 0:
            raise ValueError("The logarithm is undefined for the zero matrix.")
        prec = ctx.prec
        try:
            ctx.prec += 10
            tol = ctx.eps * 128
            I = A**0
            B = A
            n = 0
            while 1:
                B = ctx.sqrtm(B)
                n += 1
                if ctx.mnorm(B-I, 'inf') < 0.125:
                    break
            T = X = B-I
            L = X*0
            k = 1
            while 1:
                if k & 1:
                    L += T / k
                else:
                    L -= T / k
                T *= X
                if ctx.mnorm(T, 'inf') < tol:
                    break
                k += 1
                if k > ctx.prec:
                    raise ctx.NoConvergence
        finally:
            ctx.prec = prec
        L *= 2**n
        return L

    def powm(ctx, A, r):
        r"""
        Computes `A^r = \exp(A \log r)` for a matrix `A` and complex
        number `r`.

        **Examples**

        Powers and inverse powers of a matrix::

            >>> from mpmath import (mp, matrix, powm, chop, extraprec, fib,
            ...                     phi, sqrt)
            >>> mp.pretty = True
            >>> A = matrix([[4,1,4],[7,8,9],[10,2,11]])
            >>> powm(A, 2)
            [ 63.0  20.0   69.0]
            [174.0  89.0  199.0]
            [164.0  48.0  179.0]
            >>> chop(powm(powm(A, 4), 1/4.))
            [ 4.0  1.0   4.0]
            [ 7.0  8.0   9.0]
            [10.0  2.0  11.0]
            >>> powm(extraprec(20)(powm)(A, -4), -1/4.)
            [ 4.0  1.0   4.0]
            [ 7.0  8.0   9.0]
            [10.0  2.0  11.0]
            >>> chop(powm(powm(A, 1+0.5j), 1/(1+0.5j)))
            [4.000000000000001  1.0000000000000002   4.000000000000001]
            [7.000000000000001                 8.0   9.000000000000002]
            [             10.0  2.0000000000000004  11.000000000000002]
            >>> powm(extraprec(5)(powm)(A, -1.5), -1/(1.5))
            [ 4.0                1.0   4.0]
            [ 7.0  7.999999999999999   9.0]
            [10.0                2.0  11.0]

        A Fibonacci-generating matrix::

            >>> powm([[1,1],[1,0]], 10)
            [89.0  55.0]
            [55.0  34.0]
            >>> fib(10)
            55.0
            >>> powm([[1,1],[1,0]], 6.5)
            [(16.516662696425303-0.012108983738178898j)  (10.207858927108324+0.01959274725759321j)]
            [ (10.207858927108324+0.01959274725759321j)   (6.308803769316979-0.03170173099577211j)]
            >>> (phi**6.5 - (1-phi)**6.5)/sqrt(5)
            (10.207858927108326-0.019592747257593225j)
            >>> powm([[1,1],[1,0]], 6.2)
            [(14.307695300266648-0.008222855781077001j)  (8.817334648375935+0.01330486013837118j)]
            [  (8.817334648375935+0.01330486013837118j)  (5.490360651890715-0.02152771591944818j)]
            >>> (phi**6.2 - (1-phi)**6.2)/sqrt(5)
            (8.817334648375935-0.013304860138371179j)

        """
        A = ctx.matrix(A)
        r = ctx.convert(r)
        prec = ctx.prec
        try:
            ctx.prec += 10
            if ctx.isint(r):
                v = A ** int(r)
            elif ctx.isint(r*2):
                y = int(r*2)
                v = ctx.sqrtm(A) ** y
            else:
                v = ctx.expm(r*ctx.logm(A))
        finally:
            ctx.prec = prec
        v *= 1
        return v

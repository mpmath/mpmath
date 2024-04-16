import warnings

# TODO: interpret list as vectors (for multiplication)

# pickling helper
def _make_matrix(x):
    from mpmath import mp
    return mp.matrix(x)

class _matrix:
    """
    Numerical matrix.

    Specify the dimensions or the data as a nested list.
    Elements default to zero.
    Use a flat list to create a column vector easily.

    The datatype of the context (mpf for mp, mpi for iv, and float for fp) is used to store the data.

    Creating matrices
    -----------------

    Matrices in mpmath are implemented using dictionaries. Only non-zero values
    are stored, so it is cheap to represent sparse matrices.

    The most basic way to create one is to use the ``matrix`` class directly.
    You can create an empty matrix specifying the dimensions:

        >>> from mpmath import (mp, matrix, randmatrix, nprint, ones, norm,
        ...                     mnorm, inf)
        >>> mp.dps = 15
        >>> matrix(2)
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0']])
        >>> matrix(2, 3)
        matrix(
        [['0.0', '0.0', '0.0'],
         ['0.0', '0.0', '0.0']])

    Calling ``matrix`` with one dimension will create a square matrix.

    To access the dimensions of a matrix, use the ``rows`` or ``cols`` keyword:

        >>> A = matrix(3, 2)
        >>> A
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0'],
         ['0.0', '0.0']])
        >>> A.rows
        3
        >>> A.cols
        2

    You can also change the dimension of an existing matrix. This will set the
    new elements to 0. If the new dimension is smaller than before, the
    concerning elements are discarded:

        >>> A.rows = 2
        >>> A
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0']])

    Internally ``mpmathify`` is used every time an element is set. This
    is done using the syntax A[row,column], counting from 0:

        >>> A = matrix(2)
        >>> A[1,1] = 1 + 1j
        >>> A
        matrix(
        [['0.0', '0.0'],
         ['0.0', mpc(real='1.0', imag='1.0')]])

    A more comfortable way to create a matrix lets you use nested lists:

        >>> matrix([[1, 2], [3, 4]])
        matrix(
        [['1.0', '2.0'],
         ['3.0', '4.0']])

    Convenient advanced functions are available for creating various standard
    matrices, see ``zeros``, ``ones``, ``diag``, ``eye``, ``randmatrix`` and
    ``hilbert``.

    Vectors
    .......

    Vectors may also be represented by the ``matrix`` class (with rows = 1 or cols = 1).
    For vectors there are some things which make life easier. A column vector can
    be created using a flat list, a row vectors using an almost flat nested list::

        >>> matrix([1, 2, 3])
        matrix(
        [['1.0'],
         ['2.0'],
         ['3.0']])
        >>> matrix([[1, 2, 3]])
        matrix(
        [['1.0', '2.0', '3.0']])

    Optionally vectors can be accessed like lists, using only a single index::

        >>> x = matrix([1, 2, 3])
        >>> x[1]
        mpf('2.0')
        >>> x[1,0]
        mpf('2.0')

    Other
    .....

    Like you probably expected, matrices can be printed::

        >>> print randmatrix(3) # doctest:+SKIP
        [ 0.782963853573023  0.802057689719883  0.427895717335467]
        [0.0541876859348597  0.708243266653103  0.615134039977379]
        [ 0.856151514955773  0.544759264818486  0.686210904770947]

    Use ``nstr`` or ``nprint`` to specify the number of digits to print::

        >>> nprint(randmatrix(5), 3) # doctest:+SKIP
        [2.07e-1  1.66e-1  5.06e-1  1.89e-1  8.29e-1]
        [6.62e-1  6.55e-1  4.47e-1  4.82e-1  2.06e-2]
        [4.33e-1  7.75e-1  6.93e-2  2.86e-1  5.71e-1]
        [1.01e-1  2.53e-1  6.13e-1  3.32e-1  2.59e-1]
        [1.56e-1  7.27e-2  6.05e-1  6.67e-2  2.79e-1]

    As matrices are mutable, you will need to copy them sometimes::

        >>> A = matrix(2)
        >>> A
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0']])
        >>> B = A.copy()
        >>> B[0,0] = 1
        >>> B
        matrix(
        [['1.0', '0.0'],
         ['0.0', '0.0']])
        >>> A
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0']])

    Finally, it is possible to convert a matrix to a nested list. This is very useful,
    as most Python libraries involving matrices or arrays (namely NumPy or SymPy)
    support this format::

        >>> B.tolist()
        [[mpf('1.0'), mpf('0.0')], [mpf('0.0'), mpf('0.0')]]


    Matrix operations
    -----------------

    You can add and subtract matrices of compatible dimensions::

        >>> A = matrix([[1, 2], [3, 4]])
        >>> B = matrix([[-2, 4], [5, 9]])
        >>> A + B
        matrix(
        [['-1.0', '6.0'],
         ['8.0', '13.0']])
        >>> A - B
        matrix(
        [['3.0', '-2.0'],
         ['-2.0', '-5.0']])
        >>> A + ones(3)
        Traceback (most recent call last):
          ...
        ValueError: incompatible dimensions for addition

    It is possible to multiply or add matrices and scalars. In the latter case the
    operation will be done element-wise::

        >>> A * 2
        matrix(
        [['2.0', '4.0'],
         ['6.0', '8.0']])
        >>> A / 4
        matrix(
        [['0.25', '0.5'],
         ['0.75', '1.0']])
        >>> A - 1
        matrix(
        [['0.0', '1.0'],
         ['2.0', '3.0']])

    Of course you can perform matrix multiplication, if the dimensions are
    compatible, using ``@`` or ``*``. For clarity, ``@`` is
    recommended (`PEP 465 <https://www.python.org/dev/peps/pep-0465/>`), because
    the meaning of ``*`` is different in many other Python libraries such as NumPy.

        >>> A @ B
        matrix(
        [['8.0', '22.0'],
         ['14.0', '48.0']])
        >>> A * B # same as A @ B
        matrix(
        [['8.0', '22.0'],
         ['14.0', '48.0']])
        >>> matrix([[1, 2, 3]]) * matrix([[-6], [7], [-2]])
        matrix(
        [['2.0']])

    You can raise powers of square matrices::

        >>> A**2
        matrix(
        [['7.0', '10.0'],
         ['15.0', '22.0']])

    Negative powers will calculate the inverse::

        >>> A**-1
        matrix(
        [['-2.0', '1.0'],
         ['1.5', '-0.5']])
        >>> A * A**-1
        matrix(
        [['1.0', '1.0842021724855e-19'],
         ['-2.16840434497101e-19', '1.0']])



    Matrix transposition is straightforward::

        >>> A = ones(2, 3)
        >>> A
        matrix(
        [['1.0', '1.0', '1.0'],
         ['1.0', '1.0', '1.0']])
        >>> A.T
        matrix(
        [['1.0', '1.0'],
         ['1.0', '1.0'],
         ['1.0', '1.0']])

    Norms
    .....

    Sometimes you need to know how "large" a matrix or vector is. Due to their
    multidimensional nature it's not possible to compare them, but there are
    several functions to map a matrix or a vector to a positive real number, the
    so called norms.

    For vectors the p-norm is intended, usually the 1-, the 2- and the oo-norm are
    used.

        >>> x = matrix([-10, 2, 100])
        >>> norm(x, 1)
        mpf('112.0')
        >>> norm(x, 2)
        mpf('100.5186549850325')
        >>> norm(x, inf)
        mpf('100.0')

    Please note that the 2-norm is the most used one, though it is more expensive
    to calculate than the 1- or oo-norm.

    It is possible to generalize some vector norms to matrix norm::

        >>> A = matrix([[1, -1000], [100, 50]])
        >>> mnorm(A, 1)
        mpf('1050.0')
        >>> mnorm(A, inf)
        mpf('1001.0')
        >>> mnorm(A, 'F')
        mpf('1006.2310867787777')

    The last norm (the "Frobenius-norm") is an approximation for the 2-norm, which
    is hard to calculate and not available. The Frobenius-norm lacks some
    mathematical properties you might expect from a norm.
    """

    def __init__(self, *args, **kwargs):
        self._data = {}
        # LU decompostion cache, this is useful when solving the same system
        # multiple times, when calculating the inverse and when calculating the
        # determinant
        self._LU = None
        self.shape = None
        if "force_type" in kwargs:
            warnings.warn("The force_type argument was removed, it did not work"
                " properly anyway. If you want to force floating-point or"
                " interval computations, use the respective methods from `fp`"
                " or `mp` instead, e.g., `fp.matrix()` or `iv.matrix()`."
                " If you want to truncate values to integer, use .apply(int) instead.",
                DeprecationWarning)
        if isinstance(args[0], (list, tuple)):
            if not args[0]:
                self.shape = (0, 0)
            elif isinstance(args[0][0], (list, tuple)):
                # interpret nested list as matrix
                A = args[0]
                self.shape = (len(A), len(A[0]))
                for i, row in enumerate(A):
                    for j, a in enumerate(row):
                        # note: this will call __setitem__ which will call self.ctx.convert() to convert the datatype.
                        self[i, j] = a
            else:
                # interpret list as row vector
                v = args[0]
                self.shape = (len(v), 1)
                for i, e in enumerate(v):
                    self[i, 0] = e
        elif isinstance(args[0], int):
            # create empty matrix of given dimensions
            if len(args) == 1:
                self.shape = (args[0],) * 2
            else:
                if not isinstance(args[1], int):
                    raise TypeError("expected int")
                assert len(args) == 2
                self.shape = tuple(args)
        elif isinstance(args[0], _matrix):
            A = args[0]
            rows, cols = self.shape = A.shape
            for i in range(rows):
                for j in range(cols):
                    self[i, j] = A[i, j]
        elif hasattr(args[0], 'tolist'): # A is likely a numpy array
            A = args[0]
            if A.ndim == 1: # convert vector to (N, 1) matrix
                self.shape = (A.shape[0], 1)
            else:
                self.shape = A.shape
            A = self.ctx.matrix(A.tolist())
            self._data = A._data
        else:
            raise TypeError('could not interpret given arguments')

    def apply(self, f):
        """
        Return a copy of self with the function `f` applied elementwise.
        """
        rows, cols = self.shape
        new = self.ctx.matrix(rows, cols)
        for i in range(rows):
            for j in range(cols):
                new[i,j] = f(self[i,j])
        return new

    def __nstr__(self, n=None, **kwargs):
        # Build table of string representations of the elements
        res = []
        # Track per-column max lengths for pretty alignment
        maxlen = [0] * self.cols
        for i in range(self.rows):
            res.append([])
            for j in range(self.cols):
                if n:
                    string = self.ctx.nstr(self[i,j], n, **kwargs)
                else:
                    string = str(self[i,j])
                res[-1].append(string)
                maxlen[j] = max(len(string), maxlen[j])
        # Patch strings together
        rowsep = '\n'
        colsep = '  '
        for i, row in enumerate(res):
            for j, elem in enumerate(row):
                # Pad each element up to maxlen so the columns line up
                row[j] = elem.rjust(maxlen[j])
            res[i] = "[" + colsep.join(row) + "]"
        return rowsep.join(res) if self.rows or self.cols else ''

    def __str__(self):
        return self.__nstr__()

    def _toliststr(self):
        """
        Create a list string from a matrix.

        If avoid_type: avoid multiple 'mpf's.
        """
        rows, cols = self.shape
        # XXX: should be something like self.ctx._types
        typ = self.ctx.mpf
        s = '['
        for i in range(rows):
            s += '['
            for j in range(cols):
                if not isinstance(self[i,j], typ):
                    a = repr(self[i,j])
                else:
                    a = "'" + str(self[i,j]) + "'"
                s += a + ', '
            if s[-1] != '[':
                s = s[:-2]
            s += '],\n '
        if s[-1] != '[':
            s = s[:-3]
        s += ']'
        return s

    def tolist(self):
        """
        Convert the matrix to a nested list.
        """
        rows, cols = self.shape
        return [[self[i,j] for j in range(cols)] for i in range(rows)]

    def to_numpy(self):
        """
        Convert the matrix to a numpy object
        """
        import numpy as np
        return np.array(self).reshape(self.shape)

    def __repr__(self):
        if self.ctx.pretty:
            return self.__str__()
        s = 'matrix(\n'
        s += self._toliststr() + ')'
        return s

    def _get_element(self, key):
        '''
        Fast extraction of the i,j element from the matrix
            This function is for private use only because is unsafe:
                1. Does not check on the value of key it expects key to be a integer tuple (i,j)
                2. Does not check bounds
        '''
        if key in self._data:
            return self._data[key]
        else:
            return self.ctx.zero

    def _set_element(self, key, value):
        '''
        Fast assignment of the i,j element in the matrix
            This function is unsafe:
                1. Does not check on the value of key it expects key to be a integer tuple (i,j)
                2. Does not check bounds
                3. Does not check the value type
                4. Does not reset the LU cache
        '''
        if value: # only store non-zeros
            self._data[key] = value
        elif key in self._data:
            del self._data[key]


    def __getitem__(self, key):
        '''
            Getitem function for mp matrix class with slice index enabled
            it allows the following assingments
            scalar to a slice of the matrix
         B = A[:,2:6]
        '''
        n_rows, n_cols = self.shape
        # Convert vector to matrix indexing
        if isinstance(key, int) or isinstance(key,slice):
            # only sufficent for vectors
            if n_rows == 1:
                key = (0, key)
            elif n_cols == 1:
                key = (key, 0)
            else:
                raise IndexError('insufficient indices for matrix')

        row, col = key
        if isinstance(row, slice) or isinstance(col, slice):

            #Rows
            if isinstance(row, slice):
                #Check bounds
                if (row.start is None or row.start >= 0) and \
                    (row.stop is None or row.stop <= n_rows+1):
                    # Generate indices
                    rows = range(*row.indices(n_rows))
                else:
                    raise IndexError('Row index out of bounds')
            else:
                # Single row
                if row >= n_rows:
                    raise IndexError('Row index out of bounds')
                rows = [row]

            # Columns
            if isinstance(col, slice):
                # Check bounds
                if (col.start is None or col.start >= 0) and \
                    (col.stop is None or col.stop <= n_cols+1):
                    # Generate indices
                    columns = range(*col.indices(n_cols))
                else:
                    raise IndexError('Column index out of bounds')

            else:
                # Single column
                if col >= n_cols:
                    raise IndexError('Column index out of bounds')
                columns = [col]

            # Create matrix slice
            m = self.ctx.matrix(len(rows),len(columns))

            # Assign elements to the output matrix
            for i,x in enumerate(rows):
                for j,y in enumerate(columns):
                    m._set_element((i,j),self._get_element((x,y)))

            return m

        else:
            # single element extraction
            if row >= n_rows or col >= n_cols:
                raise IndexError('matrix index out of range')
            if key in self._data:
                return self._data[key]
            else:
                return self.ctx.zero

    def __setitem__(self, key, value):
        # setitem function for mp matrix class with slice index enabled
        # it allows the following assingments
        #  scalar to a slice of the matrix
        # A[:,2:6] = 2.5
        #  submatrix to matrix (the value matrix should be the same size as the slice size)
        # A[3,:] = B   where A is n x m  and B is n x 1

        n_rows, n_cols = self.shape
        # Convert vector to matrix indexing
        if isinstance(key, int) or isinstance(key,slice):
            # only sufficent for vectors
            if n_rows == 1:
                key = (0, key)
            elif n_cols == 1:
                key = (key, 0)
            else:
                raise IndexError('insufficient indices for matrix')

        row, col = key
        # Slice indexing
        if isinstance(row, slice) or isinstance(col, slice):
            # Rows
            if isinstance(row, slice):
                # Check bounds
                if (row.start is None or row.start >= 0) and \
                    (row.stop is None or row.stop <= n_rows+1):
                    # generate row indices
                    rows = range(*row.indices(n_rows))
                else:
                    raise IndexError('Row index out of bounds')
            else:
                # Single row
                rows = [row]
            # Columns
            if isinstance(col, slice):
                # Check bounds
                if (col.start is None or col.start >= 0) and \
                    (col.stop is None or col.stop <= n_cols+1):
                    # Generate column indices
                    columns = range(*col.indices(n_cols))
                else:
                    raise IndexError('Column index out of bounds')
            else:
                # Single column
                columns = [col]
            # Assign slice with a scalar
            if isinstance(value,self.ctx.matrix):
                # Assign elements to matrix if input and output dimensions match
                if len(rows) == value.rows and len(columns) == value.cols:
                    for i,x in enumerate(rows):
                        for j,y in enumerate(columns):
                            self._set_element((x,y), value._get_element((i,j)))
                else:
                    raise ValueError('Dimensions do not match')
            elif _is_ndarray(value):
                assert value.ndim == 2
                if len(rows) == value.shape[0] and len(columns) == value.shape[1]:
                    for i,x in enumerate(rows):
                        for j,y in enumerate(columns):
                            self.__set_element((x,y), value[i,j])
                else:
                    raise ValueError('Dimensions do not match')
            else:
                # Assign slice with scalars
                value = self.ctx.convert(value)
                for i in rows:
                    for j in columns:
                        self._set_element((i,j), value)
        else:
            # Single element assingment
            # Check bounds
            if row >= n_rows or col >= n_cols:
                raise IndexError('matrix index out of range')
            # Convert and store value
            value = self.ctx.convert(value)
            if value: # only store non-zeros
                self._data[key] = value
            elif key in self._data:
                del self._data[key]

        self._LU = None
        return

    def __iter__(self):
        rows, cols = self.shape
        for i in range(rows):
            for j in range(cols):
                yield self[i,j]

    def __mul__(self, other):
        rows, cols = self.shape
        if isinstance(other, self.ctx.matrix):
            # dot multiplication
            other_rows, other_cols = other.shape
            if cols != other_rows:
                raise ValueError('dimensions not compatible for multiplication')
            new = self.ctx.matrix(rows, other_cols)
            for i in range(rows):
                for j in range(other_cols):
                    new[i, j] = self.ctx.fdot((self._data[i,k], other._data[k,j])
                                              for k in range(other_rows) if (i,k) in self._data and (k,j) in other._data)
        elif _is_ndarray(other):
            import numpy as np
            new = self.ctx.matrix(np.array(self).dot(other))
        else:
            # try scalar multiplication
            new = self.ctx.matrix(rows, cols)
            for i in range(rows):
                for j in range(cols):
                    new[i, j] = other * self[i, j]
        return new

    def __matmul__(self, other):
        return self.__mul__(other)

    def __rmul__(self, other):
        # assume other is scalar and thus commutative
        if isinstance(other, self.ctx.matrix) or _is_ndarray(other):
            raise TypeError(f"other should not be type of {other.__class__}")
        return self.__mul__(other)

    def __pow__(self, other):
        # avoid cyclic import problems
        #from linalg import inverse
        if not isinstance(other, int):
            raise ValueError('only integer exponents are supported')
        rows, cols = self.shape
        if not rows == cols:
            raise ValueError('only powers of square matrices are defined')
        n = other
        if n == 0:
            return self.ctx.eye(rows)
        if n < 0:
            n = -n
            neg = True
        else:
            neg = False
        i = n
        y = 1
        z = self.copy()
        while i != 0:
            if i % 2 == 1:
                y = y * z
            z = z*z
            i = i // 2
        if neg:
            y = self.ctx.inverse(y)
        return y

    def __truediv__(self, other):
        # assume other is scalar and do element-wise divison
        assert not isinstance(other, self.ctx.matrix)
        assert not _is_ndarray(other)
        rows, cols = self.shape
        new = self.ctx.matrix(rows, cols)
        for i in range(rows):
            for j in range(cols):
                new[i,j] = self[i,j] / other
        return new

    def __add__(self, other):
        rows, cols = self.shape
        if isinstance(other, self.ctx.matrix) or _is_ndarray(other):
            if self.shape != other.shape:
                raise ValueError('incompatible dimensions for addition')
            new = self.ctx.matrix(rows, cols)
            for i in range(rows):
                for j in range(cols):
                    new[i,j] = self[i,j] + other[i,j]
            return new
        else:
            # assume other is scalar and add element-wise
            new = self.ctx.matrix(rows, cols)
            for i in range(rows):
                for j in range(cols):
                    new[i,j] += self[i,j] + other
            return new

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if ((isinstance(other, self.ctx.matrix) or _is_ndarray(other)) and
            self.shape != other.shape):
            raise ValueError('incompatible dimensions for subtraction')
        return self.__add__(other * (-1))

    def __pos__(self):
        """
        +M returns a copy of M, rounded to current working precision.
        """
        return (+1) * self

    def __neg__(self):
        return (-1) * self

    def __rsub__(self, other):
        return -self + other

    def __eq__(self, other):
        try:
            return self.shape == other.shape and self._data == other._data
        except AttributeError:
            return NotImplemented

    def __len__(self):
        rows, cols = self.shape
        if rows == 1:
            return cols
        elif cols == 1:
            return rows
        else:
            return rows # do it like numpy

    @property
    def rows(self):
        """Number of rows."""
        return self.shape[0]

    @rows.setter
    def rows(self, value):
        for key in self._data.copy():
            if key[0] >= value:
                del self._data[key]
        self.shape = (value, self.shape[1])

    @property
    def cols(self):
        """Number of columns."""
        return self.shape[1]

    @cols.setter
    def cols(self, value):
        for key in self._data.copy():
            if key[1] >= value:
                del self._data[key]
        self.shape = (self.shape[0], value)

    def transpose(self):
        rows, cols = self.shape
        new = self.ctx.matrix(cols, rows)
        for i in range(rows):
            for j in range(cols):
                new[j,i] = self[i,j]
        return new

    T = property(transpose)

    def conjugate(self):
        return self.apply(self.ctx.conj)

    def transpose_conj(self):
        return self.conjugate().transpose()

    H = property(transpose_conj)

    def copy(self):
        new = self.ctx.matrix(*self.shape)
        new._data = self._data.copy()
        new._LU = self._LU
        return new
    __copy__ = copy

    def __reduce__(self):
        return _make_matrix, (self.tolist(),)

    def __array__(self):
        from numpy import empty
        r = empty((self.rows, self.cols), object)
        for i in range(self.rows):
            for j in range(self.cols):
                r[i, j] = self[i, j]
        return r

    def column(self, n):
        rows = self.rows
        m = self.ctx.matrix(rows, 1)
        for i in range(rows):
            m[i] = self[i,n]
        return m

class MatrixMethods:

    def __init__(ctx):
        # XXX: subclass
        ctx.matrix = type('matrix', (_matrix,), {})
        ctx.matrix.ctx = ctx
        ctx.matrix.convert = ctx.convert

    def eye(ctx, n, **kwargs):
        """
        Create square identity matrix n x n.
        """
        A = ctx.matrix(n, **kwargs)
        for i in range(n):
            A[i,i] = 1
        return A

    def diag(ctx, diagonal, **kwargs):
        """
        Create square diagonal matrix using given list.

        Example:
        >>> from mpmath import diag, mp
        >>> diag([1, 2, 3])
        matrix(
        [['1.0', '0.0', '0.0'],
         ['0.0', '2.0', '0.0'],
         ['0.0', '0.0', '3.0']])
        """
        A = ctx.matrix(len(diagonal), **kwargs)
        for i in range(len(diagonal)):
            A[i,i] = diagonal[i]
        return A

    def zeros(ctx, *args, **kwargs):
        """
        Create matrix m x n filled with zeros.
        One given dimension will create square matrix n x n.

        Example:
        >>> from mpmath import zeros, mp
        >>> zeros(2)
        matrix(
        [['0.0', '0.0'],
         ['0.0', '0.0']])
        """
        if len(args) == 1:
            m = n = args[0]
        elif len(args) == 2:
            m = args[0]
            n = args[1]
        else:
            raise TypeError('zeros expected at most 2 arguments, got %i' % len(args))
        A = ctx.matrix(m, n, **kwargs)
        for i in range(m):
            for j in range(n):
                A[i,j] = 0
        return A

    def ones(ctx, *args, **kwargs):
        """
        Create matrix m x n filled with ones.
        One given dimension will create square matrix n x n.

        Example:
        >>> from mpmath import ones, mp
        >>> ones(2)
        matrix(
        [['1.0', '1.0'],
         ['1.0', '1.0']])
        """
        if len(args) == 1:
            m = n = args[0]
        elif len(args) == 2:
            m = args[0]
            n = args[1]
        else:
            raise TypeError('ones expected at most 2 arguments, got %i' % len(args))
        A = ctx.matrix(m, n, **kwargs)
        for i in range(m):
            for j in range(n):
                A[i,j] = 1
        return A

    def hilbert(ctx, m, n=None):
        """
        Create (pseudo) hilbert matrix m x n.
        One given dimension will create hilbert matrix n x n.

        The matrix is very ill-conditioned and symmetric, positive definite if
        square.
        """
        if n is None:
            n = m
        A = ctx.matrix(m, n)
        for i in range(m):
            for j in range(n):
                A[i,j] = ctx.one / (i + j + 1)
        return A

    def randmatrix(ctx, m, n=None, min=0, max=1, **kwargs):
        """
        Create a random m x n matrix.

        All values are >= min and <max.
        n defaults to m.

        Example:
        >>> from mpmath import randmatrix
        >>> randmatrix(2) # doctest:+SKIP
        matrix(
        [['0.53491598236191806', '0.57195669543302752'],
         ['0.85589992269513615', '0.82444367501382143']])
        """
        if not n:
            n = m
        A = ctx.matrix(m, n, **kwargs)
        for i in range(m):
            for j in range(n):
                A[i,j] = ctx.rand() * (max - min) + min
        return A

    def swap_row(ctx, A, i, j):
        """
        Swap row i with row j.
        """
        if i == j:
            return
        if isinstance(A, ctx.matrix):
            for k in range(A.cols):
                A[i,k], A[j,k] = A[j,k], A[i,k]
        elif isinstance(A, list):
            A[i], A[j] = A[j], A[i]
        elif _is_ndarray(A):
            A[[j,i]] = A[[i,j]]
        else:
            raise TypeError('could not interpret type')

    def extend(ctx, A, b):
        """
        Extend matrix A with column b and return result.
        """
        if not (isinstance(A, ctx.matrix) or _is_ndarray(A)):
            raise TypeError("A should be a type of ctx.matrix")
        if A.shape[0] != len(b):
            raise ValueError("Value should be equal to len(b)")
        if _is_ndarray(A):
            import numpy as np
            A = np.hstack((A, np.asarray(b)[:,None]))
        else:
            A = A.copy()
            A.cols += 1
            cols = A.shape[1]
            for i in range(A.rows):
                A[i, cols-1] = b[i]
        return A

    def norm(ctx, x, p=2):
        r"""
        Gives the entrywise `p`-norm of an iterable *x*, i.e. the vector norm
        `\left(\sum_k |x_k|^p\right)^{1/p}`, for any given `1 \le p \le \infty`.

        Special cases:

        If *x* is not iterable, this just returns ``absmax(x)``.

        ``p=1`` gives the sum of absolute values.

        ``p=2`` is the standard Euclidean vector norm.

        ``p=inf`` gives the magnitude of the largest element.

        For *x* a matrix, ``p=2`` is the Frobenius norm.
        For operator matrix norms, use :func:`~mpmath.mnorm` instead.

        You can use the string 'inf' as well as float('inf') or mpf('inf')
        to specify the infinity norm.

        **Examples**

            >>> from mpmath import matrix, norm, inf
            >>> x = matrix([-10, 2, 100])
            >>> norm(x, 1)
            mpf('112.0')
            >>> norm(x, 2)
            mpf('100.5186549850325')
            >>> norm(x, inf)
            mpf('100.0')

        """
        try:
            iter(x)
        except TypeError:
            return ctx.absmax(x)
        if _is_ndarray(x):
            x = x.ravel()

        if type(p) is not int:
            p = ctx.convert(p)
        if p == ctx.inf:
            return max(ctx.absmax(i) for i in x)
        elif p == 1:
            return ctx.fsum(x, absolute=1)
        elif p == 2:
            return ctx.sqrt(ctx.fsum(x, absolute=1, squared=1))
        elif p > 1:
            return ctx.nthroot(ctx.fsum(abs(i)**p for i in x), p)
        else:
            raise ValueError('p has to be >= 1')

    def mnorm(ctx, A, p=1):
        r"""
        Gives the matrix (operator) `p`-norm of A. Currently ``p=1`` and ``p=inf``
        are supported:

        ``p=1`` gives the 1-norm (maximal column sum)

        ``p=inf`` gives the `\infty`-norm (maximal row sum).
        You can use the string 'inf' as well as float('inf') or mpf('inf')

        ``p=2`` (not implemented) for a square matrix is the usual spectral
        matrix norm, i.e. the largest singular value.

        ``p='f'`` (or 'F', 'fro', 'Frobenius, 'frobenius') gives the
        Frobenius norm, which is the elementwise 2-norm. The Frobenius norm is an
        approximation of the spectral norm and satisfies

        .. math ::

            \frac{1}{\sqrt{\mathrm{rank}(A)}} \|A\|_F \le \|A\|_2 \le \|A\|_F

        The Frobenius norm lacks some mathematical properties that might
        be expected of a norm.

        For general elementwise `p`-norms, use :func:`~mpmath.norm` instead.

        **Examples**

            >>> from mpmath import matrix, mnorm, inf
            >>> A = matrix([[1, -1000], [100, 50]])
            >>> mnorm(A, 1)
            mpf('1050.0')
            >>> mnorm(A, inf)
            mpf('1001.0')
            >>> mnorm(A, 'F')
            mpf('1006.2310867787777')

        """
        A = ctx.matrix(A)
        if type(p) is not int:
            if type(p) is str and 'frobenius'.startswith(p.lower()):
                return ctx.norm(A, 2)
            p = ctx.convert(p)
        m, n = A.shape
        if p == 1:
            return max((ctx.fsum((A[i,j] for i in range(m)), absolute=1) for j in range(n)), default=0)
        elif p == ctx.inf:
            return max((ctx.fsum((A[i,j] for j in range(n)), absolute=1) for i in range(m)), default=0)
        else:
            raise NotImplementedError("matrix p-norm for arbitrary p")

def _is_ndarray(obj):
    '''Return whether an object is a numpy array.'''
    return hasattr(obj, '__array_interface__')

def _np_zeros(ctx, shape):
    import numpy as np
    return np.full(shape, ctx.mpf('0.0'))

def _np_eye(ctx, n):
    import numpy as np
    a = np.full((n, n), ctx.mpf('0.0'))
    a[np.diag_indices(n)] = ctx.mpf('1.0')
    return a

def _np_conj(ctx, a):
    import numpy as np
    return np.vectorize(ctx.conj)(a)

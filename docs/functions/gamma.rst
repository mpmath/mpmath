Factorials and gamma functions
------------------------------

Factorials and factorial-like sums and products are basic tools of
combinatorics and number theory.  Much like the exponential function is
fundamental to differential equations and analysis in general, the factorial
function (and its extension to complex numbers, the gamma function) is
fundamental to difference equations and functional equations.

A large selection of factorial-like functions is implemented in mpmath.  All
functions support complex arguments, and arguments may be arbitrarily large.
Results are numerical approximations, so to compute *exact* values a high
enough precision must be set manually::

    >>> from mpmath import mp, fac
    >>> mp.dps = 15
    >>> mp.pretty = True
    >>> fac(100)
    9.33262154439442e+157
    >>> print(int(_))    # most digits are wrong
    93326215443944150965646704795953882578400970373184098831012889540582227238570431295066113089288327277825849664006524270554535976289719382852181865895959724032
    >>> mp.dps = 160
    >>> fac(100)
    93326215443944152681699238856266700490715968264381621468592963895217599993229915608941463976156518286253697920827223758251185210916864000000000000000000000000.0

The gamma and polygamma functions are closely related to :doc:`zeta`.  See also
:doc:`qfunctions` for q-analogs of factorial-like functions.


Factorials
..........

.. autofunction:: mpmath.factorial
.. autofunction:: mpmath.fac2


Binomial coefficients
.....................

.. autofunction:: mpmath.binomial


Gamma function
..............

.. autofunction:: mpmath.gamma
.. autofunction:: mpmath.rgamma
.. autofunction:: mpmath.gammaprod
.. autofunction:: mpmath.loggamma


Rising and falling factorials
.............................

.. autofunction:: mpmath.rf
.. autofunction:: mpmath.ff


Beta function
.............

.. autofunction:: mpmath.beta
.. autofunction:: mpmath.betainc


Super- and hyperfactorials
..........................

.. autofunction:: mpmath.superfac
.. autofunction:: mpmath.hyperfac
.. autofunction:: mpmath.barnesg


Polygamma functions and harmonic numbers
........................................

.. autofunction:: mpmath.psi
.. autofunction:: mpmath.digamma
.. autofunction:: mpmath.harmonic

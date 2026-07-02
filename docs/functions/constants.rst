Mathematical constants
----------------------

Mpmath supports arbitrary-precision computation of various common (and less
common) mathematical constants.  These constants are implemented as lazy
objects that can evaluate to any precision.  Whenever the objects are used as
function arguments or as operands in arithmetic operations, they automagically
evaluate to the current working precision.  A lazy number can be converted to a
regular ``mpf`` using the unary ``+`` operator, or by calling it as a
function::

    >>> from mpmath import pi, mp
    >>> pi
    <pi: 3.14159~>
    >>> 2*pi
    mpf('6.2831853071795862')
    >>> +pi
    mpf('3.1415926535897931')
    >>> pi()
    mpf('3.1415926535897931')
    >>> mp.dps = 40
    >>> pi
    <pi: 3.14159~>
    >>> 2*pi
    mpf('6.28318530717958647692528676655900576839434')
    >>> +pi
    mpf('3.14159265358979323846264338327950288419717')
    >>> pi()
    mpf('3.14159265358979323846264338327950288419717')

The predefined objects ``j`` (imaginary unit), ``inf`` (positive infinity) and
``nan`` (not-a-number) are shortcuts to ``mpc`` and ``mpf`` instances with
these fixed values.

.. autofunction:: mpmath.mp.pi
.. autoattribute:: mpmath.mp.degree
.. autoattribute:: mpmath.mp.e
.. autoattribute:: mpmath.mp.phi
.. autofunction:: mpmath.mp.euler
.. autoattribute:: mpmath.mp.catalan
.. autoattribute:: mpmath.mp.apery
.. autoattribute:: mpmath.mp.khinchin
.. autoattribute:: mpmath.mp.glaisher
.. autoattribute:: mpmath.mp.mertens
.. autoattribute:: mpmath.mp.twinprime

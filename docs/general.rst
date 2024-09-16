Utility functions
===============================================

This page lists functions that perform basic operations
on numbers or aid general programming.

Conversion and printing
-----------------------

:func:`~mpmath.mpmathify` / ``convert()``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.mpmathify(x, strings=True)

:func:`~mpmath.nstr`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nstr(x, n=6, **kwargs)

:func:`~mpmath.nprint`
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nprint(x, n=6, **kwargs)

:func:`mpmath.mpf.__format__`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.mpf.__format__(s, format_spec)

:func:`mpmath.mpc.__format__`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.mpc.__format__(s, format_spec)

Arithmetic operations
---------------------

See also :func:`mpmath.sqrt`, :func:`mpmath.exp` etc., listed
in :doc:`functions/powers`

:func:`~mpmath.fadd`
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fadd

:func:`~mpmath.fsub`
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fsub

:func:`~mpmath.fneg`
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fneg

:func:`~mpmath.fmul`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fmul

:func:`~mpmath.fdiv`
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fdiv

:func:`~mpmath.fmod`
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fmod(x, y)

:func:`~mpmath.fsum`
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fsum(terms, absolute=False, squared=False)

:func:`~mpmath.fprod`
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fprod(factors)

:func:`~mpmath.fdot`
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fdot(A, B=None, conjugate=False)

Complex components
------------------

:func:`~mpmath.fabs`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fabs(x)

:func:`~mpmath.sign`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sign(x)

:func:`~mpmath.re`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.re(x)

:func:`~mpmath.im`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.im(x)

:func:`~mpmath.arg`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.arg(x)

:func:`~mpmath.conj`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.conj(x)

:func:`~mpmath.polar`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.polar(x)

:func:`~mpmath.rect`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.rect(x)

Integer and fractional parts
-----------------------------

:func:`~mpmath.floor`
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.floor(x)

:func:`~mpmath.ceil`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ceil(x)

:func:`~mpmath.nint`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nint(x)

:func:`~mpmath.frac`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.frac(x)

Tolerances and approximate comparisons
--------------------------------------

:func:`~mpmath.chop`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.chop(x, tol=None)

:func:`~mpmath.almosteq`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.almosteq(s, t, rel_eps=None, abs_eps=None)

Properties of numbers
-------------------------------------

:func:`~mpmath.isinf`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.isinf(x)

:func:`~mpmath.isnan`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.isnan(x)

:func:`~mpmath.isnormal`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.isnormal(x)

:func:`~mpmath.isfinite`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.isfinite(x)

:func:`~mpmath.isint`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.isint(x, gaussian=False)

:func:`~mpmath.ldexp`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ldexp(x, n)

:func:`~mpmath.frexp`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.frexp(x, n)

:func:`~mpmath.mag`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.mag(x)

:func:`~mpmath.nint_distance`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nint_distance(x)

.. :func:`~mpmath.absmin`
.. ^^^^^^^^^^^^^^^^^^^^^^^^
.. .. autofunction:: mpmath.absmin(x)
.. .. autofunction:: mpmath.absmax(x)

Number generation
-----------------

:func:`~mpmath.fraction`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fraction(p,q)

:func:`~mpmath.rand`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.rand()

:func:`~mpmath.arange`
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.arange(*args)

:func:`~mpmath.linspace`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.linspace(*args, **kwargs)

Precision management
--------------------

:func:`~mpmath.autoprec`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.autoprec

:func:`~mpmath.workprec`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.workprec

:func:`~mpmath.workdps`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.workdps

:func:`~mpmath.extraprec`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.extraprec

:func:`~mpmath.extradps`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.extradps

Performance and debugging
------------------------------------

:func:`~mpmath.memoize`
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.memoize

:func:`~mpmath.maxcalls`
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.maxcalls

:func:`~mpmath.monitor`
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.monitor

:func:`~mpmath.timing`
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.timing

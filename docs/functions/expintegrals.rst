Exponential integrals and error functions
-----------------------------------------

Exponential integrals give closed-form solutions to a large class of commonly occurring transcendental integrals that cannot be evaluated using elementary functions. Integrals of this type include those with an integrand of the form `t^a e^{t}` or `e^{-x^2}`, the latter giving rise to the Gaussian (or normal) probability distribution.

The most general function in this section is the incomplete gamma function, to which all others can be reduced. The incomplete gamma function, in turn, can be expressed using hypergeometric functions (see :doc:`hypergeometric`).

Incomplete gamma functions
..........................

:func:`~mpmath.gammainc`
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.gammainc(z, a=0, b=inf, regularized=False)

:func:`~mpmath.lower_gamma`
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.lower_gamma

:func:`~mpmath.upper_gamma`
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.upper_gamma

Exponential integrals
.....................

:func:`~mpmath.ei`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ei(x, **kwargs)

:func:`~mpmath.e1`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.e1(x, **kwargs)

:func:`~mpmath.expint`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.expint(*args)


Logarithmic integral
....................

:func:`~mpmath.li`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.li(x, **kwargs)


Trigonometric integrals
.......................

:func:`~mpmath.ci`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ci(x, **kwargs)

:func:`~mpmath.si`
^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.si(x, **kwargs)


Hyperbolic integrals
....................

:func:`~mpmath.chi`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.chi(x, **kwargs)

:func:`~mpmath.shi`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.shi(x, **kwargs)


Error functions
...............

:func:`~mpmath.erf`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.erf(x, **kwargs)

:func:`~mpmath.erfc`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.erfc(x, **kwargs)

:func:`~mpmath.erfi`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.erfi(x)

:func:`~mpmath.erfinv`
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.erfinv(x)

The normal distribution
....................................................

:func:`~mpmath.npdf`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.npdf(x, mu=0, sigma=1)

:func:`~mpmath.ncdf`
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ncdf(x, mu=0, sigma=1)


Fresnel integrals
......................................................

:func:`~mpmath.fresnels`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fresnels(x)

:func:`~mpmath.fresnelc`
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.fresnelc(x)

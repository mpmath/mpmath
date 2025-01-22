Exponential integrals and error functions
-----------------------------------------

Exponential integrals give closed-form solutions to a large class of commonly
occurring transcendental integrals that cannot be evaluated using elementary
functions.  Integrals of this type include those with an integrand of the form
`t^a e^{t}` or `e^{-x^2}`, the latter giving rise to the Gaussian (or normal)
probability distribution.

The most general function in this section is the incomplete gamma function, to
which all others can be reduced.  The incomplete gamma function, in turn, can
be expressed using hypergeometric functions (see :doc:`hypergeometric`).

Incomplete gamma functions
..........................

.. autofunction:: mpmath.gammainc
.. autofunction:: mpmath.lower_gamma
.. autofunction:: mpmath.upper_gamma


Exponential integrals
.....................

.. autofunction:: mpmath.ei
.. autofunction:: mpmath.e1
.. autofunction:: mpmath.expint


Logarithmic integral
....................

.. autofunction:: mpmath.li


Trigonometric integrals
.......................

.. autofunction:: mpmath.ci
.. autofunction:: mpmath.si


Hyperbolic integrals
....................

.. autofunction:: mpmath.chi
.. autofunction:: mpmath.shi


Error functions
...............

.. autofunction:: mpmath.erf
.. autofunction:: mpmath.erfc
.. autofunction:: mpmath.erfi
.. autofunction:: mpmath.erfinv


The normal distribution
.......................

.. autofunction:: mpmath.npdf
.. autofunction:: mpmath.ncdf


Fresnel integrals
.................

.. autofunction:: mpmath.fresnels
.. autofunction:: mpmath.fresnelc

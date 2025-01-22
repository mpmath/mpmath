Orthogonal polynomials
----------------------

An orthogonal polynomial sequence is a sequence of polynomials `P_0(x), P_1(x),
\ldots` of degree `0, 1, \ldots`, which are mutually orthogonal in the sense
that

.. math ::

    \int_S P_n(x) P_m(x) w(x) dx =
    \begin{cases}
    c_n \ne 0 & \text{if $m = n$} \\
    0         & \text{if $m \ne n$}
    \end{cases}

where `S` is some domain (e.g. an interval `[a,b] \in \mathbb{R}`) and `w(x)`
is a fixed *weight function*.  A sequence of orthogonal polynomials is
determined completely by `w`, `S`, and a normalization convention (e.g. `c_n =
1`).  Applications of orthogonal polynomials include function approximation and
solution of differential equations.

Orthogonal polynomials are sometimes defined using the differential equations
they satisfy (as functions of `x`) or the recurrence relations they satisfy
with respect to the order `n`.  Other ways of defining orthogonal polynomials
include differentiation formulas and generating functions.  The standard
orthogonal polynomials can also be represented as hypergeometric series (see
:doc:`hypergeometric`), more specifically using the Gauss hypergeometric
function `\,_2F_1` in most cases.  The following functions are generally
implemented using hypergeometric functions since this is computationally
efficient and easily generalizes.

For more information, see the `Wikipedia article on orthogonal polynomials
<http://en.wikipedia.org/wiki/Orthogonal_polynomials>`_.

Legendre functions
..................

.. autofunction:: mpmath.legendre
.. autofunction:: mpmath.legenp
.. autofunction:: mpmath.legenq


Chebyshev polynomials
.....................

.. autofunction:: mpmath.chebyt
.. autofunction:: mpmath.chebyu


Jacobi polynomials
..................

.. autofunction:: mpmath.jacobi


Gegenbauer polynomials
......................

.. autofunction:: mpmath.gegenbauer


Hermite polynomials
...................

.. autofunction:: mpmath.hermite


Laguerre polynomials
....................

.. autofunction:: mpmath.laguerre


Spherical harmonics
...................

.. autofunction:: mpmath.spherharm

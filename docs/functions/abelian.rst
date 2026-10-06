Abelian functions
-----------------

Riemann theta functions generalize the Jacobi theta functions from one
complex variable to several. They arose in Riemann's nineteenth-century
theory of Abelian functions and compact Riemann surfaces. Today they are
fundamental tools in complex analysis and algebraic geometry, and occur in
finite-gap and quasiperiodic solutions of integrable systems.

Theta functions are quasiperiodic. For integer vectors
:math:`m,n\in\mathbb Z^g`, the zero-characteristic function satisfies

.. math ::

    \theta(z+m+\tau n\mid\tau)
    = \exp\!\left(-\pi i n^T\tau n-2\pi i n^Tz\right)
      \theta(z\mid\tau).

Although theta functions are quasiperiodic, suitable ratios and combinations
of them give multiply-periodic Abelian functions. In genus two, such Abelian
functions are meromorphic functions of two complex variables with a period
lattice of rank four, generated in normalized coordinates by the columns of
:math:`I_2` and :math:`\tau`.


Riemann theta functions
.......................

.. autofunction:: mpmath.rtheta

.. autofunction:: mpmath.rtheta_jet

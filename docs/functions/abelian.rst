Abelian functions
-----------------

Riemann theta functions generalize the Jacobi theta functions from one
complex variable to several. They arose in Riemann's nineteenth-century
theory of Abelian functions and compact Riemann surfaces. Today they are
fundamental tools in complex analysis and algebraic geometry, and occur in
finite-gap and quasiperiodic solutions of integrable systems.

In genus :math:`g`, the argument :math:`\mathbf{z}\in\mathbb C^g` is a
column vector and the period matrix
:math:`\boldsymbol{\tau}\in\mathbb C^{g\times g}` is symmetric with
positive-definite imaginary part.

Theta functions are quasiperiodic. For integer vectors
:math:`\mathbf{m},\mathbf{n}\in\mathbb Z^g`, the zero-characteristic function
satisfies

.. math ::

    \theta(\mathbf{z}+\mathbf{m}+\boldsymbol{\tau}\,\mathbf{n}
           \mid\boldsymbol{\tau})
    = \exp\!\left(-\pi i\,\mathbf{n}^{T}\boldsymbol{\tau}\,\mathbf{n}
                  -2\pi i\,\mathbf{n}^{T}\mathbf{z}\right)
      \theta(\mathbf{z}\mid\boldsymbol{\tau}).

Although theta functions are quasiperiodic, suitable ratios and combinations
of them give multiply-periodic Abelian functions. In genus two, such Abelian
functions are meromorphic functions of two complex variables with a period
lattice of rank four, generated in normalized coordinates by the columns of
:math:`\mathbf{I}_2` and :math:`\boldsymbol{\tau}`.


Riemann theta functions
.......................

.. autofunction:: mpmath.rtheta

.. autofunction:: mpmath.rtheta_jet

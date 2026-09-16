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

The following plots show two real slices and the modulus over two real
variables for genus-two period matrices. Similar slices and surfaces are
illustrated in `DLMF section 21.4 <https://dlmf.nist.gov/21.4>`_.

.. plot::

   import matplotlib.pyplot as plt
   from mpmath import j, plot, re, rtheta

   tau = [[j, -0.5], [-0.5, j]]
   curves = [
       lambda x: re(rtheta([x, x/2], tau)),
       lambda x: re(rtheta([x, 2*x], tau)),
   ]
   fig, ax = plt.subplots()
   plot(curves, [-2, 2], axes=ax)
   ax.legend([r"$z=(x,x/2)$", r"$z=(x,2x)$"])

.. plot::

   import matplotlib.pyplot as plt
   from mpmath import j, rtheta, splot

   tau = [[j, 0.5], [0.5, j]]
   fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
   surface = lambda x, y: abs(rtheta([x, y], tau))
   splot(surface, [-1, 1], [-1, 1], points=35, keep_aspect=False,
         axes=ax, plot3d_kwargs={"cmap": "viridis"})
   ax.set_zlabel(r"$|\theta(z\mid\tau)|$")


Riemann theta jets
..................

.. autofunction:: mpmath.rtheta_jet

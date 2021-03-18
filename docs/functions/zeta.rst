Zeta functions, L-series and polylogarithms
-------------------------------------------

This section includes the Riemann zeta functions
and associated functions pertaining to analytic number theory.


Riemann and Hurwitz zeta functions
..................................................

:func:`~mpmath.zeta`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.zeta(s,a=1,derivative=0)


Dirichlet L-series
..................................................

:func:`~mpmath.altzeta`
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.altzeta(s)

:func:`~mpmath.dirichlet`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.dirichlet(s,chi,derivative=0)


Stieltjes constants
...................

:func:`~mpmath.stieltjes`
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.stieltjes(n,a=1)


Zeta function zeros
......................................

These functions are used for the study of the Riemann zeta function
in the critical strip.

:func:`~mpmath.zetazero`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.zetazero(n, verbose=False)

:func:`~mpmath.nzeros`
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nzeros(t)

:func:`~mpmath.siegelz`
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.siegelz(t)

:func:`~mpmath.siegeltheta`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.siegeltheta(t)

:func:`~mpmath.grampoint`
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.grampoint(n)

:func:`~mpmath.backlunds`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.backlunds(t)


Lerch transcendent
................................

:func:`~mpmath.lerchphi`
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.lerchphi(z,s,a)


Polylogarithms and Clausen functions
.......................................

:func:`~mpmath.polylog`
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.polylog(s,z)

:func:`~mpmath.clsin`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.clsin(s, z)

:func:`~mpmath.clcos`
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.clcos(s, z)

:func:`~mpmath.polyexp`
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.polyexp(s,z)


Zeta function variants
..........................

:func:`~mpmath.primezeta`
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.primezeta(s)

:func:`~mpmath.secondzeta`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.secondzeta(s, a=0.015, **kwargs)

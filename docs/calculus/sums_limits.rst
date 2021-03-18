Sums, products, limits and extrapolation
----------------------------------------

The functions listed here permit approximation of infinite
sums, products, and other sequence limits.
Use :func:`mpmath.fsum` and :func:`mpmath.fprod`
for summation and multiplication of finite sequences.

Summation
..........................................

:func:`~mpmath.nsum`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nsum

:func:`~mpmath.sumem`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sumem

:func:`~mpmath.sumap`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sumap

Products
...............................

:func:`~mpmath.nprod`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.nprod

Limits (``limit``)
..................

:func:`~mpmath.limit`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.limit

Extrapolation
..........................................

The following functions provide a direct interface to
extrapolation algorithms. :func:`~mpmath.nsum` and :func:`~mpmath.limit`
essentially work by calling the following functions with an increasing
number of terms until the extrapolated limit is accurate enough.

The following functions may be useful to call directly if the
precise number of terms needed to achieve a desired accuracy is
known in advance, or if one wishes to study the convergence
properties of the algorithms.


:func:`~mpmath.richardson`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.richardson

:func:`~mpmath.shanks`
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.shanks

:func:`~mpmath.levin`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.levin

:func:`~mpmath.cohen_alt`
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.cohen_alt


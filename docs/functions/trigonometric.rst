Trigonometric functions
-----------------------

Except where otherwise noted, the trigonometric functions
take a radian angle as input and the inverse trigonometric
functions return radian angles.

The ordinary trigonometric functions are single-valued
functions defined everywhere in the complex plane
(except at the poles of tan, sec, csc, and cot).
They are defined generally via the exponential function,
e.g.

.. math ::

    \cos(x) = \frac{e^{ix} + e^{-ix}}{2}.

The inverse trigonometric functions are multivalued,
thus requiring branch cuts, and are generally real-valued
only on a part of the real line. Definitions and branch cuts
are given in the documentation of each function.
The branch cut conventions used by mpmath are essentially
the same as those found in most standard mathematical software,
such as Mathematica and Python's own ``cmath`` libary.

Degree-radian conversion
...........................................................

:func:`~mpmath.degrees`
^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.degrees(x)

:func:`~mpmath.radians`
^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.radians(x)

Trigonometric functions
.......................

:func:`~mpmath.cos`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.cos(x, **kwargs)

:func:`~mpmath.sin`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sin(x, **kwargs)

:func:`~mpmath.tan`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.tan(x, **kwargs)

:func:`~mpmath.sec`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sec(x)

:func:`~mpmath.csc`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.csc(x)

:func:`~mpmath.cot`
^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.cot(x)

Trigonometric functions with modified argument
........................................................

:func:`~mpmath.cospi`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.cospi(x, **kwargs)

:func:`~mpmath.sinpi`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sinpi(x, **kwargs)

Inverse trigonometric functions
................................................

:func:`~mpmath.acos`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.acos(x, **kwargs)

:func:`~mpmath.asin`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.asin(x, **kwargs)

:func:`~mpmath.atan`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.atan(x, **kwargs)

:func:`~mpmath.atan2`
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.atan2(y, x)

:func:`~mpmath.asec`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.asec(x)

:func:`~mpmath.acsc`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.acsc(x)

:func:`~mpmath.acot`
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.acot(x)

Sinc function
.............

:func:`~mpmath.sinc`
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sinc(x)

:func:`~mpmath.sincpi`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.sincpi(x)

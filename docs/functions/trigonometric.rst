Trigonometric functions
-----------------------

Except where otherwise noted, the trigonometric functions take a radian angle
as input and the inverse trigonometric functions return radian angles.

The ordinary trigonometric functions are single-valued functions defined
everywhere in the complex plane (except at the poles of tan, sec, csc, and
cot).  They are defined generally via the exponential function, e.g.

.. math ::

    \cos(x) = \frac{e^{ix} + e^{-ix}}{2}.

The inverse trigonometric functions are multivalued, thus requiring branch
cuts, and are generally real-valued only on a part of the real line.
Definitions and branch cuts are given in the documentation of each function.
The branch cut conventions used by mpmath are essentially the same as those
found in most standard mathematical software, such as Mathematica and Python's
own ``cmath`` libary.

Degree-radian conversion
........................

.. autofunction:: mpmath.degrees
.. autofunction:: mpmath.radians

Trigonometric functions
.......................

.. autofunction:: mpmath.cos
.. autofunction:: mpmath.sin
.. autofunction:: mpmath.tan
.. autofunction:: mpmath.sec
.. autofunction:: mpmath.csc
.. autofunction:: mpmath.cot


Trigonometric functions with modified argument
..............................................

.. autofunction:: mpmath.cospi
.. autofunction:: mpmath.sinpi


Inverse trigonometric functions
...............................

.. autofunction:: mpmath.acos
.. autofunction:: mpmath.asin
.. autofunction:: mpmath.atan
.. autofunction:: mpmath.atan2
.. autofunction:: mpmath.asec
.. autofunction:: mpmath.acsc
.. autofunction:: mpmath.acot


Sinc function
.............

.. autofunction:: mpmath.sinc
.. autofunction:: mpmath.sincpi

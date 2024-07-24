Plotting
========

If `matplotlib <https://matplotlib.org/>`_ is available, the functions
:func:`~mpmath.plot` and :func:`~mpmath.cplot` can be used to plot functions
respectively as x-y graphs and in the complex plane.  Also,
:func:`~mpmath.splot` can be used to produce 3D surface plots.

Function curve plots
-----------------------

.. plot::

   from mpmath import cos, plot, sin
   plot([cos, sin], [-4, 4])

.. autofunction:: mpmath.plot

Complex function plots
-------------------------

.. plot::

   from mpmath import cplot, fp
   fp.cplot(fp.gamma, points=100000)

.. autofunction:: mpmath.cplot

3D surface plots
----------------

.. plot::

   from mpmath import cos, pi, sin, splot
   r, R = 1, 2.5
   f = lambda u, v: [r*cos(u), (R+r*sin(u))*cos(v), (R+r*sin(u))*sin(v)]
   splot(f, [0, 2*pi], [0, 2*pi])

.. autofunction:: mpmath.splot


Numerical inverse Laplace transform
-----------------------------------

One-step algorithm (``invertlaplace``)
......................................

.. autofunction:: mpmath.invertlaplace

Specific algorithms
...................

Fixed Talbot algorithm
~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: mpmath.calculus.inverselaplace.FixedTalbot
   :members:

Gaver-Stehfest algorithm
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: mpmath.calculus.inverselaplace.Stehfest
   :members:

de Hoog, Knight & Stokes algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: mpmath.calculus.inverselaplace.deHoog
   :members:

Cohen acceleration algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: mpmath.calculus.inverselaplace.Cohen
   :members:

Manual approach
...............

It is possible and sometimes beneficial to re-create some of the
functionality in ``invertlaplace``. This could be used to compute the
Laplace-space function evaluations in a different way.  For example,
the Laplace-space function evaluations could be the result of a
quadrature or sum, solution to a system of ordinary differential
equations, or possibly computed in parallel from some external library
or function call.

A trivial example showing the process (which could be implemented
using the existing interface):

>>> from mpmath import calculus, convert, exp, mp
>>> myTalbot = calculus.inverselaplace.FixedTalbot(mp)
>>> t = convert(0.25)
>>> myTalbot.calc_laplace_parameter(t)
>>> fp = lambda p: 1/(p + 1) - 1/(p + 1000)
>>> ft = lambda t: exp(-t) - exp(-1000*t)
>>> fpvec = [fp(p) for p in myTalbot.p]
>>> ft(t)-myTalbot.calc_time_domain_solution(fpvec,t,manual_prec=True)
mpf('1.928300179528890061756872185e-21')

This manual approach is also useful to look at the Laplace parameter,
order, or working precision which were computed.

>>> myTalbot.degree
34

Credit
......

The numerical inverse Laplace transform functionality was contributed
to mpmath by Kristopher L. Kuhlman in 2017. The Cohen method was contributed
to mpmath by Guillermo Navas-Palencia in 2022.

Bessel functions and related functions
--------------------------------------

The functions in this section arise as solutions to various differential
equations in physics, typically describing wavelike oscillatory behavior or a
combination of oscillation and exponential decay or growth.  Mathematically,
they are special cases of the confluent hypergeometric functions `\,_0F_1`,
`\,_1F_1` and `\,_1F_2` (see :doc:`hypergeometric`).


Bessel functions
................

.. autofunction:: mpmath.besselj
.. autofunction:: mpmath.j0
.. autofunction:: mpmath.j1
.. autofunction:: mpmath.bessely
.. autofunction:: mpmath.besseli
.. autofunction:: mpmath.besselk


Bessel function zeros
.....................

.. autofunction:: mpmath.besseljzero
.. autofunction:: mpmath.besselyzero


Hankel functions
................

.. autofunction:: mpmath.hankel1
.. autofunction:: mpmath.hankel2


Spherical Bessel functions
..........................

.. autofunction:: mpmath.spherical_jn
.. autofunction:: mpmath.spherical_yn


Kelvin functions
................

.. autofunction:: mpmath.ber
.. autofunction:: mpmath.bei
.. autofunction:: mpmath.ker
.. autofunction:: mpmath.kei


Struve functions
................

.. autofunction:: mpmath.struveh
.. autofunction:: mpmath.struvel


Anger-Weber functions
.....................

.. autofunction:: mpmath.angerj
.. autofunction:: mpmath.webere


Lommel functions
................

.. autofunction:: mpmath.lommels1
.. autofunction:: mpmath.lommels2


Airy and Scorer functions
.........................

.. autofunction:: mpmath.airyai
.. autofunction:: mpmath.airybi
.. autofunction:: mpmath.airyaizero
.. autofunction:: mpmath.airybizero
.. autofunction:: mpmath.scorergi
.. autofunction:: mpmath.scorerhi


Coulomb wave functions
......................

.. autofunction:: mpmath.coulombf
.. autofunction:: mpmath.coulombg
.. autofunction:: mpmath.coulombc


Confluent U and Whittaker functions
...................................

.. autofunction:: mpmath.hyperu(a, b, z)
.. autofunction:: mpmath.whitm(k,m,z)
.. autofunction:: mpmath.whitw(k,m,z)


Parabolic cylinder functions
............................

.. autofunction:: mpmath.pcfd
.. autofunction:: mpmath.pcfu
.. autofunction:: mpmath.pcfv
.. autofunction:: mpmath.pcfw

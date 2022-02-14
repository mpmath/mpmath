Bessel functions and related functions
--------------------------------------

The functions in this section arise as solutions to various differential equations in physics, typically describing wavelike oscillatory behavior or a combination of oscillation and exponential decay or growth. Mathematically, they are special cases of the confluent hypergeometric functions `\,_0F_1`, `\,_1F_1` and `\,_1F_2` (see :doc:`hypergeometric`).


Bessel functions
...................................................

:func:`~mpmath.besselj`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: mpmath.besselj(n,x,derivative=0)
.. autofunction:: mpmath.j0(x)
.. autofunction:: mpmath.j1(x)

:func:`~mpmath.bessely`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.bessely(n,x,derivative=0)

:func:`~mpmath.besseli`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.besseli(n,x,derivative=0)

:func:`~mpmath.besselk`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.besselk(n,x)


Bessel function zeros
...............................

:func:`~mpmath.besseljzero`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.besseljzero(v,m,derivative=0)

:func:`~mpmath.besselyzero`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.besselyzero(v,m,derivative=0)


Hankel functions
................

:func:`~mpmath.hankel1`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.hankel1(n,x)

:func:`~mpmath.hankel2`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.hankel2(n,x)


Kelvin functions
................

:func:`~mpmath.ber`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ber

:func:`~mpmath.bei`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: mpmath.bei

:func:`~mpmath.ker`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.ker

:func:`~mpmath.kei`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.kei


Struve functions
...................................................

:func:`~mpmath.struveh`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.struveh

:func:`~mpmath.struvel`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.struvel


Anger-Weber functions
...................................................

:func:`~mpmath.angerj`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.angerj

:func:`~mpmath.webere`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.webere


Lommel functions
...................................................

:func:`~mpmath.lommels1`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.lommels1

:func:`~mpmath.lommels2`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.lommels2

Airy and Scorer functions
...............................................

:func:`~mpmath.airyai`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.airyai(z, derivative=0, **kwargs)

:func:`~mpmath.airybi`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.airybi(z, derivative=0, **kwargs)

:func:`~mpmath.airyaizero`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.airyaizero(k, derivative=0)

:func:`~mpmath.airybizero`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.airybizero(k, derivative=0, complex=0)

:func:`~mpmath.scorergi`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.scorergi(z, **kwargs)

:func:`~mpmath.scorerhi`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.scorerhi(z, **kwargs)


Coulomb wave functions
...............................................

:func:`~mpmath.coulombf`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.coulombf(l,eta,z)

:func:`~mpmath.coulombg`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.coulombg(l,eta,z)

:func:`~mpmath.coulombc`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.coulombc(l,eta)

Confluent U and Whittaker functions
...................................

:func:`~mpmath.hyperu`
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.hyperu(a, b, z)

:func:`~mpmath.whitm`
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.whitm(k,m,z)

:func:`~mpmath.whitw`
^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.whitw(k,m,z)

Parabolic cylinder functions
.................................

:func:`~mpmath.pcfd`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.pcfd(n,z,**kwargs)

:func:`~mpmath.pcfu`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.pcfu(a,z,**kwargs)

:func:`~mpmath.pcfv`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.pcfv(a,z,**kwargs)

:func:`~mpmath.pcfw`
^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: mpmath.pcfw(a,z,**kwargs)

Hypergeometric functions
------------------------

The functions listed in :doc:`expintegrals`, :doc:`bessel` and
:doc:`orthogonal`, and many other functions as well, are merely particular
instances of the generalized hypergeometric function `\,_pF_q`.  The functions
listed in the following section enable efficient direct evaluation of the
underlying hypergeometric series, as well as linear combinations, limits with
respect to parameters, and analytic continuations thereof.  Extensions to
twodimensional series are also provided.  See also the basic or q-analog of the
hypergeometric series in :doc:`qfunctions`.

For convenience, most of the hypergeometric series of low order are provided as
standalone functions.  They can equivalently be evaluated using
:func:`~mpmath.hyper`.  As will be demonstrated in the respective docstrings,
all the ``hyp#f#`` functions implement analytic continuations and/or asymptotic
expansions with respect to the argument `z`, thereby permitting evaluation for
`z` anywhere in the complex plane.  Functions of higher degree can be computed
via :func:`~mpmath.hyper`, but generally only in rapidly convergent instances.

Most hypergeometric and hypergeometric-derived functions accept optional
keyword arguments to specify options for :func:`~mpmath.hypercomb` or
:func:`~mpmath.hyper`.  Some useful options are *maxprec*, *maxterms*,
*zeroprec*, *accurate_small*, *hmag*, *force_series*, *asymp_tol* and
*eliminate*.  These options give control over what to do in case of slow
convergence, extreme loss of accuracy or evaluation at zeros (these two cases
cannot generally be distinguished from each other automatically), and singular
parameter combinations.

Common hypergeometric series
............................

.. autofunction:: mpmath.hyp0f1
.. autofunction:: mpmath.hyp1f1
.. autofunction:: mpmath.hyp1f2
.. autofunction:: mpmath.hyp2f0
.. autofunction:: mpmath.hyp2f1
.. autofunction:: mpmath.hyp2f2
.. autofunction:: mpmath.hyp2f3
.. autofunction:: mpmath.hyp3f2


Generalized hypergeometric functions
....................................

.. autofunction:: mpmath.hyper
.. autofunction:: mpmath.hypercomb


Meijer G-function
.................

.. autofunction:: mpmath.meijerg


Bilateral hypergeometric series
...............................

.. autofunction:: mpmath.bihyper


Hypergeometric functions of two variables
.........................................

.. autofunction:: mpmath.hyper2d
.. autofunction:: mpmath.appellf1
.. autofunction:: mpmath.appellf2
.. autofunction:: mpmath.appellf3
.. autofunction:: mpmath.appellf4

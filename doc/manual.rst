.. -*- rest -*-

=============
Mpmath manual
=============

:Author: Fredrik Johansson
:E-mail: fredrik.johansson@gmail.com
:Updated: 2008-03-09
:Mpmath version: 0.7

.. section-numbering::

.. contents::
    :local:

About mpmath
============

Mpmath is a pure-Python library for arbitrary-precision floating-point arithmetic. It implements all the functions found in Python's ``math`` and ``cmath`` modules (``exp``, ``log``, ``sin``...), plus a few nonelementary special functions (``gamma``, ``zeta``...), and has utilities for arbitrary-precision numerical differentiation, integration, root-finding, and interval arithmetic. It supports unlimited exponents, has full support for complex numbers, and offers better performance than Python's standard ``decimal`` library.

Mpmath is lightweight (~100 KB), free (BSD license), and easy to install or include in other software due to being written entirely in Python without any external dependencies.

Installation
============

You can install the latest released version of mpmath by running::

    python easy_install.py mpmath

or, on Windows::

    C:\<pythonpath>\Scripts\easy_install.exe mpmath

Alternatively, you can manually download the latest released version of mpmath from the `mpmath website
<http://code.google.com/p/mpmath/>`_ or the `Python Package Index <http://pypi.python.org/pypi>`_. Either run the binary installer (Windows only) or extract the source archive and run::

    python setup.py install

Debian and Ubuntu users can ``apt-get`` mpmath; `package information <http://packages.debian.org/python-mpmath>`_ is available on the Debian website.

After the setup has completed, you should be able to fire up the interactive Python interpreter and do the following::

    >>> from mpmath import *
    >>> mp.dps = 50
    >>> print mpf(2) ** mpf('0.5')
    1.4142135623730950488016887242096980785696718753769
    >>> print 2*pi
    6.2831853071795864769252867665590057683943387987502

Basic arithmetic
================

Working with mpmath numbers
---------------------------

Mpmath provides two main numerical types: ``mpf`` and ``mpc``. The ``mpf`` type is analogous to Python's built-in ``float``. It holds a real number or one of the special values ``inf`` (positive infinity), ``-inf`` and ``nan`` (not-a-number, indicating an indeterminate result). You can create ``mpf`` instances from strings, integers, floats, and other ``mpf`` instances:

    >>> mpf(4)
    mpf('4.0')
    >>> mpf(2.5)
    mpf('2.5')
    >>> mpf("1.25e6")
    mpf('1250000.0')
    >>> mpf(mpf(2))
    mpf('2.0')
    >>> mpf("inf")
    mpf('+inf')

An ``mpc`` represents a complex number in rectangular form as a pair of ``mpf`` instances. It can be constructed from a Python ``complex``, a real number, or a pair of real numbers:

    >>> mpc(2,3)
    mpc(real='2.0', imag='3.0')
    >>> mpc(complex(2,3)).imag
    mpf('3.0')

You can mix ``mpf`` and ``mpc`` instances with each other and with Python numbers:

    >>> mpf(3) + 2*mpf('2.5') + 1.0
    mpf('9')
    >>> mpc(1j)**0.5
    mpc(real='0.70710678118654757', imag='0.70710678118654757')

Prettier output can be obtained by using ``str()`` or ``print``, which hide the ``mpf`` and ``mpc`` constructor signatures and suppress small rounding artifacts:

    >>> mpf("3.14159")
    mpf('3.1415899999999999')
    >>> print mpf("3.14159")
    3.14159
    >>> print mpc(1j)**0.5
    (0.707106781186548 + 0.707106781186548j)

Controlling precision
---------------------

Mpmath uses a global working precision; it does not keep track of the precision or accuracy of individual numbers. Performing an arithmetic operation or calling ``mpf()`` rounds the result to the current working precision. The working precision is controlled by a special object called ``mp``, which has the following default state:

    >>> mp
    Mpmath settings:
      mp.prec = 53                [default: 53]
      mp.dps = 15                 [default: 15]
      mp.rounding = 'nearest'     [default: 'nearest']

The term *precision* (**prec**) always refers to the arithmetic precision measured in bits. The *decimal precision* is called the **dps** (short for *decimal places*). Binary and decimal precision are related roughly according to the formula ``prec = 3.33*dps``. For example, it takes a precision of roughly 333 bits to hold an approximation of pi that is accurate to 100 decimal places (actually slightly more than 333 bits is used; see the section "Decimal issues" below).

Changing one property of the ``mp`` object automatically updates the other; usually you just want to change the ``dps`` value:

    >>> mp.dps = 100
    >>> mp.dps
    100
    >>> mp.prec
    336

When you've set the precision level, all ``mpf`` operations are carried out at that precision:

    >>> mp.dps = 50
    >>> mpf(1) / 6
    mpf('0.16666666666666666666666666666666666666666666666666656')
    >>> mp.dps = 25
    >>> mpf(2) ** mpf('0.5')
    mpf('1.414213562373095048801688713')

The precision of complex arithmetic is also controlled by the ``mp`` object:

    >>> mp.dps = 10
    >>> mpc(1,2) / 3
    mpc(real='0.3333333333321', imag='0.6666666666642')

The number of digits with which numbers are printed by default is determined by the working precision. To specify the number of digits to show without changing the working precision, use the ``nstr`` and ``nprint`` functions:

    >>> mp.dps = 15
    >>> a = mpf(1) / 6
    >>> a
    mpf('0.16666666666666666')
    >>> nstr(a, 8)
    '0.16666667'
    >>> nprint(a, 8)
    0.16666667
    >>> nstr(a, 50)
    '0.16666666666666665741480812812369549646973609924316'

Temporarily changing the precision
..................................

It is often useful to change the precision during only part of a calculation. A way to temporarily increase the precision and then restore it is as follows:

    >>> mp.prec += 2
     (...)
    >>> mp.prec -= 2

In Python 2.5, the ``with`` statement along with the mpmath functions ``workprec``, ``workdps``, ``extraprec`` and ``extradps`` can be used to temporarily change precision in a more safe manner:

    >>> from __future__ import with_statement
    >>> with workdps(20):
    ...     print mpf(1)/7
    ...     with extradps(10):
    ...         print mpf(1)/7
    ...
    0.14285714285714285714
    0.142857142857142857142857142857
    >>> mp.dps
    15

The ``with`` statement ensures that the precision gets reset when exiting the block, even in the case that an exception is raised. (The effect of the ``with`` statement can be emulated in Python 2.4 by using a ``try/finally`` block.)

The ``workprec`` family of functions can also be used as function decorators:

    >>> @workdps(6)
    ... def f():
    ...     return mpf(1)/3
    ...
    >>> f()
    mpf('0.33333331346511841')

Caveat: providing correct input
-------------------------------

Note that when creating a new ``mpf``, the value will at most be as accurate as the input. **Be careful when mixing mpmath numbers with Python floats**. When working at high precision, fractional ``mpf`` values should be created from strings or integers:

>>> mp.dps = 30
>>> mpf(10.9)   # bad
mpf('10.9000000000000003552713678800501')
>>> mpf('10.9')  # good
mpf('10.9')
>>> mpf(109) / mpf(10)   # also good
mpf('10.9')

(Binary fractions such as 0.5, 1.5, 0.75, 0.125, etc, are generally safe, however, since those can be represented exactly by Python floats.)

Magical numbers
---------------



Mathematical functions
----------------------

High-level features
===================

Numerical integration
---------------------

Numerical differentiation
-------------------------

Root-finding with the secant method
-----------------------------------

The function ``secant`` calculates a root of a given function using the secant method. A good initial guess for the location of the root is required for the method to be effective, so it is somewhat more appropriate to think of the secant method as a root-polishing method than a root-finding method.

If the rough location of the root is known, the secant method can be used to refine it to very high precision in only a few steps. If the root is a first-order root, only roughly log(prec) iterations are required. (The secant method is far less efficient for double roots.) A particularly efficient general approach is to compute the initial approximation using a machine precision solver (for example using one of SciPy's many solvers), and then refining it to high precision using mpmath's ``secant`` method.

Simple examples
...............

A simple example use of the secant method is to compute pi as the root of sin(*x*) closest to *x* = 3.

    >>> mp.dps = 30
    >>> print secant(sin, 3)
    3.14159265358979323846264338328

The secant methods can be used to find complex roots of analytic functions, although it must in that case generally be given a nonreal starting value (or else it will never leave the real line).

    >>> mp.dps = 15
    >>> print secant(lambda x: x**3 + 2*x + 1, j)
    (0.226698825758202 + 1.46771150871022j)

A nice application is to compute nontrivial roots of the Riemann zeta function with many digits (good initial values are needed for convergence):

    >>> mp.dps = 30
    >>> print secant(zeta, 0.5+14j)
    (0.5 + 14.1347251417346937904572519836j)

A useful application is to compute inverse functions, such as the Lambert W function which is the inverse of *w* exp(*w*), given the first term of the solution's asymptotic expansion as the initial value:

    >>> def lambert(x):
    ...     return secant(lambda w: w*exp(w) - x, log(1+x))
    ...
    >>> mp.dps = 15
    >>> print lambert(1)
    0.567143290409784
    >>> print lambert(1000)
    5.2496028524016

Options
.......

Strictly speaking, the secant method requires two initial values. By default, you only have to provide the first point ``x0``; ``secant`` automatically sets the second point to ``x0 + 1/4``. Manually providing also the second point can help in some cases if ``secant`` fails to converge.

By default, ``secant`` performs a maximum of 20 steps, which can be increased or decreased using the ``maxsteps`` keyword argument. You can pass ``secant`` the option ``verbose=True`` to show detailed progress.

Polynomials
-----------

Polynomial evaluation
.....................

Polynomial functions can be evaluated using ``polyval``, which takes as input a list of coefficients and the desired evaluation point. The following example evaluates ``2 + 5*x + x^3`` at ``x = 3.5``:

    >>> mp.dps = 20
    >>> polyval([2, 5, 0, 1], mpf('3.5'))
    mpf('62.375')

With ``derivative=True``, both the polynomial and its derivative are evaluated at the same point:

    >>> polyval([2, 5, 0, 1], mpf('3.5'), derivative=True)
    (mpf('62.375'), mpf('41.75'))

The point and coefficient list may contain complex numbers.

Finding roots of polynomials
............................

The function ``polyroots`` computes all *n* real or complex roots of an *n*-th degree polynomial using complex arithmetic, and returns them along with an error estimate. As a simple example, it will successfully compute the two real roots ``3*x^2 - 7*x + 2`` (which are 1/3 and 2):

    >>> roots, err = polyroots([2, -7, 3])
    >>> print err
    2.66453525910038e-16
    >>> for root in roots:
    ...     print root
    ...
    (0.333333333333333 - 9.62964972193618e-35j)
    (2.0 + 1.5395124730131e-50j)

As should be expected from the internal use of complex arithmetic, the calculated roots have small but nonzero imaginary parts.

The following example computes all the 5th roots of unity; i.e. the roots of ``x^5 - 1``:

    >>> mp.dps = 20
    >>> for a in polyroots([-1, 0, 0, 0, 0, 1])[0]:
    ...     print a
    ...
    (-0.8090169943749474241 + 0.58778525229247312917j)
    (1.0 + 0.0j)
    (0.3090169943749474241 + 0.95105651629515357212j)
    (-0.8090169943749474241 - 0.58778525229247312917j)
    (0.3090169943749474241 + -0.95105651629515357212j)

Interval arithmetic
-------------------

Technical details
=================

Doing a high-precision calculation in mpmath typically just amounts to setting the precision and entering a formula. However, some knowledge of mpmath's terminology and internal number model can be useful to avoid common errors, and is recommended for trying more advanced calculations.

Representation of numbers
-------------------------

Mpmath uses binary arithmetic. A binary floating-point number is a number of the form ``man * 2^exp`` where both ``man`` (the *mantissa*) and ``exp`` (the *exponent*) are integers. Some examples of floating-point numbers are given in the following table.

  +--------+----------+----------+
  | Number | Mantissa | Exponent |
  +========+==========+==========+
  |    3   |    3     |     0    |
  +--------+----------+----------+
  |   10   |    5     |     1    |
  +--------+----------+----------+
  |  -16   |   -1     |     4    |
  +--------+----------+----------+
  |  1.25  |    5     |    -2    |
  +--------+----------+----------+

Note that the representation as defined so far is not unique; one can always multiply the mantissa by 2 and subtract 1 from the exponent with no change in the numerical value. In mpmath, numbers are always normalized so that ``man`` is an odd number, unless it is 0; we take zero to have ``man = exp = 0``. With these conventions, every representable number has a unique representation. (Mpmath does not currently distinguish between positive and negative zero.)

Simple mathematical operations are now easy to define. Due to uniqueness, equality testing of two numbers simply amounts to separately checking equality of the mantissas and the exponents. Multiplication of nonzero numbers is straightforward: ``(m*2^e) * (n*2^f) = (m*n) * 2^(e+f)``. Addition is a bit more involved: we first need to multiply the mantissa of one of the operands by a suitable power of 2 to obtain equal exponents.

More technically, mpmath represents a floating-point number as a 4-tuple ``(sign, man, exp, bc)`` where `sign` is 0 or 1 (indicating positive vs negative) and the mantissa is nonnegative; ``bc`` (*bitcount*) is the size of the absolute value of the mantissa as measured in bits. Though redundant, keeping a separate sign field and explicitly keeping track of the bitcount significantly speeds up arithmetic (the bitcount, especially, is frequently needed but slow to compute from scratch due to the lack of a Python built-in function for the purpose).

The special numbers ``+inf``, ``-inf`` and ``nan`` are represented internally by a zero mantissa and a nonzero exponent.

For further details on how the arithmetic is implemented, refer to the mpmath source code. The basic arithmetic operations are found in the ``lib.py`` module; many functions there are commented extensively.

Precision and accuracy
----------------------

Contrary to popular superstition, floating-point numbers  do not come with an inherent "small uncertainty". Every binary floating-point number is an exact rational number. With arbitrary-precision integers used for the mantissa and exponent, floating-point numbers can be added, subtracted and multiplied *exactly*. In particular, integers and integer multiples of 1/2, 1/4, 1/8, 1/16, etc. can be represented, added and multiplied exactly in binary floating-point.

The reason why floating-point arithmetic is generally approximate is that we set a limit to the size of the mantissa for efficiency reasons. The maximum allowed width (bitcount) of the mantissa is called the precision or ``prec`` for short. Sums and products are exact as long as the absolute value of the mantissa is smaller than ``2^prec``. As soon as the mantissa becomes larger than this threshold, we truncate it to have at most  ``prec`` bits (the exponent is incremented accordingly to preserve the magnitude of the number), and it is this operation that typically introduces numerical errors. Division is also not generally exact; although we can add and multiply exactly by setting the precision high enough, no precision is high enough to represent for example 1/3 exactly (the same obviously applies for roots, trigonometric functions, etc).

Decimal issues
..............

Unfortunately for some applications, decimal fractions fall into the category of numbers that generally cannot be represented exactly in binary floating-point form. For example, none of the numbers ``0.1``, ``0.01``, ``0.001`` has an exact representation as a binary floating-point number. Mpmath does not fully solve this problem; users who need *exact* decimal fractions should look at the ``decimal`` module in Python's standard library.

There are a few subtle differences between binary and decimal precision. Precision and accuracy do not always correlate when translating from binary to decimal. As a simple example, the number 0.1 has a decimal precision of 1 digit but is an infinitely accurate representation of 1/10. Conversely, the number 2^-50 has a binary representation with 1 bit of precision that is infinitely accurate; the same number can actually be represented exactly as a decimal, but doing so requires 35 significant digits:

    0.00000000000000088817841970012523233890533447265625

Generally, it works out to just choose 1000 * 3.33 bits of precision in order to obtain 1000 decimal digits. In fact, mpmath will do the conversion automatically for you: you can enter a desired *dps* value and mpmath will automatically choose the appropriate *prec*. More precisely, mpmath uses the following formulas to translate between prec and dps::

  dps(prec) = max(1, int(round(int(prec) / C - 1)))

  prec(dps) = max(1, int(round((int(dps) + 1) * C)))

where ``C = log(10)/log(2)`` is the exact version of the "3.33" conversion ratio. Note that the dps is set 1 decimal digit lower than the corresponding binary precision. This margin is added to ensure that *n*-digit decimal numbers, when converted to binary, will retain all *n* digits correct when converted back to decimal.

  * The ``str`` decimal precision is roughly one digit less than the exact equivalent binary precision, to hide minor rounding errors and artifacts resulting from binary-decimal conversion

  * The ``repr`` decimal precision is roughly one digit greater to ensure that ``x == eval(repr(x))`` holds, i.e. that numbers can be converted to strings and back losslessly.

For example, the standard precision is 53 bits, which corresponds to a dps value of 15. The actual decimal precision given by 53 bits is 15.95 ~= 16.

The dps value controls the number of digits to display when printing numbers with ``str``, while the decimal precision used by ``repr`` is set two or three digits higher. For example, with 15 dps we have::

    >>> str(pi)
    '3.14159265358979'
    >>> repr(+pi)
    "mpf('3.1415926535897931')"

Rounding
--------

There are several different strategies for rounding a too large mantissa or a result that cannot at all be represented exactly in floating-point form (such as ``log(2)``). Mpmath supports the following rounding modes:

  +-----------+---------------------------------------------------------+
  | Name      | Direction                                               |
  +===========+=========================================================+
  | Floor     | Towards negative infinity                               |
  +-----------+---------------------------------------------------------+
  | Ceiling   | Towards positive infinity                               |
  +-----------+---------------------------------------------------------+
  | Down      | Towards 0                                               |
  +-----------+---------------------------------------------------------+
  | Up        | Away from 0                                             |
  +-----------+---------------------------------------------------------+
  | Nearest   | To nearest; to the nearest even number on a tie         |
  +-----------+---------------------------------------------------------+

The first four modes are called *directed* rounding schemes and are useful for implementing interval arithmetic; they are also fast. Rounding to nearest, which mpmath uses by default, is the slowest but most accurate method.

The arithmetic operations ``+``, ``-``, ``*`` and ``/`` acting on real floating-point numbers always round their results *correctly* in mpmath; that is, they are guaranteed to give exact results when possible, they always round in the intended direction, and they don't round to a number farther away than necessary. Exponentiation by an integer *n* preserves directions but may round too far if either the mantissa or *n* is very large.

Evaluation of transcendental functions (as well as square roots) is generally performed by computing an approximation with finite precision slightly higher than the target precision, and rounding the result. This gives correctly rounded results with a high probability, but can be wrong in exceptional cases.

Rounding for radix conversion is a slightly tricky business. When converting to a binary floating-point number from a decimal string, mpmath writes the number as an exact fraction and performs correct rounding division if the number is of reasonable size (roughly, larger than 10^-100 and smaller than 10^100). When converting from binary to decimal, mpmath first performs an approximate radix conversion with slightly increased precision, then truncates the resulting decimal number to remove long sequences of trailing 0's and 9's, and finally rounds to nearest, rounding up (away from zero) on a tie.

Exponent range
--------------

In hardware floating-point arithmetic, the size of the exponent is restricted to a fixed range: regular Python floats have a range between roughly 10^-300 and 10^300. Mpmath uses arbitrary precision integers for both the mantissa and the exponent, so numbers can be as large in magnitude as permitted by computer's memory. Mpmath can for example hold an approximation of a large Mersenne prime::

    >>> print mpf(2)**32582657 - 1
    1.24575026015369e+9808357

Or why not 1 googolplex::

    >>> print mpf(10) ** (10**100)
    1.0e+100000000000000000000000000000000000000000000000000
    00000000000000000000000000000000000000000000000000

Some care may be necessary when working with extremely large numbers. Although arithmetic is safe, it is for example futile to attempt to compute ``exp`` of either of the above two numbers. Mpmath does not complain when asked to perform such a calculation, but instead chugs away on the problem to the best of its ability, assuming that computer resources are infinite. In the worst case, this will be slow and allocate a huge amount of memory; if entirely impossible Python will at some point raise ``OverflowError: long int too large to convert to int``.

In some situations, it would be more convenient if mpmath would "round" extremely small numbers to 0 and extremely large numbers to ``inf``, and directly raise an exception or return ``nan`` if there is no reasonable chance of finishing a computation. This option is not available, but could be implemented in the future on demand.

Compatibility
-------------

The floating-point arithmetic provided by processors that conform to the IEEE 754 *double precision* standard has a precision of 53 bits and rounds to nearest. (Additional precision and rounding modes are usually available, but regular double precision arithmetic should be the most familiar to Python users, since the Python ``float`` type corresponds to an IEEE double with rounding to nearest on most systems.)

This corresponds roughly to a decimal accuracy of 15 digits, and is the default precision used by mpmath. Thus, under normal circumstances, mpmath should produce identical results to Python ``float`` operations. This is not always true, for the following reasons:

1) Hardware floats have a limited exponent range, as discussed above. Machine floats very close to the exponent limit may be rounded subnormally, meaning that they lose precision. Python may also raise an exception instead of rounding a ``float`` subnormally.

2) Hardware floating-point operations don't always round correctly. This is commonly the case for hardware implementations of transcendental functions like ``log`` and ``sin``, but even square roots seem to be inaccurate on some systems, and mpmath has been run on at least one modern system where Python's builtin ``float`` multiplication was inaccurate, causing mpmath's float compatibility tests to fail.

3) Mpmath may of course have bugs. (However, the basic arithmetic has been tested fairly thoroughly by now. (1) and (2) are the more common causes of discrepancies.)


Performance notes
=================

In rough numbers, Python floats are 100 times slower than raw hardware floats, and mpmath floats at standard precision are 100 times slower than Python floats. It's fortunate that a modern CPU does some 10^9 operations per second, at least leaving some 10^5 operations per second for mpmath (which is plenty for many uses). Because most time at low precision levels is spent on bookkeeping and interpreter overhead, the execution time increases sublinearly with small increments in precision. 50-digit arithmetic is essentially as fast as 15-digit arithmetic.  Asymptotically, mpmath arithmetic is as fast as Python big integer arithmetic, which is actually quite efficient up to several thousand digits (thanks to the use of Karatsuba multiplication).

Optimization tricks
-------------------

There are a few tricks that can significantly speed up mpmath code at low to medium precision (up to a few hundred digits):

  * Repeated type conversions from floats, strings and integers should be avoided.

  * Changing the rounding mode to *floor* can give a slight speedup.

  * The JIT compiler `psyco <http://psyco.sourceforge.net/>`_ fairly consistently speeds up mpmath about 2x.

  * An additional 2x gain is possible by using the low-level functions in ``mpmath.lib`` instead of ``mpf`` instances.

Here follows a simple example demonstrating some of these options.

Original algorithm (0.028 seconds)::

    x = mpf(1)
    for i in range(1000):
        x += 0.1

Preconverting the float constant (0.080 seconds)::

    x = mpf(1)
    one_tenth = mpf(0.1)
    for i in range(1000):
        x += one_tenth

With psyco (0.0036 seconds)::

    import psyco; psyco.full()
    x = mpf(1)
    one_tenth = mpf(0.1)
    for i in range(1000):
        x += one_tenth

With psyco and low-level functions (0.0017 seconds)::

    import psyco; psyco.full()
    x = from_int(1)
    one_tenth = from_float(0.1)
    for i in range(1000):
        x = fadd(x, one_tenth, 53, round_nearest)

The last version is 16.5 times faster than the first. Not all calculations can be sped up the same way, of course, or doing so may just be inconvenient.

Using the right tool
--------------------

Many calculations can be done with ordinary floating-point arithmetic, and only in special cases require multiprecision arithmetic (for example to avoid overflows in corner cases). In these situations, it may be possible to write code that uses fast regular floats by default, and automatically (or manually) falls backs to mpmath only when needed. Python's dynamic namespaces and ability to compile code on the fly are helpful. Here is a simple (probably not failsafe) example::

    import math
    import mpmath

    def evalmath(expr):
        try:
            r = eval(expr, math.__dict__)
        except OverflowError:
            r = eval(expr, mpmath.__dict__)
            try:
                r = float(r)
            except OverflowError:
                pass
        return r

    >>> evalmath('sin(3)')
    0.14112000805986721
    >>>
    >>> evalmath('exp(10000)')
    mpf('8.8068182256629216e+4342')
    >>>
    >>> evalmath('exp(10000) / exp(10000)')
    1.0

.. -*- rest -*-

=============
Mpmath manual
=============

:Author: Fredrik Johansson <fredrik.johansson@gmail.com>
:Updated: 2008-03-06
:Mpmath version: 0.7

.. section-numbering::

.. contents::
    :local:

About mpmath
============

Mpmath is a pure-Python library for arbitrary-precision floating-point arithmetic. It implements all the functions found in Python's ``math`` and ``cmath`` modules (``exp``, ``log``, ``sin``...), plus a few nonelementary special functions (``gamma``, ``zeta``...), and has utilities for arbitrary-precision numerical differentiation, integration, root-finding, and interval arithmetic. It has extensive support for complex numbers and is much faster (typically 10 or 100 times) than Python's standard ``decimal`` library.

Mpmath is lightweight (~100 KB), free (BSD license), and easy to install or include in other software due to being written in pure Python without any external dependencies.

Installation
============

You can install the latest released version of mpmath by running::

    python easy_install.py mpmath

or, on Windows::

    C:\<pythonpath>\Scripts\easy_install.exe mpmath

Alternatively, you can manually download the latest released version of mpmath from the `mpmath website
<http://code.google.com/p/mpmath/>`_ or the `Python Package Index <http://pypi.python.org/pypi>`_. To install mpmath, either run the binary installer (Windows only) or extract the source archive and run::

    python setup.py install

Debian users can ``apt-get`` mpmath; `package information <http://packages.debian.org/python-mpmath>`_ is available on the Debian website.

After the setup has completed, you should be able to fire up the interactive Python interpreter and do the following::

    >>> from mpmath import *
    >>> setdps(50)      # set working precision to 50 decimals
    >>> print mpf(2) ** mpf('0.5')    # mpf is an arbitrary-precision float type
    1.4142135623730950488016887242096980785696718753769
    >>> print 2*pi
    6.2831853071795864769252867665590057683943387987502

Basic principles
================

Doing a high-precision calculation in mpmath typically just amounts to setting the precision and entering a formula. However, some knowledge of mpmath's terminology and internal number model can be very useful to avoid common errors, and is recommended for trying more advanced calculations.

If you just want to try out mpmath's features, feel free to skip ahead to the next section.

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

More technically, mpmath represents a floating-point number as a 3-tuple ``(man, exp, bc)`` where ``bc`` (*bitcount*) is the size of the absolute value of the mantissa as measured in bits. Though redundant, storing the bitcount significantly speeds up arithmetic, since the bitcount is frequently needed but slow to compute from scratch due to the lack of a Python built-in function for the purpose.

Mpmath numbers can also hold the special values ``inf`` (positive infinity), ``-inf`` and ``nan`` (not-a-number, indicating an invalid result). These numbers use a nonnumerical internal representation.

For further details on how the arithmetic is implemented, refer to the mpmath source code. The basic arithmetic operations are found in the ``lib.py`` module; many functions there are commented extensively.

Precision and accuracy
----------------------

Contrary to popular superstition, floating-point numbers  do not come with an inherent "small uncertainty". Every binary floating-point number is an exact rational number. With arbitrary-precision integers used for the mantissa and exponent, floating-point numbers can be added, subtracted and multiplied *exactly*. In particular, integers and integer multiples of 1/2, 1/4, 1/8, 1/16, etc. can be represented, added and multiplied exactly in binary floating-point.

The reason why floating-point arithmetic is generally approximate is that we set a limit to the size of the mantissa for efficiency reasons. The maximum allowed width (bitcount) of the mantissa is called the precision or ``prec`` for short. Sums and products are exact as long as the absolute value of the mantissa is smaller than ``2^prec``. As soon as the mantissa becomes larger than this threshold, we truncate it to have at most  ``prec`` bits (the exponent is incremented accordingly to preserve the magnitude of the number), and it is this operation that typically introduces numerical errors. Division is also not generally exact; although we can add and multiply exactly by setting the precision high enough, no precision is high enough to represent for example 1/3 exactly.

Decimals
........

Unfortunately for some applications, decimal fractions fall into the category of numbers that generally cannot be represented exactly in binary floating-point form. For example, none of the numbers ``0.1``, ``0.01``, ``0.001`` has an exact representation as a binary floating-point number. Mpmath does not fully solve this problem; users who need *exact* decimal fractions should look at the ``decimal`` module in Python's standard library. However, mpmath can work with approximations of decimal fractions that are much better than those of standard floats. Instead of ``0.1000000000000000056``, you can have:

    0.10000000000000000000000000000000000000000028

or an approximation with any higher finite accuracy. The idea behind binary floating-point arithmetic is that one often does not need to print every value; instead, a calculation involving several steps can be performed entirely using efficient binary arithmetic, and only the final result needs to be converted to a decimal numeral that can be read by humans. If the calculation is done with precision a little higher than the target accuracy, rounding off the last few digits in the output gives a correct decimal value.

There are a few subtle differences between binary and decimal precision. In mpmath, the term *precision* (**prec**) always refers to the arithmetic precision measured in bits. The *decimal precision* is called the **dps** (short for *decimal places*). Binary and decimal precision are related roughly according to the formula ``prec = 3.33*dps``. For example, it takes a precision of roughly 333 bits to hold an approximation of pi that is accurate to 100 decimal places.

However, the meaning of "decimal precision" can depend slightly on context. Precision and accuracy are not always correlated when translating from binary to decimal. As a simple example, the number 0.1 has a decimal precision of 1 digit but is an infinitely accurate representation of 1/10. Conversely, the number 2^-50 has a binary representation with 1 bit of precision that is infinitely accurate; the same number can actually be represented exactly as a decimal, but doing so requires 35 significant digits:

    0.00000000000000088817841970012523233890533447265625

Generally, it works out to just think "I want 1000 digits, so I'll set the precision to ``1000 * 3.33 = 3330`` bits". In fact, as documented below, mpmath will do this conversion automatically for you, meaning that you can enter a desired *dps* value and mpmath will automatically choose the appropriate *prec*. More precisely, mpmath uses the following formulas to translate between prec and dps::

  dps(prec) = max(1, int(round(int(prec) / C - 1)))

  prec(dps) = max(1, int(round((int(dps) + 1) * C)))

where ``C = log(10)/log(2)`` is the exact version of the "3.33" conversion ratio. Note that the dps is set 1 decimal digit lower than the corresponding binary precision. This margin is added to ensure that *n*-digit decimal numbers, when converted to binary, will retain all *n* digits correct when converted back to decimal.

The dps value controls the number of digits to display when printing numbers with ``str``, while the decimal precision used by ``repr`` is set two digits higher. For example, with 15 dps we have::

    >>> str(pi)
    '3.14159265358979'
    >>> repr(+pi)
    "mpf('3.1415926535897931')"

In other words, the ``str`` decimal precision is roughly one digit less than the binary precision, and the ``repr`` decimal precision is roughly one digit greater. The extra precision for ``repr`` is to ensure that ``x == eval(repr(x))`` holds, i.e. that numbers can be converted to strings and back losslessly. (Note: it seems that this invariance does not hold on all precision levels, although it does in fact work at the standard precision. The conversion formula may be updated in a future version of mpmath.)

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
  | Half-down | To nearest; down if right between                       |
  +-----------+---------------------------------------------------------+
  | Half-up   | To nearest; right if right between                      |
  +-----------+---------------------------------------------------------+
  | Half-even | To nearest; to the nearest even number if right between |
  +-----------+---------------------------------------------------------+

The first four modes are called *directed* rounding schemes and are useful for implementing interval arithmetic. The three *nearby* rounding modes generally provide greater accuracy, but are on the other hand slower. Half-even rounding, which mpmath uses by default, is both the most accurate and the slowest method.

The arithmetic operations ``+``, ``-``, ``*`` and ``/`` always round their results *correctly*; that is, they are guaranteed to give exact results when possible, they always round in the intended direction, and they don't round to a number farther away than necessary. Exponentiation by an integer *n* preserves directions but may round too far if either the mantissa or *n* is very large.

Radix conversion and evaluation of transcendental functions (as well as square roots) is generally performed by computing an approximation with finite precision slightly higher than the target precision, and rounding the result. This gives correctly rounded results with a high probability, but can be wrong in bad cases.

When converting to a binary floating-point number from a decimal string, mpmath writes the number as an exact fraction and performs correct rounding division if the number is of reasonable size (roughly, larger than 10^-100 and smaller than 10^100). Similar comments apply when converting from binary to decimal: after performing an approximate radix conversion with slightly increased precision, the result is first truncated to remove long sequences of trailing 0's and 9's, and then rounded in the half-up direction to the desired number of decimal digits.

Exponent range
--------------

In hardware floating-point arithmetic, the size of the exponent is restricted to a fixed range: regular Python floats have a range between roughly 10^-300 and 10^300. Mpmath uses arbitrary precision integers for both the mantissa and the exponent, so numbers can be as large in magnitude as permitted by computer's memory. Mpmath can for example hold an approximation of a large Mersenne prime::

    >>> print (mpf(2)**32582657 - 1)
    1.24575026015369e+9808357

Or why not 1 googolplex::

    >>> print mpf(10) ** (10**100)
    1.0e+100000000000000000000000000000000000000000000000000000000000000000000000000
    00000000000000000000000000

Some care may be necessary when working with extremely large numbers. Although arithmetic is safe, it is for example futile to attempt to compute ``exp`` of either of the above two numbers. Mpmath does not complain when asked to perform such a calculation, but instead chugs away on the problem to the best of its ability, assuming that computer resources are infinite. In the worst case, this will be slow and allocate a huge amount of memory; if entirely impossible Python will at some point raise ``OverflowError: long int too large to convert to int``.

In some situations, it would be more convenient if mpmath would "round" extremely small numbers to 0 and extremely large numbers to ``inf``, and directly raise an exception or return ``nan`` if there is no reasonable chance of finishing a computation. This option is not available, but could be implemented in the future on demand.

Compatibility
-------------

The floating-point arithmetic provided by processors that conform to the IEEE 754 *double precision* standard has a precision of 53 bits and uses *half-even* rounding. (Additional precision and rounding modes are usually available, but regular double precision arithmetic should be the most familiar to Python users, since the Python ``float`` type corresponds to an IEEE double with half-even rounding on most systems.)

This corresponds roughly to a decimal accuracy of 15 digits, and is the default precision used by mpmath, which also uses half-even rounding by default. Thus, under normal circumstances, mpmath should produce identical results to Python ``float`` operations. This is not always true, for two reasons:

1) Hardware floats have a limited exponent range, as discussed above. Numbers very close to the exponent limit may be rounded subnormally, meaning that they lose precision.

2) Hardware floats don't always round correctly. (This is commonly the case for transcendental functions like ``log`` and ``sin``, but even square roots seem to be inaccurate on most systems, and mpmath has been run on at least one modern system where Python's builtin ``float`` multiplication was inaccurate, causing mpmath's comparative tests to fail.)

3) Mpmath may of course have bugs. (However, the basic arithmetic has been tested fairly thoroughly by now. (1) and (2) are the more common causes of discrepancies.)


Working with mpmath numbers
===========================

Setting precision
-----------------

Mathematical functions
----------------------

High-level functions
====================

Numerical integration
---------------------

Numerical differentiation
-------------------------

Root-finding
------------

Polynomials
-----------

Interval arithmetic
-------------------

Notes on performance
====================

In rough numbers, Python floats are 100 times slower than raw hardware floats, and mpmath floats at standard precision are 100 times slower than Python floats. It's fortunate that a modern CPU does some 10^9 operations per second, at least leaving some 10^5 operations per second for mpmath. 100,000 operations per second is fortunately plenty for many applications; mpmath also implements elementary functions like ``exp`` and ``sin`` efficiently, so that they are only slightly (3x to 5x) slower than plain arithmetic.

Because most time at low precision levels is constant overhead, the execution time increases sublinearly with small increments in precision. 50-digit arithmetic is essentially as fast as 15-digit arithmetic.  Asymptotically, mpmath arithmetic is as fast as Python big integer arithmetic, which is actually quite efficient up to 10,000 digits or so (due to the use of Karatsuba multiplication).

There are a few tricks that can speed up mpmath code at low to medium precision (up to a few hundred digits). Changing the rounding mode to *floor* gives a slight speedup, on the order 10-50%, at the cost of reduced accuracy. The JIT compiler `Psyco
<http://psyco.sourceforge.net/>`_ fairly consistently speeds up mpmath about 2x. An additional 2x gain is possible by using the low-level functions in ``mpmath.lib``.

A simple trick that can pay off in some cases is to store constants to avoid repeated type conversions. The second of the following code snippets is a whole 3x faster than the first::

    x = mpf(1)
    for i in range(1000):
        x += 0.5

    x = mpf(1)
    onehalf = mpf(0.5)
    for i in range(1000):
        x += onehalf

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

If you find that mpmath is orders of magnitude too slow for your needs, you should definitely look elsewhere, for example at the highly optimized C library `MPFR <http://www.mpfr.org/>`_, or the Python interface for the GNU Multiprecision Library, `GMPY <http://code.google.com/p/gmpy/>`_.
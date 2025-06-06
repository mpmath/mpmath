Setting up mpmath
=================

Mpmath requires at least Python 3.9.  It has been tested with CPython 3.9
through 3.14 and for PyPy 3.10

Download and installation
-------------------------

Using pip
.........

Releases are registered on PyPI, so you can install latest release
of the Mpmath with pip::

    pip install mpmath

or some specific version with::

    pip install mpmath==1.3.0

You can install also extra dependencies, e.g. `gmpy2
<https://github.com/aleaxit/gmpy>`_ support::

    pip install mpmath[gmpy2]

.. tip::

    Use :mod:`venv` to create isolated Python environment first,
    instead of installing everything system-wide.

Debian/Ubuntu
.............

Debian and Ubuntu users can install mpmath with::

    sudo apt install python3-mpmath

See `debian <http://packages.debian.org/stable/python/python3-mpmath>`_ and
`ubuntu <https://launchpad.net/ubuntu/+source/mpmath>`_ package information;
please verify that you are getting the latest version.

OpenSUSE
........

Mpmath is provided in the "Science" repository for all recent versions of
`openSUSE <https://www.opensuse.org/>`_. To add this repository to the YAST
software management tool, see
https://en.opensuse.org/SDB:Add_package_repositories

Look up https://download.opensuse.org/repositories/science/ for a list
of supported OpenSUSE versions.

Current development version
...........................

If you are a developer or like to get the latest updates as they come, be sure
to install from git::

    git clone git://github.com/mpmath/mpmath.git
    cd mpmath
    pip install -e .[develop,docs]

Checking that it works
......................

After the setup has completed, you should be able to fire up the interactive
Python interpreter and do the following::

    >>> from mpmath import mp, mpf, pi
    >>> mp.dps = 50
    >>> print(mpf(2) ** mpf('0.5'))
    1.4142135623730950488016887242096980785696718753769
    >>> print(2*pi)
    6.2831853071795864769252867665590057683943387987502

.. tip::

   :ref:`Run mpmath as a module <cli>` for interactive work::

       python -m mpmath


Using gmpy2 (optional)
----------------------

If `gmpy2 <https://github.com/aleaxit/gmpy>`_ version 2.2.0 or later is
installed on your system, mpmath will automatically detect it and transparently
use gmpy2 integers instead of Python integers.  This makes mpmath much faster,
especially at high precision (approximately above 100 digits).

To verify that mpmath uses gmpy2, check the internal variable ``BACKEND`` is
equal to 'gmpy'.

Using the gmpy2 backend can be disabled by setting the ``MPMATH_NOGMPY``
environment variable.  Note that the mode cannot be switched during runtime;
mpmath must be re-imported for this change to take effect.

Running tests
-------------

It is recommended that you run mpmath's full set of unit tests to make sure
everything works. The `pytest <https://pytest.org/>`_ is a required dependence
for testing.  The tests are located in the ``tests`` subdirectory of the mpmath
source tree.  They can be run using::

    pytest --pyargs mpmath

Developers may run tests from the source tree with::

    pytest

If any test fails, please send a detailed bug report to the `mpmath issue
tracker <https://github.com/mpmath/mpmath/issues>`_.

Compiling the documentation
---------------------------

If you downloaded the source package, the text source for these documentation
pages is included in the ``docs`` directory.  The documentation can be compiled
to pretty HTML using `Sphinx <https://www.sphinx-doc.org/>`_::

    sphinx-build --color -W --keep-going -b html docs build/sphinx/html

The create a PDF::

    sphinx-build --color -W --keep-going -b latex docs build/sphinx/latex
    make -C build/sphinx/latex all-pdf

Some additional demo scripts are available in the ``demo`` directory included
in the source package.

Mpmath under Sage
-------------------

Mpmath is a standard package in `Sage <https://sagemath.org/>`_, in version 4.1 or later of Sage.
Mpmath is preinstalled a regular Python module, and can be imported as usual within Sage::

    ----------------------------------------------------------------------
    | Sage Version 4.1, Release Date: 2009-07-09                         |
    | Type notebook() for the GUI, and license() for information.        |
    ----------------------------------------------------------------------
    sage: import mpmath
    sage: mpmath.mp.dps = 50
    sage: print(mpmath.mpf(2) ** 0.5)
    1.4142135623730950488016887242096980785696718753769

In Sage, mpmath can alternatively be imported via the interface library
``sage.libs.mpmath.all``. For example::

    sage: import sage.libs.mpmath.all as mpmath

This module provides a few extra conversion functions, including ``mpmath.call()``
which permits calling any mpmath function with Sage numbers as input, and getting 
Sage ``RealNumber`` or ``ComplexNumber`` instances
with the appropriate precision back::

    sage: w = mpmath.call(mpmath.erf, 2+3*I, prec=100)
    sage: w
    -20.829461427614568389103088452 + 8.6873182714701631444280787545*I
    sage: type(w)
    <type 'sage.rings.complex_number.ComplexNumber'>
    sage: w.prec()
    100

See the help for ``sage.libs.mpmath.all`` for further information.

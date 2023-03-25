Setting up mpmath
=================

Download and installation
-------------------------

Installer
.........

The mpmath setup files can be downloaded from the `Python Package Index <http://pypi.python.org/pypi/mpmath/>`_. Download the source package (available as both .zip and .tar.gz), extract it, open the extracted directory, and run

    ``python setup.py install``

Using pip
.........

Releases are registered on PyPI, so you can install latest release
of the Mpmath with pip

    ``pip install mpmath``

or some specific version with

    ``pip install mpmath==0.19``

Using setuptools
................

If you have `setuptools <http://pypi.python.org/pypi/setuptools>`_ installed, you can download and install mpmath in one step by running:

    ``easy_install mpmath``

or

    ``python -m easy_install mpmath``

If you have an old version of mpmath installed already, you may have to pass ``easy_install`` the ``-U`` flag to force an upgrade.

If installation fails, try deleting the following folders: .eggs, mpmath.egg-info, dist, build

Debian/Ubuntu
.............

Debian and Ubuntu users can install mpmath with

    ``sudo apt-get install python-mpmath``

See `debian <http://packages.debian.org/stable/python/python-mpmath>`_ and `ubuntu <https://launchpad.net/ubuntu/+source/mpmath>`_ package information; please verify that you are getting the latest version.

OpenSUSE
........

Mpmath is provided in the "Science" repository for all recent versions of `openSUSE <http://www.opensuse.org/en/>`_. To add this repository to the YAST software management tool, see http://en.opensuse.org/SDB:Add_package_repositories

Look up http://download.opensuse.org/repositories/science/ for a list
of supported OpenSUSE versions and use http://download.opensuse.org/repositories/science/openSUSE_11.1/
(or accordingly for your OpenSUSE version) as the repository URL for YAST.

Current development version
...........................

The git repository is https://github.com/fredrik-johansson/mpmath

Checking that it works
......................

After the setup has completed, you should be able to fire up the interactive Python interpreter and do the following::

    >>> from mpmath import *
    >>> mp.dps = 50
    >>> print(mpf(2) ** mpf('0.5'))
    1.4142135623730950488016887242096980785696718753769
    >>> print(2*pi)
    6.2831853071795864769252867665590057683943387987502

*Note: if you have are upgrading mpmath from an earlier version, you may have to manually uninstall the old version or remove the old files.*

Using gmpy (optional)
---------------------

By default, mpmath uses Python integers internally. If `gmpy <http://code.google.com/p/gmpy/>`_ version 1.03 or later is installed on your system, mpmath will automatically detect it and transparently use gmpy integers intead. This makes mpmath much faster, especially at high precision (approximately above 100 digits).

To verify that mpmath uses gmpy, check the internal variable ``BACKEND`` is not equal to 'python':

    >>> import mpmath.libmp
    >>> mpmath.libmp.BACKEND # doctest:+SKIP
    'gmpy'

The gmpy mode can be disabled by setting the MPMATH_NOGMPY environment variable. Note that the mode cannot be switched during runtime; mpmath must be re-imported for this change to take effect.

Running tests
-------------

It is recommended that you run mpmath's full set of unit tests to make sure everything works. The `py.test <https://pytest.org/>`_ is a required dependence for testing.  The tests are located in the ``tests`` subdirectory of the main mpmath directory. They can be run using::

    py.test --pyargs mpmath

Doctests can be run with::

    py.test --doctest-modules mpmath

If any test fails, please send a detailed bug report to the `mpmath issue tracker <https://github.com/fredrik-johansson/mpmath/issues>`_.

To run the tests with support for gmpy disabled, set ``MPMATH_NOGMPY`` environment variable.

To enable extra diagnostics, use, set ``MPMATH_STRICT`` environment variable.

Compiling the documentation
---------------------------

Building the documentation requires `Sphinx <http://sphinx.pocoo.org/>`_. 

Use::

    python setup.py build_sphinx -c docs -b html,latex

The create a PDF::

    make -C build/sphinx/latex all-pdf

Some additional demo scripts are available in the ``demo`` directory included in the source package.

Mpmath under Sage
-------------------

Mpmath is a standard package in `Sage <http://sagemath.org/>`_, in version 4.1 or later of Sage.
Mpmath is preinstalled a regular Python module, and can be imported as usual within Sage::

    ----------------------------------------------------------------------
    | Sage Version 4.1, Release Date: 2009-07-09                         |
    | Type notebook() for the GUI, and license() for information.        |
    ----------------------------------------------------------------------
    sage: import mpmath
    sage: mpmath.mp.dps = 50
    sage: print mpmath.mpf(2) ** 0.5
    1.4142135623730950488016887242096980785696718753769

The mpmath installation under Sage automatically use Sage integers for asymptotically fast arithmetic,
so there is no need to install GMPY::

    sage: mpmath.libmp.BACKEND
    'sage'

In Sage, mpmath can alternatively be imported via the interface library
``sage.libs.mpmath.all``. For example::

    sage: import sage.libs.mpmath.all as mpmath

This module provides a few extra conversion functions, including :func:`call`
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



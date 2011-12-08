mpmath
======

A Python library for arbitrary-precision floating-point arithmetic.

Website: http://code.google.com/p/mpmath
Main author: Fredrik Johansson <fredrik.johansson@gmail.com>

Mpmath is free software released under the New BSD License (see the
LICENSE file for details)

0. History and credits
----------------------

The following people have contributed patches and new features
to mpmath:

* Pearu Peterson <pearu.peterson@gmail.com>
* Mario Pernici <mario.pernici@mi.infn.it>
* Ondrej Certik <ondrej@certik.cz>
* Vinzent Steinberg <vinzent.steinberg@gmail.com>
* Nimish Telang <ntelang@gmail.com>
* Mike Taschuk <mtaschuk@ece.ualberta.ca>
* Case Van Horsen <casevh@gmail.com>
* Jorn Baayen <jorn.baayen@gmail.com>
* Chris Smith <smichr@gmail.com>
* Juan Arias de Reyna <arias@us.es>
* Ioannis Tziakos <itziakos@gmail.com>

For a detailed changelog, including individual contributions,
see the CHANGES file.

Fredrik's work on mpmath during summer 2008 was sponsored by Google
as part of the Google Summer of Code program.

Fredrik's work on mpmath during summer 2009 was sponsored by the
American Institute of Mathematics under the support of the National Science
Foundation Grant No. 0757627 (FRG: L-functions and Modular Forms).

Any opinions, findings, and conclusions or recommendations expressed in this
material are those of the author(s) and do not necessarily reflect the
views of the sponsors.

Credit also goes to:

* The authors of the GMP library and the Python wrapper
  gmpy, enabling mpmath to become much faster at
  high precision
* The authors of MPFR, pari/gp, MPFUN, and other arbitrary-
  precision libraries, whose documentation has been helpful
  for implementing many of the algorithms in mpmath
* Wikipedia contributors; Abramowitz & Stegun; Gradshteyn & Ryzhik;
  Wolfram Research for MathWorld and the Wolfram Functions site.
  These are the main references used for special functions
  implementations.
* George Brandl for developing the Sphinx documentation tool
  used to build mpmath's documentation

Release history:

* Version 0.17 released on February 1, 2011
* Version 0.16 released on September 24, 2010
* Version 0.15 released on June 6, 2010
* Version 0.14 released on February 5, 2010
* Version 0.13 released on August 13, 2009
* Version 0.12 released on June 9, 2009
* Version 0.11 released on January 26, 2009
* Version 0.10 released on October 15, 2008
* Version 0.9 released on August 23, 2008
* Version 0.8 released on April 20, 2008
* Version 0.7 released on March 12, 2008
* Version 0.6 released on January 13, 2008
* Version 0.5 released on November 24, 2007
* Version 0.4 released on November 3, 2007
* Version 0.3 released on October 5, 2007
* Version 0.2 released on October 2, 2007
* Version 0.1 released on September 27, 2007

1. Download & installation
--------------------------

Mpmath requires Python 2.5 or later. It has been tested
with Python 2.5, 2.6, 2.7, 3.1 and 3.2.

The latest release of mpmath can be downloaded from the mpmath
website. It should also be available in the Python Package Index at
http://pypi.python.org/pypi

To install, unpack the mpmath archive and run

``python setup.py install``

Mpmath can also be installed using

``python -m easy_install mpmath``

The latest development code is available from
http://code.google.com/p/mpmath/source/checkout or at
https://github.com/fredrik-johansson/mpmath.

See the main documentation for more detailed instructions.

2. Running tests
----------------

The unit tests in mpmath/tests/ can be run via the script
runtests.py, but it is recommended to run them with py.test
(http://codespeak.net/py/dist/index.html), especially
to generate more useful reports in case there are failures.

You may also want to check out the demo scripts in the demo
directory.

3. Documentation
----------------

Documentation in reStructuredText format is available in the
doc directory included with the source package. These files
are human-readable, but can be compiled to prettier HTML using
the build.py script (requires Sphinx, http://sphinx.pocoo.org/).

See setup.txt in the documentation for more information.

The most recent documentation is also available in HTML
format on the mpmath website:

http://mpmath.googlecode.com/svn/trunk/doc/build/index.html

4. Known problems
-----------------

Mpmath is a work in progress. Major issues include:

* Some functions may return incorrect values when given extremely
  large arguments or arguments very close to singularities.

* Directed rounding works for arithmetic operations. It is implemented
  heuristically for other operations, and their results may be off by one
  or two units in the last place (even if otherwise accurate).

* Some IEEE 754 features are not available. Inifinities and NaN are
  partially supported; denormal rounding is currently not available
  at all.

* The interface for switching precision and rounding is not finalized.
  The current method is not threadsafe.

5. Help and bug reports
-----------------------

General questions and comments can be sent to the mpmath mailinglist,
mpmath@googlegroups.com

You can also report bugs and send patches to the mpmath issue tracker,
http://code.google.com/p/mpmath/issues/list

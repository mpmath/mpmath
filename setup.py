#!/usr/bin/env python

from distutils.core import setup

setup(name='mpmath',
      description = 'Python library for arbitrary-precision floating-point arithmetic',
      version='0.9',
      url='http://mpmath.googlecode.com',
      author='Fredrik Johansson',
      author_email='fredrik.johansson@gmail.com',
      license = 'BSD',
      packages=['mpmath', 'mpmath/tests'],
      classifiers=['Topic :: Scientific/Engineering :: Mathematics']
     )

#!/usr/bin/env python

from distutils.core import setup

setup(name='mpmath',
      description = 'Python library for arbitrary-precision floating-point arithmetic',
      version='0.15',
      url='http://mpmath.googlecode.com',
      author='Fredrik Johansson',
      author_email='fredrik.johansson@gmail.com',
      license = 'BSD',
      packages=['mpmath',
                'mpmath/libmp',
                'mpmath/calculus',
                'mpmath/functions',
                'mpmath/matrices',
                'mpmath/tests'],
      classifiers=['Topic :: Scientific/Engineering :: Mathematics']
     )

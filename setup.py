#!/usr/bin/env python

from distutils.core import setup

setup(name='mpmath',
      description = 'Python library for arbitrary-precision floating-point arithmetic',
      version='1.1.0',
      url='http://mpmath.org',
      author='Fredrik Johansson',
      author_email='fredrik.johansson@gmail.com',
      license = 'BSD',
      packages=['mpmath',
                'mpmath.libmp',
                'mpmath.calculus',
                'mpmath.functions',
                'mpmath.matrices',
                'mpmath.tests'],
      classifiers=['Topic :: Scientific/Engineering :: Mathematics']
     )

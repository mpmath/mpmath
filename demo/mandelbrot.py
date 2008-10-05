"""
This script uses the cplot function in mpmath to plot the
Mandelbrot set.

By default, Python's builtin cmath module is used for
the calculations, since this is about 6x faster and
increased precision does not improve image accuracy
at typical zoom levels.

However, switching to mpmath and increasing the
working precision could be useful when exploring the
set at extremely high zoom levels.
"""

import mpmath
import cmath

USE_CMATH = True
ITERATIONS = 50
POINTS = 100000
ESCAPE_RADIUS = 8

# Full plot
RE = [-2.5, 1.5]
IM = [-1.5, 1.5]

# A pretty subplot
#RE = [-0.96, -0.80]
#IM = [-0.35, -0.2]

if USE_CMATH:
    exp, log, cplx = cmath.exp, cmath.log, complex
else:
    exp, log, cplx = mpmath.exp, mpmath.log, mpmath.mpc

def mandelbrot(z):
    z = cplx(z)
    c = z
    for i in xrange(ITERATIONS):
        zprev = z
        z = z*z + c
        if abs(z) > ESCAPE_RADIUS:
            return exp(1j*(i + 1 - log(log(abs(z)))/log(2)))
    return 0

try:
    import psyco
    psyco.full()
except ImportError:
    pass

mpmath.cplot(mandelbrot, RE, IM, points=POINTS, verbose=1)

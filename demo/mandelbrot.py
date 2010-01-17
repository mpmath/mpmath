"""
This script uses the cplot function in mpmath to plot the Mandelbrot set.
By default, the fp context is used for speed. The mp context could be used
to improve accuracy at extremely high zoom levels.
"""

import mpmath
import cmath

ctx = mpmath.fp
# ctx = mpmath.mp

ITERATIONS = 50
POINTS = 100000
ESCAPE_RADIUS = 8

# Full plot
RE = [-2.5, 1.5]
IM = [-1.5, 1.5]

# A pretty subplot
#RE = [-0.96, -0.80]
#IM = [-0.35, -0.2]

def mandelbrot(z):
    c = z
    for i in xrange(ITERATIONS):
        zprev = z
        z = z*z + c
        if abs(z) > ESCAPE_RADIUS:
            return ctx.exp(1j*(i + 1 - ctx.log(ctx.log(abs(z)))/ctx.log(2)))
    return 0

try:
    import psyco
    psyco.full()
except ImportError:
    pass

ctx.cplot(mandelbrot, RE, IM, points=POINTS, verbose=1)

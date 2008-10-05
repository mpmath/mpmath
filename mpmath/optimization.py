#!/usr/bin/env python

from mptypes import *
from calculus import diffc

class Secant():
    """
    Solver generating pairs of approximative root and error.

    Needs starting points x0 and x1 close to the root.
    x1 defaults to x0 + 0.25.

    Pro:
    * converges fast
    Contra:
    * converges slowly for multiple roots
    """
    maxsteps = 30

    def __init__(self, f, x0, **kwargs):
        if len(x0) == 1:
            self.x0 = x0[0]
            self.x1 = self.x0 + 0.25
        elif len(x0) == 2:
            self.x0 = x0[0]
            self.x1 = x0[1]
        else:
            raise ValueError('expected 2 starting points, got %i' * len(x0))
        self.f = f

    def __iter__(self):
        f = self.f
        x0 = self.x0
        x1 = self.x1
        f0 = f(x0)
        while True:
            f1 = f(x1)
            l = x1 - x0
            if l == 0:
                break
            s = (f1 - f0) / l
            if s == 0:
                break
            x0, x1 = x1, x1 - f1/s
            f0 = f1
            yield x1, abs(l)

@extraprec(10)
def findroot(f, x0, solver=Secant, tol=None, verbose=False, **kwargs):
    """
    Find a root of f using x0 as starting point or interval.

    If not abs(f(root)) < tol an exception is raised.

    Arguments:
    f : one dimensional function
    x0 : starting point, several starting points or interval (depends on solver)
    tol : the returned solution has an error smaller than this
    verbose : print additional information for each iteration if true
    solver : a generator for f and x0 returning approximative solution and error
    maxsteps : after how many steps the solver will cancel
    df : first derivative of f (used by some solvers)
    d2f : second derivative of f (used by some solvers)
    """
    if not isinstance(x0, (list, tuple)):
        x0 = (x0,)
    iterations = solver(f, x0, **kwargs)
    if 'maxsteps' in kwargs:
        maxsteps = kwargs['maxsteps']
    else:
        maxsteps = iterations.maxsteps
    if tol is None:
        tol = eps
    i = 0
    for x, error in iterations:
        if verbose:
            print 'x:', x
            print 'error:', error
        i += 1
        if error < tol or i >= maxsteps:
            break
    if abs(f(x)) > tol:
        raise ValueError('Could not find root within given tolerance.\n'
                         'Try another starting point or tweak arguments.')
    return x

def multiplicity(f, root, tol=eps, maxsteps=10, **kwargs):
    """
    Return the multiplicity of a given root of f.

    Internally, numerical derivatives are used. This is very inefficient for
    higher order derviatives. You can be specify the n-th derivative using the
    dnf keyword.
    """
    kwargs['d0f'] = f
    for i in xrange(maxsteps):
        dfstr = 'd' + str(i) + 'f'
        if dfstr in kwargs:
            df = kwargs[dfstr]
        else:
            df = lambda x: diffc(f, x, i)
        if not abs(df(root)) < tol:
            break
    return i

__docformat__ = 'plaintext'

from mpmath.mptypes import *

def polyval(coeffs, x, derivative=False):
    """
    Given coefficients [c0, c1, c2, ..., cn], evaluate
    P(x) = c0 + c1*x + c2*x**2 + ... + cn*x**n.

    If derivative=True is set, a tuple (P(x), P'(x)) is returned.
    """
    p = mpnumeric(coeffs[-1])
    q = mpf(0)
    for c in coeffs[-2::-1]:
        if derivative:
            q = p + x*q
        p = c + x*p
    if derivative:
        return p, q
    else:
        return p

def polyroots(coeffs, maxsteps=20):
    """
    Numerically locate all (complex) roots of a polynomial using the
    Durand-Kerner method.

    This function returns a tuple (roots, err) where roots is a list of
    complex numbers sorted by absolute value, and err is an estimate of
    the maximum error. The polynomial should be given as a list of
    coefficients.

        >>> nprint(polyroots([24,-14,-1,1]), 4)
        ([(2.0 + 8.968e-44j), (3.0 + 1.156e-33j), (-4.0 + 0.0j)], 5.921e-16)
        >>> nprint(polyroots([2,3,4]))
        ([(-0.375 + -0.599479j), (-0.375 + 0.599479j)], 2.22045e-16)

    """
    deg = len(coeffs) - 1
    # Must be monic
    lead = mpnumeric(coeffs[-1])
    if lead == 1:
        coeffs = map(mpnumeric, coeffs)
    else:
        coeffs = [c/lead for c in coeffs]
    f = lambda x: polyval(coeffs, x)
    roots = [mpc((0.4+0.9j)**n) for n in range(deg)]
    error = [mpf(1) for n in range(deg)]
    for step in range(maxsteps):
        if max(error).ae(0):
            break
        for i in range(deg):
            if not error[i].ae(0):
                p = roots[i]
                x = f(p)
                for j in range(deg):
                    if i != j:
                        try:
                            x /= (p-roots[j])
                        except ZeroDivisionError:
                            continue
                roots[i] = p - x
                error[i] = abs(x)
    roots.sort(key=abs)
    err = max(error)
    err = max(err, ldexp(1, -getprec()+1))
    return roots, err

__version__ = '0.9'

from mptypes import *

from functions import *
from elliptic import jacobi_theta_1, jacobi_theta_2, jacobi_theta_3, \
    jacobi_theta_4, jacobi_elliptic_sn, jacobi_elliptic_cn, jacobi_elliptic_dn

from calculus import diff, diffc, secant, polyval, polyroots, quadosc
from calculus import sumem, chebyfit, sumsh, sumrich, limit
from calculus import odeint
from quadrature import quad, quadgl, quadts, TanhSinh, GaussLegendre

from identification import pslq, identify, findpoly

from matrices import matrix, eye, diag, zeros, ones, randmatrix, mnorm_1, \
    mnorm_oo, mnorm_F, norm_p
from linalg import lu_solve, inverse, residual, qr_solve, cholesky_solve, det, \
    cond


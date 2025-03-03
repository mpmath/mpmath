'''
This script calculates the constant in Gerver's solution to the moving sofa
problem.

See Finch, S. R. "Moving Sofa Constant." §8.12 in Mathematical Constants.
Cambridge, England: Cambridge University Press, pp. 519-523, 2003.
'''

from mpmath import cos, sin, pi, quad, findroot, mp

mp.prec = 113

eqs = [lambda A, B, φ, θ: (A*(cos(θ) - cos(φ)) - 2*B*sin(φ)
                           + (θ - φ - 1)*cos(θ) - sin(θ) + cos(φ) + sin(φ)),
       lambda A, B, φ, θ: (A*(3*sin(θ) + sin(φ)) - 2*B*cos(φ)
                           + 3*(θ - φ - 1)*sin(θ) + 3*cos(θ) - sin(φ) + cos(φ)),
       lambda A, B, φ, θ: A*cos(φ) - (sin(φ) + 0.5 - 0.5*cos(φ) + B*sin(φ)),
       lambda A, B, φ, θ: ((A + pi/2 - φ - θ) - (B - (θ - φ)*(1 + A)/2
                                                 - 0.25*(θ - φ)**2))]
A, B, φ, θ = findroot(eqs, (0, 0, 0, 0))

def r(α):
    if 0 <= α < φ:
        return 0.5
    if φ <= α < θ:
        return (1 + A + α - φ)/2
    if θ <= α < pi/2 - θ:
        return A + α - φ
    return B - (pi/2 - α - φ)*(1 + A)/2 - (pi/2 - α - φ)**2/4

s = lambda α: 1 - r(α)

def u(α):
    if φ <= α < θ:
        return B - (α - φ)*(1 + A)/2 - (α - φ)**2/4
    return A + pi/2 - φ - α

def du(α):
    if φ <= α < θ:
        return -(1 + A)/2 - (α - φ)/2
    return -1

def y(α, f):
    if α > pi/2 - θ:
        i = [0, φ, θ, pi/2 - θ, α]
    elif α > θ:
        i = [0, φ, θ, α]
    elif α > φ:
        i = [0, φ, α]
    else:
        i = i = [0, α]
    return 1 - quad(lambda x: f(x)*sin(x), i)

y1 = lambda α: y(α, r)
y2 = lambda α: y(α, s)
y3 = lambda α: y2(α) - u(α)*sin(α)

S1 = quad(lambda x: y1(x)*r(x)*cos(x), [0, φ, θ, pi/2 - θ, pi/2 - φ])
S2 = quad(lambda x: y2(x)*s(x)*cos(x), [0, φ, θ])
S3 = quad(lambda x: y3(x)*(u(x)*sin(x) - du(x)*cos(x) - s(x)*cos(x)),
          [φ, θ, pi/4])

print(2*(S1 + S2 + S3))

from mpmath import *

import os

try:
    import psyco
    psyco.full()
except ImportError:
    pass

fp = open("gallery.html", "w")
fp.write("<html><head><title>Gallery of mathematical functions</title></head><body>")
fp.write("<h1>Gallery of mathematical functions</h1>")

fp.write('''<p>Generated with <a href="http://mpmath.org">mpmath</a>
and <a href="http://matplotlib.org">matplotlib</a> by
<a href="gallery.py">this script</a>.</p>

<p>The complex plots show the magnitude as brightness (0 = black, &infin; = white) and the phase as hue (positive real = red, imaginary = lime, negative real = cyan, negative imaginary = violet).</p>
''')

def plotgroup(name, fs, descs, xlim=[-5, 5], ylim=None, repoints=500,
    relim=[-5,5], imlim=[-5, 5], singularities=[], complex_only=False):
    print name, "..."
    fp.write('<h2>%s</h2>' % name)
    f1name = '%s.png' % name
    fnname = name + ' C%i.png'
    if not (complex_only or os.path.exists(f1name)):
        plot(fs, xlim=xlim, ylim=ylim, file=f1name, points=repoints, dpi=40,
            singularities=singularities)
    for k, f in enumerate(fs):
        if not os.path.exists(fnname % k):
            cplot(f, file=(fnname % k), re=relim, im=imlim, points=35000,
            dpi=40, verbose=1)
    fp.write('<table><tr>')
    if not complex_only:
        fp.write('<td><img src="%s" /></td>' % f1name)
    for k in range(len(fs)):
        fp.write('<td><img src="%s" /></td>' % (fnname % k))
    fp.write('</tr><tr align="center">')
    desc = ",<br/>".join(descs)
    desc = desc[0].upper() + desc[1:]
    if not complex_only:
        fp.write('<td>%s</td>' % desc)
    for k, d in enumerate(descs):
        fp.write('<td>Complex %s</td>' % descs[k])
    fp.write('</tr></table>')
    fp.flush()

plotgroup("Polynomials", [lambda x: x, lambda x: 0.1*x**3 - x + 1],
    ["linear function f(x) = x", "cubic, 0.1 x<sup>3</sup> - x + 1"])
plotgroup("Rational functions", [lambda x: 1/x, lambda x: (x+4)/(x**5-3j*x**3+2)],
    ["inverse 1/x", ("rational function " \
      "(x+4)/(x<sup>5</sup>-3ix<sup>3</sup>+2)")], ylim=[-5,5])

plotgroup("Nonanalytic functions", [abs, floor], ["absolute value |x|",
    "floor function &lfloor;x&rfloor;"],
    singularities=[-4,-3,-2,-1,0,1,2,3,4])
plotgroup("Roots", [sqrt, cbrt], ["square root x<sup>1/2</sup>", "cube root x<sup>1/3</sup>"])
plotgroup("Exponentiation", [exp, ln], ["exponential function exp(x)",
    "natural logarithm ln(x)"], ylim=[-5,5])

def hypexp(x):
    z = x**(x**x)
    # Avoid float overflow in matplotlib
    if abs(z) > 100:
        return sign(z)*1000
    return z

plotgroup("Hyper-exponentiation", [lambda x: x**x, hypexp],
    ["2-tetration x<sup>x</sup>", "3-tetration x<sup>(x<sup>x</sup>)</sup>"],
    ylim=[-5, 5])
plotgroup("Lambert W function", [lambertw, lambda x: lambertw(x, -1)],
    ["W<sub>0</sub>(x) (branch 0)", "W<sub>&minus;1</sub>(x) (branch &minus;1)"],
    ylim=[-5, 5])

plotgroup("Trigonometric functions 1", [sin, cos], ["sine", "cosine"])
plotgroup("Trigonometric functions 2", [tan, cot], ["tangent", "cotangent"],
    ylim=[-5, 5], singularities=[-3*pi/2,-pi,-pi/2,pi/2,pi,3*pi/2])
plotgroup("Trigonometric functions 3", [sec, csc], ["secant", "cosecant"],
    ylim=[-5, 5], singularities=[-3*pi/2,-pi,-pi/2,pi/2,pi,3*pi/2])

plotgroup("Inverse trigonometric functions 1", [asin, acos],
    ["inverse sine", "inverse cosine"])
plotgroup("Inverse trigonometric functions 2", [atan, acot],
    ["inverse tangent", "inverse cotangent"])
plotgroup("Inverse trigonometric functions 3", [asec, acsc],
    ["inverse secant", "inverse cosecant"])

plotgroup("Hyperbolic functions 1", [sinh, cosh],
    ["hyperbolic sine", "hyperbolic cosine"], ylim=[-5, 5])
plotgroup("Hyperbolic functions 2", [tanh, coth],
    ["hyperbolic tangent", "hyperbolic cotangent"], ylim=[-5, 5])
plotgroup("Hyperbolic functions 3", [sech, csch],
    ["hyperbolic secant", "hyperbolic cosecant"], ylim=[-5, 5])

plotgroup("Inverse hyperbolic functions 1", [asinh, acosh],
    ["inverse hyperbolic sine", "inverse hyperbolic cosine"])
plotgroup("Inverse hyperbolic functions 2", [atanh, acoth],
    ["inverse hyperbolic tangent", "inverse hyperbolic cotangent"])
plotgroup("Inverse hyperbolic functions 3", [asech, acsch],
    ["inverse hyperbolic secant", "inverse hyperbolic cosecant"])

plotgroup("Error function", [lambda x: exp(-x**2), erf],
    ["Gaussian, exp(&minus;x<sup>2</sup>)", "error function erf(x)"])
plotgroup("Fresnel integrals", [fresnels, fresnelc],
    ["Fresnel integral S(x)", "Fresnel integral C(x)"])
plotgroup("Exponential integrals", [ei, li],
    ["exponential integral Ei(x)", "Logarithmic integral li(x)"], ylim=[-5, 5])
plotgroup("Trigonometric integrals", [si, ci],
    ["sine integral Si(x)", "cosine integral Ci(x)"])
plotgroup("Hyperbolic integrals", [shi, chi],
    ["hyperbolic sine integral Shi(x)", "hyperbolic cosine integral Chi(x)"],
    ylim=[-10, 10])

plotgroup("Airy functions", [airyai, airybi], ["Airy function Ai(x)",
    "Airy function Bi(x)"], xlim=[-8, 1], ylim=[-1, 1])
plotgroup("Bessel function of the first kind", [j0, j1],
    ["Bessel function J<sub>0</sub>(x)", "Bessel function J<sub>1</sub>(x)"],
    xlim=[-15, 15], relim=[-15,15], imlim=[-15,15])

plotgroup("Gamma function", [gamma], ["gamma function &Gamma;(x)"],
    ylim=[-5, 5], repoints=1500)
plotgroup("Polygamma functions", [digamma, trigamma],
    ["digamma function &psi;<sup>(0)</sup>(x)",
    "trigamma function &psi;<sup>(1)</sup>(x)"], ylim=[-20, 20])

plotgroup("Riemann zeta function", [zeta], ["Riemann zeta function &zeta;(x)"],
    ylim=[-5, 5], relim=[-20, 20], imlim=[-30, 30])

plotgroup("Jacobi theta functions 1-2", [lambda x: jtheta(1,x,0.5), lambda x: jtheta(2,x,0.5)],
    ["&theta;<sub>1</sub>(x,1/2)", "&theta;<sub>2</sub>(x,1/2)"])
plotgroup("Jacobi theta functions 3-4", [lambda x: jtheta(3,x,0.5), lambda x: jtheta(4,x,0.5)],
    ["&theta;<sub>3</sub>(x,1/2)", "&theta;<sub>4</sub>(x,1/2)"])

plotgroup("Derivative of Jacobi theta functions 1-2", [lambda x: djtheta(1,x,0.5), lambda x: djtheta(2,x,0.5)],
    ["&theta;'<sub>1</sub>(x,1/2)", "&theta;'<sub>2</sub>(x,1/2)"])
plotgroup("Derivative of Jacobi theta functions 3-4", [lambda x: djtheta(3,x,0.5), lambda x: djtheta(4,x,0.5)],
    ["&theta;'<sub>3</sub>(x,1/2)", "&theta;'<sub>4</sub>(x,1/2)"])

plotgroup("Jacobi theta functions 1-2, variable nome",
    [lambda q: jtheta(1,(1+j)/3,q), lambda q: jtheta(2,(1+j)/3,q)],
    ["&theta;<sub>1</sub>((1+i)/3,q)", "&theta;<sub>2</sub>((1+i)/3,q)"],
    xlim=[-1,1], ylim=[-1, 1], relim=[-1,1], imlim=[-1,1], complex_only=True)

plotgroup("Jacobi theta functions 3-4, variable nome",
    [lambda q: jtheta(3,(1+j)/3,q), lambda q: jtheta(4,(1+j)/3,q)],
    ["&theta;<sub>3</sub>((1+i)/3,q)", "&theta;<sub>4</sub>((1+i)/3,q)"],
    xlim=[-1,1], ylim=[-1, 1], relim=[-1,1], imlim=[-1,1], complex_only=True)

plotgroup("Derivative of Jacobi theta functions 1-2, variable nome",
    [lambda q: djtheta(1,(1+j)/3,q), lambda q: djtheta(2,(1+j)/3,q)],
    ["&theta;<sub>1</sub>((1+i)/3,q)", "&theta;<sub>2</sub>((1+i)/3,q)"],
    xlim=[-1,1], ylim=[-1, 1], relim=[-1,1], imlim=[-1,1], complex_only=True)

plotgroup("Derivative of Jacobi theta functions 3-4, variable nome",
    [lambda q: djtheta(3,(1+j)/3,q), lambda q: djtheta(4,(1+j)/3,q)],
    ["&theta;'<sub>3</sub>((1+i)/3,q)", "&theta;'<sub>4</sub>((1+i)/3,q)"],
    xlim=[-1,1], ylim=[-1, 1], relim=[-1,1], imlim=[-1,1], complex_only=True)


fp.write("</body></html>")
fp.close()

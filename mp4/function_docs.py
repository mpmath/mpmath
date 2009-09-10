"""
Extended docstrings for functions.py
"""


pi = r"""
`\pi`, roughly equal to 3.141592654, represents the area of the unit
circle, the half-period of trigonometric functions, and many other
things in mathematics.

Mpmath can evaluate `\pi` to arbitrary precision::

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +pi
    3.1415926535897932384626433832795028841971693993751

This shows digits 99991-100000 of `\pi`::

    >>> mp.dps = 100000
    >>> str(pi)[-10:]
    '5549362464'

**Possible issues**

:data:`pi` always rounds to the nearest floating-point
number when used. This means that exact mathematical identities
involving `\pi` will generally not be preserved in floating-point
arithmetic. In particular, multiples of :data:`pi` (except for
the trivial case ``0*pi``) are *not* the exact roots of
:func:`sin`, but differ roughly by the current epsilon::

    >>> mp.dps = 15
    >>> sin(pi)
    1.22464679914735e-16

One solution is to use the :func:`sinpi` function instead::

    >>> sinpi(1)
    0.0

See the documentation of trigonometric functions for additional
details.
"""

degree = r"""
Represents one degree of angle, `1^{\circ} = \pi/180`, or
about 0.01745329. This constant may be evaluated to arbitrary
precision::

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +degree
    0.017453292519943295769236907684886127134428718885417

The :data:`degree` object is convenient for conversion
to radians::

    >>> sin(30 * degree)
    0.5
    >>> asin(0.5) / degree
    30.0
"""

e = r"""
The transcendental number `e` = 2.718281828... is the base of the
natural logarithm (:func:`ln`) and of the exponential function
(:func:`exp`).

Mpmath can be evaluate `e` to arbitrary precision::

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +e
    2.7182818284590452353602874713526624977572470937

This shows digits 99991-100000 of `e`::

    >>> mp.dps = 100000
    >>> str(e)[-10:]
    '2100427165'

**Possible issues**

:data:`e` always rounds to the nearest floating-point number
when used, and mathematical identities involving `e` may not
hold in floating-point arithmetic. For example, ``ln(e)``
might not evaluate exactly to 1.

In particular, don't use ``e**x`` to compute the exponential
function. Use ``exp(x)`` instead; this is both faster and more
accurate.
"""

phi = r"""
Represents the golden ratio `\phi = (1+\sqrt 5)/2`,
approximately equal to 1.6180339887. To high precision,
its value is::

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +phi
    1.6180339887498948482045868343656381177203091798058

Formulas for the golden ratio include the following::

    >>> (1+sqrt(5))/2
    1.6180339887498948482045868343656381177203091798058
    >>> findroot(lambda x: x**2-x-1, 1)
    1.6180339887498948482045868343656381177203091798058
    >>> limit(lambda n: fib(n+1)/fib(n), inf)
    1.6180339887498948482045868343656381177203091798058
"""

euler = r"""
Euler's constant or the Euler-Mascheroni constant `\gamma`
= 0.57721566... is a number of central importance to
number theory and special functions. It is defined as the limit

.. math ::

    \gamma = \lim_{n\to\infty} H_n - \log n

where `H_n = 1 + \frac{1}{2} + \ldots + \frac{1}{n}` is a harmonic
number (see :func:`harmonic`).

Evaluation of `\gamma` is supported at arbitrary precision::

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +euler
    0.57721566490153286060651209008240243104215933593992

We can also compute `\gamma` directly from the definition,
although this is less efficient::

    >>> limit(lambda n: harmonic(n)-log(n), inf)
    0.57721566490153286060651209008240243104215933593992

This shows digits 9991-10000 of `\gamma`::

    >>> mp.dps = 10000
    >>> str(euler)[-10:]
    '4679858165'

Integrals, series, and representations for `\gamma` in terms of
special functions include the following (there are many others)::

    >>> mp.dps = 25
    >>> -quad(lambda x: exp(-x)*log(x), [0,inf])
    0.5772156649015328606065121
    >>> quad(lambda x,y: (x-1)/(1-x*y)/log(x*y), [0,1], [0,1])
    0.5772156649015328606065121
    >>> nsum(lambda k: 1/k-log(1+1/k), [1,inf])
    0.5772156649015328606065121
    >>> nsum(lambda k: (-1)**k*zeta(k)/k, [2,inf])
    0.5772156649015328606065121
    >>> -diff(gamma, 1)
    0.5772156649015328606065121
    >>> limit(lambda x: 1/x-gamma(x), 0)
    0.5772156649015328606065121
    >>> limit(lambda x: zeta(x)-1/(x-1), 1)
    0.5772156649015328606065121
    >>> (log(2*pi*nprod(lambda n:
    ...     exp(-2+2/n)*(1+2/n)**n, [1,inf]))-3)/2
    0.5772156649015328606065121

For generalizations of the identities `\gamma = -\Gamma'(1)`
and `\gamma = \lim_{x\to1} \zeta(x)-1/(x-1)`, see
:func:`psi` and :func:`stieltjes` respectively.
"""

catalan = r"""
Catalan's constant `K` = 0.91596559... is given by the infinite
series

.. math ::

    K = \sum_{k=0}^{\infty} \frac{(-1)^k}{(2k+1)^2}.

Mpmath can evaluate it to arbitrary precision::

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +catalan
    0.91596559417721901505460351493238411077414937428167

One can also compute `K` directly from the definition, although
this is significantly less efficient::

    >>> nsum(lambda k: (-1)**k/(2*k+1)**2, [0, inf])
    0.91596559417721901505460351493238411077414937428167

This shows digits 9991-10000 of `K`::

    >>> mp.dps = 10000
    >>> str(catalan)[-10:]
    '9537871503'

Catalan's constant has numerous integral representations::

    >>> mp.dps = 50
    >>> quad(lambda x: -log(x)/(1+x**2), [0, 1])
    0.91596559417721901505460351493238411077414937428167
    >>> quad(lambda x: atan(x)/x, [0, 1])
    0.91596559417721901505460351493238411077414937428167
    >>> quad(lambda x: ellipk(x**2)/2, [0, 1])
    0.91596559417721901505460351493238411077414937428167
    >>> quad(lambda x,y: 1/(1+(x*y)**2), [0, 1], [0, 1])
    0.91596559417721901505460351493238411077414937428167

As well as series representations::

    >>> pi*log(sqrt(3)+2)/8 + 3*nsum(lambda n:
    ...  (fac(n)/(2*n+1))**2/fac(2*n), [0, inf])/8
    0.91596559417721901505460351493238411077414937428167
    >>> 1-nsum(lambda n: n*zeta(2*n+1)/16**n, [1,inf])
    0.91596559417721901505460351493238411077414937428167
"""

khinchin = r"""
Khinchin's constant `K` = 2.68542... is a number that
appears in the theory of continued fractions. Mpmath can evaluate
it to arbitrary precision::

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +khinchin
    2.6854520010653064453097148354817956938203822939945

An integral representation is::

    >>> I = quad(lambda x: log((1-x**2)/sincpi(x))/x/(1+x), [0, 1])
    >>> 2*exp(1/log(2)*I)
    2.6854520010653064453097148354817956938203822939945

The computation of ``khinchin`` is based on an efficient
implementation of the following series::

    >>> f = lambda n: (zeta(2*n)-1)/n*sum((-1)**(k+1)/mpf(k)
    ...     for k in range(1,2*n))
    >>> exp(nsum(f, [1,inf])/log(2))
    2.6854520010653064453097148354817956938203822939945
"""

glaisher = r"""
Glaisher's constant `A`, also known as the Glaisher-Kinkelin
constant, is a number approximately equal to 1.282427129 that
sometimes appears in formulas related to gamma and zeta functions.
It is also related to the Barnes G-function (see :func:`barnesg`).

The constant is defined  as `A = \exp(1/12-\zeta'(-1))` where
`\zeta'(s)` denotes the derivative of the Riemann zeta function
(see :func:`zeta`).

Mpmath can evaluate Glaisher's constant to arbitrary precision:

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +glaisher
    1.282427129100622636875342568869791727767688927325

We can verify that the value computed by :data:`glaisher` is
correct using mpmath's facilities for numerical
differentiation and arbitrary evaluation of the zeta function:

    >>> exp(mpf(1)/12 - diff(zeta, -1))
    1.282427129100622636875342568869791727767688927325

Here is an example of an integral that can be evaluated in
terms of Glaisher's constant:

    >>> mp.dps = 15
    >>> quad(lambda x: log(gamma(x)), [1, 1.5])
    -0.0428537406502909
    >>> -0.5 - 7*log(2)/24 + log(pi)/4 + 3*log(glaisher)/2
    -0.042853740650291

Mpmath computes Glaisher's constant by applying Euler-Maclaurin
summation to a slowly convergent series. The implementation is
reasonably efficient up to about 10,000 digits. See the source
code for additional details.

References:
http://mathworld.wolfram.com/Glaisher-KinkelinConstant.html
"""

apery = r"""
Represents Apery's constant, which is the irrational number
approximately equal to 1.2020569 given by

.. math ::

    \zeta(3) = \sum_{k=1}^\infty\frac{1}{k^3}.

The calculation is based on an efficient hypergeometric
series. To 50 decimal places, the value is given by::

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +apery
    1.2020569031595942853997381615114499907649862923405

Other ways to evaluate Apery's constant using mpmath
include::

    >>> zeta(3)
    1.2020569031595942853997381615114499907649862923405
    >>> -diff(trigamma, 1)/2
    1.2020569031595942853997381615114499907649862923405
    >>> 8*nsum(lambda k: 1/(2*k+1)**3, [0,inf])/7
    1.2020569031595942853997381615114499907649862923405
    >>> f = lambda k: 2/k**3/(exp(2*pi*k)-1)
    >>> 7*pi**3/180 - nsum(f, [1,inf])
    1.2020569031595942853997381615114499907649862923405

This shows digits 9991-10000 of Apery's constant::

    >>> mp.dps = 10000
    >>> str(apery)[-10:]
    '3189504235'
"""

mertens = r"""
Represents the Mertens or Meissel-Mertens constant, which is the
prime number analog of Euler's constant:

.. math ::

    B_1 = \lim_{N\to\infty}
        \left(\sum_{p_k \le N} \frac{1}{p_k} - \log \log N \right)

Here `p_k` denotes the `k`-th prime number. Other names for this
constant include the Hadamard-de la Vallee-Poussin constant or
the prime reciprocal constant.

The following gives the Mertens constant to 50 digits::

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +mertens
    0.2614972128476427837554268386086958590515666482612

References:
http://mathworld.wolfram.com/MertensConstant.html
"""

twinprime = r"""
Represents the twin prime constant, which is the factor `C_2`
featuring in the Hardy-Littlewood conjecture for the growth of the
twin prime counting function,

.. math ::

    \pi_2(n) \sim 2 C_2 \frac{n}{\log^2 n}.

It is given by the product over primes

.. math ::

    C_2 = \prod_{p\ge3} \frac{p(p-2)}{(p-1)^2} \approx 0.66016

Computing `C_2` to 50 digits::

    >>> from mpmath import *
    >>> mp.dps = 50; mp.pretty = True
    >>> +twinprime
    0.66016181584686957392781211001455577843262336028473

References:
http://mathworld.wolfram.com/TwinPrimesConstant.html
"""

ln = r"""
Computes the natural logarithm of `x`, `\ln x`.
See :func:`log` for additional documentation."""

sqrt = r"""
``sqrt(x)`` gives the principal square root of `x`, `\sqrt x`.
For positive real numbers, the principal root is simply the
positive square root. For arbitrary complex numbers, the principal
square root is defined to satisfy `\sqrt x = \exp(\log(x)/2)`.
The function thus has a branch cut along the negative half real axis.

For all mpmath numbers ``x``, calling ``sqrt(x)`` is equivalent to
performing ``x**0.5``.

**Examples**

Basic examples and limits::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> sqrt(10)
    3.16227766016838
    >>> sqrt(100)
    10.0
    >>> sqrt(-4)
    (0.0 + 2.0j)
    >>> sqrt(1+1j)
    (1.09868411346781 + 0.455089860562227j)
    >>> sqrt(inf)
    +inf

Square root evaluation is fast at huge precision::

    >>> mp.dps = 50000
    >>> a = sqrt(3)
    >>> str(a)[-10:]
    '9329332814'

:func:`sqrt` supports interval arguments::

    >>> mp.dps = 15
    >>> sqrt(mpi(16, 100))
    [4.0, 10.0]
    >>> sqrt(mpi(2))
    [1.4142135623730949234, 1.4142135623730951455]
    >>> sqrt(mpi(2)) ** 2
    [1.9999999999999995559, 2.0000000000000004441]

"""

cbrt = r"""
``cbrt(x)`` computes the cube root of `x`, `x^{1/3}`. This
function is faster and more accurate than raising to a floating-point
fraction::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = False
    >>> 125**(mpf(1)/3)
    mpf('4.9999999999999991')
    >>> cbrt(125)
    mpf('5.0')

Every nonzero complex number has three cube roots. This function
returns the cube root defined by `\exp(\log(x)/3)` where the
principal branch of the natural logarithm is used. Note that this
does not give a real cube root for negative real numbers::

    >>> mp.pretty = True
    >>> cbrt(-1)
    (0.5 + 0.866025403784439j)
"""

exp = r"""
Computes the exponential function,

.. math ::

    \exp(x) = e^x = \sum_{k=0}^{\infty} \frac{x^k}{k!}.

For complex numbers, the exponential function also satisfies

.. math ::

    \exp(x+yi) = e^x (\cos y + i \sin y).

**Basic examples**

Some values of the exponential function::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> exp(0)
    1.0
    >>> exp(1)
    2.718281828459045235360287
    >>> exp(-1)
    0.3678794411714423215955238
    >>> exp(inf)
    +inf
    >>> exp(-inf)
    0.0

Arguments can be arbitrarily large::

    >>> exp(10000)
    8.806818225662921587261496e+4342
    >>> exp(-10000)
    1.135483865314736098540939e-4343

Evaluation is supported for interval arguments::

    >>> exp(mpi(-inf,0))
    [0.0, 1.0]
    >>> exp(mpi(0,1))
    [1.0, 2.71828182845904523536028749558]

The exponential function can be evaluated efficiently to arbitrary
precision::

    >>> mp.dps = 10000
    >>> exp(pi)  #doctest: +ELLIPSIS
    23.140692632779269005729...8984304016040616

**Functional properties**

Numerical verification of Euler's identity for the complex
exponential function::

    >>> mp.dps = 15
    >>> exp(j*pi)+1
    (0.0 + 1.22464679914735e-16j)
    >>> chop(exp(j*pi)+1)
    0.0

This recovers the coefficients (reciprocal factorials) in the
Maclaurin series expansion of exp::

    >>> nprint(taylor(exp, 0, 5))
    [1.0, 1.0, 0.5, 0.166667, 4.16667e-2, 8.33333e-3]

The exponential function is its own derivative and antiderivative::

    >>> exp(pi)
    23.1406926327793
    >>> diff(exp, pi)
    23.1406926327793
    >>> quad(exp, [-inf, pi])
    23.1406926327793

The exponential function can be evaluated using various methods,
including direct summation of the series, limits, and solving
the defining differential equation::

    >>> nsum(lambda k: pi**k/fac(k), [0,inf])
    23.1406926327793
    >>> limit(lambda k: (1+pi/k)**k, inf)
    23.1406926327793
    >>> odefun(lambda t, x: x, 0, 1)(pi)
    23.1406926327793
"""

cosh = r"""
Computes the hyperbolic cosine of `x`,
`\cosh(x) = (e^x + e^{-x})/2`. Values and limits include::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> cosh(0)
    1.0
    >>> cosh(1)
    1.543080634815243778477906
    >>> cosh(-inf), cosh(+inf)
    (+inf, +inf)

The hyperbolic cosine is an even, convex function with
a global minimum at `x = 0`, having a Maclaurin series
that starts::

    >>> nprint(chop(taylor(cosh, 0, 5)))
    [1.0, 0.0, 0.5, 0.0, 4.16667e-2, 0.0]

Generalized to complex numbers, the hyperbolic cosine is
equivalent to a cosine with the argument rotated
in the imaginary direction, or `\cosh x = \cos ix`::

    >>> cosh(2+3j)
    (-3.724545504915322565473971 + 0.5118225699873846088344638j)
    >>> cos(3-2j)
    (-3.724545504915322565473971 + 0.5118225699873846088344638j)
"""

sinh = r"""
Computes the hyperbolic sine of `x`,
`\sinh(x) = (e^x - e^{-x})/2`. Values and limits include::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> sinh(0)
    0.0
    >>> sinh(1)
    1.175201193643801456882382
    >>> sinh(-inf), sinh(+inf)
    (-inf, +inf)

The hyperbolic sine is an odd function, with a Maclaurin
series that starts::

    >>> nprint(chop(taylor(sinh, 0, 5)))
    [0.0, 1.0, 0.0, 0.166667, 0.0, 8.33333e-3]

Generalized to complex numbers, the hyperbolic sine is
essentially a sine with a rotation `i` applied to
the argument; more precisely, `\sinh x = -i \sin ix`::

    >>> sinh(2+3j)
    (-3.590564589985779952012565 + 0.5309210862485198052670401j)
    >>> j*sin(3-2j)
    (-3.590564589985779952012565 + 0.5309210862485198052670401j)
"""

tanh = r"""
Computes the hyperbolic tangent of `x`,
`\tanh(x) = \sinh(x)/\cosh(x)`. Values and limits include::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> tanh(0)
    0.0
    >>> tanh(1)
    0.7615941559557648881194583
    >>> tanh(-inf), tanh(inf)
    (-1.0, 1.0)

The hyperbolic tangent is an odd, sigmoidal function, similar
to the inverse tangent and error function. Its Maclaurin
series is::

    >>> nprint(chop(taylor(tanh, 0, 5)))
    [0.0, 1.0, 0.0, -0.333333, 0.0, 0.133333]

Generalized to complex numbers, the hyperbolic tangent is
essentially a tangent with a rotation `i` applied to
the argument; more precisely, `\tanh x = -i \tan ix`::

    >>> tanh(2+3j)
    (0.9653858790221331242784803 - 0.009884375038322493720314034j)
    >>> j*tan(3-2j)
    (0.9653858790221331242784803 - 0.009884375038322493720314034j)
"""

cos = r"""
Computes the cosine of `x`, `\cos(x)`.

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> cos(pi/3)
    0.5
    >>> cos(100000001)
    -0.9802850113244713353133243
    >>> cos(2+3j)
    (-4.189625690968807230132555 - 9.109227893755336597979197j)
    >>> cos(inf)
    nan
    >>> nprint(chop(taylor(cos, 0, 6)))
    [1.0, 0.0, -0.5, 0.0, 4.16667e-2, 0.0, -1.38889e-3]
    >>> cos(mpi(0,1))
    [0.540302305868139717400936602301, 1.0]
    >>> cos(mpi(0,2))
    [-0.41614683654714238699756823214, 1.0]
"""

sin = r"""
Computes the sine of `x`, `\sin(x)`.

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> sin(pi/3)
    0.8660254037844386467637232
    >>> sin(100000001)
    0.1975887055794968911438743
    >>> sin(2+3j)
    (9.1544991469114295734673 - 4.168906959966564350754813j)
    >>> sin(inf)
    nan
    >>> nprint(chop(taylor(sin, 0, 6)))
    [0.0, 1.0, 0.0, -0.166667, 0.0, 8.33333e-3, 0.0]
    >>> sin(mpi(0,1))
    [0.0, 0.841470984807896506652502331201]
    >>> sin(mpi(0,2))
    [0.0, 1.0]
"""

tan = r"""
Computes the tangent of `x`, `\tan(x) = \frac{\sin(x)}{\cos(x)}`.
The tangent function is singular at `x = (n+1/2)\pi`, but
``tan(x)`` always returns a finite result since `(n+1/2)\pi`
cannot be represented exactly using floating-point arithmetic.

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> tan(pi/3)
    1.732050807568877293527446
    >>> tan(100000001)
    -0.2015625081449864533091058
    >>> tan(2+3j)
    (-0.003764025641504248292751221 + 1.003238627353609801446359j)
    >>> tan(inf)
    nan
    >>> nprint(chop(taylor(tan, 0, 6)))
    [0.0, 1.0, 0.0, 0.333333, 0.0, 0.133333, 0.0]
    >>> tan(mpi(0,1))
    [0.0, 1.55740772465490223050697482944]
    >>> tan(mpi(0,2))  # Interval includes a singularity
    [-inf, +inf]
"""

sec = r"""
Computes the secant of `x`, `\mathrm{sec}(x) = \frac{1}{\cos(x)}`.
The secant function is singular at `x = (n+1/2)\pi`, but
``sec(x)`` always returns a finite result since `(n+1/2)\pi`
cannot be represented exactly using floating-point arithmetic.

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> sec(pi/3)
    2.0
    >>> sec(10000001)
    -1.184723164360392819100265
    >>> sec(2+3j)
    (-0.04167496441114427004834991 + 0.0906111371962375965296612j)
    >>> sec(inf)
    nan
    >>> nprint(chop(taylor(sec, 0, 6)))
    [1.0, 0.0, 0.5, 0.0, 0.208333, 0.0, 8.47222e-2]
    >>> sec(mpi(0,1))
    [1.0, 1.85081571768092561791175324143]
    >>> sec(mpi(0,2))  # Interval includes a singularity
    [-inf, +inf]
"""

csc = r"""
Computes the cosecant of `x`, `\mathrm{csc}(x) = \frac{1}{\sin(x)}`.
This cosecant function is singular at `x = n \pi`, but with the
exception of the point `x = 0`, ``csc(x)`` returns a finite result
since `n \pi` cannot be represented exactly using floating-point
arithmetic.

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> csc(pi/3)
    1.154700538379251529018298
    >>> csc(10000001)
    -1.864910497503629858938891
    >>> csc(2+3j)
    (0.09047320975320743980579048 + 0.04120098628857412646300981j)
    >>> csc(inf)
    nan
    >>> csc(mpi(0,1))  # Interval includes a singularity
    [1.18839510577812121626159945235, +inf]
    >>> csc(mpi(0,2))
    [1.0, +inf]
"""

cot = r"""
Computes the cotangent of `x`,
`\mathrm{cot}(x) = \frac{1}{\tan(x)} = \frac{\cos(x)}{\sin(x)}`.
This cotangent function is singular at `x = n \pi`, but with the
exception of the point `x = 0`, ``cot(x)`` returns a finite result
since `n \pi` cannot be represented exactly using floating-point
arithmetic.

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> cot(pi/3)
    0.5773502691896257645091488
    >>> cot(10000001)
    1.574131876209625656003562
    >>> cot(2+3j)
    (-0.003739710376336956660117409 - 0.9967577965693583104609688j)
    >>> cot(inf)
    nan
    >>> cot(mpi(0,1))  # Interval includes a singularity
    [0.642092615934330703006419986575, +inf]
    >>> cot(mpi(1,2))
    [-inf, +inf]
"""

acos = r"""
Computes the inverse cosine or arccosine of `x`, `\cos^{-1}(x)`.
Since `-1 \le \cos(x) \le 1` for real `x`, the inverse
cosine is real-valued only for `-1 \le x \le 1`. On this interval,
:func:`acos` is defined to be a monotonically decreasing
function assuming values between `+\pi` and `0`.

Basic values are::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> acos(-1)
    3.141592653589793238462643
    >>> acos(0)
    1.570796326794896619231322
    >>> acos(1)
    0.0
    >>> nprint(chop(taylor(acos, 0, 6)))
    [1.5708, -1.0, 0.0, -0.166667, 0.0, -7.5e-2, 0.0]

:func:`acos` is defined so as to be a proper inverse function of
`\cos(\theta)` for `0 \le \theta < \pi`.
We have `\cos(\cos^{-1}(x)) = x` for all `x`, but
`\cos^{-1}(\cos(x)) = x` only for `0 \le \Re[x] < \pi`::

    >>> for x in [1, 10, -1, 2+3j, 10+3j]:
    ...     print cos(acos(x)), acos(cos(x))
    ...
    1.0 1.0
    (10.0 + 0.0j) 2.566370614359172953850574
    -1.0 1.0
    (2.0 + 3.0j) (2.0 + 3.0j)
    (10.0 + 3.0j) (2.566370614359172953850574 - 3.0j)

The inverse cosine has two branch points: `x = \pm 1`. :func:`acos`
places the branch cuts along the line segments `(-\infty, -1)` and
`(+1, +\infty)`. In general,

.. math ::

    \cos^{-1}(x) = \frac{\pi}{2} + i \log\left(ix + \sqrt{1-x^2} \right)

where the principal-branch log and square root are implied.
"""

asin = r"""
Computes the inverse sine or arcsine of `x`, `\sin^{-1}(x)`.
Since `-1 \le \sin(x) \le 1` for real `x`, the inverse
sine is real-valued only for `-1 \le x \le 1`.
On this interval, it is defined to be a monotonically increasing
function assuming values between `-\pi/2` and `\pi/2`.

Basic values are::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> asin(-1)
    -1.570796326794896619231322
    >>> asin(0)
    0.0
    >>> asin(1)
    1.570796326794896619231322
    >>> nprint(chop(taylor(asin, 0, 6)))
    [0.0, 1.0, 0.0, 0.166667, 0.0, 7.5e-2, 0.0]

:func:`asin` is defined so as to be a proper inverse function of
`\sin(\theta)` for `-\pi/2 < \theta < \pi/2`.
We have `\sin(\sin^{-1}(x)) = x` for all `x`, but
`\sin^{-1}(\sin(x)) = x` only for `-\pi/2 < \Re[x] < \pi/2`::

    >>> for x in [1, 10, -1, 1+3j, -2+3j]:
    ...     print chop(sin(asin(x))), asin(sin(x))
    ...
    1.0 1.0
    10.0 -0.5752220392306202846120698
    -1.0 -1.0
    (1.0 + 3.0j) (1.0 + 3.0j)
    (-2.0 + 3.0j) (-1.141592653589793238462643 - 3.0j)

The inverse sine has two branch points: `x = \pm 1`. :func:`asin`
places the branch cuts along the line segments `(-\infty, -1)` and
`(+1, +\infty)`. In general,

.. math ::

    \sin^{-1}(x) = -i \log\left(ix + \sqrt{1-x^2} \right)

where the principal-branch log and square root are implied.
"""

atan = r"""
Computes the inverse tangent or arctangent of `x`, `\tan^{-1}(x)`.
This is a real-valued function for all real `x`, with range
`(-\pi/2, \pi/2)`.

Basic values are::

    >>> from mpmath import *
    >>> mp.dps = 25
    >>> atan(-inf); mp.pretty = True
    -1.570796326794896619231322
    >>> atan(-1)
    -0.7853981633974483096156609
    >>> atan(0)
    0.0
    >>> atan(1)
    0.7853981633974483096156609
    >>> atan(inf)
    1.570796326794896619231322
    >>> nprint(chop(taylor(atan, 0, 6)))
    [0.0, 1.0, 0.0, -0.333333, 0.0, 0.2, 0.0]

The inverse tangent is often used to compute angles. However,
the atan2 function is often better for this as it preserves sign
(see :func:`atan2`).

:func:`atan` is defined so as to be a proper inverse function of
`\tan(\theta)` for `-\pi/2 < \theta < \pi/2`.
We have `\tan(\tan^{-1}(x)) = x` for all `x`, but
`\tan^{-1}(\tan(x)) = x` only for `-\pi/2 < \Re[x] < \pi/2`::

    >>> mp.dps = 25
    >>> for x in [1, 10, -1, 1+3j, -2+3j]:
    ...     print tan(atan(x)), atan(tan(x))
    ...
    1.0 1.0
    10.0 0.5752220392306202846120698
    -1.0 -1.0
    (1.0 + 3.0j) (1.000000000000000000000001 + 3.0j)
    (-2.0 + 3.0j) (1.141592653589793238462644 + 3.0j)

The inverse tangent has two branch points: `x = \pm i`. :func:`atan`
places the branch cuts along the line segments `(-i \infty, -i)` and
`(+i, +i \infty)`. In general,

.. math ::

    \tan^{-1}(x) = \frac{i}{2}\left(\log(1-ix)-\log(1+ix)\right)

where the principal-branch log is implied.
"""

acot = r"""Computes the inverse cotangent of `x`,
`\mathrm{cot}^{-1}(x) = \tan^{-1}(1/x)`."""

asec = r"""Computes the inverse secant of `x`,
`\mathrm{sec}^{-1}(x) = \cos^{-1}(1/x)`."""

acsc_ = r"""Computes the inverse cosecant of `x`,
`\mathrm{csc}^{-1}(x) = \sin^{-1}(1/x)`."""

sinpi = r"""
Computes `\sin(\pi x)`, more accurately than the expression
``sin(pi*x)``::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> sinpi(10**10), sin(pi*(10**10))
    (0.0, -2.23936276195592e-6)
    >>> sinpi(10**10+0.5), sin(pi*(10**10+0.5))
    (1.0, 0.999999999998721)
"""

cospi = r"""
Computes `\cos(\pi x)`, more accurately than the expression
``cos(pi*x)``::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> cospi(10**10), cos(pi*(10**10))
    (1.0, 0.999999999997493)
    >>> cospi(10**10+0.5), cos(pi*(10**10+0.5))
    (0.0, 1.59960492420134e-6)
"""

sinc = r"""
``sinc(x)`` computes the unnormalized sinc function, defined as

.. math ::

    \mathrm{sinc}(x) = \begin{cases}
        \sin(x)/x, & \mbox{if } x \ne 0 \\
        1,         & \mbox{if } x = 0.
    \end{cases}

See :func:`sincpi` for the normalized sinc function.

Simple values and limits include::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> sinc(0)
    1.0
    >>> sinc(1)
    0.841470984807897
    >>> sinc(inf)
    0.0

The integral of the sinc function is the sine integral Si::

    >>> quad(sinc, [0, 1])
    0.946083070367183
    >>> si(1)
    0.946083070367183
"""

sincpi = r"""
``sincpi(x)`` computes the normalized sinc function, defined as

.. math ::

    \mathrm{sinc}_{\pi}(x) = \begin{cases}
        \sin(\pi x)/(\pi x), & \mbox{if } x \ne 0 \\
        1,                   & \mbox{if } x = 0.
    \end{cases}

Equivalently, we have
`\mathrm{sinc}_{\pi}(x) = \mathrm{sinc}(\pi x)`.

The normalization entails that the function integrates
to unity over the entire real line::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> quadosc(sincpi, [-inf, inf], period=2.0)
    1.0

Like, :func:`sinpi`, :func:`sincpi` is evaluated accurately
at its roots::

    >>> sincpi(10)
    0.0
"""

floor = r"""
Computes the floor of `x`, `\lfloor x \rfloor`, defined as
the largest integer less than or equal to `x`::

    >>> from mpmath import *
    >>> mp.pretty = False
    >>> floor(3.5)
    mpf('3.0')

Note: :func:`floor` returns a floating-point number, not a
Python ``int``. If `\lfloor x \rfloor` is too large to be
represented exactly at the present working precision, the
result will be rounded, not necessarily in the floor
direction."""

ceil = r"""
Computes the ceiling of `x`, `\lceil x \rceil`, defined as
the smallest integer greater than or equal to `x`::

    >>> from mpmath import *
    >>> mp.pretty = False
    >>> ceil(3.5)
    mpf('4.0')

Note: :func:`ceil` returns a floating-point number, not a
Python ``int``. If `\lceil x \rceil` is too large to be
represented exactly at the present working precision, the
result will be rounded, not necessarily in the ceiling
direction."""


expm1 = r"""
Computes `e^x - 1`, accurately for small `x`.

Unlike the expression ``exp(x) - 1``, ``expm1(x)`` does not suffer from
potentially catastrophic cancellation::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> exp(1e-10)-1; print expm1(1e-10)
    1.00000008274037e-10
    1.00000000005e-10
    >>> exp(1e-20)-1; print expm1(1e-20)
    0.0
    1.0e-20
    >>> 1/(exp(1e-20)-1)
    Traceback (most recent call last):
      ...
    ZeroDivisionError
    >>> 1/expm1(1e-20)
    1.0e+20

Evaluation works for extremely tiny values::

    >>> expm1(0)
    0.0
    >>> expm1('1e-10000000')
    1.0e-10000000

"""

powm1 = r"""
Computes `x^y - 1`, accurately when `x^y` is very close to 1.

This avoids potentially catastrophic cancellation::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> power(0.99999995, 1e-10) - 1
    0.0
    >>> powm1(0.99999995, 1e-10)
    -5.00000012791934e-18

Powers exactly equal to 1, and only those powers, yield 0 exactly::

    >>> powm1(-j, 4)
    (0.0 + 0.0j)
    >>> powm1(3, 0)
    0.0
    >>> powm1(fadd(-1, 1e-100, exact=True), 4)
    -4.0e-100

Evaluation works for extremely tiny `y`::

    >>> powm1(2, '1e-100000')
    6.93147180559945e-100001
    >>> powm1(j, '1e-1000')
    (-1.23370055013617e-2000 + 1.5707963267949e-1000j)

"""

root = r"""
``root(z, n, k=0)`` computes an `n`-th root of `z`, i.e. returns a number
`r` that (up to possible approximation error) satisfies `r^n = z`.
(``nthroot`` is available as an alias for ``root``.)

Every complex number `z \ne 0` has `n` distinct `n`-th roots, which are
equidistant points on a circle with radius `|z|^{1/n}`, centered around the
origin. A specific root may be selected using the optional index
`k`. The roots are indexed counterclockwise, starting with `k = 0` for the root
closest to the positive real half-axis.

The `k = 0` root is the so-called principal `n`-th root, often denoted by
`\sqrt[n]{z}` or `z^{1/n}`, and also given by `\exp(\log(z) / n)`. If `z` is
a positive real number, the principal root is just the unique positive
`n`-th root of `z`. Under some circumstances, non-principal real roots exist:
for positive real `z`, `n` even, there is a negative root given by `k = n/2`;
for negative real `z`, `n` odd, there is a negative root given by `k = (n-1)/2`.

To obtain all roots with a simple expression, use
``[root(z,n,k) for k in range(n)]``.

An important special case, ``root(1, n, k)`` returns the `k`-th `n`-th root of
unity, `\zeta_k = e^{2 \pi i k / n}`. Alternatively, :func:`unitroots`
provides a slightly more convenient way to obtain the roots of unity,
including the option to compute only the primitive roots of unity.

Both `k` and `n` should be integers; `k` outside of ``range(n)`` will be
reduced modulo `n`. If `n` is negative, `x^{-1/n} = 1/x^{1/n}` (or
the equivalent reciprocal for a non-principal root with `k \ne 0`) is computed.

:func:`root` is implemented to use Newton's method for small
`n`. At high precision, this makes `x^{1/n}` not much more
expensive than the regular exponentiation, `x^n`. For very large
`n`, :func:`nthroot` falls back to use the exponential function.

**Examples**

:func:`nthroot`/:func:`root` is faster and more accurate than raising to a
floating-point fraction::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = False
    >>> 16807 ** (mpf(1)/5)
    mpf('7.0000000000000009')
    >>> root(16807, 5)
    mpf('7.0')
    >>> nthroot(16807, 5)    # Alias
    mpf('7.0')

A high-precision root::

    >>> mp.dps = 50; mp.pretty = True
    >>> nthroot(10, 5)
    1.584893192461113485202101373391507013269442133825
    >>> nthroot(10, 5) ** 5
    10.0

Computing principal and non-principal square and cube roots::

    >>> mp.dps = 15
    >>> root(10, 2)
    3.16227766016838
    >>> root(10, 2, 1)
    -3.16227766016838
    >>> root(-10, 3)
    (1.07721734501594 + 1.86579517236206j)
    >>> root(-10, 3, 1)
    -2.15443469003188
    >>> root(-10, 3, 2)
    (1.07721734501594 - 1.86579517236206j)

All the 7th roots of a complex number::

    >>> for r in [root(3+4j, 7, k) for k in range(7)]:
    ...     print r, r**7
    ...
    (1.24747270589553 + 0.166227124177353j) (3.0 + 4.0j)
    (0.647824911301003 + 1.07895435170559j) (3.0 + 4.0j)
    (-0.439648254723098 + 1.17920694574172j) (3.0 + 4.0j)
    (-1.19605731775069 + 0.391492658196305j) (3.0 + 4.0j)
    (-1.05181082538903 - 0.691023585965793j) (3.0 + 4.0j)
    (-0.115529328478668 - 1.25318497558335j) (3.0 + 4.0j)
    (0.907748109144957 - 0.871672518271819j) (3.0 + 4.0j)

Cube roots of unity::

    >>> for k in range(3): print root(1, 3, k)
    ...
    1.0
    (-0.5 + 0.866025403784439j)
    (-0.5 - 0.866025403784439j)

Some exact high order roots::

    >>> root(75**210, 105)
    5625.0
    >>> root(1, 128, 96)
    (0.0 - 1.0j)
    >>> root(4**128, 128, 96)
    (0.0 - 4.0j)

"""

unitroots = r"""
``unitroots(n)`` returns `\zeta_0, \zeta_1, \ldots, \zeta_{n-1}`,
all the distinct `n`-th roots of unity, as a list. If the option
*primitive=True* is passed, only the primitive roots are returned.

Every `n`-th root of unity satisfies `(\zeta_k)^n = 1`. There are `n` distinct
roots for each `n` (`\zeta_k` and `\zeta_j` are the same when
`k = j \pmod n`), which form a regular polygon with vertices on the unit
circle. They are ordered counterclockwise with increasing `k`, starting
with `\zeta_0 = 1`.

**Examples**

The roots of unity up to `n = 4`::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> nprint(unitroots(1))
    [1.0]
    >>> nprint(unitroots(2))
    [1.0, -1.0]
    >>> nprint(unitroots(3))
    [1.0, (-0.5 + 0.866025j), (-0.5 - 0.866025j)]
    >>> nprint(unitroots(4))
    [1.0, (0.0 + 1.0j), -1.0, (0.0 - 1.0j)]

Roots of unity form a geometric series that sums to 0::

    >>> mp.dps = 50
    >>> chop(fsum(unitroots(25)))
    0.0

Primitive roots up to `n = 4`::

    >>> mp.dps = 15
    >>> nprint(unitroots(1, primitive=True))
    [1.0]
    >>> nprint(unitroots(2, primitive=True))
    [-1.0]
    >>> nprint(unitroots(3, primitive=True))
    [(-0.5 + 0.866025j), (-0.5 - 0.866025j)]
    >>> nprint(unitroots(4, primitive=True))
    [(0.0 + 1.0j), (0.0 - 1.0j)]

There are only four primitive 12th roots::

    >>> nprint(unitroots(12, primitive=True))
    [(0.866025 + 0.5j), (-0.866025 + 0.5j), (-0.866025 - 0.5j), (0.866025 - 0.5j)]

The `n`-th roots of unity form a group, the cyclic group of order `n`.
Any primitive root `r` is a generator for this group, meaning that
`r^0, r^1, \ldots, r^{n-1}` gives the whole set of unit roots (in
some permuted order)::

    >>> for r in unitroots(6): print r
    ...
    1.0
    (0.5 + 0.866025403784439j)
    (-0.5 + 0.866025403784439j)
    -1.0
    (-0.5 - 0.866025403784439j)
    (0.5 - 0.866025403784439j)
    >>> r = unitroots(6, primitive=True)[1]
    >>> for k in range(6): print chop(r**k)
    ...
    1.0
    (0.5 - 0.866025403784439j)
    (-0.5 - 0.866025403784439j)
    -1.0
    (-0.5 + 0.866025403784438j)
    (0.5 + 0.866025403784438j)

The number of primitive roots equals the Euler totient function `\phi(n)`::

    >>> [len(unitroots(n, primitive=True)) for n in range(1,20)]
    [1, 1, 2, 2, 4, 2, 6, 4, 6, 4, 10, 4, 12, 6, 8, 8, 16, 6, 18]

"""


log = r"""
Computes the base-`b` logarithm of `x`, `\log_b(x)`. If `b` is
unspecified, :func:`log` computes the natural (base `e`) logarithm
and is equivalent to :func:`ln`. In general, the base `b` logarithm
is defined in terms of the natural logarithm as
`\log_b(x) = \ln(x)/\ln(b)`.

By convention, we take `\log(0) = -\infty`.

The natural logarithm is real if `x > 0` and complex if `x < 0` or if
`x` is complex. The principal branch of the complex logarithm is
used, meaning that `\Im(\ln(x)) = -\pi < \arg(x) \le \pi`.

**Examples**

Some basic values and limits::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> log(1)
    0.0
    >>> log(2)
    0.693147180559945
    >>> log(1000,10)
    3.0
    >>> log(4, 16)
    0.5
    >>> log(j)
    (0.0 + 1.5707963267949j)
    >>> log(-1)
    (0.0 + 3.14159265358979j)
    >>> log(0)
    -inf
    >>> log(inf)
    +inf

The natural logarithm is the antiderivative of `1/x`::

    >>> quad(lambda x: 1/x, [1, 5])
    1.6094379124341
    >>> log(5)
    1.6094379124341
    >>> diff(log, 10)
    0.1

The Taylor series expansion of the natural logarithm around
`x = 1` has coefficients `(-1)^{n+1}/n`::

    >>> nprint(taylor(log, 1, 7))
    [0.0, 1.0, -0.5, 0.333333, -0.25, 0.2, -0.166667, 0.142857]

:func:`log` supports arbitrary precision evaluation::

    >>> mp.dps = 50
    >>> log(pi)
    1.1447298858494001741434273513530587116472948129153
    >>> log(pi, pi**3)
    0.33333333333333333333333333333333333333333333333333
    >>> mp.dps = 25
    >>> log(3+4j)
    (1.609437912434100374600759 + 0.9272952180016122324285125j)
"""

power = r"""
Converts `x` and `y` to mpmath numbers and evaluates
`x^y = \exp(y \log(x))`::

    >>> from mpmath import *
    >>> mp.dps = 30; mp.pretty = True
    >>> power(2, 0.5)
    1.41421356237309504880168872421

This shows the leading few digits of a large Mersenne prime
(performing the exact calculation ``2**43112609-1`` and
displaying the result in Python would be very slow)::

    >>> power(2, 43112609)-1
    3.16470269330255923143453723949e+12978188
"""

modf = r"""
Converts `x` and `y` to mpmath numbers and returns `x \mod y`.
For mpmath numbers, this is equivalent to ``x % y``.

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> modf(100, pi)
    2.61062773871641

You can use :func:`modf` to compute fractional parts of numbers::

    >>> modf(10.25, 1)
    0.25

"""

radians = r"""
Converts the degree angle `x` to radians::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> radians(60)
    1.0471975511966
"""

degrees = r"""
Converts the radian angle `x` to a degree angle::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> degrees(pi/3)
    60.0
"""

atan2 = r"""
Computes the two-argument arctangent, `\mathrm{atan2}(y, x)`,
giving the signed angle between the positive `x`-axis and the
point `(x, y)` in the 2D plane. This function is defined for
real `x` and `y` only.

The two-argument arctangent essentially computes
`\mathrm{atan}(y/x)`, but accounts for the signs of both
`x` and `y` to give the angle for the correct quadrant. The
following examples illustrate the difference::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> atan2(1,1), atan(1/1.)
    (0.785398163397448, 0.785398163397448)
    >>> atan2(1,-1), atan(1/-1.)
    (2.35619449019234, -0.785398163397448)
    >>> atan2(-1,1), atan(-1/1.)
    (-0.785398163397448, -0.785398163397448)
    >>> atan2(-1,-1), atan(-1/-1.)
    (-2.35619449019234, 0.785398163397448)

The angle convention is the same as that used for the complex
argument; see :func:`arg`.
"""

fibonacci = r"""
``fibonacci(n)`` computes the `n`-th Fibonacci number, `F(n)`. The
Fibonacci numbers are defined by the recurrence `F(n) = F(n-1) + F(n-2)`
with the initial values `F(0) = 0`, `F(1) = 1`. :func:`fibonacci`
extends this definition to arbitrary real and complex arguments
using the formula

.. math ::

  F(z) = \frac{\phi^z - \cos(\pi z) \phi^{-z}}{\sqrt 5}

where `\phi` is the golden ratio. :func:`fibonacci` also uses this
continuous formula to compute `F(n)` for extremely large `n`, where
calculating the exact integer would be wasteful.

For convenience, :func:`fib` is available as an alias for
:func:`fibonacci`.

**Basic examples**

Some small Fibonacci numbers are::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for i in range(10):
    ...     print fibonacci(i),
    ...
    0.0 1.0 1.0 2.0 3.0 5.0 8.0 13.0 21.0 34.0

    >>> fibonacci(50)
    12586269025.0

The recurrence for `F(n)` extends backwards to negative `n`::

    >>> for i in range(10):
    ...     print fibonacci(-i),
    ...
    0.0 1.0 -1.0 2.0 -3.0 5.0 -8.0 13.0 -21.0 34.0

Large Fibonacci numbers will be computed approximately unless
the precision is set high enough::

    >>> fib(200)
    2.8057117299251e+41
    >>> mp.dps = 45
    >>> fib(200)
    280571172992510140037611932413038677189525.0

:func:`fibonacci` can compute approximate Fibonacci numbers
of stupendous size::

    >>> mp.dps = 15
    >>> fibonacci(10**25)
    3.49052338550226e+2089876402499787337692720

**Real and complex arguments**

The extended Fibonacci function is an analytic function. The
property `F(z) = F(z-1) + F(z-2)` holds for arbitrary `z`::

    >>> mp.dps = 15
    >>> fib(pi)
    2.1170270579161
    >>> fib(pi-1) + fib(pi-2)
    2.1170270579161
    >>> fib(3+4j)
    (-5248.51130728372 - 14195.962288353j)
    >>> fib(2+4j) + fib(1+4j)
    (-5248.51130728372 - 14195.962288353j)

The Fibonacci function has infinitely many roots on the
negative half-real axis. The first root is at 0, the second is
close to -0.18, and then there are infinitely many roots that
asymptotically approach `-n+1/2`::

    >>> findroot(fib, -0.2)
    -0.183802359692956
    >>> findroot(fib, -2)
    -1.57077646820395
    >>> findroot(fib, -17)
    -16.4999999596115
    >>> findroot(fib, -24)
    -23.5000000000479

**Mathematical relationships**

For large `n`, `F(n+1)/F(n)` approaches the golden ratio::

    >>> mp.dps = 50
    >>> fibonacci(101)/fibonacci(100)
    1.6180339887498948482045868343656381177203127439638
    >>> +phi
    1.6180339887498948482045868343656381177203091798058

The sum of reciprocal Fibonacci numbers converges to an irrational
number for which no closed form expression is known::

    >>> mp.dps = 15
    >>> nsum(lambda n: 1/fib(n), [1, inf])
    3.35988566624318

Amazingly, however, the sum of odd-index reciprocal Fibonacci
numbers can be expressed in terms of a Jacobi theta function::

    >>> nsum(lambda n: 1/fib(2*n+1), [0, inf])
    1.82451515740692
    >>> sqrt(5)*jtheta(2,0,(3-sqrt(5))/2)**2/4
    1.82451515740692

Some related sums can be done in closed form::

    >>> nsum(lambda k: 1/(1+fib(2*k+1)), [0, inf])
    1.11803398874989
    >>> phi - 0.5
    1.11803398874989
    >>> f = lambda k:(-1)**(k+1) / sum(fib(n)**2 for n in range(1,k+1))
    >>> nsum(f, [1, inf])
    0.618033988749895
    >>> phi-1
    0.618033988749895

**References**

1. http://mathworld.wolfram.com/FibonacciNumber.html
"""

zeta = r"""
``zeta(s)`` computes the Riemann zeta function, `\zeta(s)`.
The Riemann zeta function is defined for `\Re(s) > 1` by

.. math ::

  \zeta(s) = 1+\frac{1}{2^s}+\frac{1}{3^s}+\frac{1}{4^s}+\ldots

and for `\Re(s) \le 1` by analytic continuation. It has a pole
at `s = 1`.

**Examples**

Some exact values of the zeta function are::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> zeta(2)
    1.64493406684823
    >>> pi**2 / 6
    1.64493406684823
    >>> zeta(0)
    -0.5
    >>> zeta(-1)
    -0.0833333333333333
    >>> zeta(-2)
    0.0

:func:`zeta` supports arbitrary precision evaluation and
complex arguments::

    >>> mp.dps = 50
    >>> zeta(pi)
    1.1762417383825827588721504519380520911697389900217
    >>> zeta(1+2j)  # doctest: +NORMALIZE_WHITESPACE
    (0.5981655697623817367034568491742186771747764868876 -
    0.35185474521784529049653859679690026505229177886045j)

The Riemann zeta function has so-called nontrivial zeros on
the critical line `s = 1/2 + it`::

    >>> mp.dps = 15
    >>> findroot(zeta, 0.5+14j)
    (0.5 + 14.1347251417347j)
    >>> findroot(zeta, 0.5+21j)
    (0.5 + 21.0220396387716j)
    >>> findroot(zeta, 0.5+25j)
    (0.5 + 25.0108575801457j)

For investigation of the zeta function zeros, the Riemann-Siegel
Z-function is often more convenient than working with the Riemann
zeta function directly (see :func:`siegelz`).

For large positive `s`, `\zeta(s)` rapidly approaches 1::

    >>> zeta(30)
    1.00000000093133
    >>> zeta(100)
    1.0
    >>> zeta(inf)
    1.0

The following series converges and in fact has a simple
closed form value::

    >>> nsum(lambda k: zeta(k)-1, [2, inf])
    1.0

**Algorithm**

The primary algorithm is Borwein's algorithm for the Dirichlet
eta function. Three separate implementations are used: for general
real arguments, general complex arguments, and for integers. The
reflection formula is applied to arguments in the negative
half-plane. For very large real arguments, either direct
summation or the Euler prime product is used.

It should be noted that computation of `\zeta(s)` gets very slow
when `s` is far away from the real axis.

**References**

1. http://mathworld.wolfram.com/RiemannZetaFunction.html

2. http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P155.pdf
"""

altzeta = r"""
Gives the Dirichlet eta function, `\eta(s)`, also known as the
alternating zeta function. This function is defined in analogy
with the Riemann zeta function as providing the sum of the
alternating series

.. math ::

    \eta(s) = \sum_{k=0}^{\infty} \frac{(-1)^k}{k^s}
        = 1-\frac{1}{2^s}+\frac{1}{3^s}-\frac{1}{4^s}+\ldots

The eta function, unlike the Riemann zeta function, is an entire
function, having a finite value for all complex `s`. The special case
`\eta(1) = \log(2)` gives the value of the alternating harmonic series.

The alternating zeta function may expressed using the Riemann zeta function
as `\eta(s) = (1 - 2^{1-s}) \zeta(s)`. It can also be expressed
in terms of the Hurwitz zeta function (:func:`hurwitz`), for example using
:func:`dirichlet` (see documentation for that function).

**Examples**

Some special values are::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> altzeta(1)
    0.693147180559945
    >>> altzeta(0)
    0.5
    >>> altzeta(-1)
    0.25
    >>> altzeta(-2)
    0.0

An example of a sum that can be computed more accurately and
efficiently via :func:`altzeta` than via numerical summation::

    >>> sum(-(-1)**n / n**2.5 for n in range(1, 100))
    0.86720495150398402
    >>> altzeta(2.5)
    0.867199889012184

At positive even integers, the Dirichlet eta function
evaluates to a rational multiple of a power of `\pi`::

    >>> altzeta(2)
    0.822467033424113
    >>> pi**2/12
    0.822467033424113

Like the Riemann zeta function, `\eta(s)`, approaches 1
as `s` approaches positive infinity, although it does
so from below rather than from above::

    >>> altzeta(30)
    0.999999999068682
    >>> altzeta(inf)
    1.0
    >>> mp.pretty = False
    >>> altzeta(1000, rounding='d')
    mpf('0.99999999999999989')
    >>> altzeta(1000, rounding='u')
    mpf('1.0')

**References**

1. http://mathworld.wolfram.com/DirichletEtaFunction.html

2. http://en.wikipedia.org/wiki/Dirichlet_eta_function
"""

factorial = r"""
Computes the factorial, `x!`. For integers `n \ge 0`, we have
`n! = 1 \cdot 2 \cdots (n-1) \cdot n` and more generally the factorial
is defined for real or complex `x` by `x! = \Gamma(x+1)`.

**Examples**

Basic values and limits::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for k in range(6):
    ...     print k, fac(k)
    ...
    0 1.0
    1 1.0
    2 2.0
    3 6.0
    4 24.0
    5 120.0
    >>> fac(inf)
    +inf
    >>> fac(0.5), sqrt(pi)/2
    (0.886226925452758, 0.886226925452758)

For large positive `x`, `x!` can be approximated by
Stirling's formula::

    >>> x = 10**10
    >>> fac(x)
    2.32579620567308e+95657055186
    >>> sqrt(2*pi*x)*(x/e)**x
    2.32579597597705e+95657055186

:func:`fac` supports evaluation for astronomically large values::

    >>> fac(10**30)
    6.22311232304258e+29565705518096748172348871081098

Reciprocal factorials appear in the Taylor series of the
exponential function (among many other contexts)::

    >>> nsum(lambda k: 1/fac(k), [0, inf]), exp(1)
    (2.71828182845905, 2.71828182845905)
    >>> nsum(lambda k: pi**k/fac(k), [0, inf]), exp(pi)
    (23.1406926327793, 23.1406926327793)

"""

gamma = r"""
Computes the gamma function, `\Gamma(x)`. The gamma function is a
shifted version of the ordinary factorial, satisfying
`\Gamma(n) = (n-1)!` for integers `n > 0`. More generally, it
is defined by

.. math ::

    \Gamma(x) = \int_0^{\infty} t^{x-1} e^{-t}\, dt

for any real or complex `x` with `\Re(x) > 0` and for `\Re(x) < 0`
by analytic continuation.

**Examples**

Basic values and limits::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for k in range(1, 6):
    ...     print k, gamma(k)
    ...
    1 1.0
    2 1.0
    3 2.0
    4 6.0
    5 24.0
    >>> gamma(inf)
    +inf
    >>> gamma(0)
    Traceback (most recent call last):
      ...
    ValueError: gamma function pole

The gamma function of a half-integer is a rational multiple of
`\sqrt{\pi}`::

    >>> gamma(0.5), sqrt(pi)
    (1.77245385090552, 1.77245385090552)
    >>> gamma(1.5), sqrt(pi)/2
    (0.886226925452758, 0.886226925452758)

We can check the integral definition::

    >>> gamma(3.5)
    3.32335097044784
    >>> quad(lambda t: t**2.5*exp(-t), [0,inf])
    3.32335097044784

:func:`gamma` supports arbitrary-precision evaluation and
complex arguments::

    >>> mp.dps = 50
    >>> gamma(sqrt(3))
    0.91510229697308632046045539308226554038315280564184
    >>> mp.dps = 25
    >>> gamma(2j)
    (0.009902440080927490985955066 - 0.07595200133501806872408048j)

Arguments can also be large. Note that the gamma function grows
very quickly::

    >>> mp.dps = 15
    >>> gamma(10**20)
    1.9328495143101e+1956570551809674817225

"""

psi = r"""
Gives the polygamma function of order `m` of `z`, `\psi^{(m)}(z)`.
Special cases are known as the *digamma function* (`\psi^{(0)}(z)`),
the *trigamma function* (`\psi^{(1)}(z)`), etc. The polygamma
functions are defined as the logarithmic derivatives of the gamma
function:

.. math ::

    \psi^{(m)}(z) = \left(\frac{d}{dz}\right)^{m+1} \log \Gamma(z)

In particular, `\psi^{(0)}(z) = \Gamma'(z)/\Gamma(z)`. In the
present implementation of :func:`psi`, the order `m` must be a
nonnegative integer, while the argument `z` may be an arbitrary
complex number (with exception for the polygamma function's poles
at `z = 0, -1, -2, \ldots`).

**Examples**

For various rational arguments, the polygamma function reduces to
a combination of standard mathematical constants::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> psi(0, 1), -euler
    (-0.5772156649015328606065121, -0.5772156649015328606065121)
    >>> psi(1, '1/4'), pi**2+8*catalan
    (17.19732915450711073927132, 17.19732915450711073927132)
    >>> psi(2, '1/2'), -14*apery
    (-16.82879664423431999559633, -16.82879664423431999559633)

The polygamma functions are derivatives of each other::

    >>> diff(lambda x: psi(3, x), pi), psi(4, pi)
    (-0.1105749312578862734526952, -0.1105749312578862734526952)
    >>> quad(lambda x: psi(4, x), [2, 3]), psi(3,3)-psi(3,2)
    (-0.375, -0.375)

The digamma function diverges logarithmically as `z \to \infty`,
while higher orders tend to zero::

    >>> psi(0,inf), psi(1,inf), psi(2,inf)
    (+inf, 0.0, 0.0)

Evaluation for a complex argument::

    >>> psi(2, -1-2j)
    (0.03902435405364952654838445 + 0.1574325240413029954685366j)

Evaluation is supported for large orders `m` and/or large
arguments `z`::

    >>> psi(3, 10**100)
    2.0e-300
    >>> psi(250, 10**30+10**20*j)
    (-1.293142504363642687204865e-7010 + 3.232856260909107391513108e-7018j)

**Application to infinite series**

Any infinite series where the summand is a rational function of
the index `k` can be evaluated in closed form in terms of polygamma
functions of the roots and poles of the summand::

    >>> a = sqrt(2)
    >>> b = sqrt(3)
    >>> nsum(lambda k: 1/((k+a)**2*(k+b)), [0, inf])
    0.4049668927517857061917531
    >>> (psi(0,a)-psi(0,b)-a*psi(1,a)+b*psi(1,a))/(a-b)**2
    0.4049668927517857061917531

This follows from the series representation (`m > 0`)

.. math ::

    \psi^{(m)}(z) = (-1)^{m+1} m! \sum_{k=0}^{\infty}
        \frac{1}{(z+k)^{m+1}}.

Since the roots of a polynomial may be complex, it is sometimes
necessary to use the complex polygamma function to evaluate
an entirely real-valued sum::

    >>> nsum(lambda k: 1/(k**2-2*k+3), [0, inf])
    1.694361433907061256154665
    >>> nprint(polyroots([1,-2,3]))
    [(1.0 - 1.41421j), (1.0 + 1.41421j)]
    >>> r1 = 1-sqrt(2)*j
    >>> r2 = r1.conjugate()
    >>> (psi(0,-r2)-psi(0,-r1))/(r1-r2)
    (1.694361433907061256154665 + 0.0j)

"""

harmonic = r"""
If `n` is an integer, ``harmonic(n)`` gives a floating-point
approximation of the `n`-th harmonic number `H(n)`, defined as

.. math ::

    H(n) = 1 + \frac{1}{2} + \frac{1}{3} + \ldots + \frac{1}{n}

The firrst few harmonic numbers are::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(8):
    ...     print n, harmonic(n)
    ...
    0 0.0
    1 1.0
    2 1.5
    3 1.83333333333333
    4 2.08333333333333
    5 2.28333333333333
    6 2.45
    7 2.59285714285714

The infinite harmonic series `1 + 1/2 + 1/3 + \ldots` diverges::

    >>> harmonic(inf)
    +inf

:func:`harmonic` is evaluated using the digamma function rather
than by summing the harmonic series term by term. It can therefore
be computed quickly for arbitrarily large `n`, and even for
nonintegral arguments::

    >>> harmonic(10**100)
    230.835724964306
    >>> harmonic(0.5)
    0.613705638880109
    >>> harmonic(3+4j)
    (2.24757548223494 + 0.850502209186044j)

:func:`harmonic` supports arbitrary precision evaluation::

    >>> mp.dps = 50
    >>> harmonic(11)
    3.0198773448773448773448773448773448773448773448773
    >>> harmonic(pi)
    1.8727388590273302654363491032336134987519132374152

The harmonic series diverges, but at a glacial pace. It is possible
to calculate the exact number of terms required before the sum
exceeds a given amount, say 100::

    >>> mp.dps = 50
    >>> v = 10**findroot(lambda x: harmonic(10**x) - 100, 10)
    >>> v
    15092688622113788323693563264538101449859496.864101
    >>> v = int(ceil(v))
    >>> print v
    15092688622113788323693563264538101449859497
    >>> harmonic(v-1)
    99.999999999999999999999999999999999999999999942747
    >>> harmonic(v)
    100.000000000000000000000000000000000000000000009

"""

bernoulli = r"""
Computes the nth Bernoulli number, `B_n`, for any integer `n \ge 0`.

The Bernoulli numbers are rational numbers, but this function
returns a floating-point approximation. To obtain an exact
fraction, use :func:`bernfrac` instead.

**Examples**

Numerical values of the first few Bernoulli numbers::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(15):
    ...     print n, bernoulli(n)
    ...
    0 1.0
    1 -0.5
    2 0.166666666666667
    3 0.0
    4 -0.0333333333333333
    5 0.0
    6 0.0238095238095238
    7 0.0
    8 -0.0333333333333333
    9 0.0
    10 0.0757575757575758
    11 0.0
    12 -0.253113553113553
    13 0.0
    14 1.16666666666667

Bernoulli numbers can be approximated with arbitrary precision::

    >>> mp.dps = 50
    >>> bernoulli(100)
    -2.8382249570693706959264156336481764738284680928013e+78

Arbitrarily large `n` are supported::

    >>> mp.dps = 15
    >>> bernoulli(10**20 + 2)
    3.09136296657021e+1876752564973863312327

The Bernoulli numbers are related to the Riemann zeta function
at integer arguments::

    >>> -bernoulli(8) * (2*pi)**8 / (2*fac(8))
    1.00407735619794
    >>> zeta(8)
    1.00407735619794

**Algorithm**

For small `n` (`n < 3000`) :func:`bernoulli` uses a recurrence
formula due to Ramanujan. All results in this range are cached,
so sequential computation of small Bernoulli numbers is
guaranteed to be fast.

For larger `n`, `B_n` is evaluated in terms of the Riemann zeta
function.
"""

stieltjes = r"""
For a nonnegative integer `n`, ``stieltjes(n)`` computes the
`n`-th Stieltjes constant `\gamma_n`, defined as the
`n`-th coefficient in the Laurent series expansion of the
Riemann zeta function around the pole at `s = 1`. That is,
we have:

.. math ::

  \zeta(s) = \frac{1}{s-1} \sum_{n=0}^{\infty}
      \frac{(-1)^n}{n!} \gamma_n (s-1)^n

More generally, ``stieltjes(n, a)`` gives the corresponding
coefficient `\gamma_n(a)` for the Hurwitz zeta function
`\zeta(s,a)` (with `\gamma_n = \gamma_n(1)`).

**Examples**

The zeroth Stieltjes constant is just Euler's constant `\gamma`::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> stieltjes(0)
    0.577215664901533

Some more values are::

    >>> stieltjes(1)
    -0.0728158454836767
    >>> stieltjes(10)
    0.000205332814909065
    >>> stieltjes(30)
    0.00355772885557316
    >>> stieltjes(1000)
    -1.57095384420474e+486
    >>> stieltjes(2000)
    2.680424678918e+1109
    >>> stieltjes(1, 2.5)
    -0.23747539175716

An alternative way to compute `\gamma_1`::

    >>> diff(extradps(25)(lambda x: 1/(x-1) - zeta(x)), 1)
    -0.0728158454836767

:func:`stieltjes` supports arbitrary precision evaluation::

    >>> mp.dps = 50
    >>> stieltjes(2)
    -0.0096903631928723184845303860352125293590658061013408

**Algorithm**

:func:`stieltjes` numerically evaluates the integral in
the following representation due to Ainsworth, Howell and
Coffey [1], [2]:

.. math ::

  \gamma_n(a) = \frac{\log^n a}{2a} - \frac{\log^{n+1}(a)}{n+1} +
      \frac{2}{a} \Re \int_0^{\infty}
      \frac{(x/a-i)\log^n(a-ix)}{(1+x^2/a^2)(e^{2\pi x}-1)} dx.

For some reference values with `a = 1`, see e.g. [4].

**References**

1. O. R. Ainsworth & L. W. Howell, "An integral representation of
   the generalized Euler-Mascheroni constants", NASA Technical
   Paper 2456 (1985),
   http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19850014994_1985014994.pdf

2. M. W. Coffey, "The Stieltjes constants, their relation to the
   `\eta_j` coefficients, and representation of the Hurwitz
   zeta function", 	arXiv:0706.0343v1 http://arxiv.org/abs/0706.0343

3. http://mathworld.wolfram.com/StieltjesConstants.html

4. http://pi.lacim.uqam.ca/piDATA/stieltjesgamma.txt

"""

gammaprod = r"""
Given iterables `a` and `b`, ``gammaprod(a, b)`` computes the
product / quotient of gamma functions:

.. math ::

    \frac{\Gamma(a_0) \Gamma(a_1) \cdots \Gamma(a_p)}
         {\Gamma(b_0) \Gamma(b_1) \cdots \Gamma(b_q)}

Unlike direct calls to :func:`gamma`, :func:`gammaprod` considers
the entire product as a limit and evaluates this limit properly if
any of the numerator or denominator arguments are nonpositive
integers such that poles of the gamma function are encountered.
That is, :func:`gammaprod` evaluates

.. math ::

    \lim_{\epsilon \to 0}
    \frac{\Gamma(a_0+\epsilon) \Gamma(a_1+\epsilon) \cdots
        \Gamma(a_p+\epsilon)}
         {\Gamma(b_0+\epsilon) \Gamma(b_1+\epsilon) \cdots
        \Gamma(b_q+\epsilon)}

In particular:

* If there are equally many poles in the numerator and the
  denominator, the limit is a rational number times the remaining,
  regular part of the product.

* If there are more poles in the numerator, :func:`gammaprod`
  returns ``+inf``.

* If there are more poles in the denominator, :func:`gammaprod`
  returns 0.

**Examples**

The reciprocal gamma function `1/\Gamma(x)` evaluated at `x = 0`::

    >>> from mpmath import *
    >>> mp.dps = 15
    >>> gammaprod([], [0])
    0.0

A limit::

    >>> gammaprod([-4], [-3])
    -0.25
    >>> limit(lambda x: gamma(x-1)/gamma(x), -3, direction=1)
    -0.25
    >>> limit(lambda x: gamma(x-1)/gamma(x), -3, direction=-1)
    -0.25

"""

beta = r"""
Computes the beta function,
`B(x,y) = \Gamma(x) \Gamma(y) / \Gamma(x+y)`.
The beta function is also commonly defined by the integral
representation

.. math ::

    B(x,y) = \int_0^1 t^{x-1} (1-t)^{y-1} \, dt

**Examples**

For integer and half-integer arguments where all three gamma
functions are finite, the beta function becomes either rational
number or a rational multiple of `\pi`::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> beta(5, 2)
    0.0333333333333333
    >>> beta(1.5, 2)
    0.266666666666667
    >>> 16*beta(2.5, 1.5)
    3.14159265358979

Where appropriate, :func:`beta` evaluates limits. A pole
of the beta function is taken to result in ``+inf``::

    >>> beta(-0.5, 0.5)
    0.0
    >>> beta(-3, 3)
    -0.333333333333333
    >>> beta(-2, 3)
    +inf
    >>> beta(inf, 1)
    0.0
    >>> beta(inf, 0)
    nan

:func:`beta` supports complex numbers and arbitrary precision
evaluation::

    >>> beta(1, 2+j)
    (0.4 - 0.2j)
    >>> mp.dps = 25
    >>> beta(j,0.5)
    (1.079424249270925780135675 - 1.410032405664160838288752j)
    >>> mp.dps = 50
    >>> beta(pi, e)
    0.037890298781212201348153837138927165984170287886464

Various integrals can be computed by means of the
beta function::

    >>> mp.dps = 15
    >>> quad(lambda t: t**2.5*(1-t)**2, [0, 1])
    0.0230880230880231
    >>> beta(3.5, 3)
    0.0230880230880231
    >>> quad(lambda t: sin(t)**4 * sqrt(cos(t)), [0, pi/2])
    0.319504062596158
    >>> beta(2.5, 0.75)/2
    0.319504062596158

"""

betainc = r"""
``betainc(a, b, x1=0, x2=1, regularized=False)`` gives the generalized
incomplete beta function,

.. math ::

    I_{x_1}^{x_2}(a,b) = \int_{x_1}^{x_2} t^{a-1} (1-t)^{b-1} dt.

When `x_1 = 0, x_2 = 1`, this reduces to the ordinary (complete)
beta function `B(a,b)`; see :func:`beta`.

With the keyword argument ``regularized=True``, :func:`betainc`
computes the regularized incomplete beta function
`I_{x_1}^{x_2}(a,b) / B(a,b)`. This is the cumulative distribution of the
beta distribution with parameters `a`, `b`.

Note: implementations of the incomplete beta function in some other
software uses a different argument order. For example, Mathematica uses the
reversed argument order ``Beta[x1,x2,a,b]``. For the equivalent of SciPy's
three-argument incomplete beta integral (implicitly with `x1 = 0`), use
``betainc(a,b,0,x2,regularized=True)``.

**Examples**

Verifying that :func:`betainc` computes the integral in the
definition::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> x,y,a,b = 3, 4, 0, 6
    >>> betainc(x, y, a, b)
    -4010.4
    >>> quad(lambda t: t**(x-1) * (1-t)**(y-1), [a, b])
    -4010.4

The arguments may be arbitrary complex numbers::

    >>> betainc(0.75, 1-4j, 0, 2+3j)
    (0.2241657956955709603655887 + 0.3619619242700451992411724j)

With regularization::

    >>> betainc(1, 2, 0, 0.25, regularized=True)
    0.4375
    >>> betainc(pi, e, 0, 1, regularized=True)   # Complete
    1.0

The beta integral satisfies some simple argument transformation
symmetries::

    >>> mp.dps = 15
    >>> betainc(2,3,4,5), -betainc(2,3,5,4), betainc(3,2,1-5,1-4)
    (56.0833333333333, 56.0833333333333, 56.0833333333333)

The beta integral can often be evaluated analytically. For integer and
rational arguments, the incomplete beta function typically reduces to a
simple algebraic-logarithmic expression::

    >>> mp.dps = 25
    >>> identify(chop(betainc(0, 0, 3, 4)))
    '-(log((9/8)))'
    >>> identify(betainc(2, 3, 4, 5))
    '(673/12)'
    >>> identify(betainc(1.5, 1, 1, 2))
    '((-12+sqrt(1152))/18)'

"""

binomial = r"""
Computes the binomial coefficient

.. math ::

    {n \choose k} = \frac{n!}{k!(n-k)!}.

The binomial coefficient gives the number of ways that `k` items
can be chosen from a set of `n` items. More generally, the binomial
coefficient is a well-defined function of arbitrary real or
complex `n` and `k`, via the gamma function.

**Examples**

Generate Pascal's triangle::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(5):
    ...     nprint([binomial(n,k) for k in range(n+1)])
    ...
    [1.0]
    [1.0, 1.0]
    [1.0, 2.0, 1.0]
    [1.0, 3.0, 3.0, 1.0]
    [1.0, 4.0, 6.0, 4.0, 1.0]

There is 1 way to select 0 items from the empty set, and 0 ways to
select 1 item from the empty set::

    >>> binomial(0, 0)
    1.0
    >>> binomial(0, 1)
    0.0

:func:`binomial` supports large arguments::

    >>> binomial(10**20, 10**20-5)
    8.33333333333333e+97
    >>> binomial(10**20, 10**10)
    2.60784095465201e+104342944813

Nonintegral binomial coefficients find use in series
expansions::

    >>> nprint(taylor(lambda x: (1+x)**0.25, 0, 4))
    [1.0, 0.25, -9.375e-2, 5.46875e-2, -3.75977e-2]
    >>> nprint([binomial(0.25, k) for k in range(5)])
    [1.0, 0.25, -9.375e-2, 5.46875e-2, -3.75977e-2]

An integral representation::

    >>> n, k = 5, 3
    >>> f = lambda t: exp(-j*k*t)*(1+exp(j*t))**n
    >>> chop(quad(f, [-pi,pi])/(2*pi))
    10.0
    >>> binomial(n,k)
    10.0

"""

rf = r"""
Computes the rising factorial or Pochhammer symbol,

.. math ::

    x^{(n)} = x (x+1) \cdots (x+n-1) = \frac{\Gamma(x+n)}{\Gamma(x)}

where the rightmost expression is valid for nonintegral `n`.

**Examples**

For integral `n`, the rising factorial is a polynomial::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(5):
    ...     nprint(taylor(lambda x: rf(x,n), 0, n))
    ...
    [1.0]
    [0.0, 1.0]
    [0.0, 1.0, 1.0]
    [0.0, 2.0, 3.0, 1.0]
    [0.0, 6.0, 11.0, 6.0, 1.0]

Evaluation is supported for arbitrary arguments::

    >>> rf(2+3j, 5.5)
    (-7202.03920483347 - 3777.58810701527j)
"""

ff = r"""
Computes the falling factorial,

.. math ::

    (x)_n = x (x-1) \cdots (x-n+1) = \frac{\Gamma(x+1)}{\Gamma(x-n+1)}

where the rightmost expression is valid for nonintegral `n`.

**Examples**

For integral `n`, the falling factorial is a polynomial::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(5):
    ...     nprint(taylor(lambda x: ff(x,n), 0, n))
    ...
    [1.0]
    [0.0, 1.0]
    [0.0, -1.0, 1.0]
    [0.0, 2.0, -3.0, 1.0]
    [0.0, -6.0, 11.0, -6.0, 1.0]

Evaluation is supported for arbitrary arguments::

    >>> ff(2+3j, 5.5)
    (-720.41085888203 + 316.101124983878j)
"""

fac2 = r"""
Computes the double factorial `x!!`, defined for integers
`x > 0` by

.. math ::

    x!! = \begin{cases}
        1 \cdot 3 \cdots (x-2) \cdot x & x \;\mathrm{odd} \\
        2 \cdot 4 \cdots (x-2) \cdot x & x \;\mathrm{even}
    \end{cases}

and more generally by [1]

.. math ::

    x!! = 2^{x/2} \left(\frac{\pi}{2}\right)^{(\cos(\pi x)-1)/4}
          \Gamma\left(\frac{x}{2}+1\right).

**Examples**

The integer sequence of double factorials begins::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> nprint([fac2(n) for n in range(10)])
    [1.0, 1.0, 2.0, 3.0, 8.0, 15.0, 48.0, 105.0, 384.0, 945.0]

For large `x`, double factorials follow a Stirling-like asymptotic
approximation::

    >>> x = mpf(10000)
    >>> fac2(x)
    5.97272691416282e+17830
    >>> sqrt(pi)*x**((x+1)/2)*exp(-x/2)
    5.97262736954392e+17830

The recurrence formula `x!! = x (x-2)!!` can be reversed to
define the double factorial of negative odd integers (but
not negative even integers)::

    >>> fac2(-1), fac2(-3), fac2(-5), fac2(-7)
    (1.0, -1.0, 0.333333333333333, -0.0666666666666667)
    >>> fac2(-2)
    Traceback (most recent call last):
      ...
    ValueError: gamma function pole

With the exception of the poles at negative even integers,
:func:`fac2` supports evaluation for arbitrary complex arguments.
The recurrence formula is valid generally::

    >>> fac2(pi+2j)
    (-1.3697207890154e-12 + 3.93665300979176e-12j)
    >>> (pi+2j)*fac2(pi-2+2j)
    (-1.3697207890154e-12 + 3.93665300979176e-12j)

Double factorials should not be confused with nested factorials,
which are immensely larger::

    >>> fac(fac(20))
    5.13805976125208e+43675043585825292774
    >>> fac2(20)
    3715891200.0

Double factorials appear, among other things, in series expansions
of Gaussian functions and the error function. Infinite series
include::

    >>> nsum(lambda k: 1/fac2(k), [0, inf])
    3.05940740534258
    >>> sqrt(e)*(1+sqrt(pi/2)*erf(sqrt(2)/2))
    3.05940740534258
    >>> nsum(lambda k: 2**k/fac2(2*k-1), [1, inf])
    4.06015693855741
    >>> e * erf(1) * sqrt(pi)
    4.06015693855741

A beautiful Ramanujan sum::

    >>> nsum(lambda k: (-1)**k*(fac2(2*k-1)/fac2(2*k))**3, [0,inf])
    0.90917279454693
    >>> (gamma('9/8')/gamma('5/4')/gamma('7/8'))**2
    0.90917279454693

**References**

1. http://functions.wolfram.com/GammaBetaErf/Factorial2/27/01/0002/

2. http://mathworld.wolfram.com/DoubleFactorial.html

"""

hyper = r"""
Evaluates the generalized hypergeometric function

.. math ::

    \,_pF_q(a_1,\ldots,a_p; b_1,\ldots,b_q; z) =
    \sum_{n=0}^\infty \frac{(a_1)_n (a_2)_n \ldots (a_p)_n}
       {(b_1)_n(b_2)_n\ldots(b_q)_n} \frac{z^n}{n!}

where `(x)_n` denotes the rising factorial (see :func:`rf`).

The parameters lists ``a_s`` and ``b_s`` may contain integers,
real numbers, complex numbers, as well as exact fractions given in
the form of tuples `(p, q)`. :func:`hyper` is optimized to handle
integers and fractions more efficiently than arbitrary
floating-point parameters (since rational parameters are by
far the most common).

**Examples**

Verifying that :func:`hyper` gives the sum in the definition, by
comparison with :func:`nsum`::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> a,b,c,d = 2,3,4,5
    >>> x = 0.25
    >>> hyper([a,b],[c,d],x)
    1.078903941164934876086237
    >>> fn = lambda n: rf(a,n)*rf(b,n)/rf(c,n)/rf(d,n)*x**n/fac(n)
    >>> nsum(fn, [0, inf])
    1.078903941164934876086237

The parameters can be any combination of integers, fractions,
floats and complex numbers::

    >>> a, b, c, d, e = 1, (-1,2), pi, 3+4j, (2,3)
    >>> x = 0.2j
    >>> hyper([a,b],[c,d,e],x)
    (0.9923571616434024810831887 - 0.005753848733883879742993122j)
    >>> b, e = -0.5, mpf(2)/3
    >>> fn = lambda n: rf(a,n)*rf(b,n)/rf(c,n)/rf(d,n)/rf(e,n)*x**n/fac(n)
    >>> nsum(fn, [0, inf])
    (0.9923571616434024810831887 - 0.005753848733883879742993122j)

The `\,_0F_0` and `\,_1F_0` series are just elementary functions::

    >>> a, z = sqrt(2), pi
    >>> hyper([],[],z)
    23.14069263277926900572909
    >>> exp(z)
    23.14069263277926900572909
    >>> hyper([a],[],z)
    (-0.09069132879922920160334114 + 0.3283224323946162083579656j)
    >>> (1-z)**(-a)
    (-0.09069132879922920160334114 + 0.3283224323946162083579656j)

If any `a_k` coefficient is a nonpositive integer, the series terminates
into a finite polynomial::

    >>> hyper([1,1,1,-3],[2,5],1)
    0.7904761904761904761904762
    >>> identify(_)
    '(83/105)'

If any `b_k` is a nonpositive integer, the function is undefined (unless the
series terminates before the division by zero occurs)::

    >>> hyper([1,1,1,-3],[-2,5],1)
    Traceback (most recent call last):
      ...
    ZeroDivisionError
    >>> hyper([1,1,1,-1],[-2,5],1)
    1.1

Generally, the radius of convergence of the hypergeometric series is either 0
(if `p > q+1`), 1 (if `p = q+1`), or `\infty`. A divergent series will evaluate
if an analytic continuation or regularization is implemented (see the list of
:func:`hypXfX` functions); otherwise an exception is raised::

    >>> hyper([1,1], [1], 2)           # divergent case of 2F1, implemented
    -1.0
    >>> hyper([1,1,1,1], [1,1,1], 0.5) # convergent case of 4F3
    2.0
    >>> hyper([1,1,1,1], [1,1,1], 2)   # divergent case of 4F3, not implemented
    Traceback (most recent call last):
      ...
    NoConvergence
    >>> hyper([1,1], [], 0.5)          # regularization of 2F0, implemented
    (1.340965419580146562086448 + 0.8503366631752726568782447j)
    >>> hyper([1,1,1,1], [1], 0.5)     # regularization of 4F1, not implemented
    Traceback (most recent call last):
      ...
    NoConvergence

"""

hypercomb = r"""
Computes a weighted combination of hypergeometric functions

.. math ::

    \sum_{r=1}^N \left[ \prod_{k=1}^{l_r} {w_{r,k}}^{c_{r,k}}
    \frac{\prod_{k=1}^{m_r} \Gamma(\alpha_{r,k})}{\prod_{k=1}^{n_r}
    \Gamma(\beta_{r,k})}
    \,_{p_r}F_{q_r}(a_{r,1},\ldots,a_{r,p}; b_{r,1},    
    \ldots, b_{r,q}; z_r)\right].

Typically the parameters are linear combinations of a small set of base
parameters; :func:`hypercomb` permits computing a correct value in
the case that some of the `\alpha`, `\beta`, `b` turn out to be
nonpositive integers, or if division by zero occurs for some `w^c`,
assuming that there are opposing singularities that cancel out.
The limit is computed by evaluating the function with the base
parameters perturbed, at a higher working precision.

The first argument should be a function that takes the perturbable
base parameters ``params`` as input and returns `N` tuples
``(w, c, alpha, beta, a, b, z)``, where the coefficients ``w``, ``c``,
gamma factors ``alpha``, ``beta``, and hypergeometric coefficients
``a``, ``b`` each should be lists of numbers, and ``z`` should be a single
number.

**Examples**

The following evaluates

.. math ::

    (a-1) \frac{\Gamma(a-3)}{\Gamma(a-4)} \,_1F_1(a,a-1,z) = e^z(a-4)(a+z-1)

with `a=1, z=3`. There is a zero factor, two gamma function poles, and
the 1F1 function is singular; all singularities cancel out to give a finite
value::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> hypercomb(lambda a: [([a-1],[1],[a-3],[a-4],[a],[a-1],3)], [1])
    -180.769832308689
    >>> -9*exp(3)
    -180.769832308689

"""

hyp0f1 = r"""
Gives the hypergeometric function `\,_0F_1`, sometimes known as the
confluent limit function, defined as

.. math ::

    \,_0F_1(a,z) = \sum_{k=0}^{\infty} \frac{1}{(a)_k} \frac{z^k}{k!}.

This function satisfies the differential equation `z f''(z) + a f'(z) = f(z)`,
and is related to the Bessel function of the first kind (see :func:`besselj`).

``hyp0f1(a,z)`` is equivalent to ``hyper([],[a],z)``; see documentation for
:func:`hyper` for more information.

**Examples**

Evaluation for arbitrary arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> hyp0f1(2, 0.25)
    1.130318207984970054415392
    >>> hyp0f1((1,2), 1234567)
    6.27287187546220705604627e+964
    >>> hyp0f1(3+4j, 1000000j)
    (3.905169561300910030267132e+606 + 3.807708544441684513934213e+606j)

Evaluation is supported for arbitrarily large values of `z`,
using asymptotic expansions::

    >>> hyp0f1(1, 10**50)
    2.131705322874965310390701e+8685889638065036553022565
    >>> hyp0f1(1, -10**50)
    1.115945364792025420300208e-13

Verifying the differential equation::

    >>> a = 2.5
    >>> f = lambda z: hyp0f1(a,z)
    >>> for z in [0, 10, 3+4j]:
    ...     chop(z*diff(f,z,2) + a*diff(f,z) - f(z))
    ...
    0.0
    0.0
    0.0

"""

hyp1f1 = r"""
Gives the confluent hypergeometric function of the first kind,

.. math ::

    \,_1F_1(a,b,z) = \sum_{k=0}^{\infty} \frac{(a)_k}{(b)_k} \frac{z^k}{k!},

also known as Kummer's function and sometimes denoted by `M(a,b,z)`. This
function gives one solution to the confluent (Kummer's) differential equation

.. math ::

    z f''(z) + (b-z) f'(z) - af(z) = 0.

A second solution is given by the `U` function; see :func:`hyperu`.
Solutions are also given in an alternate form by the Whittaker
functions (:func:`whitm`, :func:`whitw`).

``hyp1f1(a,b,z)`` is equivalent
to ``hyper([a],[b],z)``; see documentation for :func:`hyper` for more
information.

**Examples**

Evaluation for real and complex values of the argument `z`, with
fixed parameters `a = 2, b = -1/3`::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> hyp1f1(2, (-1,3), 3.25)
    -2815.956856924817275640248
    >>> hyp1f1(2, (-1,3), -3.25)
    -1.145036502407444445553107
    >>> hyp1f1(2, (-1,3), 1000)
    -8.021799872770764149793693e+441
    >>> hyp1f1(2, (-1,3), -1000)
    0.000003131987633006813594535331
    >>> hyp1f1(2, (-1,3), 100+100j)
    (-3.189190365227034385898282e+48 - 1.106169926814270418999315e+49j)

Parameters may be complex::

    >>> hyp1f1(2+3j, -1+j, 10j)
    (261.8977905181045142673351 + 160.8930312845682213562172j)

Arbitrarily large values of `z` are supported::

    >>> hyp1f1(3, 4, 10**20)
    3.890569218254486878220752e+43429448190325182745
    >>> hyp1f1(3, 4, -10**20)
    6.0e-60
    >>> hyp1f1(3, 4, 10**20*j)
    (-1.935753855797342532571597e-20 - 2.291911213325184901239155e-20j)

Verifying the differential equation::

    >>> a, b = 1.5, 2
    >>> f = lambda z: hyp1f1(a,b,z)
    >>> for z in [0, -10, 3, 3+4j]:
    ...     chop(z*diff(f,z,2) + (b-z)*diff(f,z) - a*f(z))
    ...
    0.0
    0.0
    0.0
    0.0

An integral representation::

    >>> a, b = 1.5, 3
    >>> z = 1.5
    >>> hyp1f1(a,b,z)
    2.269381460919952778587441
    >>> g = lambda t: exp(z*t)*t**(a-1)*(1-t)**(b-a-1)
    >>> gammaprod([b],[a,b-a])*quad(g, [0,1])
    2.269381460919952778587441


"""

hyp1f2 = r"""
Gives the hypergeometric function `\,_1F_2(a_1,a_2;b_1,b_2; z)`.
The call ``hyp1f2(a1,b1,b2,z)`` is equivalent to
``hyper([a1],[b1,b2],z)``.

Evaluation works for complex and arbitrarily large arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> a, b, c = 1.5, (-1,3), 2.25
    >>> hyp1f2(a, b, c, 10**20)
    -1.159388148811981535941434e+8685889639
    >>> hyp1f2(a, b, c, -10**20)
    -12.60262607892655945795907
    >>> hyp1f2(a, b, c, 10**20*j)
    (4.237220401382240876065501e+6141851464 - 2.950930337531768015892987e+6141851464j)
    >>> hyp1f2(2+3j, -2j, 0.5j, 10-20j)
    (135881.9905586966432662004 - 86681.95885418079535738828j)

"""

hyp2f2 = r"""
Gives the hypergeometric function `\,_2F_2(a_1,a_2;b_1,b_2; z)`.
The call ``hyp2f2(a1,a2,b1,b2,z)`` is equivalent to
``hyper([a1,a2],[b1,b2],z)``.

Evaluation works for complex and arbitrarily large arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> a, b, c, d = 1.5, (-1,3), 2.25, 4
    >>> hyp2f2(a, b, c, d, 10**20)
    -5.275758229007902299823821e+43429448190325182663
    >>> hyp2f2(a, b, c, d, -10**20)
    2561445.079983207701073448
    >>> hyp2f2(a, b, c, d, 10**20*j)
    (2218276.509664121194836667 - 1280722.539991603850462856j)
    >>> hyp2f2(2+3j, -2j, 0.5j, 4j, 10-20j)
    (80500.68321405666957342788 - 20346.82752982813540993502j)

"""

hyp2f3 = r"""
Gives the hypergeometric function `\,_2F_3(a_1,a_2;b_1,b_2,b_3; z)`.
The call ``hyp2f3(a1,a2,b1,b2,b3,z)`` is equivalent to
``hyper([a1,a2],[b1,b2,b3],z)``.

Evaluation works for arbitrarily large arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> a1,a2,b1,b2,b3 = 1.5, (-1,3), 2.25, 4, (1,5)
    >>> hyp2f3(a1,a2,b1,b2,b3,10**20)
    -4.169178177065714963568963e+8685889590
    >>> hyp2f3(a1,a2,b1,b2,b3,-10**20)
    7064472.587757755088178629
    >>> hyp2f3(a1,a2,b1,b2,b3,10**20*j)
    (-5.163368465314934589818543e+6141851415 + 1.783578125755972803440364e+6141851416j)
    >>> hyp2f3(2+3j, -2j, 0.5j, 4j, -1-j, 10-20j)
    (-2280.938956687033150740228 + 13620.97336609573659199632j)
    >>> hyp2f3(2+3j, -2j, 0.5j, 4j, -1-j, 10000000-20000000j)
    (4.849835186175096516193e+3504 - 3.365981529122220091353633e+3504j)

"""

hyp2f1 = r"""
Gives the Gauss hypergeometric function `\,_2F_1` (often simply referred to as
*the* hypergeometric function), defined for `|z| < 1` as

.. math ::

    \,_2F_1(a,b,c,z) = \sum_{k=0}^{\infty}
        \frac{(a)_k (b)_k}{(c)_k} \frac{z^k}{k!}.

and for `|z| \ge 1` by analytic continuation, with a branch cut on `(1, \infty)`
when necessary.

Special cases of this function include many of the orthogonal polynomials as
well as the incomplete beta function and other functions. Properties of the
Gauss hypergeometric function are documented comprehensively in many references,
for example Abramowitz & Stegun, section 15.

The implementation supports the analytic continuation as well as evaluation
close to the unit circle where `|z| \approx 1`. The syntax ``hyp2f1(a,b,c,z)``
is equivalent to ``hyper([a,b],[c],z)``.

**Examples**

Evaluation with `z` inside, outside and on the unit circle, for
fixed parameters::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> hyp2f1(2, (1,2), 4, 0.75)
    1.303703703703703703703704
    >>> hyp2f1(2, (1,2), 4, -1.75)
    0.7431290566046919177853916
    >>> hyp2f1(2, (1,2), 4, 1.75)
    (1.418075801749271137026239 - 1.114976146679907015775102j)
    >>> hyp2f1(2, (1,2), 4, 1)
    1.6
    >>> hyp2f1(2, (1,2), 4, -1)
    0.8235498012182875315037882
    >>> hyp2f1(2, (1,2), 4, j)
    (0.9144026291433065674259078 + 0.2050415770437884900574923j)
    >>> hyp2f1(2, (1,2), 4, 2+j)
    (0.9274013540258103029011549 + 0.7455257875808100868984496j)
    >>> hyp2f1(2, (1,2), 4, 0.25j)
    (0.9931169055799728251931672 + 0.06154836525312066938147793j)

Evaluation with complex parameter values::

    >>> hyp2f1(1+j, 0.75, 10j, 1+5j)
    (0.8834833319713479923389638 + 0.7053886880648105068343509j)

Evaluation with `z = 1`::

    >>> hyp2f1(-2.5, 3.5, 1.5, 1)
    0.0
    >>> hyp2f1(-2.5, 3, 4, 1)
    0.06926406926406926406926407
    >>> hyp2f1(2, 3, 4, 1)
    +inf

Evaluation for huge arguments::

    >>> hyp2f1((-1,3), 1.75, 4, '1e100')
    (7.883714220959876246415651e+32 + 1.365499358305579597618785e+33j)
    >>> hyp2f1((-1,3), 1.75, 4, '1e1000000')
    (7.883714220959876246415651e+333332 + 1.365499358305579597618785e+333333j)
    >>> hyp2f1((-1,3), 1.75, 4, '1e1000000j')
    (1.365499358305579597618785e+333333 - 7.883714220959876246415651e+333332j)

An integral representation::

    >>> a,b,c,z = -0.5, 1, 2.5, 0.25
    >>> g = lambda t: t**(b-1) * (1-t)**(c-b-1) * (1-t*z)**(-a)
    >>> gammaprod([c],[b,c-b]) * quad(g, [0,1])
    0.9480458814362824478852618
    >>> hyp2f1(a,b,c,z)
    0.9480458814362824478852618

Verifying the hypergeometric differential equation::

    >>> f = lambda z: hyp2f1(a,b,c,z)
    >>> chop(z*(1-z)*diff(f,z,2) + (c-(a+b+1)*z)*diff(f,z) - a*b*f(z))
    0.0

"""

hyperu = r"""
Gives the Tricomi confluent hypergeometric function `U`, also known as
the Kummer or confluent hypergeometric function of the second kind. This
function gives a second linearly independent solution to the confluent
hypergeometric differential equation (the first is provided by `\,_1F_1`  --
see :func:`hyp1f1`).

**Examples**

Evaluation for arbitrary complex arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> hyperu(2,3,4)
    0.0625
    >>> hyperu(0.25, 5, 1000)
    0.1779949416140579573763523
    >>> hyperu(0.25, 5, -1000)
    (0.1256256609322773150118907 - 0.1256256609322773150118907j)

The `U` function may be singular at `z = 0`::

    >>> hyperu(1.5, 2, 0)
    +inf
    >>> hyperu(1.5, -2, 0)
    0.1719434921288400112603671

Verifying the differential equation::

    >>> a, b = 1.5, 2
    >>> f = lambda z: hyperu(a,b,z)
    >>> for z in [-10, 3, 3+4j]:
    ...     chop(z*diff(f,z,2) + (b-z)*diff(f,z) - a*f(z))
    ...
    0.0
    0.0
    0.0

An integral representation::

    >>> a,b,z = 2, 3.5, 4.25
    >>> hyperu(a,b,z)
    0.06674960718150520648014567
    >>> quad(lambda t: exp(-z*t)*t**(a-1)*(1+t)**(b-a-1),[0,inf]) / gamma(a)
    0.06674960718150520648014567


[1] http://www.math.ucla.edu/~cbm/aands/page_504.htm
"""

hyp2f0 = r"""
Gives the hypergeometric function `\,_2F_0`, defined formally by the
series

.. math ::

    \,_2F_0(a,b;;z) = \sum_{n=0}^{\infty} (a)_n (b)_n \frac{z^n}{n!}.

This series usually does not converge. For small enough `z`, it can be viewed
as an asymptotic series that may be summed directly with an appropriate
truncation. When this is not the case, :func:`hyp2f0` gives a regularized sum,
or equivalently, it uses a representation in terms of the
hypergeometric U function [1]. The series also converges when either `a` or `b`
is a nonpositive integer, as it then terminates into a polynomial
after `-a` or `-b` terms.

**Examples**

Evaluation is supported for arbitrary complex arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> hyp2f0((2,3), 1.25, -100)
    0.07095851870980052763312791
    >>> hyp2f0((2,3), 1.25, 100)
    (-0.03254379032170590665041131 + 0.07269254613282301012735797j)
    >>> hyp2f0(-0.75, 1-j, 4j)
    (-0.3579987031082732264862155 - 3.052951783922142735255881j)

Even with real arguments, the regularized value of 2F0 is often complex-valued,
but the imaginary part decreases exponentially as `z \to 0`. In the following
example, the first call uses complex evaluation while the second has a small
enough `z` to evaluate using the direct series and thus the returned value
is strictly real (this should be taken to indicate that the imaginary
part is less than ``eps``)::

    >>> mp.dps = 15
    >>> hyp2f0(1.5, 0.5, 0.05)
    (1.04166637647907 + 8.34584913683906e-8j)
    >>> hyp2f0(1.5, 0.5, 0.0005)
    1.00037535207621

The imaginary part can be retrieved by increasing the working precision::

    >>> mp.dps = 80
    >>> nprint(hyp2f0(1.5, 0.5, 0.009).imag)
    1.23828e-46

In the polynomial case (the series terminating), 2F0 can evaluate exactly::

    >>> mp.dps = 15
    >>> hyp2f0(-6,-6,2)
    291793.0
    >>> identify(hyp2f0(-2,1,0.25))
    '(5/8)'

The coefficients of the polynomials can be recovered using Taylor expansion::

    >>> nprint(taylor(lambda x: hyp2f0(-3,0.5,x), 0, 10))
    [1.0, -1.5, 2.25, -1.875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    >>> nprint(taylor(lambda x: hyp2f0(-4,0.5,x), 0, 10))
    [1.0, -2.0, 4.5, -7.5, 6.5625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


[1] http://www.math.ucla.edu/~cbm/aands/page_504.htm
"""


gammainc = r"""
``gammainc(z, a=0, b=inf)`` computes the (generalized) incomplete
gamma function with integration limits `[a, b]`:

.. math ::

  \Gamma(z,a,b) = \int_a^b t^{z-1} e^{-t} \, dt

The generalized incomplete gamma function reduces to the
following special cases when one or both endpoints are fixed:

* `\Gamma(z,0,\infty)` is the standard ("complete")
  gamma function, `\Gamma(z)` (available directly
  as the mpmath function :func:`gamma`)
* `\Gamma(z,a,\infty)` is the "upper" incomplete gamma
  function, `\Gamma(z,a)`
* `\Gamma(z,0,b)` is the "lower" incomplete gamma
  function, `\gamma(z,b)`.

Of course, we have
`\Gamma(z,0,x) + \Gamma(z,x,\infty) = \Gamma(z)`
for all `z` and `x`.

Note however that some authors reverse the order of the
arguments when defining the lower and upper incomplete
gamma function, so one should be careful to get the correct
definition.

If also given the keyword argument ``regularized=True``,
:func:`gammainc` computes the "regularized" incomplete gamma
function

.. math ::

  P(z,a,b) = \frac{\Gamma(z,a,b)}{\Gamma(z)}.

**Examples**

We can compare with numerical quadrature to verify that
:func:`gammainc` computes the integral in the definition::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> gammainc(2+3j, 4, 10)
    (0.00977212668627705160602312 - 0.0770637306312989892451977j)
    >>> quad(lambda t: t**(2+3j-1) * exp(-t), [4, 10])
    (0.00977212668627705160602312 - 0.0770637306312989892451977j)

Argument symmetries follow directly from the integral definition::

    >>> gammainc(3, 4, 5) + gammainc(3, 5, 4)
    0.0
    >>> gammainc(3,0,2) + gammainc(3,2,4); gammainc(3,0,4)
    1.523793388892911312363331
    1.523793388892911312363331
    >>> findroot(lambda z: gammainc(2,z,3), 1)
    3.0

Evaluation for arbitrarily large arguments::

    >>> gammainc(10, 100)
    4.083660630910611272288592e-26
    >>> gammainc(10, 10000000000000000)
    5.290402449901174752972486e-4342944819032375
    >>> gammainc(3+4j, 1000000+1000000j)
    (-1.257913707524362408877881e-434284 + 2.556691003883483531962095e-434284j)

Evaluation of a generalized incomplete gamma function automatically chooses
the representation that gives a more accurate result, depending on which
parameter is larger::

    >>> gammainc(10000000, 3) - gammainc(10000000, 2)   # Bad
    0.0
    >>> gammainc(10000000, 2, 3)   # Good
    1.755146243738946045873491e+4771204
    >>> gammainc(2, 0, 100000001) - gammainc(2, 0, 100000000)   # Bad
    0.0
    >>> gammainc(2, 100000000, 100000001)   # Good
    4.078258353474186729184421e-43429441

The incomplete gamma functions satisfy simple recurrence
relations::

    >>> mp.dps = 25
    >>> z, a = mpf(3.5), mpf(2)
    >>> gammainc(z+1, a); z*gammainc(z,a) + a**z*exp(-a)
    10.60130296933533459267329
    10.60130296933533459267329
    >>> gammainc(z+1,0,a); z*gammainc(z,0,a) - a**z*exp(-a)
    1.030425427232114336470932
    1.030425427232114336470932

Evaluation at integers and poles::

    >>> gammainc(-3, -4, -5)
    (-0.2214577048967798566234192 + 0.0j)
    >>> gammainc(-3, 0, 5)
    +inf

If `z` is an integer, the recurrence reduces the incomplete gamma
function to `P(a) \exp(-a) + Q(b) \exp(-b)` where `P` and
`Q` are polynomials::

    >>> gammainc(1, 2); exp(-2)
    0.1353352832366126918939995
    0.1353352832366126918939995
    >>> mp.dps = 50
    >>> identify(gammainc(6, 1, 2), ['exp(-1)', 'exp(-2)'])
    '(326*exp(-1) + (-872)*exp(-2))'

The incomplete gamma functions reduce to functions such as
the exponential integral Ei and the error function for special
arguments::

    >>> mp.dps = 25
    >>> gammainc(0, 4); -ei(-4)
    0.00377935240984890647887486
    0.00377935240984890647887486
    >>> gammainc(0.5, 0, 2); sqrt(pi)*erf(sqrt(2))
    1.691806732945198336509541
    1.691806732945198336509541

"""

erf = r"""
Computes the error function, `\mathrm{erf}(x)`. The error
function is the normalized antiderivative of the Gaussian function
`\exp(-t^2)`. More precisely,

.. math::

  \mathrm{erf}(x) = \frac{2}{\sqrt \pi} \int_0^x \exp(-t^2) \,dt

**Basic examples**

Simple values and limits include::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> erf(0)
    0.0
    >>> erf(1)
    0.842700792949715
    >>> erf(-1)
    -0.842700792949715
    >>> erf(inf)
    1.0
    >>> erf(-inf)
    -1.0

For large real `x`, `\mathrm{erf}(x)` approaches 1 very
rapidly::

    >>> erf(3)
    0.999977909503001
    >>> erf(5)
    0.999999999998463

The error function is an odd function::

    >>> nprint(chop(taylor(erf, 0, 5)))
    [0.0, 1.12838, 0.0, -0.376126, 0.0, 0.112838]

:func:`erf` implements arbitrary-precision evaluation and
supports complex numbers::

    >>> mp.dps = 50
    >>> erf(0.5)
    0.52049987781304653768274665389196452873645157575796
    >>> mp.dps = 25
    >>> erf(1+j)
    (1.316151281697947644880271 + 0.1904534692378346862841089j)

Evaluation is supported for large arguments::

    >>> mp.dps = 25
    >>> erf('1e1000')
    1.0
    >>> erf('-1e1000')
    -1.0
    >>> erf('1e-1000')
    1.128379167095512573896159e-1000
    >>> erf('1e7j')
    (0.0 + 8.593897639029319267398803e+43429448190317j)
    >>> erf('1e7+1e7j')
    (0.9999999858172446172631323 + 3.728805278735270407053139e-8j)

**Related functions**

See also :func:`erfc`, which is more accurate for large `x`,
and :func:`erfi` which gives the antiderivative of
`\exp(t^2)`.

The Fresnel integrals :func:`fresnels` and :func:`fresnelc`
are also related to the error function.
"""

erfc = r"""
Computes the complementary error function,
`\mathrm{erfc}(x) = 1-\mathrm{erf}(x)`.
This function avoids cancellation that occurs when naively
computing the complementary error function as ``1-erf(x)``::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> 1 - erf(10)
    0.0
    >>> erfc(10)
    2.08848758376254e-45

:func:`erfc` works accurately even for ludicrously large
arguments::

    >>> erfc(10**10)
    4.3504398860243e-43429448190325182776

Complex arguments are supported::

    >>> erfc(500+50j)
    (1.19739830969552e-107492 + 1.46072418957528e-107491j)

"""


erfi = r"""
Computes the imaginary error function, `\mathrm{erfi}(x)`.
The imaginary error function is defined in analogy with the
error function, but with a positive sign in the integrand:

.. math ::

  \mathrm{erfi}(x) = \frac{2}{\sqrt \pi} \int_0^x \exp(t^2) \,dt

Whereas the error function rapidly converges to 1 as `x` grows,
the imaginary error function rapidly diverges to infinity.
The functions are related as
`\mathrm{erfi}(x) = -i\,\mathrm{erf}(ix)` for all complex
numbers `x`.

**Examples**

Basic values and limits::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> erfi(0)
    0.0
    >>> erfi(1)
    1.65042575879754
    >>> erfi(-1)
    -1.65042575879754
    >>> erfi(inf)
    +inf
    >>> erfi(-inf)
    -inf

Note the symmetry between erf and erfi::

    >>> erfi(3j)
    (0.0 + 0.999977909503001j)
    >>> erf(3)
    0.999977909503001
    >>> erf(1+2j)
    (-0.536643565778565 - 5.04914370344703j)
    >>> erfi(2+1j)
    (-5.04914370344703 - 0.536643565778565j)

Large arguments are supported::

    >>> erfi(1000)
    1.71130938718796e+434291
    >>> erfi(10**10)
    7.3167287567024e+43429448190325182754
    >>> erfi(-10**10)
    -7.3167287567024e+43429448190325182754
    >>> erfi(1000-500j)
    (2.49895233563961e+325717 + 2.6846779342253e+325717j)
    >>> erfi(100000j)
    (0.0 + 1.0j)
    >>> erfi(-100000j)
    (0.0 - 1.0j)


"""

erfinv = r"""
Computes the inverse error function, satisfying

.. math ::

    \mathrm{erf}(\mathrm{erfinv}(x)) =
    \mathrm{erfinv}(\mathrm{erf}(x)) = x.

This function is defined only for `-1 \le x \le 1`.

**Examples**

Special values include::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> erfinv(0)
    0.0
    >>> erfinv(1)
    +inf
    >>> erfinv(-1)
    -inf

The domain is limited to the standard interval::

    >>> erfinv(2)
    Traceback (most recent call last):
      ...
    ValueError: erfinv(x) is defined only for -1 <= x <= 1

It is simple to check that :func:`erfinv` computes inverse values of
:func:`erf` as promised::

    >>> erf(erfinv(0.75))
    0.75
    >>> erf(erfinv(-0.995))
    -0.995

:func:`erfinv` supports arbitrary-precision evaluation::

    >>> mp.dps = 50
    >>> x = erf(2)
    >>> x
    0.99532226501895273416206925636725292861089179704006
    >>> erfinv(x)
    2.0

A definite integral involving the inverse error function::

    >>> mp.dps = 15
    >>> quad(erfinv, [0, 1])
    0.564189583547756
    >>> 1/sqrt(pi)
    0.564189583547756

The inverse error function can be used to generate random numbers
with a Gaussian distribution (although this is a relatively
inefficient algorithm)::

    >>> nprint([erfinv(2*rand()-1) for n in range(6)]) # doctest: +SKIP
    [-0.586747, 1.10233, -0.376796, 0.926037, -0.708142, -0.732012]

"""

npdf = r"""
``npdf(x, mu=0, sigma=1)`` evaluates the probability density
function of a normal distribution with mean value `\mu`
and variance `\sigma^2`.

Elementary properties of the probability distribution can
be verified using numerical integration::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> quad(npdf, [-inf, inf])
    1.0
    >>> quad(lambda x: npdf(x, 3), [3, inf])
    0.5
    >>> quad(lambda x: npdf(x, 3, 2), [3, inf])
    0.5

See also :func:`ncdf`, which gives the cumulative
distribution.
"""

ncdf = r"""
``ncdf(x, mu=0, sigma=1)`` evaluates the cumulative distribution
function of a normal distribution with mean value `\mu`
and variance `\sigma^2`.

See also :func:`npdf`, which gives the probability density.

Elementary properties include::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> ncdf(pi, mu=pi)
    0.5
    >>> ncdf(-inf)
    0.0
    >>> ncdf(+inf)
    1.0

The cumulative distribution is the integral of the density
function having identical mu and sigma::

    >>> mp.dps = 15
    >>> diff(ncdf, 2)
    0.053990966513188
    >>> npdf(2)
    0.053990966513188
    >>> diff(lambda x: ncdf(x, 1, 0.5), 0)
    0.107981933026376
    >>> npdf(0, 1, 0.5)
    0.107981933026376
"""

expint = r"""
:func:`expint(n,z)` gives the generalized exponential integral
or En-function,

.. math ::

    \mathrm{E}_n(z) = \int_1^{\infty} \frac{e^{-zt}}{t^n} dt,

where `n` and `z` may both be complex numbers. The case with `n = 1` is
also given by :func:`e1`.

**Examples**

Evaluation at real and complex arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> expint(1, 6.25)
    0.0002704758872637179088496194
    >>> expint(-3, 2+3j)
    (0.00299658467335472929656159 + 0.06100816202125885450319632j)
    >>> expint(2+3j, 4-5j)
    (0.001803529474663565056945248 - 0.002235061547756185403349091j)

At negative integer values of `n`, `E_n(z)` reduces to a
rational-exponential function::

    >>> f = lambda n, z: fac(n)*sum(z**k/fac(k-1) for k in range(1,n+2))/\
    ...     exp(z)/z**(n+2)
    >>> n = 3
    >>> z = 1/pi
    >>> expint(-n,z)
    584.2604820613019908668219
    >>> f(n,z)
    584.2604820613019908668219
    >>> n = 5
    >>> expint(-n,z)
    115366.5762594725451811138
    >>> f(n,z)
    115366.5762594725451811138
"""

e1 = r"""
Computes the exponential integral `\mathrm{E}_1(z)`, given by

.. math ::

    \mathrm{E}_1(z) = \int_z^{\infty} \frac{e^{-t}}{t} dt.

This is equivalent to :func:`expint` with `n = 1`.

**Examples**

Two ways to evaluate this function::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> e1(6.25)
    0.0002704758872637179088496194
    >>> expint(1,6.25)
    0.0002704758872637179088496194

The E1-function is essentially the same as the Ei-function (:func:`ei`)
with negated argument, except for an imaginary branch cut term::

    >>> e1(2.5)
    0.02491491787026973549562801
    >>> -ei(-2.5)
    0.02491491787026973549562801
    >>> e1(-2.5)
    (-7.073765894578600711923552 - 3.141592653589793238462643j)
    >>> -ei(2.5)
    -7.073765894578600711923552

"""

ei = r"""
Computes the exponential integral or Ei-function, `\mathrm{Ei}(x)`.
The exponential integral is defined as

.. math ::

  \mathrm{Ei}(x) = \int_{-\infty\,}^x \frac{e^t}{t} \, dt.

When the integration range includes `t = 0`, the exponential
integral is interpreted as providing the Cauchy principal value.

For real `x`, the Ei-function behaves roughly like
`\mathrm{Ei}(x) \approx \exp(x) + \log(|x|)`.

The Ei-function is related to the more general family of exponential
integral functions denoted by `E_n`, which are available as :func:`expint`.

**Basic examples**

Some basic values and limits are::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> ei(0)
    -inf
    >>> ei(1)
    1.89511781635594
    >>> ei(inf)
    +inf
    >>> ei(-inf)
    0.0

For `x < 0`, the defining integral can be evaluated
numerically as a reference::

    >>> ei(-4)
    -0.00377935240984891
    >>> quad(lambda t: exp(t)/t, [-inf, -4])
    -0.00377935240984891

:func:`ei` supports complex arguments and arbitrary
precision evaluation::

    >>> mp.dps = 50
    >>> ei(pi)
    10.928374389331410348638445906907535171566338835056
    >>> mp.dps = 25
    >>> ei(3+4j)
    (-4.154091651642689822535359 + 4.294418620024357476985535j)

**Related functions**

The exponential integral is closely related to the logarithmic
integral. See :func:`li` for additional information.

The exponential integral is related to the hyperbolic
and trigonometric integrals (see :func:`chi`, :func:`shi`,
:func:`ci`, :func:`si`) similarly to how the ordinary
exponential function is related to the hyperbolic and
trigonometric functions::

    >>> mp.dps = 15
    >>> ei(3)
    9.93383257062542
    >>> chi(3) + shi(3)
    9.93383257062542
    >>> ci(3j) - j*si(3j) - pi*j/2
    (9.93383257062542 + 0.0j)

Beware that logarithmic corrections, as in the last example
above, are required to obtain the correct branch in general.
For details, see [1].

The exponential integral is also a special case of the
hypergeometric function `\,_2F_2`::

    >>> z = 0.6
    >>> z*hyper([1,1],[2,2],z) + (ln(z)-ln(1/z))/2 + euler
    0.769881289937359
    >>> ei(z)
    0.769881289937359

**References**

1. Relations between Ei and other functions:
   http://functions.wolfram.com/GammaBetaErf/ExpIntegralEi/27/01/

2. Abramowitz & Stegun, section 5:
   http://www.math.sfu.ca/~cbm/aands/page_228.htm

3. Asymptotic expansion for Ei:
   http://mathworld.wolfram.com/En-Function.html
"""

li = r"""
Computes the logarithmic integral or li-function
`\mathrm{li}(x)`, defined by

.. math ::

    \mathrm{li}(x) = \int_0^x \frac{1}{\log t} \, dt

The logarithmic integral has a singularity at `x = 1`.

Note that there is a second logarithmic integral, the Li
function, defined by

.. math ::

    \mathrm{Li}(x) = \int_2^x \frac{1}{\log t} \, dt

This "offset logarithmic integral" can be computed via
:func:`li` using the simple identity
`\mathrm{Li}(x) = \mathrm{li}(x) - \mathrm{li}(2)`.

The logarithmic integral should also not be confused with
the polylogarithm (also denoted by Li), which is implemented
as :func:`polylog`.

**Examples**

Some basic values and limits::

    >>> from mpmath import *
    >>> mp.dps = 30; mp.pretty = True
    >>> li(0)
    0.0
    >>> li(1)
    -inf
    >>> li(1)
    -inf
    >>> li(2)
    1.04516378011749278484458888919
    >>> findroot(li, 2)
    1.45136923488338105028396848589
    >>> li(inf)
    +inf

The logarithmic integral can be evaluated for arbitrary
complex arguments::

    >>> mp.dps = 20
    >>> li(3+4j)
    (3.1343755504645775265 + 2.6769247817778742392j)

The logarithmic integral is related to the exponential integral::

    >>> ei(log(3))
    2.1635885946671919729
    >>> li(3)
    2.1635885946671919729

The logarithmic integral grows like `O(x/\log(x))`::

    >>> mp.dps = 15
    >>> x = 10**100
    >>> x/log(x)
    4.34294481903252e+97
    >>> li(x)
    4.3619719871407e+97

The prime number theorem states that the number of primes less
than `x` is asymptotic to `\mathrm{li}(x)`. For example,
it is known that there are exactly 1,925,320,391,606,803,968,923
prime numbers less than `10^{23}` [1]. The logarithmic integral
provides a very accurate estimate::

    >>> li(2) + li(10**23)
    1.92532039161405e+21

A definite integral is::

    >>> quad(li, [0, 1])
    -0.693147180559945
    >>> -ln(2)
    -0.693147180559945

**References**

1. http://mathworld.wolfram.com/PrimeCountingFunction.html

2. http://mathworld.wolfram.com/LogarithmicIntegral.html

"""

ci = r"""
Computes the cosine integral,

.. math ::

    \mathrm{Ci}(x) = -\int_x^{\infty} \frac{\cos t}{t}\,dt
    = \gamma + \log x + \int_0^x \frac{\cos t - 1}{t}\,dt

**Examples**

Some values and limits::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> ci(0)
    -inf
    >>> ci(1)
    0.3374039229009681346626462
    >>> ci(pi)
    0.07366791204642548599010096
    >>> ci(inf)
    0.0
    >>> ci(-inf)
    (0.0 + 3.141592653589793238462643j)
    >>> ci(2+3j)
    (1.408292501520849518759125 - 2.983617742029605093121118j)

The cosine integral behaves roughly like the sinc function
(see :func:`sinc`) for large real `x`::

    >>> ci(10**10)
    -4.875060251748226537857298e-11
    >>> sinc(10**10)
    -4.875060250875106915277943e-11
    >>> chop(limit(ci, inf))
    0.0

It has infinitely many roots on the positive real axis::

    >>> findroot(ci, 1)
    0.6165054856207162337971104
    >>> findroot(ci, 2)
    3.384180422551186426397851

We can evaluate the defining integral as a reference::

    >>> mp.dps = 15
    >>> -quadosc(lambda t: cos(t)/t, [5, inf], omega=1)
    -0.190029749656644
    >>> ci(5)
    -0.190029749656644

Some infinite series can be evaluated using the
cosine integral::

    >>> nsum(lambda k: (-1)**k/(fac(2*k)*(2*k)), [1,inf])
    -0.239811742000565
    >>> ci(1) - euler
    -0.239811742000565

"""

si = r"""
Computes the sine integral,

.. math ::

    \mathrm{Si}(x) = \int_0^x \frac{\sin t}{t}\,dt.

The sine integral is thus the antiderivative of the sinc
function (see :func:`sinc`).

**Examples**

Some values and limits::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> si(0)
    0.0
    >>> si(1)
    0.9460830703671830149413533
    >>> si(-1)
    -0.9460830703671830149413533
    >>> si(pi)
    1.851937051982466170361053
    >>> si(inf)
    1.570796326794896619231322
    >>> si(-inf)
    -1.570796326794896619231322
    >>> si(2+3j)
    (4.547513889562289219853204 + 1.399196580646054789459839j)

The sine integral approaches `\pi/2` for large real `x`::

    >>> si(10**10)
    1.570796326707584656968511
    >>> pi/2
    1.570796326794896619231322

We can evaluate the defining integral as a reference::

    >>> mp.dps = 15
    >>> quad(sinc, [0, 5])
    1.54993124494467
    >>> si(5)
    1.54993124494467

Some infinite series can be evaluated using the
sine integral::

    >>> nsum(lambda k: (-1)**k/(fac(2*k+1)*(2*k+1)), [0,inf])
    0.946083070367183
    >>> si(1)
    0.946083070367183

"""

chi = r"""
Computes the hyperbolic cosine integral, defined
in analogy with the cosine integral (see :func:`ci`) as

.. math ::

    \mathrm{Chi}(x) = -\int_x^{\infty} \frac{\cosh t}{t}\,dt
    = \gamma + \log x + \int_0^x \frac{\cosh t - 1}{t}\,dt

Some values and limits::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> chi(0)
    -inf
    >>> chi(1)
    0.8378669409802082408946786
    >>> chi(inf)
    +inf
    >>> findroot(chi, 0.5)
    0.5238225713898644064509583
    >>> chi(2+3j)
    (-0.1683628683277204662429321 + 2.625115880451325002151688j)

"""

shi = r"""
Computes the hyperbolic sine integral, defined
in analogy with the sine integral (see :func:`si`) as

.. math ::

    \mathrm{Shi}(x) = \int_0^x \frac{\sinh t}{t}\,dt.

Some values and limits::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> shi(0)
    0.0
    >>> shi(1)
    1.057250875375728514571842
    >>> shi(-1)
    -1.057250875375728514571842
    >>> shi(inf)
    +inf
    >>> shi(2+3j)
    (-0.1931890762719198291678095 + 2.645432555362369624818525j)

"""

fresnels = r"""
Computes the Fresnel sine integral

.. math ::

    S(x) = \int_0^x \sin\left(\frac{\pi t^2}{2}\right) \,dt

Note that some sources define this function
without the normalization factor `\pi/2`.

**Examples**

Some basic values and limits::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> fresnels(0)
    0.0
    >>> fresnels(inf)
    0.5
    >>> fresnels(-inf)
    -0.5
    >>> fresnels(1)
    0.4382591473903547660767567
    >>> fresnels(1+2j)
    (36.72546488399143842838788 + 15.58775110440458732748279j)

Comparing with the definition::

    >>> fresnels(3)
    0.4963129989673750360976123
    >>> quad(lambda t: sin(pi*t**2/2), [0,3])
    0.4963129989673750360976123
"""

fresnelc = r"""
Computes the Fresnel cosine integral

.. math ::

    C(x) = \int_0^x \cos\left(\frac{\pi t^2}{2}\right) \,dt

Note that some sources define this function
without the normalization factor `\pi/2`.

**Examples**

Some basic values and limits::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> fresnelc(0)
    0.0
    >>> fresnelc(inf)
    0.5
    >>> fresnelc(-inf)
    -0.5
    >>> fresnelc(1)
    0.7798934003768228294742064
    >>> fresnelc(1+2j)
    (16.08787137412548041729489 - 36.22568799288165021578758j)

Comparing with the definition::

    >>> fresnelc(3)
    0.6057207892976856295561611
    >>> quad(lambda t: cos(pi*t**2/2), [0,3])
    0.6057207892976856295561611
"""

airyai = r"""
Computes the Airy function `\mathrm{Ai}(x)`, which is
a solution of the Airy differential equation `y''-xy=0`.
The Ai-function behaves roughly like a slowly decaying
sine wave for `x < 0` and like a decreasing exponential for
`x > 0`.

Limits and values include::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> airyai(0), 1/(3**(2/3.)*gamma(2/3.))
    (0.355028053887817, 0.355028053887817)
    >>> airyai(1)
    0.135292416312881
    >>> airyai(-1)
    0.535560883292352
    >>> airyai(inf)
    0.0
    >>> airyai(-inf)
    0.0

Evaluation is supported for large arguments::

    >>> airyai(-100)
    0.176753393239553
    >>> airyai(100)
    2.63448215208818e-291
    >>> airyai(50+50j)
    (-5.31790195707456e-68 - 1.16358800377071e-67j)
    >>> airyai(-50+50j)
    (1.04124253736317e+158 + 3.3475255449236e+157j)

Huge arguments are also fine::

    >>> airyai(10**10)
    1.16223597829874e-289529654602171
    >>> airyai(-10**10)
    0.000173620644815282
    >>> airyai(10**10*(1+j))
    (5.71150868372136e-186339621747698 + 1.86724550696231e-186339621747697j)

The first negative root of the Ai function is::

    >>> findroot(airyai, -2)
    -2.33810741045977

We can verify the differential equation::

    >>> for x in [-3.4, 0, 2.5, 1+2j]:
    ...     print abs(diff(airyai, x, 2) - x*airyai(x)) < eps
    ...
    True
    True
    True
    True

The Taylor series expansion around `x = 0` starts with
the following coefficients (note that every third term
is zero)::

    >>> nprint(chop(taylor(airyai, 0, 5)))
    [0.355028, -0.258819, 0.0, 5.91713e-2, -2.15683e-2, 0.0]

The Airy functions are a special case of Bessel functions.
For `x < 0`, we have::

    >>> x = 3
    >>> airyai(-x)
    -0.378814293677658
    >>> p = 2*(x**1.5)/3
    >>> sqrt(x)*(besselj(1/3.,p) + besselj(-1/3.,p))/3
    -0.378814293677658

"""

airybi = r"""
Computes the Airy function `\mathrm{Bi}(x)`, which is
a solution of the Airy differential equation `y''-xy=0`.
The Bi-function behaves roughly like a slowly decaying
sine wave for `x < 0` and like an increasing exponential
for `x > 0`.

Limits and values include::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> airybi(0), 1/(3**(1/6.)*gamma(2/3.))
    (0.614926627446001, 0.614926627446001)
    >>> airybi(1)
    1.20742359495287
    >>> airybi(-1)
    0.103997389496945
    >>> airybi(inf)
    +inf
    >>> airybi(-inf)
    0.0

Evaluation is supported for large arguments::

    >>> airybi(-100)
    0.0242738876801601
    >>> airybi(100)
    6.0412239966702e+288
    >>> airybi(50+50j)
    (-5.32207626732144e+63 + 1.47845029116524e+65j)
    >>> airybi(-50+50j)
    (-3.3475255449236e+157 + 1.04124253736317e+158j)

Huge arguments are also fine::

    >>> mp.dps = 15
    >>> airybi(10**10)
    1.36938578794354e+289529654602165
    >>> airybi(-10**10)
    0.00177565614169293
    >>> airybi(10**10*(1+j))
    (-6.5599559310962e+186339621747689 - 6.82246272698136e+186339621747690j)

The first negative root of the Bi function is::

    >>> findroot(airybi, -1)
    -1.17371322270913

We can verify the differential equation::

    >>> for x in [-3.4, 0, 2.5, 1+2j]:
    ...     print abs(diff(airybi, x, 2) - x*airybi(x)) < eps
    ...
    True
    True
    True
    True

The Taylor series expansion around `x = 0` starts with
the following coefficients (note that every third term
is zero)::

    >>> nprint(chop(taylor(airybi, 0, 5)))
    [0.614927, 0.448288, 0.0, 0.102488, 3.73574e-2, 0.0]

The Airy functions are a special case of Bessel functions.
For `x < 0`, we have::

    >>> x = 3
    >>> airybi(-x)
    -0.198289626374927
    >>> p = 2*(x**1.5)/3
    >>> sqrt(x/3)*(besselj(-1/3.,p) - besselj(1/3.,p))
    -0.198289626374926
"""

ellipk = r"""
Evaluates the complete elliptic integral of the first kind,
`K(m)`, defined by

.. math ::

    K(m) = \int_0^{\pi/2} \frac{1}{\sqrt{1-m \sin^2 t}} dt.

Note that the argument is the parameter `m = k^2`,
not the modulus `k` which is sometimes used.

Alternatively, in terms of a hypergeometric function,
we have:

.. math ::

    K(m) = \frac{\pi}{2} \,_2F_1(1/2, 1/2, 1, m)

**Examples**

Values and limits include::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> ellipk(0)
    1.570796326794896619231322
    >>> ellipk(inf)
    (0.0 + 0.0j)
    >>> ellipk(-inf)
    0.0
    >>> ellipk(1)
    +inf
    >>> ellipk(-1)
    1.31102877714605990523242
    >>> ellipk(2)
    (1.31102877714605990523242 - 1.31102877714605990523242j)

Verifying the defining integral and hypergeometric
representation::

    >>> ellipk(0.5)
    1.85407467730137191843385
    >>> quad(lambda t: (1-0.5*sin(t)**2)**-0.5, [0, pi/2])
    1.85407467730137191843385
    >>> pi/2*hyp2f1(0.5,0.5,1,0.5)
    1.85407467730137191843385

Evaluation is supported for arbitrary complex `m`::

    >>> ellipk(3+4j)
    (0.9111955638049650086562171 + 0.6313342832413452438845091j)

A definite integral::

    >>> quad(ellipk, [0, 1])
    2.0
"""

ellipe = r"""
Evaluates the complete elliptic integral of the second kind,
`E(m)`, defined by

.. math ::

    E(m) = \int_0^{\pi/2} \sqrt{1-m \sin^2 t} dt.

Note that the argument is the parameter `m = k^2`,
not the modulus `k` which is sometimes used.

Alternatively, in terms of a hypergeometric function,
we have:

.. math ::

    E(m) = \frac{\pi}{2} \,_2F_1(1/2, -1/2, 1, m)

**Examples**

Basic values and limits::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> ellipe(0)
    1.570796326794896619231322
    >>> ellipe(1)
    1.0
    >>> ellipe(-1)
    1.910098894513856008952381
    >>> ellipe(2)
    (0.5990701173677961037199612 + 0.5990701173677961037199612j)
    >>> ellipe(inf)
    (0.0 + +infj)
    >>> ellipe(-inf)
    +inf

Verifying the defining integral and hypergeometric
representation::

    >>> ellipe(0.5)
    1.350643881047675502520175
    >>> quad(lambda t: sqrt(1-0.5*sin(t)**2), [0, pi/2])
    1.350643881047675502520175
    >>> pi/2*hyp2f1(0.5,-0.5,1,0.5)
    1.350643881047675502520175

Evaluation is supported for arbitrary complex `m`::

    >>> ellipe(0.5+0.25j)
    (1.360868682163129682716687 - 0.1238733442561786843557315j)
    >>> ellipe(3+4j)
    (1.499553520933346954333612 - 1.577879007912758274533309j)

A definite integral::

    >>> quad(ellipe, [0,1])
    1.333333333333333333333333

"""

agm = r"""
``agm(a, b)`` computes the arithmetic-geometric mean of `a` and
`b`, defined as the limit of the following iteration:

.. math ::

    a_0 = a

    b_0 = b

    a_{n+1} = \frac{a_n+b_n}{2}

    b_{n+1} = \sqrt{a_n b_n}

This function can be called with a single argument, computing
`\mathrm{agm}(a,1) = \mathrm{agm}(1,a)`.

**Examples**

It is a well-known theorem that the geometric mean of
two distinct positive numbers is less than the arithmetic
mean. It follows that the arithmetic-geometric mean lies
between the two means::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> a = mpf(3)
    >>> b = mpf(4)
    >>> sqrt(a*b)
    3.46410161513775
    >>> agm(a,b)
    3.48202767635957
    >>> (a+b)/2
    3.5

The arithmetic-geometric mean is scale-invariant::

    >>> agm(10*e, 10*pi)
    29.261085515723
    >>> 10*agm(e, pi)
    29.261085515723

As an order-of-magnitude estimate, `\mathrm{agm}(1,x) \approx x`
for large `x`::

    >>> agm(10**10)
    643448704.760133
    >>> agm(10**50)
    1.34814309345871e+48

For tiny `x`, `\mathrm{agm}(1,x) \approx -\pi/(2 \log(x/4))`::

    >>> agm('0.01')
    0.262166887202249
    >>> -pi/2/log('0.0025')
    0.262172347753122

The arithmetic-geometric mean can also be computed for complex
numbers::

    >>> agm(3, 2+j)
    (2.51055133276184 + 0.547394054060638j)

The AGM iteration converges very quickly (each step doubles
the number of correct digits), so :func:`agm` supports efficient
high-precision evaluation::

    >>> mp.dps = 10000
    >>> a = agm(1,2)
    >>> str(a)[-10:]
    '1679581912'

**Mathematical relations**

The arithmetic-geometric mean may be used to evaluate the
following two parametric definite integrals:

.. math ::

  I_1 = \int_0^{\infty}
    \frac{1}{\sqrt{(x^2+a^2)(x^2+b^2)}} \,dx

  I_2 = \int_0^{\pi/2}
    \frac{1}{\sqrt{a^2 \cos^2(x) + b^2 \sin^2(x)}} \,dx

We have::

    >>> mp.dps = 15
    >>> a = 3
    >>> b = 4
    >>> f1 = lambda x: ((x**2+a**2)*(x**2+b**2))**-0.5
    >>> f2 = lambda x: ((a*cos(x))**2 + (b*sin(x))**2)**-0.5
    >>> quad(f1, [0, inf])
    0.451115405388492
    >>> quad(f2, [0, pi/2])
    0.451115405388492
    >>> pi/(2*agm(a,b))
    0.451115405388492

A formula for `\Gamma(1/4)`::

    >>> gamma(0.25)
    3.62560990822191
    >>> sqrt(2*sqrt(2*pi**3)/agm(1,sqrt(2)))
    3.62560990822191

**Possible issues**

The branch cut chosen for complex `a` and `b` is somewhat
arbitrary.

"""

gegenbauer = r"""
Evaluates the Gegenbauer polynomial, or ultraspherical polynomial,

.. math ::

    C_n^{(a)}(z) = {n+2a-1 \choose n} \,_2F_1\left(-n, n+2a;
        a+\frac{1}{2}; \frac{1}{2}(1-z)\right).

When `n` is a nonnegative integer, this formula gives a polynomial
in `z` of degree `n`, but all parameters are permitted to be
complex numbers. With `a = 1/2`, the Gegenbauer polynomial
reduces to a Legendre polynomial.

**Examples**

Evaluation for arbitrary arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> gegenbauer(3, 0.5, -10)
    -2485.0
    >>> gegenbauer(1000, 10, 100)
    3.012757178975667428359374e+2322
    >>> gegenbauer(2+3j, -0.75, -1000j)
    (-5038991.358609026523401901 + 9414549.285447104177860806j)

Evaluation at negative integer orders::

    >>> gegenbauer(-4, 2, 1.75)
    -1.0
    >>> chop(gegenbauer(-4, 3, 1.75))
    0.0
    >>> chop(gegenbauer(-4, 2j, 1.75))
    0.0
    >>> gegenbauer(-7, 0.5, 3)
    8989.0

The Gegenbauer polynomials solve the differential equation::

    >>> n, a = 4.5, 1+2j
    >>> f = lambda z: gegenbauer(n, a, z)
    >>> for z in [0, 0.75, -0.5j]:
    ...     chop((1-z**2)*diff(f,z,2) - (2*a+1)*z*diff(f,z) + n*(n+2*a)*f(z))
    ...
    0.0
    0.0
    0.0

The Gegenbauer polynomials have generating function
`(1-2zt+t^2)^{-a}`::

    >>> a, z = 2.5, 1
    >>> taylor(lambda t: (1-2*z*t+t**2)**(-a), 0, 3)
    [1.0, 5.0, 15.0, 35.0]
    >>> [gegenbauer(n,a,z) for n in range(4)]
    [1.0, 5.0, 15.0, 35.0]

The Gegenbauer polynomials are orthogonal on `[-1, 1]` with respect
to the weight `(1-z^2)^{a-\frac{1}{2}}`::

    >>> a, n, m = 2.5, 4, 5
    >>> Cn = lambda z: gegenbauer(n, a, z)
    >>> Cm = lambda z: gegenbauer(m, a, z)
    >>> chop(quad(lambda z: Cn(z)*Cm(z)*(1-z**2)*(a-0.5), [-1, 1]))
    0.0
"""

laguerre = r"""
Gives the generalized (associated) Laguerre polynomial, defined by

.. math ::

    L_n^a(z) = \frac{\Gamma(n+b+1)}{\Gamma(b+1) \Gamma(n+1)}
        \,_1F_1(-n, a+1, z).

With `a = 0` and `n` a nonnegative integer, this reduces to an ordinary
Laguerre polynomial, the sequence of which begins
`L_0(z) = 1, L_1(z) = 1-z, L_2(z) = z^2-2z+1, \ldots`.

The Laguerre polynomials are orthogonal with respect to the weight
`z^a e^{-z}` on `[0, \infty)`.

**Examples**

Evaluation for arbitrary arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> laguerre(5, 0, 0.25)
    0.03726399739583333333333333
    >>> laguerre(1+j, 0.5, 2+3j)
    (4.474921610704496808379097 - 11.02058050372068958069241j)
    >>> laguerre(2, 0, 10000)
    49980001.0
    >>> laguerre(2.5, 0, 10000)
    -9.327764910194842158583189e+4328

The first few Laguerre polynomials, normalized to have integer
coefficients::

    >>> for n in range(7):
    ...     chop(taylor(lambda z: fac(n)*laguerre(n, 0, z), 0, n))
    ...
    [1.0]
    [1.0, -1.0]
    [2.0, -4.0, 1.0]
    [6.0, -18.0, 9.0, -1.0]
    [24.0, -96.0, 72.0, -16.0, 1.0]
    [120.0, -600.0, 600.0, -200.0, 25.0, -1.0]
    [720.0, -4320.0, 5400.0, -2400.0, 450.0, -36.0, 1.0]

Verifying orthogonality::

    >>> Lm = lambda t: laguerre(m,a,t)
    >>> Ln = lambda t: laguerre(n,a,t)
    >>> a, n, m = 2.5, 2, 3
    >>> chop(quad(lambda t: exp(-t)*t**a*Lm(t)*Ln(t), [0,inf]))
    0.0


"""

hermite = r"""
Evaluates the Hermite polynomial `H_n(z)`, which may be defined using
the recurrence

.. math ::

    H_0(z) = 1

    H_1(z) = 2z

    H_{n+1} = 2z H_n(z) - 2n H_{n-1}(z).

The Hermite polynomials are orthogonal on `(-\infty, \infty)` with
respect to the weight `e^{-z^2}`. More generally, allowing arbitrary complex
values of `n`, the Hermite function `H_n(z)` is defined as

.. math ::

    H_n(z) = (2z)^n \,_2F_0\left(-\frac{n}{2}, \frac{1-n}{2},
        -\frac{1}{z^2}\right)

for `\Re{z} > 0`, or generally

.. math ::

    H_n(z) = 2^n \sqrt{\pi} \left(
        \frac{1}{\Gamma\left(\frac{1-n}{2}\right)}
        \,_1F_1\left(-\frac{n}{2}, \frac{1}{2}, z^2\right) -
        \frac{2z}{\Gamma\left(-\frac{n}{2}\right)}
        \,_1F_1\left(\frac{1-n}{2}, \frac{3}{2}, z^2\right)
    \right).

**Examples**

Evaluation for arbitrary arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> hermite(0, 10)
    1.0
    >>> hermite(1, 10); hermite(2, 10)
    20.0
    398.0
    >>> hermite(10000, 2)
    4.950440066552087387515653e+19334
    >>> hermite(3, -10**8)
    -7999999999999998800000000.0
    >>> hermite(-3, -10**8)
    1.675159751729877682920301e+4342944819032534
    >>> hermite(2+3j, -1+2j)
    (-0.076521306029935133894219 - 0.1084662449961914580276007j)

Coefficients of the first few Hermite polynomials are::

    >>> for n in range(7):
    ...     chop(taylor(lambda z: hermite(n, z), 0, n))
    ...
    [1.0]
    [0.0, 2.0]
    [-2.0, 0.0, 4.0]
    [0.0, -12.0, 0.0, 8.0]
    [12.0, 0.0, -48.0, 0.0, 16.0]
    [0.0, 120.0, 0.0, -160.0, 0.0, 32.0]
    [-120.0, 0.0, 720.0, 0.0, -480.0, 0.0, 64.0]

Values at `z = 0`::

    >>> for n in range(-5, 9):
    ...     hermite(n, 0)
    ...
    0.02769459142039868792653387
    0.08333333333333333333333333
    0.2215567313631895034122709
    0.5
    0.8862269254527580136490837
    1.0
    0.0
    -2.0
    0.0
    12.0
    0.0
    -120.0
    0.0
    1680.0

Hermite functions satisfy the differential equation::

    >>> n = 4
    >>> f = lambda z: hermite(n, z)
    >>> z = 1.5
    >>> chop(diff(f,z,2) - 2*z*diff(f,z) + 2*n*f(z))
    0.0

Verifying orthogonality::

    >>> chop(quad(lambda t: hermite(2,t)*hermite(4,t)*exp(-t**2), [-inf,inf]))
    0.0

"""

jacobi = r"""
``jacobi(n, a, b, x)`` evaluates the Jacobi polynomial
`P_n^{(a,b)}(x)`. The Jacobi polynomials are a special
case of the hypergeometric function `\,_2F_1` given by:

.. math ::

    P_n^{(a,b)}(x) = {n+a \choose n}
      \,_2F_1\left(-n,1+a+b+n,a+1,\frac{1-x}{2}\right).

Note that this definition generalizes to nonintegral values
of `n`. When `n` is an integer, the hypergeometric series
terminates after a finite number of terms, giving
a polynomial in `x`.

**Evaluation of Jacobi polynomials**

A special evaluation is `P_n^{(a,b)}(1) = {n+a \choose n}`::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> jacobi(4, 0.5, 0.25, 1)
    2.4609375
    >>> binomial(4+0.5, 4)
    2.4609375

A Jacobi polynomial of degree `n` is equal to its
Taylor polynomial of degree `n`. The explicit
coefficients of Jacobi polynomials can therefore
be recovered easily using :func:`taylor`::

    >>> for n in range(5):
    ...     nprint(taylor(lambda x: jacobi(n,1,2,x), 0, n))
    ...
    [1.0]
    [-0.5, 2.5]
    [-0.75, -1.5, 5.25]
    [0.5, -3.5, -3.5, 10.5]
    [0.625, 2.5, -11.25, -7.5, 20.625]

For nonintegral `n`, the Jacobi "polynomial" is no longer
a polynomial::

    >>> nprint(taylor(lambda x: jacobi(0.5,1,2,x), 0, 4))
    [0.309983, 1.84119, -1.26933, 1.26699, -1.34808]

**Orthogonality**

The Jacobi polynomials are orthogonal on the interval
`[-1, 1]` with respect to the weight function
`w(x) = (1-x)^a (1+x)^b`. That is,
`w(x) P_n^{(a,b)}(x) P_m^{(a,b)}(x)` integrates to
zero if `m \ne n` and to a nonzero number if `m = n`.

The orthogonality is easy to verify using numerical
quadrature::

    >>> P = jacobi
    >>> f = lambda x: (1-x)**a * (1+x)**b * P(m,a,b,x) * P(n,a,b,x)
    >>> a = 2
    >>> b = 3
    >>> m, n = 3, 4
    >>> chop(quad(f, [-1, 1]), 1)
    0.0
    >>> m, n = 4, 4
    >>> quad(f, [-1, 1])
    1.9047619047619

**Differential equation**

The Jacobi polynomials are solutions of the differential
equation

.. math ::

  (1-x^2) y'' + (b-a-(a+b+2)x) y' + n (n+a+b+1) y = 0.

We can verify that :func:`jacobi` approximately satisfies
this equation::

    >>> from mpmath import *
    >>> mp.dps = 15
    >>> a = 2.5
    >>> b = 4
    >>> n = 3
    >>> y = lambda x: jacobi(n,a,b,x)
    >>> x = pi
    >>> A0 = n*(n+a+b+1)*y(x)
    >>> A1 = (b-a-(a+b+2)*x)*diff(y,x)
    >>> A2 = (1-x**2)*diff(y,x,2)
    >>> nprint(A2 + A1 + A0, 1)
    4.0e-12

The difference of order `10^{-12}` is as close to zero as
it could be at 15-digit working precision, since the terms
are large::

    >>> A0, A1, A2
    (26560.2328981879, -21503.7641037294, -5056.46879445852)

"""

legendre = r"""
``legendre(n, x)`` evaluates the Legendre polynomial `P_n(x)`.
The Legendre polynomials are given by the formula

.. math ::

    P_n(x) = \frac{1}{2^n n!} \frac{d^n}{dx^n} (x^2 -1)^n.

Alternatively, they can be computed recursively using

.. math ::

    P_0(x) = 1

    P_1(x) = x

    (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x).

A third definition is in terms of the hypergeometric function
`\,_2F_1`, whereby they can be generalized to arbitrary `n`:

.. math ::

    P_n(x) = \,_2F_1\left(-n, n+1, 1, \frac{1-x}{2}\right)

**Basic evaluation**

The Legendre polynomials assume fixed values at the points
`x = -1` and `x = 1`::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> nprint([legendre(n, 1) for n in range(6)])
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    >>> nprint([legendre(n, -1) for n in range(6)])
    [1.0, -1.0, 1.0, -1.0, 1.0, -1.0]

The coefficients of Legendre polynomials can be recovered
using degree-`n` Taylor expansion::

    >>> for n in range(5):
    ...     nprint(chop(taylor(lambda x: legendre(n, x), 0, n)))
    ...
    [1.0]
    [0.0, 1.0]
    [-0.5, 0.0, 1.5]
    [0.0, -1.5, 0.0, 2.5]
    [0.375, 0.0, -3.75, 0.0, 4.375]

The roots of Legendre polynomials are located symmetrically
on the interval `[-1, 1]`::

    >>> for n in range(5):
    ...     nprint(polyroots(taylor(lambda x: legendre(n, x), 0, n)[::-1]))
    ...
    []
    [0.0]
    [-0.57735, 0.57735]
    [-0.774597, 0.0, 0.774597]
    [-0.861136, -0.339981, 0.339981, 0.861136]

An example of an evaluation for arbitrary `n`::

    >>> legendre(0.75, 2+4j)
    (1.94952805264875 + 2.1071073099422j)

**Orthogonality**

The Legendre polynomials are orthogonal on `[-1, 1]` with respect
to the trivial weight `w(x) = 1`. That is, `P_m(x) P_n(x)`
integrates to zero if `m \ne n` and to `2/(2n+1)` if `m = n`::

    >>> m, n = 3, 4
    >>> quad(lambda x: legendre(m,x)*legendre(n,x), [-1, 1])
    0.0
    >>> m, n = 4, 4
    >>> quad(lambda x: legendre(m,x)*legendre(n,x), [-1, 1])
    0.222222222222222

**Differential equation**

The Legendre polynomials satisfy the differential equation

.. math ::

    ((1-x^2) y')' + n(n+1) y' = 0.

We can verify this numerically::

    >>> n = 3.6
    >>> x = 0.73
    >>> P = legendre
    >>> A = diff(lambda t: (1-t**2)*diff(lambda u: P(n,u), t), x)
    >>> B = n*(n+1)*P(n,x)
    >>> nprint(A+B,1)
    9.0e-16

"""


legenp = r"""
Calculates the (associated) Legendre function of the first kind of
degree *n* and order *m*, `P_n^m(z)`. Taking `m = 0` gives the ordinary
Legendre function of the second kind, `P_n(z)`. The parameters may be
complex numbers.

In terms of the Gauss hypergeometric function, the (associated) Legendre
function is defined as

.. math ::

    P_n^m(z) = \frac{1}{\Gamma(1-m)} \frac{(1+z)^{m/2}}{(1-z)^{m/2}}
        \,_2F_1\left(-n, n+1, 1-m, \frac{1-z}{2}\right).

With *type=3* instead of *type=2*, the alternative
definition

.. math ::

    \hat{P}_n^m(z) = \frac{1}{\Gamma(1-m)} \frac{(z+1)^{m/2}}{(z-1)^{m/2}}
        \,_2F_1\left(-n, n+1, 1-m, \frac{1-z}{2}\right).

is used. These functions correspond respectively to ``LegendreP[n,m,2,z]``
and ``LegendreP[n,m,3,z]`` in Mathematica.

The general solution of the (associated) Legendre differential equation

.. math ::

    (1-z^2) f''(z) - 2zf'(z) + \left(n(n+1)-\frac{m^2}{1-z^2}\right)f(z) = 0

is given by `C_1 P_n^m(z) + C_2 Q_n^m(z)` for arbitrary constants
`C_1`, `C_2`, where `Q_n^m(z)` is a Legendre function of the
second kind as implemented by :func:`legenq`.

**Examples**

Evaluation for arbitrary parameters and arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> legenp(2, 0, 10); legendre(2, 10)
    149.5
    149.5
    >>> legenp(-2, 0.5, 2.5)
    (1.972260393822275434196053 - 1.972260393822275434196053j)
    >>> legenp(2+3j, 1-j, -0.5+4j)
    (-3.335677248386698208736542 - 5.663270217461022307645625j)
    >>> chop(legenp(3, 2, -1.5, type=2))
    28.125
    >>> chop(legenp(3, 2, -1.5, type=3))
    -28.125

Verifying the associated Legendre differential equation::

    >>> n, m = 2, -0.5
    >>> C1, C2 = 1, -3
    >>> f = lambda z: C1*legenp(n,m,z) + C2*legenq(n,m,z)
    >>> deq = lambda z: (1-z**2)*diff(f,z,2) - 2*z*diff(f,z) + \
    ...     (n*(n+1)-m**2/(1-z**2))*f(z)
    >>> for z in [0, 2, -1.5, 0.5+2j]:
    ...     chop(deq(mpmathify(z)))
    ...
    0.0
    0.0
    0.0
    0.0
"""

legenq = r"""
Calculates the (associated) Legendre function of the second kind of
degree *n* and order *m*, `Q_n^m(z)`. Taking `m = 0` gives the ordinary
Legendre function of the second kind, `Q_n(z)`. The parameters may
complex numbers.

The Legendre functions of the second kind give a second set of
solutions to the (associated) Legendre differential equation.
(See :func:`legenp`.)
Unlike the Legendre functions of the first kind, they are not
polynomials of `z` for integer `n`, `m` but rational or logarithmic
functions with poles at `z = \pm 1`.

There are various ways to define Legendre functions of
the second kind, giving rise to different complex structure.
A version can be selected using the *type* keyword argument.
The *type=2* and *type=3* functions are given respectively by

.. math ::

    Q_n^m(z) = \frac{\pi}{2 \sin(\pi m)}
        \left( \cos(\pi m) P_n^m(z) -
        \frac{\Gamma(1+m+n)}{\Gamma(1-m+n)} P_n^{-m}(z)\right)

    \hat{Q}_n^m(z) = \frac{\pi}{2 \sin(\pi m)} e^{\pi i m}
        \left( \hat{P}_n^m(z) -
        \frac{\Gamma(1+m+n)}{\Gamma(1-m+n)} \hat{P}_n^{-m}(z)\right)

where `P` and `\hat{P}` are the *type=2* and *type=3* Legendre functions
of the first kind. The formulas above should be understood as limits
when `m` is an integer.

These functions correspond to ``LegendreQ[n,m,2,z]`` (or ``LegendreQ[n,m,z]``)
and ``LegendreQ[n,m,3,z]`` in Mathematica. The *type=3* function
is essentially the same as the function defined in
Abramowitz & Stegun (eq. 8.1.3) but with `(z+1)^{m/2}(z-1)^{m/2}` instead
of `(z^2-1)^{m/2}`, giving slightly different branches.

**Examples**

Evaluation for arbitrary parameters and arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> legenq(2, 0, 0.5)
    -0.8186632680417568557122028
    >>> legenq(-1.5, -2, 2.5)
    (0.6655964618250228714288277 + 0.3937692045497259717762649j)
    >>> legenq(2-j, 3+4j, -6+5j)
    (-10001.95256487468541686564 - 6011.691337610097577791134j)

Different versions of the function::

    >>> legenq(2, 1, 0.5)
    0.7298060598018049369381857
    >>> legenq(2, 1, 1.5)
    (-7.902916572420817192300921 + 0.1998650072605976600724502j)
    >>> legenq(2, 1, 0.5, type=3)
    (2.040524284763495081918338 - 0.7298060598018049369381857j)
    >>> chop(legenq(2, 1, 1.5, type=3))
    -0.1998650072605976600724502

"""

chebyt = r"""
``chebyt(n, x)`` evaluates the Chebyshev polynomial of the first
kind `T_n(x)`, defined by the identity

.. math ::

    T_n(\cos x) = \cos(n x).

The Chebyshev polynomials of the first kind are a special
case of the Jacobi polynomials, and by extension of the
hypergeometric function `\,_2F_1`. They can thus also be
evaluated for nonintegral `n`.

**Basic evaluation**

The coefficients of the `n`-th polynomial can be recovered
using using degree-`n` Taylor expansion::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(5):
    ...     nprint(chop(taylor(lambda x: chebyt(n, x), 0, n)))
    ...
    [1.0]
    [0.0, 1.0]
    [-1.0, 0.0, 2.0]
    [0.0, -3.0, 0.0, 4.0]
    [1.0, 0.0, -8.0, 0.0, 8.0]

**Orthogonality**

The Chebyshev polynomials of the first kind are orthogonal
on the interval `[-1, 1]` with respect to the weight
function `w(x) = 1/\sqrt{1-x^2}`::

    >>> f = lambda x: chebyt(m,x)*chebyt(n,x)/sqrt(1-x**2)
    >>> m, n = 3, 4
    >>> nprint(quad(f, [-1, 1]),1)
    0.0
    >>> m, n = 4, 4
    >>> quad(f, [-1, 1])
    1.57079632596448

"""

chebyu = r"""
``chebyu(n, x)`` evaluates the Chebyshev polynomial of the second
kind `U_n(x)`, defined by the identity

.. math ::

    U_n(\cos x) = \frac{\sin((n+1)x)}{\sin(x)}.

The Chebyshev polynomials of the second kind are a special
case of the Jacobi polynomials, and by extension of the
hypergeometric function `\,_2F_1`. They can thus also be
evaluated for nonintegral `n`.

**Basic evaluation**

The coefficients of the `n`-th polynomial can be recovered
using using degree-`n` Taylor expansion::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(5):
    ...     nprint(chop(taylor(lambda x: chebyu(n, x), 0, n)))
    ...
    [1.0]
    [0.0, 2.0]
    [-1.0, 0.0, 4.0]
    [0.0, -4.0, 0.0, 8.0]
    [1.0, 0.0, -12.0, 0.0, 16.0]

**Orthogonality**

The Chebyshev polynomials of the second kind are orthogonal
on the interval `[-1, 1]` with respect to the weight
function `w(x) = \sqrt{1-x^2}`::

    >>> f = lambda x: chebyu(m,x)*chebyu(n,x)*sqrt(1-x**2)
    >>> m, n = 3, 4
    >>> quad(f, [-1, 1])
    0.0
    >>> m, n = 4, 4
    >>> quad(f, [-1, 1])
    1.5707963267949
"""

besselj = r"""
``besselj(n, x, derivative=0)`` gives the Bessel function of the first kind
`J_n(x)`. Bessel functions of the first kind are defined as
solutions of the differential equation

.. math ::

    x^2 y'' + x y' + (x^2 - n^2) y = 0

which appears, among other things, when solving the radial
part of Laplace's equation in cylindrical coordinates. This
equation has two solutions for given `n`, where the
`J_n`-function is the solution that is nonsingular at `x = 0`.
For positive integer `n`, `J_n(x)` behaves roughly like a sine
(odd `n`) or cosine (even `n`) multiplied by a magnitude factor
that decays slowly as `x \to \pm\infty`.

Generally, `J_n` is a special case of the hypergeometric
function `\,_0F_1`:

.. math ::

    J_n(x) = \frac{x^n}{2^n \Gamma(n+1)}
             \,_0F_1\left(n+1,-\frac{x^2}{4}\right)

With *derivative* = `m \ne 0`, the `m`-th derivative

.. math ::

    \frac{d^m}{dx^m} J_n(x)

is computed.

**Examples**

Evaluation is supported for arbitrary arguments, and at
arbitrary precision::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> besselj(2, 1000)
    -0.024777229528606
    >>> besselj(4, 0.75)
    0.000801070086542314
    >>> besselj(2, 1000j)
    (-2.48071721019185e+432 + 6.41567059811949e-437j)
    >>> mp.dps = 25
    >>> besselj(0.75j, 3+4j)
    (-2.778118364828153309919653 - 1.5863603889018621585533j)
    >>> mp.dps = 50
    >>> besselj(1, pi)
    0.28461534317975275734531059968613140570981118184947

Arguments may be large::

    >>> mp.dps = 25
    >>> besselj(0, 10000)
    -0.007096160353388801477265164
    >>> besselj(0, 10**10)
    0.000002175591750246891726859055
    >>> besselj(2, 10**100)
    7.337048736538615712436929e-51
    >>> besselj(2, 10**5*j)
    (-3.540725411970948860173735e+43426 + 4.4949812409615803110051e-43433j)

The Bessel functions of the first kind satisfy simple
symmetries around `x = 0`::

    >>> mp.dps = 15
    >>> nprint([besselj(n,0) for n in range(5)])
    [1.0, 0.0, 0.0, 0.0, 0.0]
    >>> nprint([besselj(n,pi) for n in range(5)])
    [-0.304242, 0.284615, 0.485434, 0.333458, 0.151425]
    >>> nprint([besselj(n,-pi) for n in range(5)])
    [-0.304242, -0.284615, 0.485434, -0.333458, 0.151425]

Roots of Bessel functions are often used::

    >>> nprint([findroot(j0, k) for k in [2, 5, 8, 11, 14]])
    [2.40483, 5.52008, 8.65373, 11.7915, 14.9309]
    >>> nprint([findroot(j1, k) for k in [3, 7, 10, 13, 16]])
    [3.83171, 7.01559, 10.1735, 13.3237, 16.4706]

The roots are not periodic, but the distance between successive
roots asymptotically approaches `2 \pi`. Bessel functions of
the first kind have the following normalization::

    >>> quadosc(j0, [0, inf], period=2*pi)
    1.0
    >>> quadosc(j1, [0, inf], period=2*pi)
    1.0

For `n = 1/2` or `n = -1/2`, the Bessel function reduces to a
trigonometric function::

    >>> x = 10
    >>> besselj(0.5, x), sqrt(2/(pi*x))*sin(x)
    (-0.13726373575505, -0.13726373575505)
    >>> besselj(-0.5, x), sqrt(2/(pi*x))*cos(x)
    (-0.211708866331398, -0.211708866331398)

Derivatives of any order can be computed (negative orders
correspond to integration)::

    >>> mp.dps = 25
    >>> besselj(0, 7.5, 1)
    -0.1352484275797055051822405
    >>> diff(lambda x: besselj(0,x), 7.5)
    -0.1352484275797055051822405
    >>> besselj(0, 7.5, 10)
    -0.1377811164763244890135677
    >>> diff(lambda x: besselj(0,x), 7.5, 10)
    -0.1377811164763244890135677
    >>> besselj(0,7.5,-1) - besselj(0,3.5,-1)
    -0.1241343240399987693521378
    >>> quad(j0, [3.5, 7.5])
    -0.1241343240399987693521378

Differentiation with a noninteger order gives the fractional derivative
in the sense of the Riemann-Liouville differintegral, as computed by
:func:`differint`::

    >>> mp.dps = 15
    >>> besselj(1, 3.5, 0.75)
    -0.385977722939384
    >>> differint(lambda x: besselj(1, x), 3.5, 0.75)
    -0.385977722939384

"""

besseli = r"""
``besseli(n, x, derivative=0)`` gives the modified Bessel function of the
first kind,

.. math ::

    I_n(x) = i^{-n} J_n(ix).

With *derivative* = `m \ne 0`, the `m`-th derivative

.. math ::

    \frac{d^m}{dx^m} I_n(x)

is computed.

**Examples**

Some values of `I_n(x)`::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> besseli(0,0)
    1.0
    >>> besseli(1,0)
    0.0
    >>> besseli(0,1)
    1.266065877752008335598245
    >>> besseli(3.5, 2+3j)
    (-0.2904369752642538144289025 - 0.4469098397654815837307006j)

Arguments may be large::

    >>> besseli(2, 1000)
    2.480717210191852440616782e+432
    >>> besseli(2, 10**10)
    4.299602851624027900335391e+4342944813
    >>> besseli(2, 6000+10000j)
    (-2.114650753239580827144204e+2603 + 4.385040221241629041351886e+2602j)

For integers `n`, the following integral representation holds::

    >>> mp.dps = 15
    >>> n = 3
    >>> x = 2.3
    >>> quad(lambda t: exp(x*cos(t))*cos(n*t), [0,pi])/pi
    0.349223221159309
    >>> besseli(n,x)
    0.349223221159309

Derivatives and antiderivatives of any order can be computed::

    >>> mp.dps = 25
    >>> besseli(2, 7.5, 1)
    195.8229038931399062565883
    >>> diff(lambda x: besseli(2,x), 7.5)
    195.8229038931399062565883
    >>> besseli(2, 7.5, 10)
    153.3296508971734525525176
    >>> diff(lambda x: besseli(2,x), 7.5, 10)
    153.3296508971734525525176
    >>> besseli(2,7.5,-1) - besseli(2,3.5,-1)
    202.5043900051930141956876
    >>> quad(lambda x: besseli(2,x), [3.5, 7.5])
    202.5043900051930141956876

"""

bessely = r"""
``bessely(n, x, derivative=0)`` gives the Bessel function of the second kind,

.. math ::

    Y_n(x) = \frac{J_n(x) \cos(\pi n) - J_{-n}(x)}{\sin(\pi n)}.

For `n` an integer, this formula should be understood as a
limit. With *derivative* = `m \ne 0`, the `m`-th derivative

.. math ::

    \frac{d^m}{dx^m} Y_n(x)

is computed.

**Examples**

Some values of `Y_n(x)`::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> bessely(0,0), bessely(1,0), bessely(2,0)
    (-inf, -inf, -inf)
    >>> bessely(1, pi)
    0.3588729167767189594679827
    >>> bessely(0.5, 3+4j)
    (9.242861436961450520325216 - 3.085042824915332562522402j)

Arguments may be large::

    >>> bessely(0, 10000)
    0.00364780555898660588668872
    >>> bessely(2.5, 10**50)
    -4.8952500412050989295774e-26
    >>> bessely(2.5, -10**50)
    (0.0 + 4.8952500412050989295774e-26j)

Derivatives and antiderivatives of any order can be computed::

    >>> bessely(2, 3.5, 1)
    0.3842618820422660066089231
    >>> diff(lambda x: bessely(2, x), 3.5)
    0.3842618820422660066089231
    >>> bessely(0.5, 3.5, 1)
    -0.2066598304156764337900417
    >>> diff(lambda x: bessely(0.5, x), 3.5)
    -0.2066598304156764337900417
    >>> diff(lambda x: bessely(2, x), 0.5, 10)
    -208173867409.5547350101511
    >>> bessely(2, 0.5, 10)
    -208173867409.5547350101511
    >>> bessely(2, 100.5, 100)
    0.02668487547301372334849043
    >>> quad(lambda x: bessely(2,x), [1,3])
    -1.377046859093181969213262
    >>> bessely(2,3,-1) - bessely(2,1,-1)
    -1.377046859093181969213262

"""

besselk = r"""
``besselk(n, x)`` gives the modified Bessel function of the
second kind,

.. math ::

    K_n(x) = \frac{\pi}{2} \frac{I_{-n}(x)-I_{n}(x)}{\sin(\pi n)}

For `n` an integer, this formula should be understood as a
limit.

**Examples**

Evaluation is supported for arbitrary complex arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> besselk(0,1)
    0.4210244382407083333356274
    >>> besselk(0, -1)
    (0.4210244382407083333356274 - 3.97746326050642263725661j)
    >>> besselk(3.5, 2+3j)
    (-0.02090732889633760668464128 + 0.2464022641351420167819697j)
    >>> besselk(2+3j, 0.5)
    (0.9615816021726349402626083 + 0.1918250181801757416908224j)

Arguments may be large::

    >>> besselk(0, 100)
    4.656628229175902018939005e-45
    >>> besselk(1, 10**6)
    4.131967049321725588398296e-434298
    >>> besselk(1, 10**6*j)
    (0.001140348428252385844876706 - 0.0005200017201681152909000961j)
    >>> besselk(4.5, fmul(10**50, j, exact=True))
    (1.561034538142413947789221e-26 + 1.243554598118700063281496e-25j)

The point `x = 0` is a singularity (logarithmic if `n = 0`)::

    >>> besselk(0,0)
    +inf
    >>> besselk(1,0)
    +inf
    >>> for n in range(-4, 5):
    ...     print besselk(n, '1e-1000')
    ...
    4.8e+4001
    8.0e+3000
    2.0e+2000
    1.0e+1000
    2302.701024509704096466802
    1.0e+1000
    2.0e+2000
    8.0e+3000
    4.8e+4001

"""

hankel1 = r"""
``hankel1(n,x)`` computes the Hankel function of the first kind,
which is the complex combination of Bessel functions given by

.. math ::

    H_n^{(1)}(x) = J_n(x) + i Y_n(x).

**Examples**

The Hankel function is generally complex-valued::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> hankel1(2, pi)
    (0.4854339326315091097054957 - 0.0999007139290278787734903j)
    >>> hankel1(3.5, pi)
    (0.2340002029630507922628888 - 0.6419643823412927142424049j)
"""

hankel2 = r"""
``hankel2(n,x)`` computes the Hankel function of the second kind,
which is the complex combination of Bessel functions given by

.. math ::

    H_n^{(2)}(x) = J_n(x) - i Y_n(x).

**Examples**

The Hankel function is generally complex-valued::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> hankel2(2, pi)
    (0.4854339326315091097054957 + 0.0999007139290278787734903j)
    >>> hankel2(3.5, pi)
    (0.2340002029630507922628888 + 0.6419643823412927142424049j)
"""

lambertw = r"""
The Lambert W function `W(z)` is defined as the inverse function
of `w \exp(w)`. In other words, the value of `W(z)` is such that
`z = W(z) \exp(W(z))` for any complex number `z`.

The Lambert W function is a multivalued function with infinitely
many branches. Each branch gives a separate solution of the
equation `w \exp(w)`. All branches are supported by
:func:`lambertw`:

* ``lambertw(z)`` gives the principal solution (branch 0)

* ``lambertw(z, k)`` gives the solution on branch `k`

The Lambert W function has two partially real branches: the
principal branch (`k = 0`) is real for real `z > -1/e`, and the
`k = -1` branch is real for `-1/e < z < 0`. All branches except
`k = 0` have a logarithmic singularity at `z = 0`.

**Basic examples**

The Lambert W function is the inverse of `w \exp(w)`::

    >>> from mpmath import *
    >>> mp.dps = 35; mp.pretty = True
    >>> w = lambertw(1)
    >>> w
    0.56714329040978387299996866221035555
    >>> w*exp(w)
    1.0

Any branch gives a valid inverse::

    >>> w = lambertw(1, k=3)
    >>> w    # doctest: +NORMALIZE_WHITESPACE
    (-2.8535817554090378072068187234910812 +
      17.113535539412145912607826671159289j)
    >>> w*exp(w)
    (1.0 + 3.5075477124212226194278700785075126e-36j)

**Applications to equation-solving**

The Lambert W function may be used to solve various kinds of
equations, such as finding the value of the infinite power
tower `z^{z^{z^{\ldots}}}`::

    >>> def tower(z, n):
    ...     if n == 0:
    ...         return z
    ...     return z ** tower(z, n-1)
    ...
    >>> tower(0.5, 100)
    0.641185744504986
    >>> mp.dps = 50
    >>> -lambertw(-log(0.5))/log(0.5)
    0.6411857445049859844862004821148236665628209571911

**Properties**

The Lambert W function grows roughly like the natural logarithm
for large arguments::

    >>> mp.dps = 15
    >>> lambertw(1000)
    5.2496028524016
    >>> log(1000)
    6.90775527898214
    >>> lambertw(10**100)
    224.843106445119
    >>> log(10**100)
    230.258509299405

The principal branch of the Lambert W function has a rational
Taylor series expansion around `z = 0`::

    >>> nprint(taylor(lambertw, 0, 6), 10)
    [0.0, 1.0, -1.0, 1.5, -2.666666667, 5.208333333, -10.8]

Some special values and limits are::

    >>> mp.dps = 15
    >>> lambertw(0)
    0.0
    >>> lambertw(1)
    0.567143290409784
    >>> lambertw(e)
    1.0
    >>> lambertw(inf)
    +inf
    >>> lambertw(0, k=-1)
    -inf
    >>> lambertw(0, k=3)
    -inf
    >>> lambertw(inf, k=3)
    (+inf + 18.8495559215388j)

The `k = 0` and `k = -1` branches join at `z = -1/e` where
`W(z) = -1` for both branches. Since `-1/e` can only be represented
approximately with mpmath numbers, evaluating the Lambert W function
at this point only gives `-1` approximately::

    >>> mp.dps = 25
    >>> lambertw(-1/e, 0)
    -0.999999999999837133022867
    >>> lambertw(-1/e, -1)
    -1.00000000000016286697718

If `-1/e` happens to round in the negative direction, there might be
a small imaginary part::

    >>> mp.dps = 15
    >>> lambertw(-1/e)
    (-1.0 + 8.22007971511612e-9j)

**Possible issues**

The evaluation can become inaccurate very close to the branch point
at `-1/e`. In some corner cases, :func:`lambertw` might currently
fail to converge, or can end up on the wrong branch.

**Algorithm**

Halley's iteration is used to invert `w \exp(w)`, using a first-order
asymptotic approximation (`O(\log(w))` or `O(w)`) as the initial
estimate.

The definition, implementation and choice of branches is based
on Corless et al, "On the Lambert W function", Adv. Comp. Math. 5
(1996) 329-359, available online here:
http://www.apmaths.uwo.ca/~djeffrey/Offprints/W-adv-cm.pdf

TODO: use a series expansion when extremely close to the branch point
at `-1/e` and make sure that the proper branch is chosen there
"""

barnesg = r"""
Evaluates the Barnes G-function, which generalizes the
superfactorial (:func:`superfac`) and by extension also the
hyperfactorial (:func:`hyperfac`) to the complex numbers
in an analogous way to how the gamma function generalizes
the ordinary factorial.

The Barnes G-function may be defined in terms of a Weierstrass
product:

.. math ::

    G(z+1) = (2\pi)^{z/2} e^{-[z(z+1)+\gamma z^2]/2}
    \prod_{n=1}^\infty
    \left[\left(1+\frac{z}{n}\right)^ne^{-z+z^2/(2n)}\right]

For positive integers `n`, we have have relation to superfactorials
`G(n) = \mathrm{sf}(n-2) = 0! \cdot 1! \cdots (n-2)!`.

**Examples**

Some elementary values and limits of the Barnes G-function::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> barnesg(1), barnesg(2), barnesg(3)
    (1.0, 1.0, 1.0)
    >>> barnesg(4)
    2.0
    >>> barnesg(5)
    12.0
    >>> barnesg(6)
    288.0
    >>> barnesg(7)
    34560.0
    >>> barnesg(8)
    24883200.0
    >>> barnesg(inf)
    +inf
    >>> barnesg(0), barnesg(-1), barnesg(-2)
    (0.0, 0.0, 0.0)

Closed-form values are known for some rational arguments::

    >>> barnesg('1/2')
    0.603244281209446
    >>> sqrt(exp(0.25+log(2)/12)/sqrt(pi)/glaisher**3)
    0.603244281209446
    >>> barnesg('1/4')
    0.29375596533861
    >>> nthroot(exp('3/8')/exp(catalan/pi)/
    ...      gamma(0.25)**3/sqrt(glaisher)**9, 4)
    0.29375596533861

The Barnes G-function satisfies the functional equation
`G(z+1) = \Gamma(z) G(z)`::

    >>> z = pi
    >>> barnesg(z+1)
    2.39292119327948
    >>> gamma(z)*barnesg(z)
    2.39292119327948

The asymptotic growth rate of the Barnes G-function is related to
the Glaisher-Kinkelin constant::

    >>> limit(lambda n: barnesg(n+1)/(n**(n**2/2-mpf(1)/12)*
    ...     (2*pi)**(n/2)*exp(-3*n**2/4)), inf)
    0.847536694177301
    >>> exp('1/12')/glaisher
    0.847536694177301

The Barnes G-function can be differentiated in closed form::

    >>> z = 3
    >>> diff(barnesg, z)
    0.264507203401607
    >>> barnesg(z)*((z-1)*psi(0,z)-z+(log(2*pi)+1)/2)
    0.264507203401607

Evaluation is supported for arbitrary arguments and at arbitrary
precision::

    >>> barnesg(6.5)
    2548.7457695685
    >>> barnesg(-pi)
    0.00535976768353037
    >>> barnesg(3+4j)
    (-0.000676375932234244 - 4.42236140124728e-5j)
    >>> mp.dps = 50
    >>> barnesg(1/sqrt(2))
    0.81305501090451340843586085064413533788206204124732
    >>> q = barnesg(10j)
    >>> q.real
    0.000000000021852360840356557241543036724799812371995850552234
    >>> q.imag
    -0.00000000000070035335320062304849020654215545839053210041457588
    >>> mp.dps = 15
    >>> barnesg(100)
    3.10361006263698e+6626
    >>> barnesg(-101)
    0.0
    >>> barnesg(-10.5)
    5.94463017605008e+25
    >>> barnesg(-10000.5)
    -6.14322868174828e+167480422
    >>> barnesg(1000j)
    (5.21133054865546e-1173597 + 4.27461836811016e-1173597j)
    >>> barnesg(-1000+1000j)
    (2.43114569750291e+1026623 + 2.24851410674842e+1026623j)


**References**

1. Whittaker & Watson, *A Course of Modern Analysis*,
   Cambridge University Press, 4th edition (1927), p.264
2. http://en.wikipedia.org/wiki/Barnes_G-function
3. http://mathworld.wolfram.com/BarnesG-Function.html

"""

superfac = r"""
Computes the superfactorial, defined as the product of
consecutive factorials

.. math ::

    \mathrm{sf}(n) = \prod_{k=1}^n k!

For general complex `z`, `\mathrm{sf}(z)` is defined
in terms of the Barnes G-function (see :func:`barnesg`).

**Examples**

The first few superfactorials are (OEIS A000178)::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(10):
    ...     print n, superfac(n)
    ...
    0 1.0
    1 1.0
    2 2.0
    3 12.0
    4 288.0
    5 34560.0
    6 24883200.0
    7 125411328000.0
    8 5.05658474496e+15
    9 1.83493347225108e+21

Superfactorials grow very rapidly::

    >>> superfac(1000)
    3.24570818422368e+1177245
    >>> superfac(10**10)
    2.61398543581249e+467427913956904067453

Evaluation is supported for arbitrary arguments::

    >>> mp.dps = 25
    >>> superfac(pi)
    17.20051550121297985285333
    >>> superfac(2+3j)
    (-0.005915485633199789627466468 + 0.008156449464604044948738263j)
    >>> diff(superfac, 1)
    0.2645072034016070205673056

**References**

1. http://www.research.att.com/~njas/sequences/A000178

"""


hyperfac = r"""
Computes the hyperfactorial, defined for integers as the product

.. math ::

    H(n) = \prod_{k=1}^n k^k.


The hyperfactorial satisfies the recurrence formula `H(z) = z^z H(z-1)`.
It can be defined more generally in terms of the Barnes G-function (see
:func:`barnesg`) and the gamma function by the formula

.. math ::

    H(z) = \frac{\Gamma(z+1)^z}{G(z)}.

The extension to complex numbers can also be done via
the integral representation

.. math ::

    H(z) = (2\pi)^{-z/2} \exp \left[
        {z+1 \choose 2} + \int_0^z \log(t!)\,dt
        \right].

**Examples**

The rapidly-growing sequence of hyperfactorials begins
(OEIS A002109)::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(10):
    ...     print n, hyperfac(n)
    ...
    0 1.0
    1 1.0
    2 4.0
    3 108.0
    4 27648.0
    5 86400000.0
    6 4031078400000.0
    7 3.3197663987712e+18
    8 5.56964379417266e+25
    9 2.15779412229419e+34

Some even larger hyperfactorials are::

    >>> hyperfac(1000)
    5.46458120882585e+1392926
    >>> hyperfac(10**10)
    4.60408207642219e+489142638002418704309

The hyperfactorial can be evaluated for arbitrary arguments::

    >>> hyperfac(0.5)
    0.880449235173423
    >>> diff(hyperfac, 1)
    0.581061466795327
    >>> hyperfac(pi)
    205.211134637462
    >>> hyperfac(-10+1j)
    (3.01144471378225e+46 - 2.45285242480185e+46j)

The recurrence property of the hyperfactorial holds
generally::

    >>> z = 3-4*j
    >>> hyperfac(z)
    (-4.49795891462086e-7 - 6.33262283196162e-7j)
    >>> z**z * hyperfac(z-1)
    (-4.49795891462086e-7 - 6.33262283196162e-7j)
    >>> z = mpf(-0.6)
    >>> chop(z**z * hyperfac(z-1))
    1.28170142849352
    >>> hyperfac(z)
    1.28170142849352

The hyperfactorial may also be computed using the integral
definition::

    >>> z = 2.5
    >>> hyperfac(z)
    15.9842119922237
    >>> (2*pi)**(-z/2)*exp(binomial(z+1,2) +
    ...     quad(lambda t: loggamma(t+1), [0, z]))
    15.9842119922237

:func:`hyperfac` supports arbitrary-precision evaluation::

    >>> mp.dps = 50
    >>> hyperfac(10)
    215779412229418562091680268288000000000000000.0
    >>> hyperfac(1/sqrt(2))
    0.89404818005227001975423476035729076375705084390942

**References**

1. http://www.research.att.com/~njas/sequences/A002109
2. http://mathworld.wolfram.com/Hyperfactorial.html

"""


loggamma = r"""
Computes the log-gamma function. Unlike `\ln(\Gamma(z))`, which
has infinitely many complex branch cuts, the log-gamma function
only has a single branch cut along the negative half-axis.
The functions are identical only on (and very close to) the positive
half-axis; elsewhere they differ by `2 n \pi i` (the real parts
agree)::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> loggamma(13.2), log(gamma(13.2))
    (20.494004194566, 20.494004194566)
    >>> loggamma(3+4j)
    (-1.75662678460378 + 4.74266443803466j)
    >>> log(gamma(3+4j))
    (-1.75662678460378 - 1.54052086914493j)

Note: this is a placeholder implementation. It is slower than
:func:`gamma`, and is in particular *not* faster than :func:`gamma`
for large arguments.
"""

siegeltheta = r"""
Computes the Riemann-Siegel theta function,

.. math ::

    \theta(t) = \frac{
    \log\Gamma\left(\frac{1+2it}{4}\right) -
    \log\Gamma\left(\frac{1-2it}{4}\right)
    }{2i} - \frac{\log \pi}{2} t.

The Riemann-Siegel theta function is important in
providing the phase factor for the Z-function
(see :func:`siegelz`). Evaluation is supported for real and
complex arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> siegeltheta(0)
    0.0
    >>> siegeltheta(inf)
    +inf
    >>> siegeltheta(-inf)
    -inf
    >>> siegeltheta(1)
    -1.767547952812290388302216
    >>> siegeltheta(10+0.25j)
    (-3.068638039426838572528867 + 0.05804937947429712998395177j)

The Riemann-Siegel theta function has odd symmetry around `t = 0`,
two local extreme points and three real roots including 0 (located
symmetrically)::

    >>> nprint(chop(taylor(siegeltheta, 0, 5)))
    [0.0, -2.68609, 0.0, 2.69433, 0.0, -6.40218]
    >>> findroot(diffun(siegeltheta), 7)
    6.28983598883690277966509
    >>> findroot(siegeltheta, 20)
    17.84559954041086081682634

For large `t`, there is a famous asymptotic formula
for `\theta(t)`, to first order given by::

    >>> t = mpf(10**6)
    >>> siegeltheta(t)
    5488816.353078403444882823
    >>> -t*log(2*pi/t)/2-t/2
    5488816.745777464310273645
"""

grampoint = r"""
Gives the `n`-th Gram point `g_n`, defined as the solution
to the equation `\theta(g_n) = \pi n` where `\theta(t)`
is the Riemann-Siegel theta function (:func:`siegeltheta`).

The first few Gram points are::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> grampoint(0)
    17.84559954041086081682634
    >>> grampoint(1)
    23.17028270124630927899664
    >>> grampoint(2)
    27.67018221781633796093849
    >>> grampoint(3)
    31.71797995476405317955149

Checking the definition::

    >>> siegeltheta(grampoint(3))
    9.42477796076937971538793
    >>> 3*pi
    9.42477796076937971538793

A large Gram point::

    >>> grampoint(10**10)
    3293531632.728335454561153

Gram points are useful when studying the Z-function
(:func:`siegelz`). See the documentation of that function
for additional examples.

:func:`grampoint` can solve the defining equation for
nonintegral `n`. There is a fixed point where `g(x) = x`::

    >>> findroot(lambda x: grampoint(x) - x, 10000)
    9146.698193171459265866198

**References**

1. http://mathworld.wolfram.com/GramPoint.html

"""

siegelz = r"""
Computes the Z-function, also known as the Riemann-Siegel Z function,

.. math ::

    Z(t) = e^{i \theta(t)} \zeta(1/2+it)

where `\zeta(s)` is the Riemann zeta function (:func:`zeta`)
and where `\theta(t)` denotes the Riemann-Siegel theta function
(see :func:`siegeltheta`).

Evaluation is supported for real and complex arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> siegelz(1)
    -0.7363054628673177346778998
    >>> siegelz(3+4j)
    (-0.1852895764366314976003936 - 0.2773099198055652246992479j)

The Z-function has a Maclaurin expansion::

    >>> nprint(chop(taylor(siegelz, 0, 4)))
    [-1.46035, 0.0, 2.73588, 0.0, -8.39357]

The Z-function `Z(t)` is equal to `\pm |\zeta(s)|` on the
critical line `s = 1/2+it` (i.e. for real arguments `t`
to `Z`).  Its zeros coincide with those of the Riemann zeta
function::

    >>> findroot(siegelz, 14)
    14.13472514173469379045725
    >>> findroot(siegelz, 20)
    21.02203963877155499262848
    >>> findroot(zeta, 0.5+14j)
    (0.5 + 14.13472514173469379045725j)
    >>> findroot(zeta, 0.5+20j)
    (0.5 + 21.02203963877155499262848j)

Since the Z-function is real-valued on the critical line
(and unlike `|\zeta(s)|` analytic), it is useful for
investigating the zeros of the Riemann zeta function.
For example, one can use a root-finding algorithm based
on sign changes::

    >>> findroot(siegelz, [100, 200], solver='bisect')
    176.4414342977104188888926

To locate roots, Gram points `g_n` which can be computed
by :func:`grampoint` are useful. If `(-1)^n Z(g_n)` is
positive for two consecutive `n`, then `Z(t)` must have
a zero between those points::

    >>> g10 = grampoint(10)
    >>> g11 = grampoint(11)
    >>> (-1)**10 * siegelz(g10) > 0
    True
    >>> (-1)**11 * siegelz(g11) > 0
    True
    >>> findroot(siegelz, [g10, g11], solver='bisect')
    56.44624769706339480436776
    >>> g10, g11
    (54.67523744685325626632663, 57.54516517954725443703014)

"""

zetazero = r"""
Returns the `n`-th nontrivial zero of the Riemann zeta function.
The zero is computed using :func:`findroot`, using a table lookup
for the initial point.

The zeros are located on the critical line with real part 1/2::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> zetazero(1)
    (0.5 + 14.13472514173469379045725j)
    >>> zetazero(2)
    (0.5 + 21.02203963877155499262848j)
    >>> zetazero(20)
    (0.5 + 77.14484006887480537268266j)

Negative indices give the conjugate zeros (`n = 0` is undefined)::

    >>> zetazero(-1)
    (0.5 - 14.13472514173469379045725j)

The default table only provides `n` up to 100. For larger `n` up to
100,000,  :func:`zetazero` will automatically download a table
(1.8 MB) from the website of Andrew Odlyzko [1]. This requires a
fast connection to the internet. Alternatively, you can supply the
url to a custom table. The table should be a file listing the
imaginary parts as float literals, separated by line breaks.

1. http://www.dtc.umn.edu/~odlyzko/zeta_tables/
"""

riemannr = r"""
Evaluates the Riemann R function, a smooth approximation of the
prime counting function `\pi(x)` (see :func:`primepi`). The Riemann
R function gives a fast numerical approximation useful e.g. to
roughly estimate the number of primes in a given interval.

The Riemann R function is computed using the rapidly convergent Gram
series,

.. math ::

    R(x) = 1 + \sum_{k=1}^{\infty}
        \frac{\log^k x}{k k! \zeta(k+1)}.

From the Gram series, one sees that the Riemann R function is a
well-defined analytic function (except for a branch cut along
the negative real half-axis); it can be evaluated for arbitrary
real or complex arguments.

The Riemann R function gives a very accurate approximation
of the prime counting function. For example, it is wrong by at
most 2 for `x < 1000`, and for `x = 10^9` differs from the exact
value of `\pi(x)` by 79, or less than two parts in a million.
It is about 10 times more accurate than the logarithmic integral
estimate (see :func:`li`), which however is even faster to evaluate.
It is orders of magnitude more accurate than the extremely
fast `x/\log x` estimate.

**Examples**

For small arguments, the Riemann R function almost exactly
gives the prime counting function if rounded to the nearest
integer::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> primepi(50), riemannr(50)
    (15, 14.9757023241462)
    >>> max(abs(primepi(n)-round(riemannr(n))) for n in range(100))
    1.0
    >>> max(abs(primepi(n)-round(riemannr(n))) for n in range(300))
    2.0

The Riemann R function can be evaluated for arguments far too large
for exact determination of `\pi(x)` to be computationally
feasible with any presently known algorithm::

    >>> riemannr(10**30)
    1.46923988977204e+28
    >>> riemannr(10**100)
    4.3619719871407e+97
    >>> riemannr(10**1000)
    4.3448325764012e+996

A comparison of the Riemann R function and logarithmic integral estimates
for `\pi(x)` using exact values of `\pi(10^n)` up to `n = 9`.
The fractional error is shown in parentheses::

    >>> exact = [4,25,168,1229,9592,78498,664579,5761455,50847534]
    >>> for n, p in enumerate(exact):
    ...     n += 1
    ...     r, l = riemannr(10**n), li(10**n)
    ...     rerr, lerr = nstr((r-p)/p,3), nstr((l-p)/p,3)
    ...     print "%i %i %s(%s) %s(%s)" % (n, p, r, rerr, l, lerr)
    ...
    1 4 4.56458314100509(1.41e-1) 6.1655995047873(5.41e-1)
    2 25 25.6616332669242(2.65e-2) 30.1261415840796(2.05e-1)
    3 168 168.359446281167(2.14e-3) 177.609657990152(5.72e-2)
    4 1229 1226.93121834343(-1.68e-3) 1246.13721589939(1.39e-2)
    5 9592 9587.43173884197(-4.76e-4) 9629.8090010508(3.94e-3)
    6 78498 78527.3994291277(3.75e-4) 78627.5491594622(1.65e-3)
    7 664579 664667.447564748(1.33e-4) 664918.405048569(5.11e-4)
    8 5761455 5761551.86732017(1.68e-5) 5762209.37544803(1.31e-4)
    9 50847534 50847455.4277214(-1.55e-6) 50849234.9570018(3.35e-5)

The derivative of the Riemann R function gives the approximate
probability for a number of magnitude `x` to be prime::

    >>> diff(riemannr, 1000)
    0.141903028110784
    >>> mpf(primepi(1050) - primepi(950)) / 100
    0.15

Evaluation is supported for arbitrary arguments and at arbitrary
precision::

    >>> mp.dps = 30
    >>> riemannr(7.5)
    3.72934743264966261918857135136
    >>> riemannr(-4+2j)
    (-0.551002208155486427591793957644 + 2.16966398138119450043195899746j)

"""

primepi = r"""
Evaluates the prime counting function, `\pi(x)`, which gives
the number of primes less than or equal to `x`. The argument
`x` may be fractional.

The prime counting function is very expensive to evaluate
precisely for large `x`, and the present implementation is
not optimized in any way. For numerical approximation of the
prime counting function, it is better to use :func:`primepi2`
or :func:`riemannr`.

Some values of the prime counting function::

    >>> from mpmath import *
    >>> [primepi(k) for k in range(20)]
    [0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 8]
    >>> primepi(3.5)
    2
    >>> primepi(100000)
    9592

"""

primepi2 = r"""
Returns an interval (as an ``mpi`` instance) providing bounds
for the value of the prime counting function `\pi(x)`. For small
`x`, :func:`primepi2` returns an exact interval based on
the output of :func:`primepi`. For `x > 2656`, a loose interval
based on Schoenfeld's inequality

.. math ::

    |\pi(x) - \mathrm{li}(x)| < \frac{\sqrt x \log x}{8 \pi}

is returned. This estimate is rigorous assuming the truth of
the Riemann hypothesis, and can be computed very quickly.

**Examples**

Exact values of the prime counting function for small `x`::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> primepi2(10)
    [4.0, 4.0]
    >>> primepi2(100)
    [25.0, 25.0]
    >>> primepi2(1000)
    [168.0, 168.0]

Loose intervals are generated for moderately large `x`:

    >>> primepi2(10000), primepi(10000)
    ([1209.0, 1283.0], 1229)
    >>> primepi2(50000), primepi(50000)
    ([5070.0, 5263.0], 5133)

As `x` increases, the absolute error gets worse while the relative
error improves. The exact value of `\pi(10^{23})` is
1925320391606803968923, and :func:`primepi2` gives 9 significant
digits::

    >>> p = primepi2(10**23)
    >>> p
    [1.9253203909477020467e+21, 1.925320392280406229e+21]
    >>> p.delta / p.a
    6.9219865355293e-10

A more precise, nonrigorous estimate for `\pi(x)` can be
obtained using the Riemann R function (:func:`riemannr`).
For large enough `x`, the value returned by :func:`primepi2`
essentially amounts to a small perturbation of the value returned by
:func:`riemannr`::

    >>> primepi2(10**100)
    [4.3619719871407024816e+97, 4.3619719871407032404e+97]
    >>> riemannr(10**100)
    4.3619719871407e+97
"""

primezeta = r"""
Computes the prime zeta function, which is defined
in analogy with the Riemann zeta function (:func:`zeta`)
as

.. math ::

    P(s) = \sum_p \frac{1}{p^s}

where the sum is taken over all prime numbers `p`. Although
this sum only converges for `\mathrm{Re}(s) > 1`, the
function is defined by analytic continuation in the
half-plane `\mathrm{Re}(s) > 0`.

**Examples**

Arbitrary-precision evaluation for real and complex arguments is
supported::

    >>> from mpmath import *
    >>> mp.dps = 30; mp.pretty = True
    >>> primezeta(2)
    0.452247420041065498506543364832
    >>> primezeta(pi)
    0.15483752698840284272036497397
    >>> mp.dps = 50
    >>> primezeta(3)
    0.17476263929944353642311331466570670097541212192615
    >>> mp.dps = 20
    >>> primezeta(3+4j)
    (-0.12085382601645763295 - 0.013370403397787023602j)

The prime zeta function has a logarithmic pole at `s = 1`,
with residue equal to the difference of the Mertens and
Euler constants::

    >>> primezeta(1)
    +inf
    >>> extradps(25)(lambda x: primezeta(1+x)+log(x))(+eps)
    -0.31571845205389007685
    >>> mertens-euler
    -0.31571845205389007685

The analytic continuation to `0 < \mathrm{Re}(s) \le 1`
is implemented. In this strip the function exhibits
very complex behavior; on the unit interval, it has poles at
`1/n` for every squarefree integer `n`::

    >>> primezeta(0.5)         # Pole at s = 1/2
    (-inf + 3.1415926535897932385j)
    >>> primezeta(0.25)
    (-1.0416106801757269036 + 0.52359877559829887308j)
    >>> primezeta(0.5+10j)
    (0.54892423556409790529 + 0.45626803423487934264j)

Although evaluation works in principle for any `\mathrm{Re}(s) > 0`,
it should be noted that the evaluation time increases exponentially
as `s` approaches the imaginary axis.

For large `\mathrm{Re}(s)`, `P(s)` is asymptotic to `2^{-s}`::

    >>> primezeta(inf)
    0.0
    >>> primezeta(10), mpf(2)**-10
    (0.00099360357443698021786, 0.0009765625)
    >>> primezeta(1000)
    9.3326361850321887899e-302
    >>> primezeta(1000+1000j)
    (-3.8565440833654995949e-302 - 8.4985390447553234305e-302j)

**References**

Carl-Erik Froberg, "On the prime zeta function",
BIT 8 (1968), pp. 187-202.

"""

bernpoly = r"""
Evaluates the Bernoulli polynomial `B_n(z)`.

The first few Bernoulli polynomials are::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(6):
    ...     nprint(chop(taylor(lambda x: bernpoly(n,x), 0, n)))
    ...
    [1.0]
    [-0.5, 1.0]
    [0.166667, -1.0, 1.0]
    [0.0, 0.5, -1.5, 1.0]
    [-3.33333e-2, 0.0, 1.0, -2.0, 1.0]
    [0.0, -0.166667, 0.0, 1.66667, -2.5, 1.0]

At `z = 0`, the Bernoulli polynomial evaluates to a
Bernoulli number (see :func:`bernoulli`)::

    >>> bernpoly(12, 0), bernoulli(12)
    (-0.253113553113553, -0.253113553113553)
    >>> bernpoly(13, 0), bernoulli(13)
    (0.0, 0.0)

"""

polylog = r"""
Computes the polylogarithm, defined by the sum

.. math ::

    \mathrm{Li}_s(z) = \sum_{k=1}^{\infty} \frac{z^k}{k^s}.

This series is convergent only for `|z| < 1`, so elsewhere
the analytic continuation is implied.

The polylogarithm should not be confused with the logarithmic
integral (also denoted by Li or li), which is implemented
as :func:`li`.

**Examples**

The polylogarithm satisfies a huge number of functional identities.
A sample of polylogarithm evaluations is shown below::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> polylog(1,0.5), log(2)
    (0.693147180559945, 0.693147180559945)
    >>> polylog(2,0.5), (pi**2-6*log(2)**2)/12
    (0.582240526465012, 0.582240526465012)
    >>> polylog(2,-phi), -log(phi)**2-pi**2/10
    (-1.21852526068613, -1.21852526068613)
    >>> polylog(3,0.5), 7*zeta(3)/8-pi**2*log(2)/12+log(2)**3/6
    (0.53721319360804, 0.53721319360804)

:func:`polylog` can evaluate the analytic continuation of the
polylogarithm when `s` is an integer::

    >>> polylog(2, 10)
    (0.536301287357863 - 7.23378441241546j)
    >>> polylog(2, -10)
    -4.1982778868581
    >>> polylog(2, 10j)
    (-3.05968879432873 + 3.71678149306807j)
    >>> polylog(-2, 10)
    -0.150891632373114
    >>> polylog(-2, -10)
    0.067618332081142
    >>> polylog(-2, 10j)
    (0.0384353698579347 + 0.0912451798066779j)

Some more examples, with arguments on the unit circle (note that
the series definition cannot be used for computation here)::

    >>> polylog(2,j)
    (-0.205616758356028 + 0.915965594177219j)
    >>> j*catalan-pi**2/48
    (-0.205616758356028 + 0.915965594177219j)
    >>> polylog(3,exp(2*pi*j/3))
    (-0.534247512515375 + 0.765587078525922j)
    >>> -4*zeta(3)/9 + 2*j*pi**3/81
    (-0.534247512515375 + 0.765587078525921j)

Polylogarithms of different order are related by integration
and differentiation::

    >>> s, z = 3, 0.5
    >>> polylog(s+1, z)
    0.517479061673899
    >>> quad(lambda t: polylog(s,t)/t, [0, z])
    0.517479061673899
    >>> z*diff(lambda t: polylog(s+2,t), z)
    0.517479061673899

Taylor series expansions around `z = 0` are::

    >>> for n in range(-3, 4):
    ...     nprint(taylor(lambda x: polylog(n,x), 0, 5))
    ...
    [0.0, 1.0, 8.0, 27.0, 64.0, 125.0]
    [0.0, 1.0, 4.0, 9.0, 16.0, 25.0]
    [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    [0.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    [0.0, 1.0, 0.5, 0.333333, 0.25, 0.2]
    [0.0, 1.0, 0.25, 0.111111, 6.25e-2, 4.0e-2]
    [0.0, 1.0, 0.125, 3.7037e-2, 1.5625e-2, 8.0e-3]

The series defining the polylogarithm is simultaneously
a Taylor series and an L-series. For certain values of `z`, the
polylogarithm reduces to a pure zeta function::

    >>> polylog(pi, 1), zeta(pi)
    (1.17624173838258, 1.17624173838258)
    >>> polylog(pi, -1), -altzeta(pi)
    (-0.909670702980385, -0.909670702980385)

Evaluation for arbitrary, nonintegral `s` is supported
for `z` within the unit circle:

    >>> polylog(3+4j, 0.25)
    (0.24258605789446 - 0.00222938275488344j)
    >>> nsum(lambda k: 0.25**k / k**(3+4j), [1,inf])
    (0.24258605789446 - 0.00222938275488344j)

It is also currently supported outside of the unit circle for `z`
not too large in magnitude::

    >>> polylog(1+j, 20+40j)
    (-7.1421172179728 - 3.92726697721369j)
    >>> polylog(1+j, 200+400j)
    Traceback (most recent call last):
      ...
    NotImplementedError: polylog for arbitrary s and z

**References**

1. Richard Crandall, "Note on fast polylogarithm computation"
   http://people.reed.edu/~crandall/papers/Polylog.pdf
2. http://en.wikipedia.org/wiki/Polylogarithm
3. http://mathworld.wolfram.com/Polylogarithm.html

"""

bell = r"""
For `n` a nonnegative integer, ``bell(n,x)`` evaluates the Bell
polynomial `B_n(x)`, the first few of which are

.. math ::

    B_0(x) = 1

    B_1(x) = x

    B_2(x) = x^2+x

    B_3(x) = x^3+3x^2+x

If `x = 1` or :func:`bell` is called with only one argument, it
gives the `n`-th Bell number `B_n`, which is the number of
partitions of a set with `n` elements. By setting the precision to
at least `\log_{10} B_n` digits, :func:`bell` provides fast
calculation of exact Bell numbers.

In general, :func:`bell` computes

.. math ::

    B_n(x) = e^{-x} \left(\mathrm{sinc}(\pi n) + E_n(x)\right)

where `E_n(x)` is the generalized exponential function implemented
by :func:`polyexp`. This is an extension of Dobinski's formula [1],
where the modification is the sinc term ensuring that `B_n(x)` is
continuous in `n`; :func:`bell` can thus be evaluated,
differentiated, etc for arbitrary complex arguments.

**Examples**

Simple evaluations::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> bell(0, 2.5)
    1.0
    >>> bell(1, 2.5)
    2.5
    >>> bell(2, 2.5)
    8.75

Evaluation for arbitrary complex arguments::

    >>> bell(5.75+1j, 2-3j)
    (-10767.71345136587098445143 - 15449.55065599872579097221j)

The first few Bell polynomials::

    >>> for k in range(7):
    ...     nprint(taylor(lambda x: bell(k,x), 0, k))
    ...
    [1.0]
    [0.0, 1.0]
    [0.0, 1.0, 1.0]
    [0.0, 1.0, 3.0, 1.0]
    [0.0, 1.0, 7.0, 6.0, 1.0]
    [0.0, 1.0, 15.0, 25.0, 10.0, 1.0]
    [0.0, 1.0, 31.0, 90.0, 65.0, 15.0, 1.0]

The first few Bell numbers and complementary Bell numbers::

    >>> [int(bell(k)) for k in range(10)]
    [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147]
    >>> [int(bell(k,-1)) for k in range(10)]
    [1, -1, 0, 1, 1, -2, -9, -9, 50, 267]

Large Bell numbers::

    >>> mp.dps = 50
    >>> bell(50)
    185724268771078270438257767181908917499221852770.0
    >>> bell(50,-1)
    -29113173035759403920216141265491160286912.0

Some even larger values::

    >>> mp.dps = 25
    >>> bell(1000,-1)
    -1.237132026969293954162816e+1869
    >>> bell(1000)
    2.989901335682408421480422e+1927
    >>> bell(1000,2)
    6.591553486811969380442171e+1987
    >>> bell(1000,100.5)
    9.101014101401543575679639e+2529

A determinant identity satisfied by Bell numbers::

    >>> mp.dps = 15
    >>> N = 8
    >>> det([[bell(k+j) for j in range(N)] for k in range(N)])
    125411328000.0
    >>> superfac(N-1)
    125411328000.0

**References**

1. http://mathworld.wolfram.com/DobinskisFormula.html

"""

polyexp = r"""
Evaluates the polyexponential function, defined for arbitrary
complex `s`, `z` by the series

.. math ::

    E_s(z) = \sum_{k=1}^{\infty} \frac{k^s}{k!} z^k.

`E_s(z)` is constructed from the exponential function analogously
to how the polylogarithm is constructed from the ordinary
logarithm; as a function of `s` (with `z` fixed), `E_s` is an L-series
It is an entire function of both `s` and `z`.

The polyexponential function provides a generalization of the
Bell polynomials `B_n(x)` (see :func:`bell`) to noninteger orders `n`.
In terms of the Bell polynomials,

.. math ::

    E_s(z) = e^z B_s(z) - \mathrm{sinc}(\pi s).

Note that `B_n(x)` and `e^{-x} E_n(x)` are identical if `n`
is a nonzero integer, but not otherwise. In particular, they differ
at `n = 0`.

**Examples**

Evaluating a series::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> nsum(lambda k: sqrt(k)/fac(k), [1,inf])
    2.101755547733791780315904
    >>> polyexp(0.5,1)
    2.101755547733791780315904

Evaluation for arbitrary arguments::

    >>> polyexp(-3-4j, 2.5+2j)
    (2.351660261190434618268706 + 1.202966666673054671364215j)

Evaluation is accurate for tiny function values::

    >>> polyexp(4, -100)
    3.499471750566824369520223e-36

If `n` is a nonpositive integer, `E_n` reduces to a special
instance of the hypergeometric function `\,_pF_q`::

    >>> n = 3
    >>> x = pi
    >>> polyexp(-n,x)
    4.042192318847986561771779
    >>> x*hyper([1]*(n+1), [2]*(n+1), x)
    4.042192318847986561771779

"""

cyclotomic = r"""
Evaluates the cyclotomic polynomial `\Phi_n(x)`, defined by

.. math ::

    \Phi_n(x) = \prod_{\zeta} (x - \zeta)

where `\zeta` ranges over all primitive `n`-th roots of unity
(see :func:`unitroots`). An equivalent representation, used
for computation, is

.. math ::

    \Phi_n(x) = \prod_{d\mid n}(x^d-1)^{\mu(n/d)} = \Phi_n(x)

where `\mu(m)` denotes the Moebius function. The cyclotomic
polynomials are integer polynomials, the first of which can be
written explicitly as

.. math ::

    \Phi_0(x) = 1

    \Phi_1(x) = x - 1

    \Phi_2(x) = x + 1

    \Phi_3(x) = x^3 + x^2 + 1

    \Phi_4(x) = x^2 + 1

    \Phi_5(x) = x^4 + x^3 + x^2 + x + 1

    \Phi_6(x) = x^2 - x + 1

**Examples**

The coefficients of low-order cyclotomic polynomials can be recovered
using Taylor expansion::

    >>> from mpmath import *
    >>> mp.dps = 15; mp.pretty = True
    >>> for n in range(9):
    ...     p = chop(taylor(lambda x: cyclotomic(n,x), 0, 10))
    ...     print n,; nprint(p[:10+1-p[::-1].index(1)])
    ...
    0 [1.0]
    1 [-1.0, 1.0]
    2 [1.0, 1.0]
    3 [1.0, 1.0, 1.0]
    4 [1.0, 0.0, 1.0]
    5 [1.0, 1.0, 1.0, 1.0, 1.0]
    6 [1.0, -1.0, 1.0]
    7 [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    8 [1.0, 0.0, 0.0, 0.0, 1.0]

The definition as a product over primitive roots may be checked
by computing the product explicitly (for a real argument, this
method will generally introduce numerical noise in the imaginary
part)::

    >>> mp.dps = 25
    >>> z = 3+4j
    >>> cyclotomic(10, z)
    (-419.0 - 360.0j)
    >>> fprod(z-r for r in unitroots(10, primitive=True))
    (-419.0 - 360.0j)
    >>> z = 3
    >>> cyclotomic(10, z)
    61.0
    >>> fprod(z-r for r in unitroots(10, primitive=True))
    (61.0 - 3.146045605088568607055454e-25j)

Up to permutation, the roots of a given cyclotomic polynomial
can be checked to agree with the list of primitive roots::

    >>> p = taylor(lambda x: cyclotomic(6,x), 0, 6)[:3]
    >>> for r in polyroots(p[::-1]):
    ...     print r
    ...
    (0.5 - 0.8660254037844386467637232j)
    (0.5 + 0.8660254037844386467637232j)
    >>>
    >>> for r in unitroots(6, primitive=True):
    ...     print r
    ...
    (0.5 + 0.8660254037844386467637232j)
    (0.5 - 0.8660254037844386467637232j)

"""

meijerg = r"""
Evaluates the Meijer G-function, defined as

.. math ::

    G^{m,n}_{p,q} \left( \left. \begin{matrix}
         a_1, \dots, a_n ; a_{n+1} \dots a_p \\
         b_1, \dots, b_m ; b_{m+1} \dots b_q
    \end{matrix}\; \right| \; z ; r \right) =
    \frac{1}{2 \pi i} \int_L
    \frac{\prod_{j=1}^m \Gamma(b_j+s) \prod_{j=1}^n\Gamma(1-a_j-s)}
         {\prod_{j=n+1}^{p}\Gamma(a_j+s) \prod_{j=m+1}^q \Gamma(1-b_j-s)}
         z^{-s/r} ds

for an appropriate choice of the contour `L` (see references).

There are `p` elements `a_j`.
The argument *a_s* should be a pair of lists, the first containing the
`n` elements `a_1, \ldots, a_n` and the second containing
the `p-n` elements `a_{n+1}, \ldots a_p`.

There are `q` elements `b_j`.
The argument *b_s* should be a pair of lists, the first containing the
`m` elements `b_1, \ldots, b_m` and the second containing
the `q-m` elements `b_{m+1}, \ldots b_q`.

The implicit tuple `(m, n, p, q)` constitutes the order or degree of the
Meijer G-function, and is determined by the lengths of the coefficient
vectors. Confusingly, the indices in this tuple appear in a different order
from the coefficients, but this notation is standard. The many examples
given below should hopefully clear up any potential confusion.

**Algorithm**

The Meijer G-function is evaluated by rewriting it as a combination of `m`
hypergeometric series each of degree `\,_pF_{q-1}`. Currently,
evaluation will only work if `\,_pF_{q-1}` can be evaluated via
:func:`hyper` for the given `z`, i.e. if either `z` is small enough
for direct convergence or if the hypergeometric function in question is
implemented with the appropriate analytic continuation / asymptotic
expansion for large `z`.

Keyword arguments are forwarded to :func:`hypercomb`. If :func:`meijerg`
fails to return an accurate value, passing ``check_cancellation=True``
may help.

**Examples**

Many standard functions are special cases of the Meijer G-function
(possibly rescaled and/or with branch cut corrections). We define
some test parameters::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> a = mpf(0.75)
    >>> b = mpf(1.5)
    >>> z = mpf(2.25)

The exponential function:
`e^z = G^{1,0}_{0,1} \left( \left. \begin{matrix} - \\ 0 \end{matrix} \;
\right| \; -z \right)`

    >>> meijerg([[],[]], [[0],[]], -z)
    9.487735836358525720550369
    >>> exp(z)
    9.487735836358525720550369

The natural logarithm:
`\log(1+z) = G^{1,2}_{2,2} \left( \left. \begin{matrix} 1, 1 \\ 1, 0
\end{matrix} \; \right| \; -z \right)`

    >>> meijerg([[1,1],[]], [[1],[0]], z)
    1.178654996341646117219023
    >>> log(1+z)
    1.178654996341646117219023

A rational function:
`\frac{z}{z+1} = G^{1,2}_{2,2} \left( \left. \begin{matrix} 1, 1 \\ 1, 1
\end{matrix} \; \right| \; z \right)`

    >>> meijerg([[1,1],[]], [[1],[1]], z)
    0.6923076923076923076923077
    >>> z/(z+1)
    0.6923076923076923076923077

The sine and cosine functions:

`\frac{1}{\sqrt \pi} \sin(2 \sqrt z) = G^{1,0}_{0,2} \left( \left. \begin{matrix}
- \\ \frac{1}{2}, 0 \end{matrix} \; \right| \; z \right)`

`\frac{1}{\sqrt \pi} \cos(2 \sqrt z) = G^{1,0}_{0,2} \left( \left. \begin{matrix}
- \\ 0, \frac{1}{2} \end{matrix} \; \right| \; z \right)`

    >>> meijerg([[],[]], [[0.5],[0]], (z/2)**2)
    0.4389807929218676682296453
    >>> sin(z)/sqrt(pi)
    0.4389807929218676682296453
    >>> meijerg([[],[]], [[0],[0.5]], (z/2)**2)
    -0.3544090145996275423331762
    >>> cos(z)/sqrt(pi)
    -0.3544090145996275423331762

Bessel functions:

`J_a(2 \sqrt z) = G^{1,0}_{0,2} \left( \left.
\begin{matrix} - \\ \frac{a}{2}, -\frac{a}{2}
\end{matrix} \; \right| \; z \right)`

`Y_a(2 \sqrt z) = G^{2,0}_{1,3} \left( \left.
\begin{matrix} \frac{-a-1}{2} \\ \frac{a}{2}, -\frac{a}{2}, \frac{-a-1}{2}
\end{matrix} \; \right| \; z \right)`

`(-z)^{a/2} z^{-a/2} I_a(2 \sqrt z) = G^{1,0}_{0,2} \left( \left.
\begin{matrix} - \\ \frac{a}{2}, -\frac{a}{2}
\end{matrix} \; \right| \; -z \right)`

`2 K_a(2 \sqrt z) = G^{2,0}_{0,2} \left( \left.
\begin{matrix} - \\ \frac{a}{2}, -\frac{a}{2}
\end{matrix} \; \right| \; z \right)`

As the example with the Bessel *I* function shows, a branch
factor is required for some arguments when inverting the square root.

    >>> meijerg([[],[]], [[a/2],[-a/2]], (z/2)**2)
    0.5059425789597154858527264
    >>> besselj(a,z)
    0.5059425789597154858527264
    >>> meijerg([[],[(-a-1)/2]], [[a/2,-a/2],[(-a-1)/2]], (z/2)**2)
    0.1853868950066556941442559
    >>> bessely(a, z)
    0.1853868950066556941442559
    >>> meijerg([[],[]], [[a/2],[-a/2]], -(z/2)**2)
    (0.8685913322427653875717476 + 2.096964974460199200551738j)
    >>> (-z)**(a/2) / z**(a/2) * besseli(a, z)
    (0.8685913322427653875717476 + 2.096964974460199200551738j)
    >>> 0.5*meijerg([[],[]], [[a/2,-a/2],[]], (z/2)**2)
    0.09334163695597828403796071
    >>> besselk(a,z)
    0.09334163695597828403796071

Error functions:

`\sqrt{\pi} z^{2(a-1)} \mathrm{erfc}(z) = G^{2,0}_{1,2} \left( \left.
\begin{matrix} a \\ a-1, a-\frac{1}{2}
\end{matrix} \; \right| \; z, \frac{1}{2} \right)`

    >>> meijerg([[],[a]], [[a-1,a-0.5],[]], z, 0.5)
    0.00172839843123091957468712
    >>> sqrt(pi) * z**(2*a-2) * erfc(z)
    0.00172839843123091957468712

A Meijer G-function of higher degree, (1,1,2,3):

    >>> meijerg([[a],[b]], [[a],[b,a-1]], z)
    1.55984467443050210115617
    >>> sin((b-a)*pi)/pi*(exp(z)-1)*z**(a-1)
    1.55984467443050210115617

A Meijer G-function of still higher degree, (4,1,2,4), that can
be expanded as a messy combination of exponential integrals:

    >>> meijerg([[a],[2*b-a]], [[b,a,b-0.5,-1-a+2*b],[]], z)
    0.3323667133658557271898061
    >>> chop(4**(a-b+1)*sqrt(pi)*gamma(2*b-2*a)*z**a*\
    ...     expint(2*b-2*a, -2*sqrt(-z))*expint(2*b-2*a, 2*sqrt(-z)))
    0.3323667133658557271898061

**References**

1. http://en.wikipedia.org/wiki/Meijer_G-function

2. http://mathworld.wolfram.com/MeijerG-Function.html

3. http://functions.wolfram.com/HypergeometricFunctions/MeijerG/

4. http://functions.wolfram.com/HypergeometricFunctions/MeijerG1/

"""

clsin = r"""
Computes the Clausen sine function, defined formally by the series

.. math ::

    \mathrm{Cl}_s(z) = \sum_{k=1}^{\infty} \frac{\sin(kz)}{k^s}.

The special case `\mathrm{Cl}_2(z)` (i.e. ``clsin(2,z)``) is the classical
"Clausen function". More generally, the Clausen function is defined for
complex `s` and `z`, even when the series does not converge. The
Clausen function is related to the polylogarithm (:func:`polylog`) as

.. math ::

    \mathrm{Cl}_s(z) = \frac{1}{2i}\left(\mathrm{Li}_s\left(e^{iz}\right) -
                       \mathrm{Li}_s\left(e^{-iz}\right)\right)

    = \mathrm{Im}\left[\mathrm{Li}_s(e^{iz})\right] \quad (s, z \in \mathbb{R}),

and this representation can be taken to provide the analytic continuation of the
series. The complementary function :func:`clcos` gives the corresponding
cosine sum.

**Examples**

Evaluation for arbitrarily chosen `s` and `z`::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> s, z = 3, 4
    >>> clsin(s, z); nsum(lambda k: sin(z*k)/k**s, [1,inf])
    -0.6533010136329338746275795
    -0.6533010136329338746275795

Using `z + \pi` instead of `z` gives an alternating series::

    >>> clsin(s, z+pi)
    0.8860032351260589402871624
    >>> nsum(lambda k: (-1)**k*sin(z*k)/k**s, [1,inf])
    0.8860032351260589402871624

With `s = 1`, the sum can be expressed in closed form
using elementary functions::

    >>> z = 1 + sqrt(3)
    >>> clsin(1, z)
    0.2047709230104579724675985
    >>> chop((log(1-exp(-j*z)) - log(1-exp(j*z)))/(2*j))
    0.2047709230104579724675985
    >>> nsum(lambda k: sin(k*z)/k, [1,inf])
    0.2047709230104579724675985

The classical Clausen function `\mathrm{Cl}_2(\theta)` gives the
value of the integral `\int_0^{\theta} -\ln(2\sin(x/2)) dx` for
`0 < \theta < 2 \pi`::

    >>> cl2 = lambda t: clsin(2, t)
    >>> cl2(3.5)
    -0.2465045302347694216534255
    >>> -quad(lambda x: ln(2*sin(0.5*x)), [0, 3.5])
    -0.2465045302347694216534255

This function is symmetric about `\theta = \pi` with zeros and extreme
points::

    >>> cl2(0); cl2(pi/3); chop(cl2(pi)); cl2(5*pi/3); chop(cl2(2*pi))
    0.0
    1.014941606409653625021203
    0.0
    -1.014941606409653625021203
    0.0

Catalan's constant is a special value::

    >>> cl2(pi/2)
    0.9159655941772190150546035
    >>> +catalan
    0.9159655941772190150546035

The Clausen sine function can be expressed in closed form when
`s` is an odd integer (becoming zero when `s` < 0)::

    >>> z = 1 + sqrt(2)
    >>> clsin(1, z); (pi-z)/2
    0.3636895456083490948304773
    0.3636895456083490948304773
    >>> clsin(3, z); pi**2/6*z - pi*z**2/4 + z**3/12
    0.5661751584451144991707161
    0.5661751584451144991707161
    >>> clsin(-1, z)
    0.0
    >>> clsin(-3, z)
    0.0

It can also be expressed in closed form for even integer `s \le 0`,
providing a finite sum for series such as
`\sin(z) + \sin(2z) + \sin(3z) + \ldots`::

    >>> z = 1 + sqrt(2)
    >>> clsin(0, z)
    0.1903105029507513881275865
    >>> cot(z/2)/2
    0.1903105029507513881275865
    >>> clsin(-2, z)
    -0.1089406163841548817581392
    >>> -cot(z/2)*csc(z/2)**2/4
    -0.1089406163841548817581392

Evaluation for complex `s`, `z` in a nonconvergent case::

    >>> s, z = -1-j, 1+2j
    >>> clsin(s, z)
    (-0.593079480117379002516034 + 0.9038644233367868273362446j)
    >>> extraprec(20)(nsum)(lambda k: sin(k*z)/k**s, [1,inf])
    (-0.593079480117379002516034 + 0.9038644233367868273362446j)

"""

clcos = r"""
Computes the Clausen cosine function, defined formally by the series

.. math ::

    \mathrm{\widetilde{Cl}}_s(z) = \sum_{k=1}^{\infty} \frac{\cos(kz)}{k^s}.

This function is complementary to the Clausen sine function
:func:`clsin`. In terms of the polylogarithm,

.. math ::

    \mathrm{\widetilde{Cl}}_s(z) =
        \frac{1}{2}\left(\mathrm{Li}_s\left(e^{iz}\right) +
        \mathrm{Li}_s\left(e^{-iz}\right)\right)

    = \mathrm{Re}\left[\mathrm{Li}_s(e^{iz})\right] \quad (s, z \in \mathbb{R}).

**Examples**

Evaluation for arbitrarily chosen `s` and `z`::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> s, z = 3, 4
    >>> clcos(s, z); nsum(lambda k: cos(z*k)/k**s, [1,inf])
    -0.6518926267198991308332759
    -0.6518926267198991308332759

Using `z + \pi` instead of `z` gives an alternating series::

    >>> s, z = 3, 0.5
    >>> clcos(s, z+pi)
    -0.8155530586502260817855618
    >>> nsum(lambda k: (-1)**k*cos(z*k)/k**s, [1,inf])
    -0.8155530586502260817855618

With `s = 1`, the sum can be expressed in closed form
using elementary functions::

    >>> z = 1 + sqrt(3)
    >>> clcos(1, z)
    -0.6720334373369714849797918
    >>> chop(-0.5*(log(1-exp(j*z))+log(1-exp(-j*z))))
    -0.6720334373369714849797918
    >>> -log(abs(2*sin(0.5*z)))    # Equivalent to above when z is real
    -0.6720334373369714849797918
    >>> nsum(lambda k: cos(k*z)/k, [1,inf])
    -0.6720334373369714849797918

It can also be expressed in closed form when `s` is an even integer.
For example,

    >>> clcos(2,z)
    -0.7805359025135583118863007
    >>> pi**2/6 - pi*z/2 + z**2/4
    -0.7805359025135583118863007

The case `s = 0` gives the renormalized sum of
`\cos(z) + \cos(2z) + \cos(3z) + \ldots` (which happens to be the same for
any value of `z`)::

    >>> clcos(0, z)
    -0.5
    >>> nsum(lambda k: cos(k*z), [1,inf])
    -0.5

Also the sums

.. math ::

    \cos(z) + 2\cos(2z) + 3\cos(3z) + \ldots

and

.. math ::

    \cos(z) + 2^n \cos(2z) + 3^n \cos(3z) + \ldots

for higher integer powers `n = -s` can be done in closed form. They are zero
when `n` is positive and even (`s` negative and even)::

    >>> clcos(-1, z); 1/(2*cos(z)-2)
    -0.2607829375240542480694126
    -0.2607829375240542480694126
    >>> clcos(-3, z); (2+cos(z))*csc(z/2)**4/8
    0.1472635054979944390848006
    0.1472635054979944390848006
    >>> clcos(-2, z); clcos(-4, z); clcos(-6, z)
    0.0
    0.0
    0.0

With `z = \pi`, the series reduces to that of the Riemann zeta function
(more generally, if `z = p \pi/q`, it is a finite sum over Hurwitz zeta
function values)::

    >>> clcos(2.5, 0); zeta(2.5)
    1.34148725725091717975677
    1.34148725725091717975677
    >>> clcos(2.5, pi); -altzeta(2.5)
    -0.8671998890121841381913472
    -0.8671998890121841381913472

Evaluation for complex `s`, `z` in a nonconvergent case::

    >>> s, z = -1-j, 1+2j
    >>> clcos(s, z)
    (0.9407430121562251476136807 + 0.715826296033590204557054j)
    >>> extraprec(15)(nsum)(lambda k: cos(k*z)/k**s, [1,inf])
    (0.9407430121562251476136807 + 0.7158262960335902045570541j)

"""

whitm = r"""
Evaluates the Whittaker function `M(k,m,z)`, which gives a solution
to the Whittaker differential equation

.. math ::

    \frac{d^2f}{dz^2} + \left(-\frac{1}{4}+\frac{k}{z}+
      \frac{(\frac{1}{4}-m^2)}{z^2}\right) f = 0.

A second solution is given by :func:`whitw`.

The Whittaker functions are defined in Abramowitz & Stegun, section 13.1.
They are alternate forms of the confluent hypergeometric functions
`\,_1F_1` and `U`:

.. math ::

    M(k,m,z) = e^{-\frac{1}{2}z} z^{\frac{1}{2}+m}
        \,_1F_1(\tfrac{1}{2}+m-k, 1+2m, z)

    W(k,m,z) = e^{-\frac{1}{2}z} z^{\frac{1}{2}+m}
        U(\tfrac{1}{2}+m-k, 1+2m, z).

**Examples**

Evaluation for arbitrary real and complex arguments is supported::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> whitm(1, 1, 1)
    0.7302596799460411820509668
    >>> whitm(1, 1, -1)
    (0.0 - 1.417977827655098025684246j)
    >>> whitm(j, j/2, 2+3j)
    (3.245477713363581112736478 - 0.822879187542699127327782j)
    >>> whitm(2, 3, 100000)
    4.303985255686378497193063e+21707

Evaluation at zero::

    >>> whitm(1,-1,0); whitm(1,-0.5,0); whitm(1,0,0)
    +inf
    nan
    0.0

We can verify that :func:`whitm` numerically satisfies the
differential equation for arbitrarily chosen values::

    >>> k = mpf(0.25)
    >>> m = mpf(1.5)
    >>> f = lambda z: whitm(k,m,z)
    >>> for z in [-1, 2.5, 3, 1+2j]:
    ...     chop(diff(f,z,2) + (-0.25 + k/z + (0.25-m**2)/z**2)*f(z))
    ...
    0.0
    0.0
    0.0
    0.0

An integral involving both :func:`whitm` and :func:`whitmw`,
verifying evaluation along the real axis::

    >>> quad(lambda x: exp(-x)*whitm(3,2,x)*whitw(1,-2,x), [0,inf])
    3.438869842576800225207341
    >>> 128/(21*sqrt(pi))
    3.438869842576800225207341

"""

whitw = r"""
Evaluates the Whittaker function `W(k,m,z)`, which gives a second
solution to the Whittaker differential equation. (See :func:`whitm`.)

**Examples**

Evaluation for arbitrary real and complex arguments is supported::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> whitw(1, 1, 1)
    1.19532063107581155661012
    >>> whitw(1, 1, -1)
    (-0.9424875979222187313924639 - 0.2607738054097702293308689j)
    >>> whitw(j, j/2, 2+3j)
    (0.1782899315111033879430369 - 0.01609578360403649340169406j)
    >>> whitw(2, 3, 100000)
    1.887705114889527446891274e-21705
    >>> whitw(-1, -1, 100)
    1.905250692824046162462058e-24

Evaluation at zero::

    >>> for m in [-1, -0.5, 0, 0.5, 1]:
    ...     whitw(1, m, 0)
    ...
    +inf
    nan
    0.0
    nan
    +inf

We can verify that :func:`whitw` numerically satisfies the
differential equation for arbitrarily chosen values::

    >>> k = mpf(0.25)
    >>> m = mpf(1.5)
    >>> f = lambda z: whitw(k,m,z)
    >>> for z in [-1, 2.5, 3, 1+2j]:
    ...     chop(diff(f,z,2) + (-0.25 + k/z + (0.25-m**2)/z**2)*f(z))
    ...
    0.0
    0.0
    0.0
    0.0

"""

ber = r"""
Computes the Kelvin function ber, which for real arguments gives the real part
of the Bessel J function of a rotated argument

.. math ::

    J_n\left(x e^{3\pi i/4}\right) = \mathrm{ber}_n(x) + i \mathrm{bei}_n(x).

The imaginary part is given by :func:`bei`.

**Examples**

Verifying the defining relation::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> n, x = 2, 3.5
    >>> ber(n,x)
    1.442338852571888752631129
    >>> bei(n,x)
    -0.948359035324558320217678
    >>> besselj(n, x*root(1,8,3))
    (1.442338852571888752631129 - 0.948359035324558320217678j)

The ber and bei functions are also defined by analytic continuation
for complex arguments::

    >>> ber(1+j, 2+3j)
    (4.675445984756614424069563 - 15.84901771719130765656316j)
    >>> bei(1+j, 2+3j)
    (15.83886679193707699364398 + 4.684053288183046528703611j)

"""

bei = r"""
Computes the Kelvin function bei, which for real arguments gives the
imaginary part of the Bessel J function of a rotated argument.
See :func:`ber`.
"""

ker = r"""
Computes the Kelvin function ker, which for real arguments gives the real part
of the (rescaled) Bessel K function of a rotated argument

.. math ::

    e^{-\pi i/2} K_n\left(x e^{3\pi i/4}\right) = \mathrm{ker}_n(x) + i \mathrm{kei}_n(x).

The imaginary part is given by :func:`kei`.

**Examples**

Verifying the defining relation::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> n, x = 2, 4.5
    >>> ker(n,x)
    0.02542895201906369640249801
    >>> kei(n,x)
    -0.02074960467222823237055351
    >>> exp(-n*pi*j/2) * besselk(n, x*root(1,8,1))
    (0.02542895201906369640249801 - 0.02074960467222823237055351j)

The ker and kei functions are also defined by analytic continuation
for complex arguments::

    >>> ker(1+j, 3+4j)
    (1.586084268115490421090533 - 2.939717517906339193598719j)
    >>> kei(1+j, 3+4j)
    (-2.940403256319453402690132 - 1.585621643835618941044855j)

"""

kei = r"""
Computes the Kelvin function kei, which for real arguments gives the
imaginary part of the (rescaled) Bessel K function of a rotated argument.
See :func:`ker`.
"""

struveh = r"""
Gives the Struve function

.. math ::

    \,\mathbf{H}_n(z) = 
    \sum_{k=0}^\infty \frac{(-1)^k}{\Gamma(k+\frac{3}{2})
        \Gamma(k+n+\frac{3}{2})} {\left({\frac{z}{2}}\right)}^{2k+n+1}

which is a solution to the Struve differential equation

.. math ::

    z^2 f''(z) + z f'(z) + (z^2-n^2) f(z) = \frac{2 z^{n+1}}{\pi (2n-1)!!}.

**Examples**

Evaluation for arbitrary real and complex arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> struveh(0, 3.5)
    0.3608207733778295024977797
    >>> struveh(-1, 10)
    -0.255212719726956768034732
    >>> struveh(1, -100.5)
    0.5819566816797362287502246
    >>> struveh(2.5, 10000000000000)
    3153915652525200060.308937
    >>> struveh(2.5, -10000000000000)
    (0.0 - 3153915652525200060.308937j)
    >>> struveh(1+j, 1000000+4000000j)
    (-3.066421087689197632388731e+1737173 - 1.596619701076529803290973e+1737173j)

A Struve function of half-integer order is elementary; for example:

    >>> z = 3
    >>> struveh(0.5, 3)
    0.9167076867564138178671595
    >>> sqrt(2/(pi*z))*(1-cos(z))
    0.9167076867564138178671595

Numerically verifying the differential equation::

    >>> z = mpf(4.5)
    >>> n = 3
    >>> f = lambda z: struveh(n,z)
    >>> lhs = z**2*diff(f,z,2) + z*diff(f,z) + (z**2-n**2)*f(z)
    >>> rhs = 2*z**(n+1)/fac2(2*n-1)/pi
    >>> lhs
    17.40359302709875496632744
    >>> rhs
    17.40359302709875496632744

"""

struvel = r"""
Gives the modified Struve function

.. math ::

    \,\mathbf{L}_n(z) = -i e^{-n\pi i/2} \mathbf{H}_n(i z)

which solves to the modified Struve differential equation

.. math ::

    z^2 f''(z) + z f'(z) - (z^2+n^2) f(z) = \frac{2 z^{n+1}}{\pi (2n-1)!!}.

**Examples**

Evaluation for arbitrary real and complex arguments::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> struvel(0, 3.5)
    7.180846515103737996249972
    >>> struvel(-1, 10)
    2670.994904980850550721511
    >>> struvel(1, -100.5)
    1.757089288053346261497686e+42
    >>> struvel(2.5, 10000000000000)
    4.160893281017115450519948e+4342944819025
    >>> struvel(2.5, -10000000000000)
    (0.0 - 4.160893281017115450519948e+4342944819025j)
    >>> struvel(1+j, 700j)
    (-0.1721150049480079451246076 + 0.1240770953126831093464055j)
    >>> struvel(1+j, 1000000+4000000j)
    (-2.973341637511505389128708e+434290 - 5.164633059729968297147448e+434290j)

Numerically verifying the differential equation::

    >>> z = mpf(3.5)
    >>> n = 3
    >>> f = lambda z: struvel(n,z)
    >>> lhs = z**2*diff(f,z,2) + z*diff(f,z) - (z**2+n**2)*f(z)
    >>> rhs = 2*z**(n+1)/fac2(2*n-1)/pi
    >>> lhs
    6.368850306060678353018165
    >>> rhs
    6.368850306060678353018165
"""

appellf1 = r"""
Gives the Appell F1 hypergeometric function of two variables,

.. math ::

    F_1(a,b_1,b_2,c,z_1,z_2) = \sum_{m=0}^{\infty}
        \sum_{n=0}^{\infty}
        \frac{(a)_{m+n} (b_1)_m (b_2)_n}{m! \,n! \,(c)_{m+n}} z_1^m z_2^n.

This series is only generally convergent when `|z_1| < 1` and `|z_2| < 1`,
although :func:`appellf1` can evaluate the continuation
in many cases.

**Examples**

Evaluation is supported for real and complex parameters::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> appellf1(1,0,0.5,1,0.5,0.25)
    1.154700538379251529018298
    >>> appellf1(1,1+j,0.5,1,0.5,0.5j)
    (1.138403860350148085179415 + 1.510544741058517621110615j)

For some integer parameters, the F1 series reduces to a polynomial::

    >>> appellf1(2,-4,-3,1,2,5)
    -816.0
    >>> appellf1(-5,1,2,1,4,5)
    -20528.0

The analytic continuation with respect to either `z_1` or `z_2`,
and sometimes with respect to both, can be evaluated::

    >>> appellf1(2,3,4,5,100,0.5)
    (0.0006231042714165329279738662 + 0.0000005769149277148425774499857j)
    >>> appellf1('1.1', '0.3', '0.2+2j', '0.4', '0.2', 1.5+3j)
    (-0.1782604566893954897128702 + 0.002472407104546216117161499j)
    >>> appellf1(1,2,3,4,10,12)
    -0.07122993830066776374929313

For certain arguments, F1 reduces to an ordinary hypergeometric function::

    >>> appellf1(1,2,3,5,0.5,0.25)
    1.547902270302684019335555
    >>> 4*hyp2f1(1,2,5,'1/3')/3
    1.547902270302684019335555
    >>> appellf1(1,2,3,4,0,1.5)
    (-1.717202506168937502740238 - 2.792526803190927323077905j)
    >>> hyp2f1(1,3,4,1.5)
    (-1.717202506168937502740238 - 2.792526803190927323077905j)

The Appell F1 function allows for closed-form evaluation of various
integrals, such as any integral of the form
`\int x^r (x+a)^p (x+b)^q dx`::

    >>> def integral(a,b,p,q,r,x1,x2):
    ...     a,b,p,q,r,x1,x2 = map(mpmathify, [a,b,p,q,r,x1,x2])
    ...     f = lambda x: x**r * (x+a)**p * (x+b)**q
    ...     def F(x):
    ...         v = x**(r+1)/(r+1) * (a+x)**p * (b+x)**q
    ...         v *= (1+x/a)**(-p)
    ...         v *= (1+x/b)**(-q)
    ...         v *= appellf1(r+1,-p,-q,2+r,-x/a,-x/b)
    ...         return v
    ...     print "Num. quad:", quad(f, [x1,x2])
    ...     print "Appell F1:", F(x2)-F(x1)
    ...
    >>> integral('1/5','4/3','-2','3','1/2',0,1)
    Num. quad: 9.073335358785776206576981
    Appell F1: 9.073335358785776206576981
    >>> integral('3/2','4/3','-2','3','1/2',0,1)
    Num. quad: 1.092829171999626454344678
    Appell F1: 1.092829171999626454344678
    >>> integral('3/2','4/3','-2','3','1/2',12,25)
    Num. quad: 1106.323225040235116498927
    Appell F1: 1106.323225040235116498927

Also incomplete elliptic integrals fall into this category [1]::

    >>> def E(z, m):
    ...     if (pi/2).ae(z):
    ...         return ellipe(m)
    ...     return 2*round(re(z)/pi)*ellipe(m) + mpf(-1)**round(re(z)/pi)*\
    ...         sin(z)*appellf1(0.5,0.5,-0.5,1.5,sin(z)**2,m*sin(z)**2)
    ...
    >>> z, m = 1, 0.5
    >>> E(z,m); quad(lambda t: sqrt(1-m*sin(t)**2), [0,pi/4,3*pi/4,z])
    0.9273298836244400669659042
    0.9273298836244400669659042
    >>> z, m = 3, 2
    >>> E(z,m); quad(lambda t: sqrt(1-m*sin(t)**2), [0,pi/4,3*pi/4,z])
    (1.057495752337234229715836 + 1.198140234735592207439922j)
    (1.057495752337234229715836 + 1.198140234735592207439922j)

**References**

1. http://functions.wolfram.com/EllipticIntegrals/EllipticE2/26/01/
"""

hurwitz = r"""
Computes the Hurwitz zeta function

.. math ::

    \zeta(s,a) = \sum_{k=0}^\infty \frac{1}{(a+k)^s}.

With `a = 1`, this reduces to the ordinary Riemann zeta function.
Optionally, ``hurwitz(s, a, n)`` computes
the `n`-th derivative with respect to `s`,

.. math ::

    \zeta^{(n)}(s,a) = (-1)^n \sum_{k=0}^\infty \frac{\log^n(a+k)}{(a+k)^s}.

Although these series only converge for `\Re(s) > 1`, the Hurwitz
zeta function is defined through analytic continuation for arbitrary
complex `s \ne 1` (`s = 1` is a pole). The implementation uses Euler-Maclaurin
summation along with reflection formulas in some cases when `\Re(s) < 0`.

The parameter `a` is usually a rational number `a = p/q`, and may be specified
as such by passing an integer tuple `(p, q)`. Evaluation is supported for
arbitrary complex `a`, but may be slow and/or inaccurate when `\Re(s) < 0` for
nonrational `a` or when computing derivatives.

**Examples**

Some basic evaluations:

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> hurwitz(2, 3); -5./4 + pi**2/6
    0.3949340668482264364724152
    0.3949340668482264364724152
    >>> hurwitz(2, (3,4)); pi**2 - 8*catalan
    2.541879647671606498397663
    2.541879647671606498397663

For positive integer values of `s`, the Hurwitz zeta function can be evaluated
using polygamma functions::

    >>> hurwitz(4, (1,5)); psi(3, '1/5')/6
    625.5408324774542966919938
    625.5408324774542966919938

Some Riemann zeta function values::

    >>> hurwitz(2)
    1.644934066848226436472415
    >>> 1-sum((hurwitz(k)-1)/k for k in range(2,85))
    0.5772156649015328606065121
    >>> +euler
    0.5772156649015328606065121
    >>> hurwitz(inf)
    1.0

Evaluation on the critical line::

    >>> findroot(hurwitz, 0.5+14j)
    (0.5 + 14.13472514173469379045725j)

A derivative identity; evaluation for `s < 1`::

    >>> hurwitz(0, 3+4j, 1)
    (-2.675565317808456852310934 + 4.742664438034657928194889j)
    >>> loggamma(3+4j) - ln(2*pi)/2
    (-2.675565317808456852310934 + 4.742664438034657928194889j)

A high order derivative::

    >>> hurwitz(2, 1, 20)
    2432902008176640000.000242

A fourth derivative at a complex value::

    >>> hurwitz(3+4j, 5.5+2j, 4)
    (-0.140075548947797130681075 - 0.3109263360275413251313634j)

Generating a Taylor series at `s = 2` using derivatives::

    >>> for k in range(11): print hurwitz(2,1,k)/fac(k), "*", "(s-2)^%i" % k
    ...
    1.644934066848226436472415 * (s-2)^0
    -0.9375482543158437537025741 * (s-2)^1
    0.9946401171494505117104293 * (s-2)^2
    -1.000024300473840810940657 * (s-2)^3
    1.000061933072352565457512 * (s-2)^4
    -1.000006869443931806408941 * (s-2)^5
    1.000000173233769531820592 * (s-2)^6
    -0.9999999569989868493432399 * (s-2)^7
    0.9999999937218844508684206 * (s-2)^8
    -0.9999999996355013916608284 * (s-2)^9
    1.000000000004610645020747 * (s-2)^10

Evaluation at zero and for negative integer `s`::

    >>> hurwitz(0); zeta(0)
    -0.5
    -0.5
    >>> hurwitz(0, 10)
    -9.5
    >>> hurwitz(-1); zeta(-1)
    -0.08333333333333333333333333
    -0.08333333333333333333333333
    >>> hurwitz(-2); zeta(-2)
    0.0
    0.0
    >>> hurwitz(-2, (2,3)); mpf(1)/81
    0.01234567901234567901234568
    0.01234567901234567901234568

Evaluation for negative `s`, with rational and nonrational `a`::

    >>> hurwitz(-3+4j, (5,4))
    (0.2899236037682695182085988 + 0.06561206166091757973112783j)
    >>> hurwitz(-3.25, 1/pi)
    -0.0005117269627574430494396877
    >>> hurwitz(-3.5, pi(dps=30), 1)
    11.15636039044000329471003
    >>> hurwitz(-100.5, (8,3))
    -4.68162300487989766727122e+77
    >>> hurwitz(-10.5, (-8,3))
    (-0.01521913704446246609237979 + 29907.72510874248161608216j)
    >>> hurwitz(-1000.5, (-8,3))
    (1.031911949062334538202567e+1770 + 1.519555750556794218804724e+426j)
    >>> hurwitz(-1+j, 3+4j)
    (-16.32988355630802510888631 - 22.17706465801374033261383j)
    >>> hurwitz(-1+j, 3+4j, 2)
    (32.48985276392056641594055 - 51.11604466157397267043655j)
    >>> diff(lambda s: hurwitz(s, 3+4j), -1+j, 2)
    (32.48985276392056641594055 - 51.11604466157397267043655j)

"""

dirichlet = r"""
Evaluates the Dirichlet L-function

.. math ::

    L(s,\chi) = \sum_{k=1}^\infty \frac{\chi(k)}{k^s}.

where `\chi` is a periodic sequence of length `q` which should be supplied
in the form of a list `[\chi(0), \chi(1), \ldots, \chi(q-1)]`.
Strictly, `\chi` should be a Dirichlet character, but any periodic
sequence will work.

For example, ``dirichlet(s, [1])`` gives the ordinary
Riemann zeta function and ``dirichlet(s, [-1,1])`` gives
the alternating zeta function (Dirichlet eta function).

Also the derivative with respect to `s` (currently only a first
derivative) can be evaluated.

**Examples**

The ordinary Riemann zeta function::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> dirichlet(3, [1]); zeta(3)
    1.202056903159594285399738
    1.202056903159594285399738
    >>> dirichlet(1, [1])
    +inf

The alternating zeta function::

    >>> dirichlet(1, [-1,1]); ln(2)
    0.6931471805599453094172321
    0.6931471805599453094172321

The following defines the Dirichlet beta function
`\beta(s) = \sum_{k=0}^\infty \frac{(-1)^k}{(2k+1)^s}` and verifies
several values of this function::

    >>> B = lambda s, d=0: dirichlet(s, [0, 1, 0, -1], d)
    >>> B(0); 1./2
    0.5
    0.5
    >>> B(1); pi/4
    0.7853981633974483096156609
    0.7853981633974483096156609
    >>> B(2); +catalan
    0.9159655941772190150546035
    0.9159655941772190150546035
    >>> B(2,1); diff(B, 2)
    0.08158073611659279510291217
    0.08158073611659279510291217
    >>> B(-1,1); 2*catalan/pi
    0.5831218080616375602767689
    0.5831218080616375602767689
    >>> B(0,1); log(gamma(0.25)**2/(2*pi*sqrt(2)))
    0.3915943927068367764719453
    0.3915943927068367764719454
    >>> B(1,1); 0.25*pi*(euler+2*ln2+3*ln(pi)-4*ln(gamma(0.25)))
    0.1929013167969124293631898
    0.1929013167969124293631898

A custom L-series of period 3::

    >>> dirichlet(2, [2,0,1])
    0.7059715047839078092146831
    >>> 2*nsum(lambda k: (3*k)**-2, [1,inf]) + \
    ...   nsum(lambda k: (3*k+2)**-2, [0,inf])
    0.7059715047839078092146831

"""

coulombf = r"""
Calculates the regular Coulomb wave function

.. math ::

    F_l(\eta,z) = C_l(\eta) z^{l+1} e^{-iz} \,_1F_1(l+1-i\eta, 2l+2, 2iz)

where the normalization constant `C_l(\eta)` is as calculated by
:func:`coulombc`. This function solves the differential equation

.. math ::

    f''(z) + \left(1-\frac{2\eta}{z}-\frac{l(l+1)}{z^2}\right) f(z) = 0.

A second linearly independent solution is given by the irregular
Coulomb wave function `G_l(\eta,z)` (see :func:`coulombg`)
and thus the general solution is
`f(z) = C_1 F_l(\eta,z) + C_2 G_l(\eta,z)` for arbitrary
constants `C_1`, `C_2`.
Physically, the Coulomb wave functions give the radial solution
to the Schrodinger equation for a point particle in a `1/z` potential; `z` is
then the radius and `l`, `\eta` are quantum numbers.

The Coulomb wave functions with real parameters are defined
in Abramowitz & Stegun, section 14. However, all parameters are permitted
to be complex in this implementation (see references).

**Examples**

Evaluation is supported for arbitrary magnitudes of `z`::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> coulombf(2, 1.5, 3.5)
    0.4080998961088761187426445
    >>> coulombf(-2, 1.5, 3.5)
    0.7103040849492536747533465
    >>> coulombf(2, 1.5, '1e-10')
    4.143324917492256448770769e-33
    >>> coulombf(2, 1.5, 1000)
    0.4482623140325567050716179
    >>> coulombf(2, 1.5, 10**10)
    -0.066804196437694360046619

Verifying the differential equation::

    >>> l, eta, z = 2, 3, mpf(2.75)
    >>> A, B = 1, 2
    >>> f = lambda z: A*coulombf(l,eta,z) + B*coulombg(l,eta,z)
    >>> chop(diff(f,z,2) + (1-2*eta/z - l*(l+1)/z**2)*f(z))
    0.0

A Wronskian relation satisfied by the Coulomb wave functions::

    >>> l = 2
    >>> eta = 1.5
    >>> F = lambda z: coulombf(l,eta,z)
    >>> G = lambda z: coulombg(l,eta,z)
    >>> for z in [3.5, -1, 2+3j]:
    ...     chop(diff(F,z)*G(z) - F(z)*diff(G,z))
    ...
    1.0
    1.0
    1.0

Another Wronskian relation::

    >>> F = coulombf
    >>> G = coulombg
    >>> for z in [3.5, -1, 2+3j]:
    ...     chop(F(l-1,eta,z)*G(l,eta,z)-F(l,eta,z)*G(l-1,eta,z) - l/sqrt(l**2+eta**2))
    ...
    0.0
    0.0
    0.0

An integral identity connecting the regular and irregular wave functions::

    >>> l, eta, z = 4+j, 2-j, 5+2j
    >>> coulombf(l,eta,z) + j*coulombg(l,eta,z)
    (0.7997977752284033239714479 + 0.9294486669502295512503127j)
    >>> g = lambda t: exp(-t)*t**(l-j*eta)*(t+2*j*z)**(l+j*eta)
    >>> j*exp(-j*z)*z**(-l)/fac(2*l+1)/coulombc(l,eta)*quad(g, [0,inf])
    (0.7997977752284033239714479 + 0.9294486669502295512503127j)

Some test case with complex parameters, taken from Michel [2]::

    >>> mp.dps = 15
    >>> coulombf(1+0.1j, 50+50j, 100.156)
    (-1.02107292320897e+15 - 2.83675545731519e+15j)
    >>> coulombg(1+0.1j, 50+50j, 100.156)
    (2.83675545731519e+15 - 1.02107292320897e+15j)
    >>> coulombf(1e-5j, 10+1e-5j, 0.1+1e-6j)
    (4.30566371247811e-14 - 9.03347835361657e-19j)
    >>> coulombg(1e-5j, 10+1e-5j, 0.1+1e-6j)
    (778709182061.134 + 18418936.2660553j)

The following reproduces a table in Abramowitz & Stegun, at twice
the precision::

    >>> mp.dps = 10
    >>> eta = 2; z = 5
    >>> for l in [5, 4, 3, 2, 1, 0]:
    ...     print l, coulombf(l,eta,z), diff(lambda z: coulombf(l,eta,z), z)
    ...
    5 0.09079533488 0.1042553261
    4 0.2148205331 0.2029591779
    3 0.4313159311 0.320534053
    2 0.7212774133 0.3952408216
    1 0.9935056752 0.3708676452
    0 1.143337392 0.2937960375

**References**

1. I.J. Thompson & A.R. Barnett, "Coulomb and Bessel Functions of Complex
   Arguments and Order", J. Comp. Phys., vol 64, no. 2, June 1986.

2. N. Michel, "Precise Coulomb wave functions for a wide range of
   complex `l`, `\eta` and `z`", http://arxiv.org/abs/physics/0702051v1

"""

coulombg = r"""
Calculates the irregular Coulomb wave function

.. math ::

    G_l(\eta,z) = \frac{F_l(\eta,z) \cos(\chi) - F_{-l-1}(\eta,z)}{\sin(\chi)}

where `\chi = \sigma_l - \sigma_{-l-1} - (l+1/2) \pi`
and `\sigma_l(\eta) = (\ln \Gamma(1+l+i\eta)-\ln \Gamma(1+l-i\eta))/(2i)`.

See :func:`coulombf` for additional information.

**Examples**

Evaluation is supported for arbitrary magnitudes of `z`::

    >>> from mpmath import *
    >>> mp.dps = 25; mp.pretty = True
    >>> coulombg(-2, 1.5, 3.5)
    1.380011900612186346255524
    >>> coulombg(2, 1.5, 3.5)
    1.919153700722748795245926
    >>> coulombg(-2, 1.5, '1e-10')
    201126715824.7329115106793
    >>> coulombg(-2, 1.5, 1000)
    0.1802071520691149410425512
    >>> coulombg(-2, 1.5, 10**10)
    0.652103020061678070929794

The following reproduces a table in Abramowitz & Stegun,
at twice the precision::

    >>> mp.dps = 10
    >>> eta = 2; z = 5
    >>> for l in [1, 2, 3, 4, 5]:
    ...     print l, coulombg(l,eta,z), -diff(lambda z: coulombg(l,eta,z), z)
    ...
    1 1.08148276 0.6028279961
    2 1.496877075 0.5661803178
    3 2.048694714 0.7959909551
    4 3.09408669 1.731802374
    5 5.629840456 4.549343289

Evaluation close to the singularity at `z = 0`::

    >>> mp.dps = 15
    >>> coulombg(0,10,1)
    3088184933.67358
    >>> coulombg(0,10,'1e-10')
    5554866000719.8
    >>> coulombg(0,10,'1e-100')
    5554866221524.1

Evaluation with a half-integer value for `l`::

    >>> coulombg(1.5, 1, 10)
    0.852320038297334

"""

coulombc = r"""
Gives the normalizing Gamow constant for Coulomb wave functions,

.. math ::

    C_l(\eta) = 2^l \exp\left(-\pi \eta/2 + [\ln \Gamma(1+l+i\eta) +
        \ln \Gamma(1+l-i\eta)]/2 - \ln \Gamma(2l+2)\right),

where the log gamma function with continuous imaginary part
away from the negative half axis (see :func:`loggamma`) is implied.

This function is used internally for the calculation of
Coulomb wave functions, and automatically cached to make multiple
evaluations with fixed `l`, `\eta` fast.


"""

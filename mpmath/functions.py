import gammazeta
import libintmath
import rational

from libhyper import NoConvergence

class SpecialFunctions(object):
    """
    This class implements special functions using high-level code.

    Elementary and some other functions (e.g. gamma function, basecase
    hypergeometric series) are assumed to be predefined by the context as
    "builtins" or "low-level" functions.
    """
    defined_functions = {}

    # The series for the Jacobi theta functions converge for |q| < 1;
    # in the current implementation they throw a ValueError for
    # abs(q) > THETA_Q_LIM
    THETA_Q_LIM = 1 - 10**-7

    def __init__(self):
        cls = self.__class__
        for name in cls.defined_functions:
            f, wrap = cls.defined_functions[name]
            cls.wrap_specfun(name, f, wrap)

        self.mpq_1 = self._mpq(rational.mpq_1)
        self.mpq_0 = self._mpq(rational.mpq_0)
        self.mpq_1_2 = self._mpq(rational.mpq_1_2)
        self.mpq_3_2 = self._mpq(rational.mpq_3_2)
        self.mpq_1_4 = self._mpq(rational.mpq_1_4)
        self.mpq_1_16 = self._mpq(rational.mpq_1_16)
        self.mpq_3_16 = self._mpq(rational.mpq_3_16)
        self.mpq_5_2 = self._mpq(rational.mpq_5_2)
        self.mpq_3_4 = self._mpq(rational.mpq_3_4)
        self.mpq_7_4 = self._mpq(rational.mpq_7_4)
        self.mpq_5_4 = self._mpq(rational.mpq_5_4)

    # Default -- do nothing
    @classmethod
    def wrap_specfun(cls, name, f, wrap):
        setattr(cls, name, f)

    # Optional fast versions of common functions in common cases.
    # If not overridden, default (generic hypergeometric series)
    # implementations will be used
    def _besselj(ctx, n, z): raise NotImplementedError
    def _erf(ctx, z): raise NotImplementedError
    def _erfc(ctx, z): raise NotImplementedError
    def _gamma_upper_int(ctx, z, a): raise NotImplementedError
    def _expint_int(ctx, n, z): raise NotImplementedError


def defun_wrapped(f):
    SpecialFunctions.defined_functions[f.__name__] = f, True

def defun(f):
    SpecialFunctions.defined_functions[f.__name__] = f, False

def defun_static(f):
    setattr(SpecialFunctions, f.__name__, f)

defun_static(gammazeta.bernfrac)

@defun_wrapped
def cot(ctx, z): return ctx.one / ctx.tan(z)

@defun_wrapped
def sec(ctx, z): return ctx.one / ctx.cos(z)

@defun_wrapped
def csc(ctx, z): return ctx.one / ctx.sin(z)

@defun_wrapped
def coth(ctx, z): return ctx.one / ctx.tanh(z)

@defun_wrapped
def sech(ctx, z): return ctx.one / ctx.cosh(z)

@defun_wrapped
def csch(ctx, z): return ctx.one / ctx.sinh(z)

@defun_wrapped
def acot(ctx, z): return ctx.atan(ctx.one / z)

@defun_wrapped
def asec(ctx, z): return ctx.acos(ctx.one / z)

@defun_wrapped
def acsc(ctx, z): return ctx.asin(ctx.one / z)

@defun_wrapped
def acoth(ctx, z): return ctx.atanh(ctx.one / z)

@defun_wrapped
def asech(ctx, z): return ctx.acosh(ctx.one / z)

@defun_wrapped
def acsch(ctx, z): return ctx.asinh(ctx.one / z)


@defun
def sign(ctx, x):
    x = ctx.convert(x)
    if not x or ctx.isnan(x):
        return x
    if ctx.is_real_type(x):
        return ctx.mpf(cmp(x, 0))
    return x / abs(x)

@defun
def agm(ctx, a, b=1):
    if b == 1:
        return ctx.agm1(a)
    a = ctx.convert(a)
    b = ctx.convert(b)
    return ctx._agm(a, b)

@defun_wrapped
def sinc(ctx, x):
    if ctx.isinf(x):
        return 1/x
    if not x:
        return x+1
    return ctx.sin(x)/x

@defun_wrapped
def sincpi(ctx, x):
    if ctx.isinf(x):
        return 1/x
    if not x:
        return x+1
    return ctx.sinpi(x)/(ctx.pi*x)

# TODO: tests; improve implementation
@defun_wrapped
def expm1(ctx, x):
    if not x:
        return ctx.zero
    # exp(x) - 1 ~ x
    if ctx.mag(x) < -ctx.prec:
        return x + 0.5*x**2
    # TODO: accurately eval the smaller of the real/imag parts
    return sum_accurately(ctx, lambda: iter([ctx.exp(x),-1]),1)

@defun_wrapped
def powm1(ctx, x, y):
    mag = ctx.mag
    one = ctx.one
    w = x**y - one
    M = mag(w)
    # Only moderate cancellation
    if M > -8:
        return w
    # Check for the only possible exact cases
    if not w:
        if (not y) or (x in (1, -1, 1j, -1j) and ctx.isint(y)):
            return w
    x1 = x - one
    magy = mag(y)
    lnx = ctx.ln(x)
    # Small y: x^y - 1 ~ log(x)*y + O(log(x)^2 * y^2)
    if magy + mag(lnx) < -ctx.prec:
        return lnx*y + (lnx*y)**2/2
    # TODO: accurately eval the smaller of the real/imag part
    return sum_accurately(ctx, lambda: iter([x**y, -1]), 1)

@defun
def _rootof1(ctx, k, n):
    k = int(k)
    n = int(n)
    k %= n
    if not k:
        return ctx.one
    elif 2*k == n:
        return -ctx.one
    elif 4*k == n:
        return ctx.j
    elif 4*k == 3*n:
        return -ctx.j
    return ctx.exp(2*ctx.pi*k/n*ctx.j)

@defun
def root(ctx, x, n, k=0):
    n = int(n)
    x = ctx.convert(x)
    if k:
        # Special case: there is an exact real root
        if (n & 1 and 2*k == n-1) and (not ctx.im(x)) and (ctx.re(x) < 0):
            return -ctx.root(-x, n)
        # Multiply by root of unity
        prec = ctx.prec
        try:
            ctx.prec += 10
            v = ctx.root(x, n, 0) * ctx._rootof1(k, n)
        finally:
            ctx.prec = prec
        return +v
    return ctx._nthroot(x, n)

@defun
def unitroots(ctx, n, primitive=False):
    gcd = libintmath.gcd
    prec = ctx.prec
    try:
        ctx.prec += 10
        if primitive:
            v = [ctx._rootof1(k,n) for k in range(n) if gcd(k,n) == 1]
        else:
            # TODO: this can be done *much* faster
            v = [ctx._rootof1(k,n) for k in range(n)]
    finally:
        ctx.prec = prec
    return [+x for x in v]

@defun
def arg(ctx, x):
    x = ctx.convert(x)
    return ctx.atan2(x.imag, x.real)

@defun
def fabs(ctx, x):
    return abs(ctx.convert(x))

@defun
def re(ctx, x):
    return ctx.convert(x).real

@defun
def im(ctx, x):
    return ctx.convert(x).imag

@defun
def conj(ctx, x):
    return ctx.convert(x).conjugate()

@defun
def log(ctx, x, b=None):
    if b is None:
        return ctx.ln(x)
    wp = ctx.prec + 20
    return ctx.ln(x, prec=wp) / ctx.ln(b, prec=wp)

@defun
def log10(ctx, x):
    return ctx.log(x, 10)

@defun
def modf(ctx, x, y):
    return ctx.convert(x) % ctx.convert(y)

@defun
def degrees(ctx, x):
    return x / ctx.degree

@defun
def radians(ctx, x):
    return x * ctx.degree

@defun
def psi0(ctx, z):
    """Shortcut for psi(0,z) (the digamma function)"""
    return ctx.psi(0, z)

@defun
def psi1(ctx, z):
    """Shortcut for psi(1,z) (the trigamma function)"""
    return ctx.psi(1, z)

@defun
def psi2(ctx, z):
    """Shortcut for psi(2,z) (the tetragamma function)"""
    return ctx.psi(2, z)

@defun
def psi3(ctx, z):
    """Shortcut for psi(3,z) (the pentagamma function)"""
    return ctx.psi(3, z)

#polygamma = MultiPrecisionArithmetic.polygamma = psi
#digamma = MultiPrecisionArithmetic.digamma = psi0
#trigamma = MultiPrecisionArithmetic.trigamma = psi1
#tetragamma = MultiPrecisionArithmetic.tetragamma = psi2
#pentagamma = MultiPrecisionArithmetic.pentagamma = psi3

@defun
def stieltjes(ctx, n, a=1):
    n = ctx.convert(n)
    a = ctx.convert(a)
    if n < 0:
        return ctx.bad_domain("Stieltjes constants defined for n >= 0")
    if hasattr(ctx, "stieltjes_cache"):
        stieltjes_cache = ctx.stieltjes_cache
    else:
        stieltjes_cache = ctx.stieltjes_cache = {}
    if a == 1:
        if n == 0:
            return +ctx.euler
        if n in stieltjes_cache:
            prec, s = stieltjes_cache[n]
            if prec >= ctx.prec:
                return +s
    mag = 1
    def f(x):
        xa = x/a
        v = (xa-ctx.j)*ctx.ln(a-ctx.j*x)**n/(1+xa**2)/(ctx.exp(2*ctx.pi*x)-1)
        return v.real / mag
    orig = ctx.prec
    try:
        # Normalize integrand by approx. magnitude to
        # speed up quadrature (which uses absolute error)
        if n > 50:
            ctx.prec = 20
            mag = ctx.quad(f, [0,ctx.inf], maxdegree=3)
        ctx.prec = orig + 10 + int(n**0.5)
        s = ctx.quad(f, [0,ctx.inf], maxdegree=20)
        v = ctx.ln(a)**n/(2*a) - ctx.ln(a)**(n+1)/(n+1) + 2*s/a*mag
    finally:
        ctx.prec = orig
    if a == 1 and ctx.isint(n):
        stieltjes_cache[n] = (ctx.prec, v)
    return +v

@defun
def gammaprod(ctx, a, b, _infsign=False):
    a = [ctx.convert(x) for x in a]
    b = [ctx.convert(x) for x in b]
    poles_num = []
    poles_den = []
    regular_num = []
    regular_den = []
    for x in a: [regular_num, poles_num][ctx.isnpint(x)].append(x)
    for x in b: [regular_den, poles_den][ctx.isnpint(x)].append(x)
    # One more pole in numerator or denominator gives 0 or inf
    if len(poles_num) < len(poles_den): return ctx.zero
    if len(poles_num) > len(poles_den):
        # Get correct sign of infinity for x+h, h -> 0 from above
        # XXX: hack, this should be done properly
        if _infsign:
            a = [x and x*(1+ctx.eps) or x+ctx.eps for x in poles_num]
            b = [x and x*(1+ctx.eps) or x+ctx.eps for x in poles_den]
            return ctx.sign(ctx.gammaprod(a+regular_num,b+regular_den)) * ctx.inf
        else:
            return ctx.inf
    # All poles cancel
    # lim G(i)/G(j) = (-1)**(i+j) * gamma(1-j) / gamma(1-i)
    p = ctx.one
    orig = ctx.prec
    try:
        ctx.prec = orig + 15
        while poles_num:
            i = poles_num.pop()
            j = poles_den.pop()
            p *= (-1)**(i+j) * ctx.gamma(1-j) / ctx.gamma(1-i)
        for x in regular_num: p *= ctx.gamma(x)
        for x in regular_den: p /= ctx.gamma(x)
    finally:
        ctx.prec = orig
    return +p

@defun
def beta(ctx, x, y):
    x = ctx.convert(x)
    y = ctx.convert(y)
    if ctx.isinf(y):
        x, y = y, x
    if ctx.isinf(x):
        if x == ctx.inf and not y.imag:
            if y == ctx.ninf:
                return ctx.nan
            if y > 0:
                return ctx.zero
            if ctx.isint(y):
                return ctx.nan
            if y < 0:
                return ctx.sign(ctx.gamma(y)) * ctx.inf
        return ctx.nan
    return ctx.gammaprod([x, y], [x+y])

@defun_wrapped
def betainc(ctx, a, b, x1=0, x2=1, regularized=False):
    if x1 == x2:
        v = 0
    elif not x1:
        if x1 == 0 and x2 == 1:
            v = ctx.beta(a, b)
        else:
            v = x2**a * ctx.hyp2f1(a, 1-b, a+1, x2) / a
    else:
        m, d = ctx.nint_distance(a)
        if m <= 0:
            if d < -ctx.prec:
                h = +ctx.eps
                ctx.prec *= 2
                a += h
            elif d < -4:
                ctx.prec -= d
        s1 = x2**a * ctx.hyp2f1(a,1-b,a+1,x2)
        s2 = x1**a * ctx.hyp2f1(a,1-b,a+1,x1)
        v = (s1 - s2) / a
    if regularized:
        v /= ctx.beta(a,b)
    return v

@defun
def binomial(ctx, n, k):
    return ctx.gammaprod([n+1], [k+1, n-k+1])

@defun
def rf(ctx, x, n):
    return ctx.gammaprod([x+n], [x])

@defun
def ff(ctx, x, n):
    return ctx.gammaprod([x+1], [x-n+1])

@defun_wrapped
def fac2(ctx, x):
    if ctx.isinf(x):
        if x == ctx.inf:
            return x
        return ctx.nan
    return 2**(x/2)*(ctx.pi/2)**((ctx.cospi(x)-1)/4)*ctx.gamma(x/2+1)

@defun
def hyper(ctx, a_s, b_s, z, **kwargs):
    """
    Hypergeometric function, general case.
    """
    p = len(a_s)
    q = len(b_s)
    z = ctx.convert(z)
    degree = p, q
    # Handle common cases
    if p == 0:
        if q == 1:
            return ctx.hyp0f1(b_s[0], z, **kwargs)
        elif q == 0:
            return ctx.exp(z)
    elif p == 1:
        if q == 1:
            return ctx.hyp1f1(a_s[0], b_s[0], z, **kwargs)
        elif q == 2:
            return ctx.hyp1f2(a_s[0], b_s[0], b_s[1], z, **kwargs)
        elif q == 0:
            return (1-z)**(-a_s[0])
    elif p == 2:
        if q == 1:
            return ctx.hyp2f1(a_s[0], a_s[1], b_s[0], z, **kwargs)
        elif q == 2:
            return ctx.hyp2f2(a_s[0], a_s[1], b_s[0], b_s[1], z, **kwargs)
        elif q == 3:
            return ctx.hyp2f3(a_s[0], a_s[1], b_s[0], b_s[1], b_s[2], z, **kwargs)
        elif q == 0:
            return ctx.hyp2f0(a_s[0], a_s[1], z, **kwargs)
    elif p == q+1:
        return ctx._hypq1fq(p, q, a_s, b_s, z, **kwargs)
    coeffs = map(ctx.convert_maybe_rational, a_s+b_s)
    coeffs, types = zip(*coeffs)
    return ctx.hypsum(p, q, types, coeffs, z, **kwargs)

@defun
def hyp0f1(ctx, b, z, **kwargs):
    """
    Hypergeometric 0F1.
    """
    z = ctx.convert(z)
    b, btype = ctx.convert_maybe_rational(b)
    if z:
        magz = ctx.mag(z)
    else:
        magz = 0
    if magz >= 8 and not kwargs.get('force_series'):
        try:
            # http://functions.wolfram.com/HypergeometricFunctions/
            # Hypergeometric0F1/06/02/03/0004/
            # We don't need hypercomb because the only possible singularity
            # occurs when the value is undefined. However, we should perhaps
            # still check for cancellation...
            # TODO: handle the all-real case more efficiently!
            # TODO: figure out how much precision is needed (exponential growth)
            orig = ctx.prec
            try:
                ctx.prec += 12 + magz//2
                w = ctx.sqrt(-z)
                jw = ctx.j*w
                u = 1/(4*jw)
                c = ctx.mpq_1_2 - b
                E = ctx.exp(2*jw)
                H1 = (-jw)**c/E*ctx.hyp2f0(b-ctx.mpq_1_2, ctx.mpq_3_2-b, -u,
                    force_series=True)
                H2 = (jw)**c*E*ctx.hyp2f0(b-ctx.mpq_1_2, ctx.mpq_3_2-b, u,
                    force_series=True)
                v = ctx.gamma(b)/(2*ctx.sqrt(ctx.pi))*(H1 + H2)
            finally:
                ctx.prec = orig
            if ctx.is_real_type(b) and ctx.is_real_type(z):
                v = v.real
            return +v
        except NoConvergence:
            pass
    return ctx.hypsum(0, 1, (btype,), [b], z, **kwargs)

@defun
def hyp1f1(ctx, a, b, z, **kwargs):
    """
    Hypergeometric 1F1.
    """
    z = ctx.convert(z)
    a, atype = ctx.convert_maybe_rational(a)
    b, btype = ctx.convert_maybe_rational(b)
    if not z:
        return ctx.one+z
    magz = ctx.mag(z)
    if magz >= 7 and not (ctx.isint(a) and ctx.re(a) <= 0):
        if ctx.isinf(z):
            if ctx.sign(a) == ctx.sign(b) == ctx.sign(z) == 1:
                return ctx.inf
            return ctx.nan * z
        try:
            try:
                ctx.prec += magz
                sector = z.imag < 0 and z.real <= 0
                def h(a,b):
                    if sector:
                        E = ctx.exp(-ctx.j * ctx.pi*a)
                    else:
                        E = ctx.exp(ctx.j * ctx.pi*a)
                    rz = 1/z
                    T1 = ([E,z], [1,-a], [b], [b-a], [a, 1+a-b], [], -rz)
                    T2 = ([ctx.exp(z),z], [1,a-b], [b], [a], [b-a, 1-a], [], rz)
                    return T1, T2
                v = ctx.hypercomb(h, [a,b], force_series=True)
                if ctx.is_real_type(a) and ctx.is_real_type(b) and ctx.is_real_type(z):
                    v = v.real
                return +v
            except NoConvergence:
                pass
        finally:
            ctx.prec -= magz
    #print "TEST", ctx.prec, a._v, b._v, z._v
    v = ctx.hypsum(1, 1, (atype, btype), [a, b], z, **kwargs)
    #print "1F1", complex(a), complex(b), complex(z), complex(v)
    return v

def _hyp2f1_gosper(ctx,a,b,c,z,**kwargs):
    # Use Gosper's recurrence
    # See http://www.math.utexas.edu/pipermail/maxima/2006/000126.html
    _a,_b,_c,_z = a, b, c, z
    orig = ctx.prec
    maxprec = kwargs.get('maxprec', 100*orig)
    extra = 10
    while 1:
        ctx.prec = orig + extra
        #a = ctx.convert(_a)
        #b = ctx.convert(_b)
        #c = ctx.convert(_c)
        z = ctx.convert(_z)
        d = ctx.mpf(0)
        e = ctx.mpf(1)
        f = ctx.mpf(0)
        k = 0
        # Common subexpression elimination, unfortunately making
        # things a bit unreadable. The formula is quite messy to begin
        # with, though...
        abz = a*b*z
        ch = c * ctx.mpq_1_2
        c1h = (c+1) * ctx.mpq_1_2
        nz = 1-z
        g = z/nz
        abg = a*b*g
        cba = c-b-a
        z2 = z-2
        tol = -ctx.prec - 10
        nstr = ctx.nstr
        nprint = ctx.nprint
        mag = ctx.mag
        maxmag = ctx.ninf
        while 1:
            kch = k+ch
            kakbz = (k+a)*(k+b)*z / (4*(k+1)*kch*(k+c1h))
            d1 = kakbz*(e-(k+cba)*d*g)
            e1 = kakbz*(d*abg+(k+c)*e)
            ft = d*(k*(cba*z+k*z2-c)-abz)/(2*kch*nz)
            f1 = f + e - ft
            maxmag = max(maxmag, mag(f1))
            if mag(f1-f) < tol:
                break
            d, e, f = d1, e1, f1
            k += 1
        cancellation = maxmag - mag(f1)
        if cancellation < extra:
            break
        else:
            extra += cancellation
            if extra > maxprec:
                raise NoConvergence
    return f1

@defun
def hyp2f1(ctx,a,b,c,z,**kwargs):
    """
    Hypergeometric 2F1.
    """
    z = ctx.convert(z)
    a, atype = ctx.convert_maybe_rational(a)
    b, btype = ctx.convert_maybe_rational(b)
    c, ctype = ctx.convert_maybe_rational(c)

    if z == 1:
        # TODO: the following logic can be simplified
        convergent = ctx.re(c-a-b) > 0
        finite = (ctx.isint(a) and a <= 0) or (ctx.isint(b) and b <= 0)
        zerodiv = ctx.isint(c) and c <= 0 and not \
            ((ctx.isint(a) and c <= a <= 0) or (ctx.isint(b) and c <= b <= 0))
        #print "bz", a, b, c, z, convergent, finite, zerodiv
        # Gauss's theorem gives the value if convergent
        if (convergent or finite) and not zerodiv:
            return ctx.gammaprod([c, c-a-b], [c-a, c-b], _infsign=True)
        # Otherwise, there is a pole and we take the
        # sign to be that when approaching from below
        # XXX: this evaluation is not necessarily correct in all cases
        return ctx.hyp2f1(a,b,c,1-ctx.eps*2) * ctx.inf

    # Equal to 1 (first term), unless there is a subsequent
    # division by zero
    if not z:
        # Division by zero but power of z is higher than
        # first order so cancels
        if c or a == 0 or b == 0:
            return 1+z
        # Indeterminate
        return ctx.nan

    # Hit zero denominator unless numerator goes to 0 first
    if ctx.isint(c) and c <= 0:
        if (ctx.isint(a) and c <= a <= 0) or \
           (ctx.isint(b) and c <= b <= 0):
            pass
        else:
            # Pole in series
            return ctx.inf
    elif a == c:
        return (1-z)**(-b)
    elif b == c:
        return (1-z)**(-a)

    absz = abs(z)

    # Fast case: standard series converges rapidly,
    # possibly in finitely many terms
    if absz <= 0.8 or (ctx.isint(a) and a <= 0 and a >= -1000) or \
                      (ctx.isint(b) and b <= 0 and b >= -1000):
        return ctx.hypsum(2, 1, (atype, btype, ctype), [a, b, c], z, **kwargs)

    orig = ctx.prec
    try:
        ctx.prec += 10

        # Use 1/z transformation
        if absz >= 1.3:
            def h(a,b):
                t = ctx.mpq_1-c; ab = a-b; rz = 1/z
                T1 = ([-z],[-a], [c,-ab],[b,c-a], [a,t+a],[ctx.mpq_1+ab],  rz)
                T2 = ([-z],[-b], [c,ab],[a,c-b], [b,t+b],[ctx.mpq_1-ab],  rz)
                return T1, T2
            v = ctx.hypercomb(h, [a,b], **kwargs)

        # Use 1-z transformation
        elif abs(1-z) <= 0.75:
            def h(a,b):
                t = c-a-b; ca = c-a; cb = c-b; rz = 1-z
                T1 = [], [], [c,t], [ca,cb], [a,b], [1-t], rz
                T2 = [rz], [t], [c,a+b-c], [a,b], [ca,cb], [1+t], rz
                return T1, T2
            v = ctx.hypercomb(h, [a,b], **kwargs)

        # Remaining part of unit circle
        else:
            v = _hyp2f1_gosper(ctx,a,b,c,z,**kwargs)

        # TODO: for z ~= -1, we could use a series expansion around -1
        # using the fact that F(a,b,c,-1) = F(a,c-b,c,0.5)/2^a

    finally:
        ctx.prec = orig
    return +v

@defun
def _hypq1fq(ctx, p, q, a_s, b_s, z, **kwargs):
    r"""
    Evaluates 3F2, 4F3, 5F4, ...
    """
    z = ctx.convert(z)
    a_s = map(ctx.convert_maybe_rational, a_s)
    b_s = map(ctx.convert_maybe_rational, b_s)
    a_s, a_types = zip(*a_s)
    b_s, b_types = zip(*b_s)
    a_s = list(a_s)
    b_s = list(b_s)
    absz = abs(z)
    ispoly = False
    for a in a_s:
        if ctx.isint(a) and a <= 0:
            ispoly = True
            break
    # Direct summation
    if absz < 1 or ispoly:
        try:
            return ctx.hypsum(p, q, a_types+b_types, a_s+b_s, z, **kwargs)
        except NoConvergence:
            if absz > 1.1 or ispoly:
                raise
    # Use expansion at |z-1| -> 0.
    # Reference: Wolfgang Buhring, "Generalized Hypergeometric Functions at
    #   Unit Argument", Proc. Amer. Math. Soc., Vol. 114, No. 1 (Jan. 1992),
    #   pp.145-153
    # The current implementation has several problems:
    # 1. We only implement it for 3F2. The expansion coefficients are
    #    given by extremely messy nested sums in the higher degree cases
    #    (see reference). Is efficient sequential generation of the coefficients
    #    possible in the > 3F2 case?
    # 2. Although the series converges, it may do so slowly, so we need
    #    convergence acceleration. The acceleration implemented by
    #    nsum does not always help, so results returned are sometimes
    #    inaccurate! Can we do better?
    # 3. We should check conditions for convergence, and possibly
    #    do a better job of cancelling out gamma poles if possible.
    if z == 1:
        # XXX: should also check for division by zero in the
        # denominator of the series (cf. hyp2f1)
        S = ctx.re(sum(b_s)-sum(a_s))
        if S <= 0:
            #return ctx.hyper(a_s, b_s, 1-ctx.eps*2, **kwargs) * ctx.inf
            return ctx.hyper(a_s, b_s, 0.9, **kwargs) * ctx.inf
    if (p,q) == (3,2) and abs(z-1) < 0.05:   # and kwargs.get('sum1')
        #print "Using alternate summation (experimental)"
        a1,a2,a3 = a_s
        b1,b2 = b_s
        u = b1+b2-a3
        initial = ctx.gammaprod([b2-a3,b1-a3,a1,a2],[b2-a3,b1-a3,1,u])
        def term(k, _cache={0:initial}):
            u = b1+b2-a3+k
            if k in _cache:
                t = _cache[k]
            else:
                t = _cache[k-1]
                t *= (b1+k-a3-1)*(b2+k-a3-1)
                t /= k*(u-1)
                _cache[k] = t
            return t * ctx.hyp2f1(a1,a2,u,z)
        try:
            S = ctx.nsum(term, [0,ctx.inf], verbose=kwargs.get('verbose'),
                strict=True)
            return S * ctx.gammaprod([b1,b2],[a1,a2,a3])
        except ctx.NoConvergence:
            pass
    # Try to use convergence acceleration on and close to the unit circle.
    # Problem: the convergence acceleration degenerates as |z-1| -> 0,
    # except for special cases. Everywhere else, the Shanks transformation
    # is very efficient.
    if absz < 1.1 and z.real <= 1:
        def term(k, _cache={0:ctx.one}):
            k = int(k)
            if k in _cache:
                return _cache[k]
            t = _cache[k-1]
            m = k-1
            for j in xrange(p): t *= (a_s[j]+m)
            for j in xrange(q): t /= (b_s[j]+m)
            t *= z
            t /= k
            _cache[k] = t
            return t
        return ctx.nsum(term, [0,ctx.inf], verbose=kwargs.get('verbose'),
            strict=True)
    # Use 1/z transformation
    # http://functions.wolfram.com/HypergeometricFunctions/
    #   HypergeometricPFQ/06/01/05/02/0004/
    def h(*args):
        a_s = list(args[:p])
        b_s = list(args[p:])
        Ts = []
        recz = ctx.one/z
        negz = ctx.fneg(z, exact=True)
        for k in range(q+1):
            ak = a_s[k]
            C = [negz]
            Cp = [-ak]
            Gn = b_s + [ak] + [a_s[j]-ak for j in range(q+1) if j != k]
            Gd = a_s + [b_s[j]-ak for j in range(q)]
            Fn = [ak] + [ak-b_s[j]+1 for j in range(q)]
            Fd = [1-a_s[j]+ak for j in range(q+1) if j != k]
            Ts.append((C, Cp, Gn, Gd, Fn, Fd, recz))
        return Ts
    return ctx.hypercomb(h, a_s+b_s, **kwargs)

def _check_need_perturb(ctx, terms, prec, discard_known_zeros):
    perturb = recompute = False
    extraprec = 0
    discard = []
    for term_index, term in enumerate(terms):
        w_s, c_s, alpha_s, beta_s, a_s, b_s, z = term
        have_singular_nongamma_weight = False
        # Avoid division by zero in leading factors (TODO:
        # also check for near division by zero?)
        for k, w in enumerate(w_s):
            if not w:
                if ctx.re(c_s[k]) <= 0:
                    perturb = recompute = True
                    have_singular_nongamma_weight = True
        pole_count = [0, 0, 0]
        # Check for gamma and series poles and near-poles
        for data_index, data in enumerate([alpha_s, beta_s, b_s]):
            for i, x in enumerate(data):
                n, d = ctx.nint_distance(x)
                # Poles
                if n > 0:
                    continue
                if d == ctx.ninf:
                    # OK if we have a polynomial
                    # ------------------------------
                    ok = False
                    if data_index == 2:
                        for u in a_s:
                            if ctx.isnpint(u) and u >= int(n):
                                ok = True
                                break
                    if ok:
                        continue
                    pole_count[data_index] += 1
                    # ------------------------------
                    #perturb = recompute = True
                    #return perturb, recompute, extraprec
                elif d < -4:
                    extraprec += -d
                    recompute = True
        if discard_known_zeros and pole_count[1] > pole_count[0] + pole_count[2] \
            and not have_singular_nongamma_weight:
            discard.append(term_index)
        elif sum(pole_count):
            perturb = recompute = True
    return perturb, recompute, extraprec, discard

_hypercomb_msg = """
hypercomb() failed to converge to the requested %i bits of accuracy
using a working precision of %i bits. The function value may be zero or
infinite; try passing zeroprec=N or infprec=M to bound finite values between
2^(-N) and 2^M. Otherwise try a higher maxprec or maxterms.
"""

@defun
def hypercomb(ctx, function, params=[], discard_known_zeros=True, **kwargs):
    orig = ctx.prec
    sumvalue = ctx.zero
    dist = ctx.nint_distance
    ninf = ctx.ninf
    orig_params = params[:]
    verbose = kwargs.get('verbose', False)
    maxprec = kwargs.get('maxprec', ctx._default_hyper_maxprec(orig))
    kwargs['maxprec'] = maxprec   # For calls to hypsum
    zeroprec = kwargs.get('zeroprec')
    infprec = kwargs.get('infprec')
    perturbed_reference_value = None
    hextra = 0
    try:
        while 1:
            ctx.prec += 10
            if ctx.prec > maxprec:
                raise ValueError(_hypercomb_msg % (orig, ctx.prec))
            orig2 = ctx.prec
            params = orig_params[:]
            terms = function(*params)
            if verbose:
                print
                print "ENTERING hypercomb main loop"
                print "prec =", ctx.prec
                print "hextra", hextra
            perturb, recompute, extraprec, discard = \
                _check_need_perturb(ctx, terms, orig, discard_known_zeros)
            ctx.prec += extraprec
            if perturb:
                if "hmag" in kwargs:
                    hmag = kwargs["hmag"]
                elif ctx._fixed_precision:
                    hmag = int(ctx.prec*0.3)
                else:
                    hmag = orig + 10 + hextra
                h = ctx.ldexp(ctx.one, -hmag)
                ctx.prec = orig2 + 10 + hmag + 10
                for k in range(len(params)):
                    params[k] += h
                    # Heuristically ensure that the perturbations
                    # are "independent" so that two perturbations
                    # don't accidentally cancel each other out
                    # in a subtraction.
                    h += h/(k+1)
            if recompute:
                terms = function(*params)
            if discard_known_zeros:
                terms = [term for (i, term) in enumerate(terms) if i not in discard]
            if not terms:
                return ctx.zero
            evaluated_terms = []
            for term_index, term_data in enumerate(terms):
                w_s, c_s, alpha_s, beta_s, a_s, b_s, z = term_data
                if verbose:
                    print
                    print "  Evaluating term %i/%i : %iF%i" % \
                        (term_index+1, len(terms), len(a_s), len(b_s))
                    print "    powers", ctx.nstr(w_s), ctx.nstr(c_s)
                    print "    gamma", ctx.nstr(alpha_s), ctx.nstr(beta_s)
                    print "    hyper", ctx.nstr(a_s), ctx.nstr(b_s)
                    print "    z", ctx.nstr(z)
                v = ctx.hyper(a_s, b_s, z, **kwargs)
                for a in alpha_s: v *= ctx.gamma(a)
                for b in beta_s: v /= ctx.gamma(b)
                for w, c in zip(w_s, c_s): v *= ctx.power(w, c)
                if verbose:
                    print "    Value:", v
                evaluated_terms.append(v)

            if len(terms) == 1 and (not perturb):
                sumvalue = evaluated_terms[0]
                break

            if ctx._fixed_precision:
                sumvalue = ctx.fsum(evaluated_terms)
                break

            sumvalue = ctx.fsum(evaluated_terms)
            term_magnitudes = [ctx.mag(x) for x in evaluated_terms]
            max_magnitude = max(term_magnitudes)
            sum_magnitude = ctx.mag(sumvalue)
            cancellation = max_magnitude - sum_magnitude
            if verbose:
                print
                print "  Cancellation:", cancellation, "bits"
                print "  Increased precision:", ctx.prec - orig, "bits"

            precision_ok = cancellation < ctx.prec - orig

            if zeroprec is None:
                zero_ok = False
            else:
                zero_ok = max_magnitude - ctx.prec < -zeroprec
            if infprec is None:
                inf_ok = False
            else:
                inf_ok = max_magnitude > infprec

            if precision_ok and (not perturb) or ctx.isnan(cancellation):
                break
            elif precision_ok:
                if perturbed_reference_value is None:
                    hextra += 20
                    perturbed_reference_value = sumvalue
                    continue
                elif ctx.mag(sumvalue - perturbed_reference_value) <= \
                        ctx.mag(sumvalue) - orig:
                    break
                elif zero_ok:
                    sumvalue = ctx.zero
                    break
                elif inf_ok:
                    sumvalue = ctx.inf
                    break
                elif 'hmag' in kwargs:
                    break
                else:
                    hextra *= 2
                    perturbed_reference_value = sumvalue
            # Increase precision
            else:
                increment = min(max(cancellation, orig//2), max(extraprec,orig))
                ctx.prec += increment
                if verbose:
                    print "  Must start over with increased precision"
                continue
    finally:
        ctx.prec = orig
    return +sumvalue



@defun
def hyp2f2(ctx,a1,a2,b1,b2,z,**kwargs):
    a1, a1type = ctx.convert_maybe_rational(a1)
    a2, a2type = ctx.convert_maybe_rational(a2)
    b1, b1type = ctx.convert_maybe_rational(b1)
    b2, b2type = ctx.convert_maybe_rational(b2)
    z = ctx.convert(z)

    absz = abs(z)
    magz = ctx.mag(z)
    orig = ctx.prec

    # Asymptotic expansion is ~ exp(z)
    asymp_extraprec = magz

    # Asymptotic series is in terms of 3F1
    can_use_asymptotic = (not kwargs.get('force_series')) and \
        (ctx.mag(absz) > 3)

    # TODO: much of the following could be shared with 2F3 instead of
    # copypasted
    if can_use_asymptotic:
        #print "using asymp"
        try:
            try:
                ctx.prec += asymp_extraprec
                # http://functions.wolfram.com/HypergeometricFunctions/
                # Hypergeometric2F2/06/02/02/0002/
                def h(a1,a2,b1,b2):
                    X = a1+a2-b1-b2
                    A2 = a1+a2
                    B2 = b1+b2
                    c = {}
                    c[0] = ctx.one
                    c[1] = (A2-1)*X+b1*b2-a1*a2
                    s1 = 0
                    k = 0
                    tprev = 0
                    while 1:
                        if k not in c:
                            uu1 = 1-B2+2*a1+a1**2+2*a2+a2**2-A2*B2+a1*a2+b1*b2+(2*B2-3*(A2+1))*k+2*k**2
                            uu2 = (k-A2+b1-1)*(k-A2+b2-1)*(k-X-2)
                            c[k] = ctx.one/k * (uu1*c[k-1]-uu2*c[k-2])
                        t1 = c[k] * z**(-k)
                        if abs(t1) < 0.1*ctx.eps:
                            #print "Convergence :)"
                            break
                        # Quit if the series doesn't converge quickly enough
                        if k > 5 and abs(tprev) / abs(t1) < 1.5:
                            #print "No convergence :("
                            raise NoConvergence
                        s1 += t1
                        tprev = t1
                        k += 1
                    S = ctx.exp(z)*s1
                    T1 = [z,S], [X,1], [b1,b2],[a1,a2],[],[],0
                    T2 = [-z],[-a1],[b1,b2,a2-a1],[a2,b1-a1,b2-a1],[a1,a1-b1+1,a1-b2+1],[a1-a2+1],-1/z
                    T3 = [-z],[-a2],[b1,b2,a1-a2],[a1,b1-a2,b2-a2],[a2,a2-b1+1,a2-b2+1],[-a1+a2+1],-1/z
                    return T1, T2, T3
                v = ctx.hypercomb(h, [a1,a2,b1,b2])
                if sum(ctx.is_real_type(u) for u in [a1,a2,b1,b2,z]) == 5:
                    v = ctx.re(v)
                return v
            except NoConvergence:
                pass
        finally:
            ctx.prec = orig

    return ctx.hypsum(2, 2, (a1type, a2type, b1type, b2type), [a1, a2, b1, b2], z, **kwargs)



@defun
def hyp1f2(ctx,a1,b1,b2,z,**kwargs):
    a1, a1type = ctx.convert_maybe_rational(a1)
    b1, b1type = ctx.convert_maybe_rational(b1)
    b2, b2type = ctx.convert_maybe_rational(b2)
    z = ctx.convert(z)

    absz = abs(z)
    magz = ctx.mag(z)
    orig = ctx.prec

    # Asymptotic expansion is ~ exp(sqrt(z))
    asymp_extraprec = z and magz//2

    # Asymptotic series is in terms of 3F0
    can_use_asymptotic = (not kwargs.get('force_series')) and \
        (ctx.mag(absz) > 19) and \
        (ctx.sqrt(absz) > 1.5*orig) #and \
        #ctx._hyp_check_convergence([a1, a1-b1+1, a1-b2+1], [],
        #    1/absz, orig+40+asymp_extraprec)

    # TODO: much of the following could be shared with 2F3 instead of
    # copypasted
    if can_use_asymptotic:
        #print "using asymp"
        try:
            try:
                ctx.prec += asymp_extraprec
                # http://functions.wolfram.com/HypergeometricFunctions/
                # Hypergeometric1F2/06/02/03/
                def h(a1,b1,b2):
                    X = ctx.mpq_1_2*(a1-b1-b2+ctx.mpq_1_2)
                    c = {}
                    c[0] = ctx.one
                    c[1] = 2*(ctx.mpq_1_4*(3*a1+b1+b2-2)*(a1-b1-b2)+b1*b2-ctx.mpq_3_16)
                    c[2] = 2*(b1*b2+ctx.mpq_1_4*(a1-b1-b2)*(3*a1+b1+b2-2)-ctx.mpq_3_16)**2+\
                        ctx.mpq_1_16*(-16*(2*a1-3)*b1*b2 + \
                        4*(a1-b1-b2)*(-8*a1**2+11*a1+b1+b2-2)-3)
                    s1 = 0
                    s2 = 0
                    k = 0
                    tprev = 0
                    while 1:
                        if k not in c:
                            uu1 = (3*k**2+(-6*a1+2*b1+2*b2-4)*k + 3*a1**2 - \
                                (b1-b2)**2 - 2*a1*(b1+b2-2) + ctx.mpq_1_4)
                            uu2 = (k-a1+b1-b2-ctx.mpq_1_2)*(k-a1-b1+b2-ctx.mpq_1_2)*\
                                (k-a1+b1+b2-ctx.mpq_5_2)
                            c[k] = ctx.one/(2*k)*(uu1*c[k-1]-uu2*c[k-2])
                        w = c[k] * (-z)**(-0.5*k)
                        t1 = (-ctx.j)**k * ctx.mpf(2)**(-k) * w
                        t2 = ctx.j**k * ctx.mpf(2)**(-k) * w
                        if abs(t1) < 0.1*ctx.eps:
                            #print "Convergence :)"
                            break
                        # Quit if the series doesn't converge quickly enough
                        if k > 5 and abs(tprev) / abs(t1) < 1.5:
                            #print "No convergence :("
                            raise NoConvergence
                        s1 += t1
                        s2 += t2
                        tprev = t1
                        k += 1
                    S = ctx.exp(ctx.j*(ctx.pi*X+2*ctx.sqrt(-z)))*s1 + \
                        ctx.exp(-ctx.j*(ctx.pi*X+2*ctx.sqrt(-z)))*s2
                    T1 = [0.5*S, ctx.pi, -z], [1, -0.5, X], [b1, b2], [a1],\
                        [], [], 0
                    T2 = [-z], [-a1], [b1,b2],[b1-a1,b2-a1], \
                        [a1,a1-b1+1,a1-b2+1], [], 1/z
                    return T1, T2
                v = ctx.hypercomb(h, [a1,b1,b2])
                if sum(ctx.is_real_type(u) for u in [a1,b1,b2,z]) == 4:
                    v = ctx.re(v)
                return v
            except NoConvergence:
                pass
        finally:
            ctx.prec = orig

    #print "not using asymp"
    return ctx.hypsum(1, 2, (a1type, b1type, b2type), [a1, b1, b2], z, **kwargs)



@defun
def hyp2f3(ctx,a1,a2,b1,b2,b3,z,**kwargs):
    a1, a1type = ctx.convert_maybe_rational(a1)
    a2, a2type = ctx.convert_maybe_rational(a2)
    b1, b1type = ctx.convert_maybe_rational(b1)
    b2, b2type = ctx.convert_maybe_rational(b2)
    b3, b3type = ctx.convert_maybe_rational(b3)
    z = ctx.convert(z)

    #return ctx.hypsum(a1r+a2r, a1f+a2f, a1c+a2c, b1r+b2r+b3r,
    #    b1f+b2f+b3f, b1c+b2c+b3c, z)

    absz = abs(z)
    magz = ctx.mag(z)

    # Asymptotic expansion is ~ exp(sqrt(z))
    asymp_extraprec = z and magz//2
    orig = ctx.prec

    # Asymptotic series is in terms of 4F1
    # The square root below empirically provides a plausible criterion
    # for the leading series to converge
    can_use_asymptotic = (not kwargs.get('force_series')) and \
        (ctx.mag(absz) > 19) and (ctx.sqrt(absz) > 1.5*orig)

    if can_use_asymptotic:
        #print "using asymp"
        try:
            try:
                ctx.prec += asymp_extraprec
                # http://functions.wolfram.com/HypergeometricFunctions/
                # Hypergeometric2F3/06/02/03/01/0002/
                def h(a1,a2,b1,b2,b3):
                    X = ctx.mpq_1_2*(a1+a2-b1-b2-b3+ctx.mpq_1_2)
                    A2 = a1+a2
                    B3 = b1+b2+b3
                    A = a1*a2
                    B = b1*b2+b3*b2+b1*b3
                    R = b1*b2*b3
                    c = {}
                    c[0] = ctx.one
                    c[1] = 2*(B - A + ctx.mpq_1_4*(3*A2+B3-2)*(A2-B3) - ctx.mpq_3_16)
                    c[2] = ctx.mpq_1_2*c[1]**2 + ctx.mpq_1_16*(-16*(2*A2-3)*(B-A) + 32*R +\
                        4*(-8*A2**2 + 11*A2 + 8*A + B3 - 2)*(A2-B3)-3)
                    s1 = 0
                    s2 = 0
                    k = 0
                    tprev = 0
                    while 1:
                        if k not in c:
                            uu1 = (k-2*X-3)*(k-2*X-2*b1-1)*(k-2*X-2*b2-1)*\
                                (k-2*X-2*b3-1)
                            uu2 = (4*(k-1)**3 - 6*(4*X+B3)*(k-1)**2 + \
                                2*(24*X**2+12*B3*X+4*B+B3-1)*(k-1) - 32*X**3 - \
                                24*B3*X**2 - 4*B - 8*R - 4*(4*B+B3-1)*X + 2*B3-1)
                            uu3 = (5*(k-1)**2+2*(-10*X+A2-3*B3+3)*(k-1)+2*c[1])
                            c[k] = ctx.one/(2*k)*(uu1*c[k-3]-uu2*c[k-2]+uu3*c[k-1])
                        w = c[k] * (-z)**(-0.5*k)
                        t1 = (-ctx.j)**k * ctx.mpf(2)**(-k) * w
                        t2 = ctx.j**k * ctx.mpf(2)**(-k) * w
                        if abs(t1) < 0.1*ctx.eps:
                            break
                        # Quit if the series doesn't converge quickly enough
                        if k > 5 and abs(tprev) / abs(t1) < 1.5:
                            raise NoConvergence
                        s1 += t1
                        s2 += t2
                        tprev = t1
                        k += 1
                    S = ctx.exp(ctx.j*(ctx.pi*X+2*ctx.sqrt(-z)))*s1 + \
                        ctx.exp(-ctx.j*(ctx.pi*X+2*ctx.sqrt(-z)))*s2
                    T1 = [0.5*S, ctx.pi, -z], [1, -0.5, X], [b1, b2, b3], [a1, a2],\
                        [], [], 0
                    T2 = [-z], [-a1], [b1,b2,b3,a2-a1],[a2,b1-a1,b2-a1,b3-a1], \
                        [a1,a1-b1+1,a1-b2+1,a1-b3+1], [a1-a2+1], 1/z
                    T3 = [-z], [-a2], [b1,b2,b3,a1-a2],[a1,b1-a2,b2-a2,b3-a2], \
                        [a2,a2-b1+1,a2-b2+1,a2-b3+1],[-a1+a2+1], 1/z
                    return T1, T2, T3
                v = ctx.hypercomb(h, [a1,a2,b1,b2,b3])
                if sum(ctx.is_real_type(u) for u in [a1,a2,b1,b2,b3,z]) == 6:
                    v = ctx.re(v)
                return v
            except NoConvergence:
                pass
        finally:
            ctx.prec = orig

    return ctx.hypsum(2, 3, (a1type, a2type, b1type, b2type, b3type), [a1, a2, b1, b2, b3], z, **kwargs)

@defun
def hyp3f2(ctx, a1, a2, a3, b1, b2, z, **kwargs):
    return ctx._hypq1fq(3, 2, [a1,a2,a3], [b1,b2], z, **kwargs)

@defun
def hyp2f0(ctx, a, b, z, **kwargs):
    """
    Hypergeometric 2F0.
    """
    a, atype = ctx.convert_maybe_rational(a)
    b, btype = ctx.convert_maybe_rational(b)
    z = ctx.convert(z)
    # We want to try aggressively to use the asymptotic expansion,  
    # and fall back only when absolutely necessary
    try:
        return ctx.hypsum(2, 0, (atype,btype), [a,b], z, maxterms=ctx.prec)
    except NoConvergence:
        if kwargs.get('force_series'):
            raise
        pass
    def h(a, b):
        w = ctx.sinpi(b)
        rz = -1/z
        T1 = ([ctx.pi,w,rz],[1,-1,a],[],[a-b+1,b],[a],[b],rz)
        T2 = ([-ctx.pi,w,rz],[1,-1,1+a-b],[],[a,2-b],[a-b+1],[2-b],rz)
        return T1, T2
    return ctx.hypercomb(h, [a, 1+a-b], **kwargs)

@defun
def hyperu(ctx, a,b,z, **kwargs):
    a, atype = ctx.convert_maybe_rational(a)
    b, btype = ctx.convert_maybe_rational(b)
    z = ctx.convert(z)
    if not z:
        if ctx.re(b) <= 1:
            return ctx.gammaprod([1-b],[a-b+1])
        else:
            return ctx.inf + z
    bb = 1+a-b
    bb, bbtype = ctx.convert_maybe_rational(bb)
    try:
        orig = ctx.prec
        try:
            ctx.prec += 10
            v = ctx.hypsum(2, 0, (atype, bbtype), [a, bb], -1/z, maxterms=ctx.prec)
            return v / z**a
        finally:
            ctx.prec = orig
    except NoConvergence:
        pass
    def h(a,b):
        w = ctx.sinpi(b)
        T1 = ([ctx.pi,w],[1,-1],[],[a-b+1,b],[a],[b],z)
        T2 = ([-ctx.pi,w,z],[1,-1,1-b],[],[a,2-b],[a-b+1],[2-b],z)
        return T1, T2
    return ctx.hypercomb(h, [a,b], **kwargs)

@defun_wrapped
def _erf_complex(ctx, z):
    z2 = ctx.square_exp_arg(z, -1)
    #z2 = -z**2
    v = (2/ctx.sqrt(ctx.pi))*z * ctx.hyp1f1((1,2),(3,2), z2)
    if not z.real:
        v = v.imag*ctx.j
    return v

@defun_wrapped
def _erfc_complex(ctx, z):
    if ctx.re(z) > 2:
        z2 = ctx.square_exp_arg(z)
        nz2 = ctx.fneg(z2, exact=True)
        v = ctx.exp(nz2)/ctx.sqrt(ctx.pi) * ctx.hyperu((1,2),(1,2), z2)
    else:
        v = 1 - ctx._erf_complex(z)
    if not z.real:
        v = 1+v.imag*ctx.j
    return v

@defun
def erf(ctx, z):
    z = ctx.convert(z)
    if ctx.is_real_type(z):
        try:
            return ctx._erf(z.real)
        except NotImplementedError:
            pass
    if ctx.is_complex_type(z) and not z.imag:
        try:
            return type(z)(ctx._erf(z.real))
        except NotImplementedError:
            pass
    return ctx._erf_complex(z)

@defun
def erfc(ctx, z):
    z = ctx.convert(z)
    if ctx.is_real_type(z):
        try:
            return ctx._erfc(z.real)
        except NotImplementedError:
            pass
    if ctx.is_complex_type(z) and not z.imag:
        try:
            return type(z)(ctx._erfc(z.real))
        except NotImplementedError:
            pass
    return ctx._erfc_complex(z)

@defun
def square_exp_arg(ctx, z, mult=1):
    z2 = ctx.fmul(z, z, prec=ctx.prec*4+20)
    if mult != 1:
        z2 = ctx.fmul(z2, mult, exact=True)
    return z2

@defun_wrapped
def erfi(ctx, z):
    if not z:
        return z
    z2 = ctx.square_exp_arg(z)
    v = (2/ctx.sqrt(ctx.pi)*z) * ctx.hyp1f1((1,2), (3,2), z2)
    if not z.real:
        v = v.imag*ctx.j
    return v

@defun_wrapped
def erfinv(ctx, x):
    if x.imag or (x < -1) or (x > 1):
        return ctx.bad_domain("erfinv(x) is defined only for -1 <= x <= 1")
    if ctx.isnan(x): return x
    if not x: return x
    if x == 1: return ctx.inf
    if x == -1: return ctx.ninf
    if abs(x) < 0.9:
        a = 0.53728*x**3 + 0.813198*x
    else:
        # An asymptotic formula
        u = ctx.ln(2/ctx.pi/(abs(x)-1)**2)
        a = ctx.sign(x) * ctx.sqrt(u - ctx.ln(u))/ctx.sqrt(2)
    ctx.prec += 10
    return ctx.findroot(lambda t: ctx.erf(t)-x, a)

@defun_wrapped
def npdf(ctx, x, mu=0, sigma=1):
    sigma = ctx.convert(sigma)
    return ctx.exp(-(x-mu)**2/(2*sigma**2)) / (sigma*ctx.sqrt(2*ctx.pi))

@defun_wrapped
def ncdf(ctx, x, mu=0, sigma=1):
    a = (x-mu)/(sigma*ctx.sqrt(2))
    if a < 0:
        return ctx.erfc(-a)/2
    else:
        return (1+ctx.erf(a))/2
@defun
def gammainc(ctx, z, a=0, b=None, regularized=False):
    regularized = bool(regularized)
    z = ctx.convert(z)
    if a is None:
        a = ctx.zero
        lower_modified = False
    else:
        a = ctx.convert(a)
        lower_modified = a != ctx.zero
    if b is None:
        b = ctx.inf
        upper_modified = False
    else:
        b = ctx.convert(b)
        upper_modified = b != ctx.inf
    # Complete gamma function
    if not (upper_modified or lower_modified):
        if regularized:
            if ctx.re(z) < 0:
                return ctx.inf
            elif ctx.re(z) > 0:
                return ctx.one
            else:
                return ctx.nan
        return ctx.gamma(z)
    if a == b:
        return ctx.zero
    # Standardize
    if ctx.re(a) > ctx.re(b):
        return -ctx.gammainc(z, b, a, regularized)
    # Generalized gamma
    if upper_modified and lower_modified:
        return +ctx._gamma3(z, a, b, regularized)
    # Upper gamma
    elif lower_modified:
        return ctx._upper_gamma(z, a, regularized)
    # Lower gamma
    elif upper_modified:
        return ctx._lower_gamma(z, b, regularized)

@defun
def _lower_gamma(ctx, z, b, regularized=False):
    # Pole
    if ctx.isnpint(z):
        return type(z)(ctx.inf)
    G = [z] * regularized
    negb = ctx.fneg(b, exact=True)
    def h(z):
        T1 = [ctx.exp(negb), b, z], [1, z, -1], [], G, [1], [1+z], b
        return (T1,)
    return ctx.hypercomb(h, [z])

@defun
def _upper_gamma(ctx, z, a, regularized=False):
    # Fast integer case, when available
    if ctx.isint(z):
        try:
            if regularized:
                # Gamma pole
                if ctx.isnpint(z):
                    return type(z)(ctx.zero)
                orig = ctx.prec
                try:
                    ctx.prec += 10
                    return ctx._gamma_upper_int(z, a) / ctx.gamma(z)
                finally:
                    ctx.prec = orig
            else:
                return ctx._gamma_upper_int(z, a)
        except NotImplementedError:
            pass
    nega = ctx.fneg(a, exact=True)
    G = [z] * regularized
    # Use 2F0 series when possible; fall back to lower gamma representation
    try:
        def h(z):
            r = z-1
            return [([ctx.exp(nega), a], [1, r], [], G, [1, -r], [], 1/nega)]
        return ctx.hypercomb(h, [z], force_series=True)
    except NoConvergence:
        def h(z):
            T1 = [], [1, z-1], [z], G, [], [], 0
            T2 = [-ctx.exp(nega), a, z], [1, z, -1], [], G, [1], [1+z], a
            return T1, T2
        return ctx.hypercomb(h, [z])

@defun
def _gamma3(ctx, z, a, b, regularized=False):
    pole = ctx.isnpint(z)
    if regularized and pole:
        return ctx.zero
    try:
        ctx.prec += 15
        # We don't know in advance whether it's better to write as a difference
        # of lower or upper gamma functions, so try both
        T1 = ctx.gammainc(z, a, regularized=regularized)
        T2 = ctx.gammainc(z, b, regularized=regularized)
        R = T1 - T2
        if ctx.mag(R) - max(ctx.mag(T1), ctx.mag(T2)) > -10:
            return R
        if not pole:
            T1 = ctx.gammainc(z, 0, b, regularized=regularized)
            T2 = ctx.gammainc(z, 0, a, regularized=regularized)
            R = T1 - T2
            # May be ok, but should probably at least print a warning
            # about possible cancellation
            if 1: #ctx.mag(R) - max(ctx.mag(T1), ctx.mag(T2)) > -10:
                return R
    finally:
        ctx.prec -= 15
    raise NotImplementedError

@defun_wrapped
def expint(ctx, n, z):
    if ctx.isint(n) and ctx.is_real_type(z):
        try:
            return ctx._expint_int(n, z)
        except NotImplementedError:
            pass
    if ctx.isnan(n) or ctx.isnan(z):
        return z*n
    if z == ctx.inf:
        return 1/z
    if z == 0:
        # integral from 1 to infinity of t^n
        if ctx.re(n) <= 1:
            # TODO: reasonable sign of infinity
            return type(z)(ctx.inf)
        else:
            return ctx.one/(n-1)
    if n == 0:
        return ctx.exp(-z)/z
    if n == -1:
        return ctx.exp(-z)*(z+1)/z**2
    return z**(n-1) * ctx.gammainc(1-n, z)

@defun_wrapped
def li(ctx, z, offset=False):
    if offset:
        if z == 2:
            return ctx.zero
        return ctx.ei(ctx.ln(z)) - ctx.ei(ctx.ln2)
    if not z:
        return z
    if z == 1:
        return ctx.ninf
    return ctx.ei(ctx.ln(z))

@defun_wrapped
def chi(ctx, z):
    if not z:
        return ctx.ninf
    if z == ctx.inf or z == ctx.ninf:
        return ctx.inf
    z2 = ctx.square_exp_arg(z, 0.25)
    return ctx.euler + ctx.ln(z) + z2*ctx.hyp2f3(1,1,2,2,(3,2),z2)

@defun_wrapped
def shi(ctx, z):
    if z == ctx.inf:
        return z
    if z == ctx.ninf:
        return z
    z2 = ctx.square_exp_arg(z, 0.25)
    return z*ctx.hyp1f2((1,2),(3,2),(3,2),z2)

@defun_wrapped
def fresnels(ctx, z):
    if z == ctx.inf:
        return ctx.mpf(0.5)
    if z == ctx.ninf:
        return ctx.mpf(-0.5)
    #return ctx.pi*z**3/6*ctx.hypsum(1,2,('Q','Q','Q'),\
    #    [ctx.mpq_3_4,ctx.mpq_3_2,ctx.mpq_7_4],-ctx.pi**2*z**4/16)
    return ctx.pi*z**3/6*ctx.hyp1f2((3,4),(3,2),(7,4),-ctx.pi**2*z**4/16)

@defun_wrapped
def fresnelc(ctx, z):
    if z == ctx.inf:
        return ctx.mpf(0.5)
    if z == ctx.ninf:
        return ctx.mpf(-0.5)
    #return z*ctx.hypsum(1,2,('Q','Q','Q'),\
    #    [ctx.mpq_1_4,ctx.mpq_1_2,ctx.mpq_5_4],-ctx.pi**2*z**4/16)
    return z*ctx.hyp1f2((1,4),(1,2),(5,4),-ctx.pi**2*z**4/16)

@defun_wrapped
def airyai(ctx, z):
    if z == ctx.inf or z == ctx.ninf:
        return ctx.zero
    if z:
        # Account for exponential scaling
        ctx.prec += max(0, int(1.5*ctx.mag(z)))
    if z.real > 4:
        # We could still use 1F1, but it results in huge cancellation;
        # the following expansion is better
        w = z**1.5
        r = -ctx.mpf(3)/(4*w)
        v = ctx.exp(-2*w/3)/(2*ctx.sqrt(ctx.pi)*ctx.nthroot(z,4))
        v *= ctx.hyp2f0((1,6),(5,6),r)
        return v
    elif z.real > 1:
        # If not using asymptotic series:
        # cancellation: both terms are ~ 2^(z^1.5),
        # result is ~ 2^(-z^1.5), so need ~2*z^1.5 extra bits
        ctx.prec += 2*int(z.real**1.5)
    z3 = z**3 / 9
    a = ctx.hyp0f1((2,3), z3) / (ctx.cbrt(9) * ctx.gamma(ctx.mpf(2)/3))
    b = z * ctx.hyp0f1((4,3), z3) / (ctx.cbrt(3) * ctx.gamma(ctx.mpf(1)/3))
    return a - b

@defun_wrapped
def airybi(ctx, z):
    if z == ctx.inf:
        return z
    if z == ctx.ninf:
        return 1/z
    if z:
        # Account for exponential scaling
        ctx.prec += max(0, int(1.5*ctx.mag(z)))
    z3 = z**3 / 9
    rt = ctx.nthroot(3, 6)
    a = ctx.hyp0f1((2,3), z3) / (rt * ctx.gamma(ctx.mpf(2)/3))
    b = z * rt * ctx.hyp0f1((4,3), z3) / ctx.gamma(ctx.mpf(1)/3)
    return a + b

@defun_wrapped
def hermite(ctx, n, z, **kwargs):
    if not z:
        try:
            return 2**n * ctx.sqrt(ctx.pi) / ctx.gamma(0.5*(1-n))
        except ValueError:
            return 0.0*(n+z)
    if ctx.re(z) > 0 or (ctx.re(z) == 0 and ctx.im(z) > 0) or ctx.isnpint(-n):
        prec = ctx.prec
        ctx.prec = ctx.prec*4+20
        z2 = -z**(-2)
        ctx.prec = prec
        return (2*z)**n * ctx.hyp2f0(-0.5*n, -0.5*(n-1), z2, **kwargs)
    else:
        prec = ctx.prec
        ctx.prec = ctx.prec*4+20
        z2 = z**2
        ctx.prec = prec
        return ctx.hermite(n,-z) + 2**(n+2)*ctx.sqrt(ctx.pi) * (-z) / \
            ctx.gamma(-0.5*n) * ctx.hyp1f1((1-n)*0.5, 1.5, z2, **kwargs)

@defun_wrapped
def gegenbauer(ctx, n, a, z, **kwargs):
    # Special cases: a+0.5, a*2 poles
    if ctx.isnpint(a):
        return 0*(z+n)
    if ctx.isnpint(a+0.5):
        # TODO: something else is required here
        # E.g.: gegenbauer(-2, -0.5, 3) == -12
        if ctx.isnpint(n+1):
            raise NotImplementedError("Gegenbauer function with two limits")
        def h(a):
            a2 = 2*a
            T = [], [], [n+a2], [n+1, a2], [-n, n+a2], [a+0.5], 0.5*(1-z)
            return [T]
        return ctx.hypercomb(h, [a], **kwargs)
    def h(n):
        a2 = 2*a
        T = [], [], [n+a2], [n+1, a2], [-n, n+a2], [a+0.5], 0.5*(1-z)
        return [T]
    return ctx.hypercomb(h, [n], **kwargs)

@defun_wrapped
def jacobi(ctx, n, a, b, x, **kwargs):
    if not ctx.isnpint(a):
        def h(n):
            return (([], [], [a+n+1], [n+1, a+1], [-n, a+b+n+1], [a+1], (1-x)*0.5),)
        return ctx.hypercomb(h, [n], **kwargs)
    if not ctx.isint(b):
        def h(n, a):
            return (([], [], [-b], [n+1, -b-n], [-n, a+b+n+1], [b+1], (x+1)*0.5),)
        return ctx.hypercomb(h, [n, a], **kwargs)
    # XXX: determine appropriate limit
    return ctx.binomial(n+a,n) * ctx.hyp2f1(-n,1+n+a+b,a+1,(1-x)/2, **kwargs)

@defun_wrapped
def laguerre(ctx, n, a, z, **kwargs):
    # XXX: limits, poles
    #if ctx.isnpint(n):
    #    return 0*(a+z)
    def h(a):
        return (([], [], [a+n+1], [a+1, n+1], [-n], [a+1], z),)
    return ctx.hypercomb(h, [a], **kwargs)

@defun_wrapped
def legendre(ctx, n, x, **kwargs):
    if ctx.isint(n):
        n = int(n)
        # Accuracy near zeros
        if (n + (n < 0)) & 1:
            if not x:
                return x
            mag = ctx.mag(x)
            if mag < -2*ctx.prec-10:
                return x
            if mag < -5:
                ctx.prec += -mag
    return ctx.hyp2f1(-n,n+1,1,(1-x)/2, **kwargs)

@defun
def legenp(ctx, n, m, z, type=2, **kwargs):
    # Legendre function, 1st kind
    n = ctx.convert(n)
    m = ctx.convert(m)
    # Faster
    if not m:
        return ctx.legendre(n, z, **kwargs)
    # TODO: correct evaluation at singularities
    if type == 2:
        def h(n,m):
            g = m*0.5
            T = [1+z, 1-z], [g, -g], [], [1-m], [-n, n+1], [1-m], 0.5*(1-z)
            return (T,)
        return ctx.hypercomb(h, [n,m], **kwargs)
    if type == 3:
        def h(n,m):
            g = m*0.5
            T = [z+1, z-1], [g, -g], [], [1-m], [-n, n+1], [1-m], 0.5*(1-z)
            return (T,)
        return ctx.hypercomb(h, [n,m], **kwargs)
    raise ValueError("requires type=2 or type=3")

@defun
def legenq(ctx, n, m, z, type=2, **kwargs):
    # Legendre function, 2nd kind
    n = ctx.convert(n)
    m = ctx.convert(m)
    z = ctx.convert(z)
    if z in (1, -1):
        #if ctx.isint(m):
        #    return ctx.nan
        #return ctx.inf  # unsigned
        return ctx.nan
    if type == 2:
        def h(n, m):
            s = 2 * ctx.sinpi(m) / ctx.pi
            c = ctx.cospi(m)
            a = 1+z
            b = 1-z
            u = m/2
            w = (1-z)/2
            T1 = [s, c, a, b], [-1, 1, u, -u], [], [1-m], \
                [-n, n+1], [1-m], w
            T2 = [-s, a, b], [-1, -u, u], [n+m+1], [n-m+1, m+1], \
                [-n, n+1], [m+1], w
            return T1, T2
        return ctx.hypercomb(h, [n, m], **kwargs)
    if type == 3:
        # The following is faster when there only is a single series
        # Note: not valid for -1 < z < 0 (?)
        if abs(z) > 1:
            def h(n, m):
                T1 = [ctx.exp(ctx.j*m*ctx.pi), 2, ctx.pi, z, z-1, z+1], \
                     [1, -n-1, 0.5, -n-m-1, 0.5*m, 0.5*m], \
                     [n+m+1], [n+1.5], \
                     [0.5*(2+n+m), 0.5*(1+n+m)], [n+1.5], z**(-2)
                return [T1]
            return ctx.hypercomb(h, [n, m], **kwargs)
        else:
            # not valid for 1 < z < inf ?
            def h(n, m):
                s = 2 * ctx.sinpi(m) / ctx.pi
                c = ctx.exp(ctx.pi*ctx.j*m)
                a = 1+z
                b = z-1
                u = m/2
                w = (1-z)/2
                T1 = [s, c, a, b], [-1, 1, u, -u], [], [1-m], \
                    [-n, n+1], [1-m], w
                T2 = [-s, c, a, b], [-1, 1, -u, u], [n+m+1], [n-m+1, m+1], \
                    [-n, n+1], [m+1], w
                return T1, T2
            return ctx.hypercomb(h, [n, m], **kwargs)
    raise ValueError("requires type=2 or type=3")

@defun_wrapped
def chebyt(ctx, n, x, **kwargs):
    return ctx.hyp2f1(-n,n,(1,2),(1-x)/2, **kwargs)

@defun_wrapped
def chebyu(ctx, n, x, **kwargs):
    return (n+1) * ctx.hyp2f1(-n, n+2, (3,2), (1-x)/2, **kwargs)

@defun
def j0(ctx, x):
    """Computes the Bessel function `J_0(x)`. See :func:`besselj`."""
    return ctx.besselj(0, x)

@defun
def j1(ctx, x):
    """Computes the Bessel function `J_1(x)`.  See :func:`besselj`."""
    return ctx.besselj(1, x)

@defun
def besselj(ctx, n, z, derivative=0, **kwargs):
    if type(n) is int:
        n_isint = True
    else:
        n = ctx.convert(n)
        n_isint = ctx.isint(n)
        if n_isint:
            n = int(n)
    if n_isint and n < 0:
        return (-1)**n * ctx.besselj(-n, z, derivative, **kwargs)
    z = ctx.convert(z)
    M = ctx.mag(z)
    if derivative:
        d = ctx.convert(derivative)
        # TODO: the integer special-casing shouldn't be necessary.
        # However, the hypergeometric series gets inaccurate for large d
        # because of inaccurate pole cancellation at a pole far from
        # zero (needs to be fixed in hypercomb or hypsum)
        if ctx.isint(d) and d >= 0:
            d = int(d)
            orig = ctx.prec
            try:
                ctx.prec += 15
                v = ctx.fsum((-1)**k * ctx.binomial(d,k) * ctx.besselj(2*k+n-d,z)
                    for k in range(d+1))
            finally:
                ctx.prec = orig
            v *= ctx.mpf(2)**(-d)
        else:
            def h(n,d):
                r = ctx.fmul(ctx.fmul(z, z, prec=ctx.prec+M), -0.25, exact=True)
                B = [0.5*(n-d+1), 0.5*(n-d+2)]
                T = [([2,ctx.pi,z],[d-2*n,0.5,n-d],[],B,[(n+1)*0.5,(n+2)*0.5],B+[n+1],r)]
                return T
            v = ctx.hypercomb(h, [n,d], **kwargs)
    else:
    # Fast case: J_n(x), n int, appropriate magnitude for fixed-point calculation
        if (not derivative) and n_isint and abs(M) < 10 and abs(n) < 20:
            try:
                return ctx._besselj(n, z)
            except NotImplementedError:
                pass
        if not z:
            if not n:
                v = ctx.one + n+z
            elif ctx.re(n) > 0:
                v = n*z
            else:
                v = ctx.inf + z + n
        else:
            #v = 0
            orig = ctx.prec
            try:
                # XXX: workaround for accuracy in low level hypergeometric series
                # when alternating, large arguments
                ctx.prec += min(3*abs(M), ctx.prec)
                w = ctx.fmul(z, 0.5, exact=True)
                def h(n):
                    r = ctx.fneg(ctx.fmul(w, w, prec=max(0,ctx.prec+M)), exact=True)
                    return [([w], [n], [], [n+1], [], [n+1], r)]
                v = ctx.hypercomb(h, [n], **kwargs)
            finally:
                ctx.prec = orig
        v = +v
    return v

@defun
def besseli(ctx, n, z, derivative=0, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    if not z:
        if derivative:
            raise ValueError
        if not n:
            # I(0,0) = 1
            return 1+n+z
        if ctx.isint(n):
            return 0*(n+z)
        r = ctx.re(n)
        if r == 0:
            return ctx.nan*(n+z)
        elif r > 0:
            return 0*(n+z)
        else:
            return ctx.inf+(n+z)
    M = ctx.mag(z)
    if derivative:
        d = ctx.convert(derivative)
        def h(n,d):
            r = ctx.fmul(ctx.fmul(z, z, prec=ctx.prec+M), 0.25, exact=True)
            B = [0.5*(n-d+1), 0.5*(n-d+2), n+1]
            T = [([2,ctx.pi,z],[d-2*n,0.5,n-d],[n+1],B,[(n+1)*0.5,(n+2)*0.5],B,r)]
            return T
        v = ctx.hypercomb(h, [n,d], **kwargs)
    else:
        def h(n):
            w = ctx.fmul(z, 0.5, exact=True)
            r = ctx.fmul(w, w, prec=max(0,ctx.prec+M))
            return [([w], [n], [], [n+1], [], [n+1], r)]
        v = ctx.hypercomb(h, [n], **kwargs)
    return v

@defun_wrapped
def bessely(ctx, n, z, derivative=0, **kwargs):
    if not z:
        if derivative:
            # Not implemented
            raise ValueError
        if not n:
            # ~ log(z/2)
            return -ctx.inf + (n+z)
        if ctx.im(n):
            return nan * (n+z)
        r = ctx.re(n)
        q = n+0.5
        if ctx.isint(q):
            if n > 0:
                return -ctx.inf + (n+z)
            else:
                return 0 * (n+z)
        if r < 0 and int(ctx.floor(q)) % 2:
            return ctx.inf + (n+z)
        else:
            return ctx.ninf + (n+z)
    # XXX: use hypercomb
    ctx.prec += 10
    m, d = ctx.nint_distance(n)
    if d < -ctx.prec:
        h = +ctx.eps
        ctx.prec *= 2
        n += h
    elif d < 0:
        ctx.prec -= d
    # TODO: avoid cancellation for imaginary arguments
    return (ctx.besselj(n,z,derivative,**kwargs)*ctx.cospi(n) - \
        ctx.besselj(-n,z,derivative,**kwargs))/ctx.sinpi(n)

@defun_wrapped
def besselk(ctx, n, z, **kwargs):
    if not z:
        return ctx.inf
    M = ctx.mag(z)
    if M < 1:
        # Represent as limit definition
        def h(n):
            r = (z/2)**2
            T1 = [z, 2], [-n, n-1], [n], [], [], [1-n], r
            T2 = [z, 2], [n, -n-1], [-n], [], [], [1+n], r
            return T1, T2
    # We could use the limit definition always, but it leads
    # to very bad cancellation (of exponentially large terms)
    # for large real z
    # Instead represent in terms of 2F0
    else:
        ctx.prec += M
        def h(n):
            return [([ctx.pi/2, z, ctx.exp(-z)], [0.5,-0.5,1], [], [], \
                [n+0.5, 0.5-n], [], -1/(2*z))]
    return ctx.hypercomb(h, [n], **kwargs)

@defun_wrapped
def hankel1(ctx,n,x,**kwargs):
    return ctx.besselj(n,x,**kwargs) + ctx.j*ctx.bessely(n,x,**kwargs)

@defun_wrapped
def hankel2(ctx,n,x,**kwargs):
    return ctx.besselj(n,x,**kwargs) - ctx.j*ctx.bessely(n,x,**kwargs)

@defun_wrapped
def whitm(ctx,k,m,z,**kwargs):
    if z == 0:
        # M(k,m,z) = 0^(1/2+m)
        if ctx.re(m) > -0.5:
            return z
        elif ctx.re(m) < -0.5:
            return ctx.inf + z
        else:
            return ctx.nan * z
    x = ctx.fmul(-0.5, z, exact=True)
    y = 0.5+m
    return ctx.exp(x) * z**y * ctx.hyp1f1(y-k, 1+2*m, z, **kwargs)

@defun_wrapped
def whitw(ctx,k,m,z,**kwargs):
    if z == 0:
        g = abs(ctx.re(m))
        if g < 0.5:
            return z
        elif g > 0.5:
            return ctx.inf + z
        else:
            return ctx.nan * z
    x = ctx.fmul(-0.5, z, exact=True)
    y = 0.5+m
    return ctx.exp(x) * z**y * ctx.hyperu(y-k, 1+2*m, z, **kwargs)

@defun
def struveh(ctx,n,z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/StruveH/26/01/02/
    def h(n):
        return [([z/2, 0.5*ctx.sqrt(ctx.pi)], [n+1, -1], [], [n+1.5], [1], [1.5, n+1.5], -(z/2)**2)]
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def struvel(ctx,n,z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/StruveL/26/01/02/
    def h(n):
        return [([z/2, 0.5*ctx.sqrt(ctx.pi)], [n+1, -1], [], [n+1.5], [1], [1.5, n+1.5], (z/2)**2)]
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def ber(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinBer2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [ctx.cospi(0.75*n), z/2], [1, n], [], [n+1], [], [0.5, 0.5*(n+1), 0.5*n+1], r
        T2 = [-ctx.sinpi(0.75*n), z/2], [1, n+2], [], [n+2], [], [1.5, 0.5*(n+3), 0.5*n+1], r
        return T1, T2
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def bei(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinBei2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [ctx.cospi(0.75*n), z/2], [1, n+2], [], [n+2], [], [1.5, 0.5*(n+3), 0.5*n+1], r
        T2 = [ctx.sinpi(0.75*n), z/2], [1, n], [], [n+1], [], [0.5, 0.5*(n+1), 0.5*n+1], r
        return T1, T2
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def ker(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinKer2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [2, z, 4*ctx.cospi(0.25*n)], [-n-3, n, 1], [-n], [], [], [0.5, 0.5*(1+n), 0.5*(n+2)], r
        T2 = [2, z, -ctx.sinpi(0.25*n)], [-n-3, 2+n, 1], [-n-1], [], [], [1.5, 0.5*(3+n), 0.5*(n+2)], r
        T3 = [2, z, 4*ctx.cospi(0.75*n)], [n-3, -n, 1], [n], [], [], [0.5, 0.5*(1-n), 1-0.5*n], r
        T4 = [2, z, -ctx.sinpi(0.75*n)], [n-3, 2-n, 1], [n-1], [], [], [1.5, 0.5*(3-n), 1-0.5*n], r
        return T1, T2, T3, T4
    return ctx.hypercomb(h, [n], **kwargs)

@defun
def kei(ctx, n, z, **kwargs):
    n = ctx.convert(n)
    z = ctx.convert(z)
    # http://functions.wolfram.com/Bessel-TypeFunctions/KelvinKei2/26/01/02/0001/
    def h(n):
        r = -(z/4)**4
        T1 = [-ctx.cospi(0.75*n), 2, z], [1, n-3, 2-n], [n-1], [], [], [1.5, 0.5*(3-n), 1-0.5*n], r
        T2 = [-ctx.sinpi(0.75*n), 2, z], [1, n-1, -n], [n], [], [], [0.5, 0.5*(1-n), 1-0.5*n], r
        T3 = [-ctx.sinpi(0.25*n), 2, z], [1, -n-1, n], [-n], [], [], [0.5, 0.5*(n+1), 0.5*(n+2)], r
        T4 = [-ctx.cospi(0.25*n), 2, z], [1, -n-3, n+2], [-n-1], [], [], [1.5, 0.5*(n+3), 0.5*(n+2)], r
        return T1, T2, T3, T4
    return ctx.hypercomb(h, [n], **kwargs)

@defun_wrapped
def lambertw(ctx, z, k=0):
    k = int(k)
    if ctx.isnan(z):
        return z
    ctx.prec += 20
    absz = abs(z)
    # We must be extremely careful near the singularities at -1/e and 0
    u = ctx.exp(-1)
    if absz <= u:
        if not z:
            # w(0,0) = 0; for all other branches we hit the pole
            if not k:
                return z
            return ctx.ninf
        if not k:
            w = z
        # For small real z < 0, the -1 branch aves roughly like log(-z)
        elif k == -1 and not ctx.im(z) and ctx.re(z) < 0:
            w = ctx.ln(-z)
        # Use a simple asymptotic approximation.
        else:
            w = ctx.ln(z)
            # The branches are roughly logarithmic. This approximation
            # gets better for large |k|; need to check that this always
            # works for k ~= -1, 0, 1.
            if k: w += k * 2*ctx.pi*ctx.j
    elif k == 0 and ctx.im(z) and absz <= 0.7:
        # Both the W(z) ~= z and W(z) ~= ln(z) approximations break
        # down around z ~= -0.5 (converging to the wrong branch), so patch 
        # with a constant approximation (adjusted for sign)
        if abs(z+0.5) < 0.1:
            if ctx.im(z) > 0:
                w = ctx.mpc(0.7+0.7j)
            else:
                w = ctx.mpc(0.7-0.7j)
        else:
            w = z
    else:
        if z == ctx.inf:
            if k == 0:
                return z
            else:
                return z + 2*k*ctx.pi*ctx.j
        if z == ctx.ninf:
            return (-z) + (2*k+1)*ctx.pi*ctx.j
        # Simple asymptotic approximation as above
        w = ctx.ln(z)
        if k: w += k * 2*ctx.pi*ctx.j
    # Use Halley iteration to solve w*exp(w) = z
    two = ctx.mpf(2)
    weps = ctx.ldexp(ctx.eps, 15)
    for i in xrange(100):
        ew = ctx.exp(w)
        wew = w*ew
        wewz = wew-z
        wn = w - wewz/(wew+ew-(w+two)*wewz/(two*w+two))
        if abs(wn-w) < weps*abs(wn):
            return wn
        else:
            w = wn
    ctx.warn("Lambert W iteration failed to converge for %s" % z)
    return wn

@defun_wrapped
def barnesg(ctx, z):
    if ctx.isinf(z):
        if z == ctx.inf:
            return z
        return ctx.nan
    if ctx.isnan(z):
        return z
    if (not z.imag) and z.real <= 0 and ctx.isint(z.real):
        return z*0
    # Account for size (would not be needed if computing log(G))
    if abs(z) > 5:
        ctx.dps += 2*ctx.log(abs(z),2)
    # Reflection formula
    if ctx.re(z) < -ctx.dps:
        w = 1-z
        pi2 = 2*ctx.pi
        u = ctx.exp(ctx.j*pi2*w)
        v = ctx.j*ctx.pi/12 - ctx.j*ctx.pi*w**2/2 + w*ctx.ln(1-u) - \
            ctx.j*ctx.polylog(2, u)/pi2
        v = ctx.barnesg(2-z)*ctx.exp(v)/pi2**w
        if ctx.is_real_type(z):
            v = v.real
        return v
    # Estimate terms for asymptotic expansion
    # TODO: fixme, obviously
    N = ctx.dps // 2 + 5
    G = 1
    while abs(z) < N or ctx.re(z) < 1:
        G /= ctx.gamma(z)
        z += 1
    z -= 1
    s = ctx.mpf(1)/12
    s -= ctx.log(ctx.glaisher)
    s += z*ctx.log(2*ctx.pi)/2
    s += (z**2/2-ctx.mpf(1)/12)*ctx.log(z)
    s -= 3*z**2/4
    z2k = z2 = z**2
    for k in xrange(1, N+1):
        t = ctx.bernoulli(2*k+2) / (4*k*(k+1)*z2k)
        if abs(t) < ctx.eps:
            #print k, N      # check how many terms were needed
            break
        z2k *= z2
        s += t
    #if k == N:
    #    print "warning: series for barnesg failed to converge", ctx.dps
    return G*ctx.exp(s)

@defun
def superfac(ctx, z):
    return ctx.barnesg(z+2)

@defun_wrapped
def hyperfac(ctx, z):
    # XXX: estimate needed extra bits accurately
    if z == ctx.inf:
        return z
    if abs(z) > 5:
        extra = 4*int(ctx.log(abs(z),2))
    else:
        extra = 0
    ctx.prec += extra
    if not z.imag and z.real < 0 and ctx.isint(z.real):
        n = int(ctx.re(z))
        h = ctx.hyperfac(-n-1)
        if ((n+1)//2) & 1:
            h = -h
        if ctx.is_complex_type(z):
            return h + 0j
        return h
    zp1 = z+1
    # Wrong branch cut
    #v = ctx.gamma(zp1)**z
    #ctx.prec -= extra
    #return v / ctx.barnesg(zp1)
    v = ctx.exp(z*ctx.loggamma(zp1))
    ctx.prec -= extra
    return v / ctx.barnesg(zp1)

@defun_wrapped
def loggamma(ctx, z):
    a = z.real
    b = z.imag
    if not b and a > 0:
        return ctx.ln(ctx.gamma(z))
    u = ctx.arg(z)
    w = ctx.ln(ctx.gamma(z))
    if b:
        gi = -b - u/2 + a*u + b*ctx.ln(abs(z))
        n = ctx.floor((gi-w.imag)/(2*ctx.pi)+0.5) * (2*ctx.pi)
        return w + n*ctx.j
    elif a < 0:
        n = int(ctx.floor(a))
        w += (n-(n%2))*ctx.pi*ctx.j
    return w

@defun_wrapped
def siegeltheta(ctx, t):
    if t.imag:
        # XXX: cancellation occurs
        a = ctx.loggamma(0.25+0.5j*t)
        b = ctx.loggamma(0.25-0.5j*t)
        return -ctx.ln(ctx.pi)/2*t - 0.5j*(a-b)
    else:
        if ctx.isinf(t):
            return t
        return ctx.loggamma(0.25+0.5j*t).imag - ctx.ln(ctx.pi)/2*t

@defun_wrapped
def grampoint(ctx, n):
    # asymptotic expansion, from
    # http://mathworld.wolfram.com/GramPoint.html
    g = 2*ctx.pi*ctx.exp(1+ctx.lambertw((8*n+1)/(8*ctx.e)))
    return ctx.findroot(lambda t: ctx.siegeltheta(t)-ctx.pi*n, g)

@defun_wrapped
def siegelz(ctx, t):
    v = ctx.exp(ctx.j*ctx.siegeltheta(t))*ctx.zeta(0.5+ctx.j*t)
    if ctx.is_real_type(t):
        return v.real
    return v

_zeta_zeros = [
14.134725142,21.022039639,25.010857580,30.424876126,32.935061588,
37.586178159,40.918719012,43.327073281,48.005150881,49.773832478,
52.970321478,56.446247697,59.347044003,60.831778525,65.112544048,
67.079810529,69.546401711,72.067157674,75.704690699,77.144840069,
79.337375020,82.910380854,84.735492981,87.425274613,88.809111208,
92.491899271,94.651344041,95.870634228,98.831194218,101.317851006,
103.725538040,105.446623052,107.168611184,111.029535543,111.874659177,
114.320220915,116.226680321,118.790782866,121.370125002,122.946829294,
124.256818554,127.516683880,129.578704200,131.087688531,133.497737203,
134.756509753,138.116042055,139.736208952,141.123707404,143.111845808,
146.000982487,147.422765343,150.053520421,150.925257612,153.024693811,
156.112909294,157.597591818,158.849988171,161.188964138,163.030709687,
165.537069188,167.184439978,169.094515416,169.911976479,173.411536520,
174.754191523,176.441434298,178.377407776,179.916484020,182.207078484,
184.874467848,185.598783678,187.228922584,189.416158656,192.026656361,
193.079726604,195.265396680,196.876481841,198.015309676,201.264751944,
202.493594514,204.189671803,205.394697202,207.906258888,209.576509717,
211.690862595,213.347919360,214.547044783,216.169538508,219.067596349,
220.714918839,221.430705555,224.007000255,224.983324670,227.421444280,
229.337413306,231.250188700,231.987235253,233.693404179,236.524229666,
]

def _load_zeta_zeros(url):
    import urllib
    d = urllib.urlopen(url)
    L = [float(x) for x in d.readlines()]
    # Sanity check
    assert round(L[0]) == 14
    _zeta_zeros[:] = L

@defun
def zetazero(ctx, n, url='http://www.dtc.umn.edu/~odlyzko/zeta_tables/zeros1'):
    n = int(n)
    if n < 0:
        return ctx.zetazero(-n).conjugate()
    if n == 0:
        raise ValueError("n must be nonzero")
    if n > len(_zeta_zeros) and n <= 100000:
        _load_zeta_zeros(url)
    if n > len(_zeta_zeros):
        raise NotImplementedError("n too large for zetazeros")
    return ctx.mpc(0.5, ctx.findroot(ctx.siegelz, _zeta_zeros[n-1]))

@defun_wrapped
def riemannr(ctx, x):
    if x == 0:
        return ctx.zero
    # Check if a simple asymptotic estimate is accurate enough
    if abs(x) > 1000:
        a = ctx.li(x)
        b = 0.5*ctx.li(ctx.sqrt(x))
        if abs(b) < abs(a)*ctx.eps:
            return a
    if abs(x) < 0.01:
        # XXX
        ctx.prec += int(-ctx.log(abs(x),2))
    # Sum Gram's series
    s = t = ctx.one
    u = ctx.ln(x)
    k = 1
    while abs(t) > abs(s)*ctx.eps:
        t = t * u / k
        s += t / (k * ctx.zeta(k+1))
        k += 1
    return s

@defun_static
def primepi(ctx, x):
    x = int(x)
    if x < 2:
        return 0
    from gammazeta import list_primes
    return len(list_primes(x))

@defun_wrapped
def primepi2(ctx, x):
    x = int(x)
    if x < 2:
        return ctx.mpi(0,0)
    if x < 2657:
        return ctx.mpi(ctx.primepi(x))
    mid = ctx.li(x)
    # Schoenfeld's estimate for x >= 2657, assuming RH
    err = ctx.sqrt(x,rounding='u')*ctx.ln(x,rounding='u')/8/ctx.pi(rounding='d')
    a = ctx.floor((ctx.mpi(mid)-err).a, rounding='d')
    b = ctx.ceil((ctx.mpi(mid)+err).b, rounding='u')
    return ctx.mpi(a, b)

@defun_wrapped
def primezeta(ctx, s):
    if ctx.isnan(s):
        return s
    if ctx.re(s) <= 0:
        raise ValueError("prime zeta function defined only for re(s) > 0")
    if s == 1:
        return ctx.inf
    if s == 0.5:
        return ctx.mpc(ctx.ninf, ctx.pi)
    r = ctx.re(s)
    if r > ctx.prec:
        return 0.5**s
    else:
        wp = ctx.prec + int(r)
        def terms():
            orig = ctx.prec
            # zeta ~ 1+eps; need to set precision
            # to get logarithm accurately
            k = 0
            while 1:
                k += 1
                u = libintmath.moebius(k)
                if not u:
                    continue
                ctx.prec = wp
                t = u*ctx.ln(ctx.zeta(k*s))/k
                if not t:
                    return
                #print ctx.prec, ctx.nstr(t)
                ctx.prec = orig
                yield t
    return sum_accurately(ctx, terms)

@defun_wrapped
def bernpoly(ctx, n, z):
    # Slow implementation:
    #return sum(ctx.binomial(n,k)*ctx.bernoulli(k)*z**(n-k) for k in xrange(0,n+1))
    n = int(n)
    if n <= 3:
        if n < 0:
            raise ValueError("Bernoulli polynomials only defined for n >= 0")
        if n == 0: return z ** 0
        if n == 1: return z - 0.5
        if n == 2: return (6*z*(z-1)+1)/6
        if n == 3: return z*(z*(z-1.5)+0.5)
    if abs(z) > 2:
        s = t = ctx.one
        r = ctx.one/z
        for k in xrange(1,n+1):
            t = t*(n+1-k)/k*r
            if not (k > 2 and k & 1):
                u = t*ctx.bernoulli(k)
                s += u
                if abs(u) < ctx.eps:
                    break
        return z**n * s
    s = ctx.bernoulli(n)
    t = ctx.one
    for k in xrange(1,n+1):
        t = t*(n+1-k)/k * z
        m = n-k
        if not (m > 2 and m & 1):
            u = t*ctx.bernoulli(m)
            s += u
            #if abs(u) < ctx.eps:
            #    break
    return s

# TODO: this should be implemented low-level
def polylog_series(ctx, s, z):
    tol = +ctx.eps
    l = ctx.zero
    k = 1
    zk = z
    while 1:
        term = zk / k**s
        l += term
        if abs(term) < tol:
            break
        zk *= z
        k += 1
    return l

def polylog_continuation(ctx, n, z):
    if n < 0:
        return z*0
    twopij = 2j * ctx.pi
    a = -twopij**n/ctx.fac(n) * ctx.bernpoly(n, ctx.ln(z)/twopij)
    if ctx.is_real_type(z) and z < 0:
        a = a.real
    if z.imag < 0 or (z.imag == 0 and z.real >= 1):
        a -= twopij*ctx.ln(z)**(n-1)/ctx.fac(n-1)
    return a

def polylog_unitcircle(ctx, n, z):
    tol = +ctx.eps
    if n > 1:
        l = ctx.zero
        logz = ctx.ln(z)
        logmz = ctx.one
        m = 0
        while 1:
            if (n-m) != 1:
                term = ctx.zeta(n-m) * logmz / ctx.fac(m)
                if term and abs(term) < tol:
                    break
                l += term
            logmz *= logz
            m += 1
        l += ctx.ln(z)**(n-1)/ctx.fac(n-1)*(ctx.harmonic(n-1)-ctx.ln(-ctx.ln(z)))
    elif n < 1:  # else
        l = ctx.fac(-n)*(-ctx.ln(z))**(n-1)
        logz = ctx.ln(z)
        logkz = ctx.one
        k = 0
        while 1:
            b = ctx.bernoulli(k-n+1)
            if b:
                term = b*logkz/(ctx.fac(k)*(k-n+1))
                if abs(term) < tol:
                    break
                l -= term
            logkz *= logz
            k += 1
    else:
        raise ValueError
    if ctx.is_real_type(z) and z < 0:
        l = l.real
    return l

def polylog_general(ctx, s, z):
    v = ctx.zero
    u = ctx.ln(z)
    if not abs(u) < 5: # theoretically |u| < 2*pi
        raise NotImplementedError("polylog for arbitrary s and z")
    t = 1
    k = 0
    while 1:
        term = ctx.zeta(s-k) * t
        if abs(term) < ctx.eps:
            break
        v += term
        k += 1
        t *= u
        t /= k
    return ctx.gamma(1-s)*(-u)**(s-1) + v

@defun_wrapped
def polylog(ctx, s, z):
    s = ctx.convert(s)
    z = ctx.convert(z)
    if z == 1:
        return ctx.zeta(s)
    if z == -1:
        return -ctx.altzeta(s)
    if s == 0:
        return z/(1-z)
    if s == 1:
        return -ctx.ln(1-z)
    if s == -1:
        return z/(1-z)**2
    if abs(z) <= 0.75 or (not ctx.isint(s) and abs(z) < 0.9):
        return polylog_series(ctx, s, z)
    if abs(z) >= 1.4 and ctx.isint(s):
        return (-1)**(s+1)*polylog_series(ctx, s, 1/z) + polylog_continuation(ctx, s, z)
    if ctx.isint(s):
        return polylog_unitcircle(ctx, int(s), z)
    return polylog_general(ctx, s, z)

    #raise NotImplementedError("polylog for arbitrary s and z")
    # This could perhaps be used in some cases
    #from quadrature import quad
    #return quad(lambda t: t**(s-1)/(exp(t)/z-1),[0,inf])/gamma(s)

@defun
def clsin(ctx, s, z):
    if ctx.isint(s) and s < 0 and int(s) % 2 == 1:
        return z*0
    a = ctx.exp(ctx.j*z)
    if ctx.is_real_type(z) and ctx.is_real_type(s):
        return ctx.im(ctx.polylog(s,a))
    b = 1/a
    return (-0.5j)*(ctx.polylog(s,a) - ctx.polylog(s,b))

@defun_wrapped
def clcos(ctx, s, z):
    if ctx.isint(s) and s < 0 and int(s) % 2 == 0:
        return z*0
    a = ctx.exp(ctx.j*z)
    if ctx.is_real_type(z) and ctx.is_real_type(s):
        return ctx.re(ctx.polylog(s,a))
    b = 1/a
    return 0.5*(ctx.polylog(s,a) + ctx.polylog(s,b))

@defun
def hurwitz(ctx, s, a=1, derivative=0):
    d = int(derivative)
    s = ctx.convert(s)
    if s == 1:
        return ctx.inf
    abss = abs(s)
    if abss == ctx.inf:
        if ctx.re(s) == ctx.inf:
            if d == 0:
                return ctx.one
            return ctx.zero
        return s*0
    elif ctx.isnan(abss):
        return 1/s
    if ctx.re(s) > 2*ctx.prec and a == 1 and not derivative:
        return ctx.one + ctx.power(2, -s)
    prec = ctx.prec
    try:
        ctx.prec += 10
        v = ctx._hurwitz(s, a, d)
    finally:
        ctx.prec = prec
    return +v

@defun
def _hurwitz(ctx, s, a=1, d=0):
    # We strongly want to special-case rational a
    a, atype = ctx.convert_maybe_rational(a)
    prec = ctx.prec
    # TODO: implement reflection for derivatives
    res = ctx.re(s)
    negs = -s
    try:
        if res < 0 and not d: #and (ar or af):
            # Integer reflection formula
            if ctx.isnpint(s):
                n = int(res)
                if n <= 0:
                    return ctx.bernpoly(1-n, a) / (n-1)
            t = 1-s
            # We now require a to be standardized
            v = 0
            shift = 0
            b = a
            while ctx.re(b) > 1:
                b -= 1
                v -= b**negs
                shift -= 1
            while ctx.re(b) <= 0:
                v += b**negs
                b += 1
                shift += 1
            # Rational reflection formula
            if atype == 'Q' or atype == 'Z':
                try:
                    p, q = a
                except:
                    assert a == int(a)
                    p = int(a)
                    q = 1
                p += shift*q
                assert 1 <= p <= q
                g = ctx.fsum(ctx.cospi(t/2-2*k*b)*ctx._hurwitz(t,(k,q)) \
                    for k in range(1,q+1))
                g *= 2*ctx.gamma(t)/(2*ctx.pi*q)**t
                v += g
                return v
            # General reflection formula
            else:
                C1 = ctx.cospi(t/2)
                C2 = ctx.sinpi(t/2)
                # Clausen functions; could maybe use polylog directly
                if C1: C1 *= ctx.clcos(t, 2*ctx.pi*a)
                if C2: C2 *= ctx.clsin(t, 2*ctx.pi*a)
                v += 2*ctx.gamma(t)/(2*ctx.pi)**t*(C1+C2)
                return v
    except NotImplementedError:
        pass
    a = ctx.convert(a)
    tol = -prec
    # Estimate number of terms for Euler-Maclaurin summation; could be improved
    M1 = 0
    M2 = prec // 3
    N = M2
    lsum = 0
    # This speeds up the recurrence for derivatives
    if ctx.isint(s):
        s = int(s.real)
    s1 = s-1
    while 1:
        # Truncated L-series
        if d:
            l = ctx.fsum(ctx.log(n+a)**d * (n+a)**negs for n in range(M1,M2))
        else:
            l = ctx.fsum((n+a)**negs for n in range(M1,M2))
        lsum += l
        M2a = M2+a
        logM2a = ctx.log(M2a)
        logM2ad = logM2a**d
        logs = [logM2ad]
        logr = 1/logM2a
        rM2a = 1/M2a
        M2as = rM2a**s
        if d:
            tailsum = ctx.gammainc(d+1, s1*logM2a) / s1**(d+1)
        else:
            tailsum = 1/((s-1)*(M2a)**s1)
        tailsum += 0.5 * logM2ad * M2as
        U = [1]
        r = M2as
        fact = 2
        for j in range(1, N+1):
            # TODO: the following could perhaps be tidied a bit
            j2 = 2*j
            if j == 1:
                upds = [1]
            else:
                upds = [j2-2, j2-1]
            for m in upds:
                D = min(m,d+1)
                if m <= d:
                    logs.append(logs[-1] * logr)
                Un = [0]*(D+1)
                for i in xrange(D): Un[i] = (1-m-s)*U[i]
                for i in xrange(1,D+1): Un[i] += (d-(i-1))*U[i-1]
                U = Un
                r *= rM2a
            t = ctx.fdot(U, logs) * r * ctx.bernoulli(j2)/(-fact)
            tailsum += t
            if ctx.mag(t) < tol:
                return (-1)**d * (lsum + tailsum)
            fact *= (j2+1)*(j2+2)
        M1, M2 = M2, M2*2

@defun
def dirichlet(ctx, s, chi=[1], derivative=0):
    s = ctx.convert(s)
    q = len(chi)
    d = int(derivative)
    if d > 2:
        raise NotImplementedError("arbitrary order derivatives")
    prec = ctx.prec
    try:
        ctx.prec += 10
        if s == 1:
            have_pole = True
            for x in chi:
                if x and x != 1:
                    have_pole = False
                    h = +ctx.eps
                    ctx.prec *= 2*(d+1)
                    s += h
            if have_pole:
                return +ctx.inf
        z = ctx.zero
        for p in range(1,q+1):
            if chi[p%q]:
                if d == 1:
                    z += chi[p%q] * (ctx.hurwitz(s, (p,q), 1) - \
                        ctx.hurwitz(s, (p,q))*ctx.log(q))
                else:
                    z += chi[p%q] * ctx.hurwitz(s, (p,q))
        z /= q**s
    finally:
        ctx.prec = prec
    return +z

# Experimental code; could be used elsewhere
def sum_accurately(ctx, terms, check_step=1):
    orig = ctx.prec
    extra = 10
    while 1:
        ctx.prec = orig + extra
        max_term = ctx.ninf
        s = 0
        k = 0
        for term in terms():
            s += term
            if not k % check_step and term:
                abs_term = abs(term)
                abs_sum = abs(s)
                max_term = max(max_term, abs_term)
                if abs_term <= ctx.eps*abs_sum:
                    break
            k += 1
        if abs_sum:
            cancellation = int(max(0,ctx.log(max_term/abs_sum,2)))
        else:
            cancellation = ctx.prec
        if cancellation < extra:
            break
        else:
            extra += cancellation
    return s

@defun_wrapped
def bell(ctx, n, x=1):
    x = ctx.convert(x)
    if not n:
        if ctx.isnan(x):
            return x
        return type(x)(1)
    if ctx.isinf(x) or ctx.isinf(n) or ctx.isnan(x) or ctx.isnan(n):
        return x**n
    if n == 1: return x
    if n == 2: return x*(x+1)
    if x == 0: return ctx.sincpi(n)
    return _polyexp(ctx, n, x, True) / ctx.exp(x)

def _polyexp(ctx, n, x, extra=False):
    def _terms():
        if extra:
            yield ctx.sincpi(n)
        t = x
        k = 1
        while 1:
            yield k**n * t
            k += 1
            t = t*x/k
    return sum_accurately(ctx, _terms, check_step=4)

@defun_wrapped
def polyexp(ctx, s, z):
    if ctx.isinf(z) or ctx.isinf(s) or ctx.isnan(z) or ctx.isnan(s):
        return z**s
    if z == 0: return z*s
    if s == 0: return ctx.expm1(z)
    if s == 1: return ctx.exp(z)*z
    if s == 2: return ctx.exp(z)*z*(z+1)
    return _polyexp(ctx, s, z)

@defun_wrapped
def cyclotomic(ctx, n, z):
    n = int(n)
    assert n >= 0
    p = ctx.one
    if n == 0:
        return p
    if n == 1:
        return z - p
    if n == 2:
        return z + p
    # Use divisor product representation. Unfortunately, this sometimes
    # includes singularities for roots of unity, which we have to cancel out.
    # Matching zeros/poles pairwise, we have (1-z^a)/(1-z^b) ~ a/b + O(z-1).
    moebius = libintmath.moebius
    a_prod = 1
    b_prod = 1
    num_zeros = 0
    num_poles = 0
    for d in range(1,n+1):
        if not n % d:
            w = moebius(n//d)
            # Use powm1 because it is important that we get 0 only
            # if it really is exactly 0
            b = -ctx.powm1(z, d)
            if b:
                p *= b**w
            else:
                if w == 1:
                    a_prod *= d
                    num_zeros += 1
                elif w == -1:
                    b_prod *= d
                    num_poles += 1
    #print n, num_zeros, num_poles
    if num_zeros:
        if num_zeros > num_poles:
            p *= 0
        else:
            p *= a_prod
            p /= b_prod
    return p

@defun
def meijerg(ctx, a_s, b_s, z, r=1, **kwargs):
    an, ap = a_s
    bm, bq = b_s
    n = len(an)
    p = n + len(ap)
    m = len(bm)
    q = m + len(bq)
    a = an+ap
    b = bm+bq
    a = map(ctx.convert, a)
    b = map(ctx.convert, b)
    z = ctx.convert(z)
    #print "m,n,p,q =", m,n,p,q
    #assert p<q or (p==q and (m+n > p)) or (p==q and (m+n) == p and abs(z) < 1)
    def h(*args):
        a = args[:p]
        b = args[p:]
        #assert len(a) == p
        #assert len(b) == q
        terms = []
        for k in range(m):
            bases = [z]
            expts = [b[k]/r]
            gn = [b[j]-b[k] for j in range(m) if j != k]
            gn += [1-a[j]+b[k] for j in range(n)]
            gd = [a[j]-b[k] for j in range(n,p)]
            gd += [1-b[j]+b[k] for j in range(m,q)]
            hn = [1-a[j]+b[k] for j in range(p)]
            hd = [1-b[j]+b[k] for j in range(q) if j != k]
            hz = (-ctx.one)**(p-m-n) * z**(ctx.one/r)
            terms.append((bases, expts, gn, gd, hn, hd, hz))
        return terms
    return ctx.hypercomb(h, a+b, **kwargs)

@defun_wrapped
def appellf1(ctx,a,b1,b2,c,z1,z2,**kwargs):
    # Assume z1 smaller
    # We will use z1 for the outer loop
    if abs(z1) > abs(z2):
        z1, z2 = z2, z1
        b1, b2 = b2, b1
    def ok(x):
        return abs(x) < 0.99
    # Finite cases
    if ctx.isnpint(a):
        pass
    elif ctx.isnpint(b1):
        pass
    elif ctx.isnpint(b2):
        z1, z2, b1, b2 = z2, z1, b2, b1
    else:
        #print z1, z2
        # Note: ok if |z2| > 1, because
        # 2F1 implements analytic continuation
        if not ok(z1):
            u1 = (z1-z2)/(z1-1)
            if not ok(u1):
                raise ValueError("Analytic continuation not implemented")
            #print "Using analytic continuation"
            return (1-z1)**(-b1)*(1-z2)**(c-a-b2)*\
                ctx.appellf1(c-a,b1,c-b1-b2,c,u1,z2,**kwargs)
    #print "inner is", a, b2, c
    one = ctx.one
    s = 0
    t = 1
    k = 0
    while 1:
        h = ctx.hyp2f1(a,b2,c,z2,zeroprec=ctx.prec,**kwargs)
        term = t * h
        if abs(term) < ctx.eps and abs(h) > 10*ctx.eps:
            break
        s += term
        k += 1
        t = (t*a*b1*z1) / (c*k)
        c += one
        a += one
        b1 += one
    return s

@defun_wrapped
def coulombc(ctx, l, eta, _cache={}):
    if (l, eta) in _cache and _cache[l,eta][0] >= ctx.prec:
        return +_cache[l,eta][1]
    G3 = ctx.loggamma(2*l+2)
    G1 = ctx.loggamma(1+l+ctx.j*eta)
    G2 = ctx.loggamma(1+l-ctx.j*eta)
    v = 2**l * ctx.exp((-ctx.pi*eta+G1+G2)/2 - G3)
    if not (ctx.im(l) or ctx.im(eta)):
        v = ctx.re(v)
    _cache[l,eta] = (ctx.prec, v)
    return v

@defun_wrapped
def coulombf(ctx, l, eta, z, w=1, chop=True, **kwargs):
    # Regular Coulomb wave function
    # Note: w can be either 1 or -1; the other may be better in some cases
    # TODO: check that chop=True chops when and only when it should
    #ctx.prec += 10
    def h(l, eta):
        try:
            jw = ctx.j*w
            jwz = ctx.fmul(jw, z, exact=True)
            jwz2 = ctx.fmul(jwz, -2, exact=True)
            C = ctx.coulombc(l, eta)
            T1 = [C, z, ctx.exp(jwz)], [1, l+1, 1], [], [], [1+l+jw*eta], \
                [2*l+2], jwz2
        except ValueError:
            T1 = [0], [-1], [], [], [], [], 0
        return (T1,)
    v = ctx.hypercomb(h, [l,eta], **kwargs)
    if chop and (not ctx.im(l)) and (not ctx.im(eta)) and (not ctx.im(z)) and \
        (ctx.re(z) >= 0):
        v = ctx.re(v)
    return v

@defun_wrapped
def _coulomb_chi(ctx, l, eta, _cache={}):
    if (l, eta) in _cache and _cache[l,eta][0] >= ctx.prec:
        return _cache[l,eta][1]
    def terms():
        l2 = -l-1
        jeta = ctx.j*eta
        return [ctx.loggamma(1+l+jeta) * (-0.5j),
            ctx.loggamma(1+l-jeta) * (0.5j),
            ctx.loggamma(1+l2+jeta) * (0.5j),
            ctx.loggamma(1+l2-jeta) * (-0.5j),
            -(l+0.5)*ctx.pi]
    v = sum_accurately(ctx, terms, 1)
    _cache[l,eta] = (ctx.prec, v)
    return v

@defun_wrapped
def coulombg(ctx, l, eta, z, w=1, chop=True, **kwargs):
    # Irregular Coulomb wave function
    # Note: w can be either 1 or -1; the other may be better in some cases
    # TODO: check that chop=True chops when and only when it should
    if not l.imag:
        l = l.real  # XXX: for isint
    def h(l, eta):
        # Force perturbation for integers and half-integers
        if ctx.isint(l*2):
            T1 = [0], [-1], [], [], [], [], 0
            return (T1,)
        l2 = -l-1
        try:
            chi = ctx._coulomb_chi(l, eta)
            jw = ctx.j*w
            s = ctx.sin(chi); c = ctx.cos(chi)
            C1 = ctx.coulombc(l,eta)
            C2 = ctx.coulombc(l2,eta)
            u = ctx.exp(jw*z)
            x = -2*jw*z
            T1 = [s, C1, z, u, c], [-1, 1, l+1, 1, 1], [], [], \
                [1+l+jw*eta], [2*l+2], x
            T2 = [-s, C2, z, u],   [-1, 1, l2+1, 1],    [], [], \
                [1+l2+jw*eta], [2*l2+2], x
            return T1, T2
        except ValueError:
            T1 = [0], [-1], [], [], [], [], 0
            return (T1,)
    v = ctx.hypercomb(h, [l,eta], **kwargs)
    if chop and (not ctx.im(l)) and (not ctx.im(eta)) and (not ctx.im(z)) and \
        (ctx.re(z) >= 0):
        v = ctx.re(v)
    return v

@defun
def calculate_nome(ctx, k):
    k = ctx.convert(k)
    if abs(k) > ctx.one:             # range error
        raise ValueError
    if k == ctx.zero:
        return ctx.zero
    elif k == ctx.one:
        return ctx.one
    else:
        kprimesquared = ctx.one - k**2
        kprime = ctx.sqrt(kprimesquared)
        top = ctx.ellipk(kprimesquared)
        bottom = ctx.ellipk(k**2)
        argument = -ctx.pi*top/bottom
        nome = ctx.exp(argument)
        return nome

@defun
def _jacobi_theta2(ctx, z, q):
    extra1 = 10
    extra2 = 20
    # the loops below break when the fixed precision quantities
    # a and b go to zero;
    # right shifting small negative numbers by wp one obtains -1, not zero,
    # so the condition a**2 + b**2 > MIN is used to break the loops.
    MIN = 2
    if z == ctx.zero:
        if (not q.imag):
            wp = ctx.prec + extra1
            x = ctx.to_fixed(q, wp)
            x2 = (x*x) >> wp
            a = b = x2
            s = x2
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                s += a
            s = (1 << (wp+1)) + (s << 1)
            s = ctx.ldexp(s, -wp)
        else:
            wp = ctx.prec + extra1
            xre = ctx.to_fixed(q.real, wp)
            xim = ctx.to_fixed(q.imag, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp-1)
            are = bre = x2re
            aim = bim = x2im
            sre = (1<<wp) + are
            sim = aim
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                sre += are
                sim += aim
            sre = (sre << 1)
            sim = (sim << 1)
            sre = ctx.ldexp(sre, -wp)
            sim = ctx.ldexp(sim, -wp)
            s = ctx.mpc(sre, sim)
    else:
        if (not q.imag) and (not z.imag):
            wp = ctx.prec + extra1
            x = ctx.to_fixed(q, wp)
            x2 = (x*x) >> wp
            a = b = x2
            c1, s1 = ctx.cos_sin(z, prec=wp)
            cn = c1 = ctx.to_fixed(c1, wp)
            sn = s1 = ctx.to_fixed(s1, wp)
            c2 = (c1*c1 - s1*s1) >> wp
            s2 = (c1 * s1) >> (wp - 1)
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            s = c1 + ((a * cn) >> wp)
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
                s += (a * cn) >> wp
            s = (s << 1)
            s = ctx.ldexp(s, -wp)
            s *= ctx.nthroot(q, 4)
            return s
        # case z real, q complex
        elif not z.imag:
            wp = ctx.prec + extra2
            xre = ctx.to_fixed(q.real, wp)
            xim = ctx.to_fixed(q.imag, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = x2re
            aim = bim = x2im
            c1, s1 = ctx.cos_sin(z, prec=wp)
            cn = c1 = ctx.to_fixed(c1, wp)
            sn = s1 = ctx.to_fixed(s1, wp)
            c2 = (c1*c1 - s1*s1) >> wp
            s2 = (c1 * s1) >> (wp - 1)
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            sre = c1 + ((are * cn) >> wp)
            sim = ((aim * cn) >> wp)
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
                sre += ((are * cn) >> wp)
                sim += ((aim * cn) >> wp)
            sre = (sre << 1)
            sim = (sim << 1)
            sre = ctx.ldexp(sre, -wp)
            sim = ctx.ldexp(sim, -wp)
            s = ctx.mpc(sre, sim)
        #case z complex, q real
        elif not q.imag:
            wp = ctx.prec + extra2
            x = ctx.to_fixed(q, wp)
            x2 = (x*x) >> wp
            a = b = x2
            prec0 = ctx.prec
            ctx.prec = wp
            c1 = ctx.cos(z)
            s1 = ctx.sin(z)
            ctx.prec = prec0
            cnre = c1re = ctx.to_fixed(c1.real, wp)
            cnim = c1im = ctx.to_fixed(c1.imag, wp)
            snre = s1re = ctx.to_fixed(s1.real, wp)
            snim = s1im = ctx.to_fixed(s1.imag, wp)
            #c2 = (c1*c1 - s1*s1) >> wp
            c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
            c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
            #s2 = (c1 * s1) >> (wp - 1)
            s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
            s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
            #cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            sre = c1re + ((a * cnre) >> wp)
            sim = c1im + ((a * cnim) >> wp)
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
                t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
                t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
                t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                sre += ((a * cnre) >> wp)
                sim += ((a * cnim) >> wp)
            sre = (sre << 1)
            sim = (sim << 1)
            sre = ctx.ldexp(sre, -wp)
            sim = ctx.ldexp(sim, -wp)
            s = ctx.mpc(sre, sim)
        # case z and q complex
        else:
            wp = ctx.prec + extra2
            xre = ctx.to_fixed(q.real, wp)
            xim = ctx.to_fixed(q.imag, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = x2re
            aim = bim = x2im
            prec0 = ctx.prec
            ctx.prec = wp
            # cos(z), siz(z) with z complex
            c1 = ctx.cos(z)
            s1 = ctx.sin(z)
            ctx.prec = prec0
            cnre = c1re = ctx.to_fixed(c1.real, wp)
            cnim = c1im = ctx.to_fixed(c1.imag, wp)
            snre = s1re = ctx.to_fixed(s1.real, wp)
            snim = s1im = ctx.to_fixed(s1.imag, wp)
            c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
            c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
            s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
            s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            n = 1
            termre = c1re
            termim = c1im
            sre = c1re + ((are * cnre - aim * cnim) >> wp)
            sim = c1im + ((are * cnim + aim * cnre) >> wp)
            n = 3
            termre = ((are * cnre - aim * cnim) >> wp)
            termim = ((are * cnim + aim * cnre) >> wp)
            sre = c1re + ((are * cnre - aim * cnim) >> wp)
            sim = c1im + ((are * cnim + aim * cnre) >> wp)
            n = 5
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                #cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
                t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
                t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
                t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
                t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                termre = ((are * cnre - aim * cnim) >> wp)
                termim = ((aim * cnre + are * cnim) >> wp)
                sre += ((are * cnre - aim * cnim) >> wp)
                sim += ((aim * cnre + are * cnim) >> wp)
                n += 2
            sre = (sre << 1)
            sim = (sim << 1)
            sre = ctx.ldexp(sre, -wp)
            sim = ctx.ldexp(sim, -wp)
            s = ctx.mpc(sre, sim)
    s *= ctx.nthroot(q, 4)
    return s

@defun
def _djacobi_theta2(ctx, z, q, nd):
    MIN = 2
    extra1 = 10
    extra2 = 20
    if (not q.imag) and (not z.imag):
        wp = ctx.prec + extra1
        x = ctx.to_fixed(q, wp)
        x2 = (x*x) >> wp
        a = b = x2
        c1, s1 = ctx.cos_sin(z, prec=wp)
        cn = c1 = ctx.to_fixed(c1, wp)
        sn = s1 = ctx.to_fixed(s1, wp)
        c2 = (c1*c1 - s1*s1) >> wp
        s2 = (c1 * s1) >> (wp - 1)
        cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        if (nd&1):
            s = s1 + ((a * sn * 3**nd) >> wp)
        else:
            s = c1 + ((a * cn * 3**nd) >> wp)
        n = 2
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
            if nd&1:
                s += (a * sn * (2*n+1)**nd) >> wp
            else:
                s += (a * cn * (2*n+1)**nd) >> wp
            n += 1
        s = -(s << 1)
        s = ctx.ldexp(s, -wp)
        # case z real, q complex
    elif not z.imag:
        wp = ctx.prec + extra2
        xre = ctx.to_fixed(q.real, wp)
        xim = ctx.to_fixed(q.imag, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = x2re
        aim = bim = x2im
        c1, s1 = ctx.cos_sin(z, prec=wp)
        cn = c1 = ctx.to_fixed(c1, wp)
        sn = s1 = ctx.to_fixed(s1, wp)
        c2 = (c1*c1 - s1*s1) >> wp
        s2 = (c1 * s1) >> (wp - 1)
        cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        if (nd&1):
            sre = s1 + ((are * sn * 3**nd) >> wp)
            sim = ((aim * sn * 3**nd) >> wp)
        else:
            sre = c1 + ((are * cn * 3**nd) >> wp)
            sim = ((aim * cn * 3**nd) >> wp)
        n = 5
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp

            if (nd&1):
                sre += ((are * sn * n**nd) >> wp)
                sim += ((aim * sn * n**nd) >> wp)
            else:
                sre += ((are * cn * n**nd) >> wp)
                sim += ((aim * cn * n**nd) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)
    #case z complex, q real
    elif not q.imag:
        wp = ctx.prec + extra2
        x = ctx.to_fixed(q, wp)
        x2 = (x*x) >> wp
        a = b = x2
        prec0 = ctx.prec
        ctx.prec = wp
        c1 = ctx.cos(z)
        s1 = ctx.sin(z)
        ctx.prec = prec0
        cnre = c1re = ctx.to_fixed(c1.real, wp)
        cnim = c1im = ctx.to_fixed(c1.imag, wp)
        snre = s1re = ctx.to_fixed(s1.real, wp)
        snim = s1im = ctx.to_fixed(s1.imag, wp)
        #c2 = (c1*c1 - s1*s1) >> wp
        c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
        c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
        #s2 = (c1 * s1) >> (wp - 1)
        s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
        s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
        #cn, sn = (cn*c2 - sn*s2) >> wp, (sn*c2 + cn*s2) >> wp
        t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
        t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
        t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
        t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
        cnre = t1
        cnim = t2
        snre = t3
        snim = t4
        if (nd&1):
            sre = s1re + ((a * snre * 3**nd) >> wp)
            sim = s1im + ((a * snim * 3**nd) >> wp)
        else:
            sre = c1re + ((a * cnre * 3**nd) >> wp)
            sim = c1im + ((a * cnim * 3**nd) >> wp)
        n = 5
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if (nd&1):
                sre += ((a * snre * n**nd) >> wp)
                sim += ((a * snim * n**nd) >> wp)
            else:
                sre += ((a * cnre * n**nd) >> wp)
                sim += ((a * cnim * n**nd) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)
    # case z and q complex
    else:
        wp = ctx.prec + extra2
        xre = ctx.to_fixed(q.real, wp)
        xim = ctx.to_fixed(q.imag, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = x2re
        aim = bim = x2im
        prec0 = ctx.prec
        ctx.prec = wp
        # cos(2*z), sin(2*z) with z complex
        c1 = ctx.cos(z)
        s1 = ctx.sin(z)
        ctx.prec = prec0
        cnre = c1re = ctx.to_fixed(c1.real, wp)
        cnim = c1im = ctx.to_fixed(c1.imag, wp)
        snre = s1re = ctx.to_fixed(s1.real, wp)
        snim = s1im = ctx.to_fixed(s1.imag, wp)
        c2re = (c1re*c1re - c1im*c1im - s1re*s1re + s1im*s1im) >> wp
        c2im = (c1re*c1im - s1re*s1im) >> (wp - 1)
        s2re = (c1re*s1re - c1im*s1im) >> (wp - 1)
        s2im = (c1re*s1im + c1im*s1re) >> (wp - 1)
        t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
        t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
        t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
        t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
        cnre = t1
        cnim = t2
        snre = t3
        snim = t4
        if (nd&1):
            sre = s1re + (((are * snre - aim * snim) * 3**nd) >> wp)
            sim = s1im + (((are * snim + aim * snre)* 3**nd) >> wp)
        else:
            sre = c1re + (((are * cnre - aim * cnim) * 3**nd) >> wp)
            sim = c1im + (((are * cnim + aim * cnre)* 3**nd) >> wp)
        n = 5
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            #cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            t1 = (cnre*c2re - cnim*c2im - snre*s2re + snim*s2im) >> wp
            t2 = (cnre*c2im + cnim*c2re - snre*s2im - snim*s2re) >> wp
            t3 = (snre*c2re - snim*c2im + cnre*s2re - cnim*s2im) >> wp
            t4 = (snre*c2im + snim*c2re + cnre*s2im + cnim*s2re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if (nd&1):
                sre += (((are * snre - aim * snim) * n**nd) >> wp)
                sim += (((aim * snre + are * snim) * n**nd) >> wp)
            else:
                sre += (((are * cnre - aim * cnim) * n**nd) >> wp)
                sim += (((aim * cnre + are * cnim) * n**nd) >> wp)
            n += 2
        sre = -(sre << 1)
        sim = -(sim << 1)
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)
    s *= ctx.nthroot(q, 4)
    if (nd&1):
        return (-1)**(nd//2) * s
    else:
        return (-1)**(1 + nd//2) * s

@defun
def _jacobi_theta3(ctx, z, q):
    extra1 = 10
    extra2 = 20
    MIN = 2
    if z == ctx.zero:
        if not q.imag:
            wp = ctx.prec + extra1
            x = ctx.to_fixed(q, wp)
            s = x
            a = b = x
            x2 = (x*x) >> wp
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                s += a
            s = (1 << wp) + (s << 1)
            s = ctx.ldexp(s, -wp)
            return s
        else:
            wp = ctx.prec + extra1
            xre = ctx.to_fixed(q.real, wp)
            xim = ctx.to_fixed(q.imag, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            sre = are = bre = xre
            sim = aim = bim = xim
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                sre += are
                sim += aim
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = ctx.ldexp(sre, -wp)
            sim = ctx.ldexp(sim, -wp)
            s = ctx.mpc(sre, sim)
            return s
    else:
        if (not q.imag) and (not z.imag):
            s = 0
            wp = ctx.prec + extra1
            x = ctx.to_fixed(q, wp)
            a = b = x
            x2 = (x*x) >> wp
            c1, s1 = ctx.cos_sin(z*2, prec=wp)
            c1 = ctx.to_fixed(c1, wp)
            s1 = ctx.to_fixed(s1, wp)
            cn = c1
            sn = s1
            s += (a * cn) >> wp
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
                s += (a * cn) >> wp
            s = (1 << wp) + (s << 1)
            s = ctx.ldexp(s, -wp)
            return s
        # case z real, q complex
        elif not z.imag:
            wp = ctx.prec + extra2
            xre = ctx.to_fixed(q.real, wp)
            xim = ctx.to_fixed(q.imag, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = xre
            aim = bim = xim
            c1, s1 = ctx.cos_sin(z*2, prec=wp)
            c1 = ctx.to_fixed(c1, wp)
            s1 = ctx.to_fixed(s1, wp)
            cn = c1
            sn = s1
            sre = (are * cn) >> wp
            sim = (aim * cn) >> wp
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
                sre += (are * cn) >> wp
                sim += (aim * cn) >> wp
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = ctx.ldexp(sre, -wp)
            sim = ctx.ldexp(sim, -wp)
            s = ctx.mpc(sre, sim)
            return s
        #case z complex, q real
        elif not q.imag:
            wp = ctx.prec + extra2
            x = ctx.to_fixed(q, wp)
            a = b = x
            x2 = (x*x) >> wp
            prec0 = ctx.prec
            ctx.prec = wp
            c1 = ctx.cos(2*z)
            s1 = ctx.sin(2*z)
            ctx.prec = prec0
            cnre = c1re = ctx.to_fixed(c1.real, wp)
            cnim = c1im = ctx.to_fixed(c1.imag, wp)
            snre = s1re = ctx.to_fixed(s1.real, wp)
            snim = s1im = ctx.to_fixed(s1.imag, wp)
            sre = (a * cnre) >> wp
            sim = (a * cnim) >> wp
            while abs(a) > MIN:
                b = (b*x2) >> wp
                a = (a*b) >> wp
                t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
                t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
                t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
                t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                sre += (a * cnre) >> wp
                sim += (a * cnim) >> wp
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = ctx.ldexp(sre, -wp)
            sim = ctx.ldexp(sim, -wp)
            s = ctx.mpc(sre, sim)
            return s
        # case z and q complex
        else:
            wp = ctx.prec + extra2
            xre = ctx.to_fixed(q.real, wp)
            xim = ctx.to_fixed(q.imag, wp)
            x2re = (xre*xre - xim*xim) >> wp
            x2im = (xre*xim) >> (wp - 1)
            are = bre = xre
            aim = bim = xim
            prec0 = ctx.prec
            ctx.prec = wp
            # cos(2*z), sin(2*z) with z complex
            c1 = ctx.cos(2*z)
            s1 = ctx.sin(2*z)
            ctx.prec = prec0
            cnre = c1re = ctx.to_fixed(c1.real, wp)
            cnim = c1im = ctx.to_fixed(c1.imag, wp)
            snre = s1re = ctx.to_fixed(s1.real, wp)
            snim = s1im = ctx.to_fixed(s1.imag, wp)
            sre = (are * cnre - aim * cnim) >> wp
            sim = (aim * cnre + are * cnim) >> wp
            while are**2 + aim**2 > MIN:
                bre, bim = (bre * x2re - bim * x2im) >> wp, \
                           (bre * x2im + bim * x2re) >> wp
                are, aim = (are * bre - aim * bim) >> wp,   \
                           (are * bim + aim * bre) >> wp
                t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
                t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
                t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
                t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
                cnre = t1
                cnim = t2
                snre = t3
                snim = t4
                sre += (are * cnre - aim * cnim) >> wp
                sim += (aim * cnre + are * cnim) >> wp
            sre = (1 << wp) + (sre << 1)
            sim = (sim << 1)
            sre = ctx.ldexp(sre, -wp)
            sim = ctx.ldexp(sim, -wp)
            s = ctx.mpc(sre, sim)
            return s

@defun
def _djacobi_theta3(ctx, z, q, nd):
    """nd=1,2,3 order of the derivative with respect to z"""
    MIN = 2
    extra1 = 10
    extra2 = 20
    if (not q.imag) and (not z.imag):
        s = 0
        wp = ctx.prec + extra1
        x = ctx.to_fixed(q, wp)
        a = b = x
        x2 = (x*x) >> wp
        c1, s1 = ctx.cos_sin(z*2, prec=wp)
        c1 = ctx.to_fixed(c1, wp)
        s1 = ctx.to_fixed(s1, wp)
        cn = c1
        sn = s1
        if (nd&1):
            s += (a * sn) >> wp
        else:
            s += (a * cn) >> wp
        n = 2
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            if nd&1:
                s += (a * sn * n**nd) >> wp
            else:
                s += (a * cn * n**nd) >> wp
            n += 1
        s = -(s << (nd+1))
        s = ctx.ldexp(s, -wp)
    # case z real, q complex
    elif not z.imag:
        wp = ctx.prec + extra2
        xre = ctx.to_fixed(q.real, wp)
        xim = ctx.to_fixed(q.imag, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = xre
        aim = bim = xim
        c1, s1 = ctx.cos_sin(z*2, prec=wp)
        c1 = ctx.to_fixed(c1, wp)
        s1 = ctx.to_fixed(s1, wp)
        cn = c1
        sn = s1
        if (nd&1):
            sre = (are * sn) >> wp
            sim = (aim * sn) >> wp
        else:
            sre = (are * cn) >> wp
            sim = (aim * cn) >> wp
        n = 2
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            cn, sn = (cn*c1 - sn*s1) >> wp, (sn*c1 + cn*s1) >> wp
            if nd&1:
                sre += (are * sn * n**nd) >> wp
                sim += (aim * sn * n**nd) >> wp
            else:
                sre += (are * cn * n**nd) >> wp
                sim += (aim * cn * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)
    #case z complex, q real
    elif not q.imag:
        wp = ctx.prec + extra2
        x = ctx.to_fixed(q, wp)
        a = b = x
        x2 = (x*x) >> wp
        prec0 = ctx.prec
        ctx.prec = wp
        c1 = ctx.cos(2*z)
        s1 = ctx.sin(2*z)
        ctx.prec = prec0
        cnre = c1re = ctx.to_fixed(c1.real, wp)
        cnim = c1im = ctx.to_fixed(c1.imag, wp)
        snre = s1re = ctx.to_fixed(s1.real, wp)
        snim = s1im = ctx.to_fixed(s1.imag, wp)
        if (nd&1):
            sre = (a * snre) >> wp
            sim = (a * snim) >> wp
        else:
            sre = (a * cnre) >> wp
            sim = (a * cnim) >> wp
        n = 2
        while abs(a) > MIN:
            b = (b*x2) >> wp
            a = (a*b) >> wp
            t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
            t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
            t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
            t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if (nd&1):
                sre += (a * snre * n**nd) >> wp
                sim += (a * snim * n**nd) >> wp
            else:
                sre += (a * cnre * n**nd) >> wp
                sim += (a * cnim * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)
    # case z and q complex
    else:
        wp = ctx.prec + extra2
        xre = ctx.to_fixed(q.real, wp)
        xim = ctx.to_fixed(q.imag, wp)
        x2re = (xre*xre - xim*xim) >> wp
        x2im = (xre*xim) >> (wp - 1)
        are = bre = xre
        aim = bim = xim
        prec0 = ctx.prec
        ctx.prec = wp
        # cos(2*z), sin(2*z) with z complex
        c1 = ctx.cos(2*z)
        s1 = ctx.sin(2*z)
        ctx.prec = prec0
        cnre = c1re = ctx.to_fixed(c1.real, wp)
        cnim = c1im = ctx.to_fixed(c1.imag, wp)
        snre = s1re = ctx.to_fixed(s1.real, wp)
        snim = s1im = ctx.to_fixed(s1.imag, wp)
        if (nd&1):
            sre = (are * snre - aim * snim) >> wp
            sim = (aim * snre + are * snim) >> wp
        else:
            sre = (are * cnre - aim * cnim) >> wp
            sim = (aim * cnre + are * cnim) >> wp
        n = 2
        while are**2 + aim**2 > MIN:
            bre, bim = (bre * x2re - bim * x2im) >> wp, \
                       (bre * x2im + bim * x2re) >> wp
            are, aim = (are * bre - aim * bim) >> wp,   \
                       (are * bim + aim * bre) >> wp
            t1 = (cnre*c1re - cnim*c1im - snre*s1re + snim*s1im) >> wp
            t2 = (cnre*c1im + cnim*c1re - snre*s1im - snim*s1re) >> wp
            t3 = (snre*c1re - snim*c1im + cnre*s1re - cnim*s1im) >> wp
            t4 = (snre*c1im + snim*c1re + cnre*s1im + cnim*s1re) >> wp
            cnre = t1
            cnim = t2
            snre = t3
            snim = t4
            if(nd&1):
                sre += ((are * snre - aim * snim) * n**nd) >> wp
                sim += ((aim * snre + are * snim) * n**nd) >> wp
            else:
                sre += ((are * cnre - aim * cnim) * n**nd) >> wp
                sim += ((aim * cnre + are * cnim) * n**nd) >> wp
            n += 1
        sre = -(sre << (nd+1))
        sim = -(sim << (nd+1))
        sre = ctx.ldexp(sre, -wp)
        sim = ctx.ldexp(sim, -wp)
        s = ctx.mpc(sre, sim)
    if (nd&1):
        return (-1)**(nd//2) * s
    else:
        return (-1)**(1 + nd//2) * s

@defun
def _jacobi_theta2a(ctx, z, q):
    """
    case z.imag != 0
    theta(2, z, q) =
    q**1/4 * Sum(q**(n*n + n) * exp(j*(2*n + 1)*z), n=-inf, inf)
    max term for minimum (2*n+1)*log(q).real - 2* z.imag
    n0 = int(z.imag/log(q).real - 1/2)
    theta(2, z, q) =
    q**1/4 * Sum(q**(n*n + n) * exp(j*(2*n + 1)*z), n=n0, inf) +
    q**1/4 * Sum(q**(n*n + n) * exp(j*(2*n + 1)*z), n, n0-1, -inf)
    """
    n = n0 = int(z.imag/ctx.log(q).real - 1/2)
    e2 = ctx.exp(2*ctx.j*z)
    e = e0 = ctx.exp(ctx.j*(2*n+1)*z)
    a = q**(n*n + n)
    # leading term
    term = a * e
    s = term
    eps1 = ctx.eps*abs(term)
    while 1:
        n += 1
        e = e * e2
        term = q**(n*n + n) * e
        if abs(term) < eps1:
            break
        s += term
    e = e0
    e2 = ctx.exp(-2j*z)
    n = n0
    while 1:
        n -= 1
        e = e * e2
        term = q**(n*n + n) * e
        if abs(term) < eps1:
            break
        s += term
    s = s * ctx.nthroot(q, 4)
    return s

@defun
def _jacobi_theta3a(ctx, z, q):
    """
    case z.imag != 0
    theta3(z, q) = Sum(q**(n*n) * exp(j*2*n*z), n, -inf, inf)
    max term for n*abs(log(q).real) + z.imag ~= 0
    n0 = int(- z.imag/abs(log(q).real))
    """
    n = n0 = int(- z.imag/abs(ctx.log(q).real))
    e2 = ctx.exp(2*ctx.j*z)
    e = e0 = ctx.exp(ctx.j*2*n*z)
    s = term = q**(n*n) * e
    eps1 = ctx.eps*abs(term)
    while 1:
        n += 1
        e = e * e2
        term = q**(n*n) * e
        if abs(term) < eps1:
            break
        s += term
    e = e0
    e2 = ctx.exp(-2*ctx.j*z)
    n = n0
    while 1:
        n -= 1
        e = e * e2
        term = q**(n*n) * e
        if abs(term) < eps1:
            break
        s += term
    return s

@defun
def _djacobi_theta2a(ctx, z, q, nd):
    """
    case z.imag != 0
    dtheta(2, z, q, nd) =
    j* q**1/4 * Sum(q**(n*n + n) * (2*n+1)*exp(j*(2*n + 1)*z), n=-inf, inf)
    max term for (2*n0+1)*log(q).real - 2* z.imag ~= 0
    n0 = int(z.imag/log(q).real - 1/2)
    """
    n = n0 = int(z.imag/ctx.log(q).real - 1/2)
    e2 = ctx.exp(2*ctx.j*z)
    e = e0 = ctx.exp(ctx.j*(2*n + 1)*z)
    a = q**(n*n + n)
    # leading term
    term = (2*n+1)**nd * a * e
    s = term
    eps1 = ctx.eps*abs(term)
    while 1:
        n += 1
        e = e * e2
        term = (2*n+1)**nd * q**(n*n + n) * e
        if abs(term) < eps1:
            break
        s += term
    e = e0
    e2 = ctx.exp(-2*ctx.j*z)
    n = n0
    while 1:
        n -= 1
        e = e * e2
        term = (2*n+1)**nd * q**(n*n + n) * e
        if abs(term) < eps1:
            break
        s += term
    return ctx.j**nd * s * ctx.nthroot(q, 4)

@defun
def _djacobi_theta3a(ctx, z, q, nd):
    """
    case z.imag != 0
    djtheta3(z, q, nd) = (2*j)**nd *
      Sum(q**(n*n) * n**nd * exp(j*2*n*z), n, -inf, inf)
    max term for minimum n*abs(log(q).real) + z.imag
    """
    n = n0 = int(-z.imag/abs(ctx.log(q).real))
    e2 = ctx.exp(2*ctx.j*z)
    e = e0 = ctx.exp(ctx.j*2*n*z)
    a = q**(n*n) * e
    s = term = n**nd * a
    if n != 0:
        eps1 = ctx.eps*abs(term)
    else:
        eps1 = ctx.eps*abs(a)
    while 1:
        n += 1
        e = e * e2
        a = q**(n*n) * e
        term = n**nd * a
        if n != 0:
            aterm = abs(term)
        else:
            aterm = abs(a)
        if aterm < eps1:
            break
        s += term
    e = e0
    e2 = ctx.exp(-2*ctx.j*z)
    n = n0
    while 1:
        n -= 1
        e = e * e2
        a = q**(n*n) * e
        term = n**nd * a
        if n != 0:
            aterm = abs(term)
        else:
            aterm = abs(a)
        if aterm < eps1:
            break
        s += term
    return (2*ctx.j)**nd * s

@defun
def jtheta(ctx, n, z, q, derivative=0):
    if derivative:
        return ctx._djtheta(n, z, q, derivative)

    z = ctx.convert(z)
    q = ctx.convert(q)

    # Implementation note
    # If z.imag is close to zero, _jacobi_theta2 and _jacobi_theta3
    # are used,
    # which compute the series starting from n=0 using fixed precision
    # numbers;
    # otherwise  _jacobi_theta2a and _jacobi_theta3a are used, which compute
    # the series starting from n=n0, which is the largest term.

    # TODO: write _jacobi_theta2a and _jacobi_theta3a using fixed-point

    if abs(q) > ctx.THETA_Q_LIM:
        raise ValueError('abs(q) > THETA_Q_LIM = %f' % ctx.THETA_Q_LIM)

    extra = 10
    if z:
        M = ctx.mag(z)
        if M > 5 or (n == 1 and M < -5):
            extra += 2*abs(M)
    cz = 0.5
    extra2 = 50
    prec0 = ctx.prec
    try:
        ctx.prec += extra
        if n == 1:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(ctx.log(q).real):
                    ctx.dps += extra2
                    res = ctx._jacobi_theta2(z - ctx.pi/2, q)
                else:
                    ctx.dps += 10
                    res = ctx._jacobi_theta2a(z - ctx.pi/2, q)
            else:
                res = ctx._jacobi_theta2(z - ctx.pi/2, q)
        elif n == 2:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(ctx.log(q).real):
                    ctx.dps += extra2
                    res = ctx._jacobi_theta2(z, q)
                else:
                    ctx.dps += 10
                    res = ctx._jacobi_theta2a(z, q)
            else:
                res = ctx._jacobi_theta2(z, q)
        elif n == 3:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(ctx.log(q).real):
                    ctx.dps += extra2
                    res = ctx._jacobi_theta3(z, q)
                else:
                    ctx.dps += 10
                    res = ctx._jacobi_theta3a(z, q)
            else:
                res = ctx._jacobi_theta3(z, q)
        elif n == 4:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(ctx.log(q).real):
                    ctx.dps += extra2
                    res = ctx._jacobi_theta3(z, -q)
                else:
                    ctx.dps += 10
                    res = ctx._jacobi_theta3a(z, -q)
            else:
                res = ctx._jacobi_theta3(z, -q)
        else:
            raise ValueError
    finally:
        ctx.prec = prec0
    return res

@defun
def _djtheta(ctx, n, z, q, derivative=1):
    z = ctx.convert(z)
    q = ctx.convert(q)
    nd = int(derivative)

    if abs(q) > ctx.THETA_Q_LIM:
        raise ValueError('abs(q) > THETA_Q_LIM = %f' % ctx.THETA_Q_LIM)
    extra = 10 + ctx.prec * nd // 10
    if z:
        M = ctx.mag(z)
        if M > 5 or (n != 1 and M < -5):
            extra += 2*abs(M)
    cz = 0.5
    extra2 = 50
    prec0 = ctx.prec
    try:
        ctx.prec += extra
        if n == 1:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(ctx.log(q).real):
                    ctx.dps += extra2
                    res = ctx._djacobi_theta2(z - ctx.pi/2, q, nd)
                else:
                    ctx.dps += 10
                    res = ctx._djacobi_theta2a(z - ctx.pi/2, q, nd)
            else:
                res = ctx._djacobi_theta2(z - ctx.pi/2, q, nd)
        elif n == 2:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(ctx.log(q).real):
                    ctx.dps += extra2
                    res = ctx._djacobi_theta2(z, q, nd)
                else:
                    ctx.dps += 10
                    res = ctx._djacobi_theta2a(z, q, nd)
            else:
                res = ctx._djacobi_theta2(z, q, nd)
        elif n == 3:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(ctx.log(q).real):
                    ctx.dps += extra2
                    res = ctx._djacobi_theta3(z, q, nd)
                else:
                    ctx.dps += 10
                    res = ctx._djacobi_theta3a(z, q, nd)
            else:
                res = ctx._djacobi_theta3(z, q, nd)
        elif n == 4:
            if abs(z.imag) != 0:
                if abs(z.imag) < cz * abs(ctx.log(q).real):
                    ctx.dps += extra2
                    res = ctx._djacobi_theta3(z, -q, nd)
                else:
                    ctx.dps += 10
                    res = ctx._djacobi_theta3a(z, -q, nd)
            else:
                res = ctx._djacobi_theta3(z, -q, nd)
        else:
            raise ValueError
    finally:
        ctx.prec = prec0
    return +res

@defun
def jsn(ctx, u, m):
    if abs(m) < ctx.eps:
        return ctx.sin(u)
    elif m == ctx.one:
        return ctx.tanh(u)
    else:
        extra = 10
    try:
        ctx.prec += extra
        q = ctx.calculate_nome(ctx.sqrt(m))
        v3 = ctx.jtheta(3, 0, q)
        v2 = ctx.jtheta(2, 0, q)        # mathworld says v4
        arg1 = u / (v3*v3)
        v1 = ctx.jtheta(1, arg1, q)
        v4 = ctx.jtheta(4, arg1, q)
        sn = (v3/v2)*(v1/v4)
    finally:
        ctx.prec -= extra
    return sn

@defun
def jcn(ctx, u, m):
    if abs(m) < ctx.eps:
        return ctx.cos(u)
    elif m == ctx.one:
        return ctx.sech(u)
    else:
        extra = 10
    try:
        ctx.prec += extra
        q = ctx.calculate_nome(ctx.sqrt(m))
        v3 = ctx.jtheta(3, 0, q)
        v2 = ctx.jtheta(2, 0, q)
        v04 = ctx.jtheta(4, 0, q)
        arg1 = u / (v3*v3)
        v1 = ctx.jtheta(2, arg1, q)
        v4 = ctx.jtheta(4, arg1, q)
        cn = (v04/v2)*(v1/v4)
    finally:
        ctx.prec -= extra
    return +cn

@defun
def jdn(ctx, u, m):
    if m == ctx.zero:
        return ctx.one
    elif m == ctx.one:
        return ctx.sech(u)
    else:
        extra = 10
    try:
        ctx.prec += extra
        q = ctx.calculate_nome(ctx.sqrt(m))
        v3 = ctx.jtheta(3, 0, q)
        v2 = ctx.jtheta(2, 0, q)
        v04 = ctx.jtheta(4, 0, q)
        arg1 = u / (v3*v3)
        v1 = ctx.jtheta(3, arg1, q)
        v4 = ctx.jtheta(4, arg1, q)
        cn = (v04/v3)*(v1/v4)
    finally:
        ctx.prec -= extra
    return +cn

if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    try:
        import psyco; psyco.full()
    except ImportError:
        pass
    import sys
    from timeit import default_timer as clock
    filter = []
    for i, arg in enumerate(sys.argv):
        if 'functions.py' in arg:
            filter = sys.argv[i+1:]
            break
    import doctest
    globs = globals().copy()
    for obj in globs: #sorted(globs.keys()):
        if filter:
            if not sum([pat in obj for pat in filter]):
                continue
        print obj,
        t1 = clock()
        doctest.run_docstring_examples(globs[obj], {})
        t2 = clock()
        print round(t2-t1, 3)


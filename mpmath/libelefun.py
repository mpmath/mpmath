"""
This module implements computation of elementary transcendental
functions (powers, logarithms, trigonometric and hyperbolic
functions, inverse trigonometric and hyperbolic) for real
floating-point numbers.

For complex and interval implementations of the same functions,
see libmpc and libmpi.

"""


from libmpf import *

def mpf_pow(s, t, prec, rnd=round_fast):
    """Compute s**t. Raises ComplexResult if s is negative and t is
    fractional."""
    ssign, sman, sexp, sbc = s
    tsign, tman, texp, tbc = t
    if ssign and texp < 0:
        raise ComplexResult("negative number raised to a fractional power")
    if texp >= 0:
        return mpf_pow_int(s, (-1)**tsign * (tman<<texp), prec, rnd)
    # s**(n/2) = sqrt(s)**n
    if texp == -1:
        if tman == 1:
            if tsign:
                return mpf_div(fone, mpf_sqrt(s, prec+10, reciprocal_rnd[rnd]), prec, rnd)
            return mpf_sqrt(s, prec, rnd)
        else:
            if tsign:
                return mpf_pow_int(mpf_sqrt(s, prec+10, reciprocal_rnd[rnd]), -tman, prec, rnd)
            return mpf_pow_int(mpf_sqrt(s, prec+10, rnd), tman, prec, rnd)
    # General formula: s**t = exp(t*log(s))
    # TODO: handle rnd direction of the logarithm carefully
    c = mpf_log(s, prec+10, rnd)
    return mpf_exp(mpf_mul(t, c), prec, rnd)


def int_pow_fixed(y, n, prec):
    """n-th power of a fixed point number with precision prec

       Returns the power in the form man, exp, 
       man * 2**exp ~= y**n
    """
    if n == 2:
        return (y*y), 0
    bc = bitcount(y)
    exp = 0
    workprec = 2 * (prec + 4*bitcount(n) + 4)
    _, pm, pe, pbc = fone
    while 1:
        if n & 1:
            pm = pm*y
            pe = pe+exp
            pbc += bc - 2
            pbc = pbc + bctable[int(pm >> pbc)]
            if pbc > workprec:
                pm = pm >> (pbc-workprec)
                pe += pbc - workprec
                pbc = workprec
            n -= 1
            if not n:
                break
        y = y*y
        exp = exp+exp
        bc = bc + bc - 2
        bc = bc + bctable[int(y >> bc)]
        if bc > workprec:
            y = y >> (bc-workprec)
            exp += bc - workprec
            bc = workprec
        n = n // 2
    return pm, pe

# froot(s, n, prec, rnd) computes the real n-th root of a
# positive mpf tuple s.
# To compute the root we start from a 50-bit estimate for r 
# generated with ordinary floating-point arithmetic, and then refine
# the value to full accuracy using the iteration

#            1  /                     y       \
#   r     = --- | (n-1)  * r   +  ----------  |
#    n+1     n  \           n     r_n**(n-1)  /

# which is simply Newton's method applied to the equation r**n = y.
# With giant_steps(start, prec+extra) = [p0,...,pm, prec+extra]
# and y = man * 2**-shift  one has
# (man * 2**exp)**(1/n) = 
# y**(1/n) * 2**(start-prec/n) * 2**(p0-start) * ... * 2**(prec+extra-pm) *
# 2**((exp+shift-(n-1)*prec)/n -extra))
# The last factor is accounted for in the last line of froot.

def nthroot_fixed(y, n, prec, exp1):
    start = 50
    try:
        y1 = rshift(y, prec - n*start)
        r = MP_BASE(y1**(1.0/n))
    except OverflowError:
        y1 = from_int(y1, start)
        fn = from_int(n)
        fn = mpf_rdiv_int(1, fn, start)
        r = mpf_pow(y1, fn, start)
        r = to_int(r)
    extra = 10
    extra1 = n
    prevp = start
    for p in giant_steps(start, prec+extra):
        pm, pe = int_pow_fixed(r, n-1, prevp)
        r2 = rshift(pm, (n-1)*prevp - p - pe - extra1)
        B = lshift(y, 2*p-prec+extra1)//r2
        r = (B + (n-1) * lshift(r, p-prevp))//n
        prevp = p
    return r

def mpf_nthroot(s, n, prec, rnd=round_fast):
    """nth-root of a positive number

    Use the Newton method when faster, otherwise use x**(1/n)
    """
    sign, man, exp, bc = s
    if sign:
        raise ComplexResult("nth root of a negative number")
    if not man:
        if s == fnan:
            return fnan
        if s == fzero:
            if n > 0:
                return fzero
            if n == 0:
                return fone
            return finf
        # Infinity
        if not n:
            return fnan
        if n < 0:
            return fzero
        return finf
    flag_inverse = False
    if n < 2:
        if n == 0:
            return fone
        if n == 1:
            return mpf_pos(s, prec, rnd)
        if n == -1:
            return mpf_div(fone, s, prec, rnd)
        # n < 0
        rnd = reciprocal_rnd[rnd]
        flag_inverse = True
        extra_inverse = 5
        prec += extra_inverse
        n = -n
    if n > 20 and (n >= 20000 or prec < int(233 + 28.3 * n**0.62)):
        prec2 = prec + 10
        fn = from_int(n)
        nth = mpf_rdiv_int(1, fn, prec2)
        r = mpf_pow(s, nth, prec2, rnd)
        s = normalize(r[0], r[1], r[2], r[3], prec, rnd)
        if flag_inverse:
            return mpf_div(fone, s, prec-extra_inverse, rnd)
        else:
            return s
    # Convert to a fixed-point number with prec2 bits.
    prec2 = prec + 2*n - (prec%n)
    # a few tests indicate that
    # for 10 < n < 10**4 a bit more precision is needed
    if n > 10:
        prec2 += prec2//10
        prec2 = prec2 - prec2%n
    # Mantissa may have more bits than we need. Trim it down.
    shift = bc - prec2
    # Adjust exponents to make prec2 and exp+shift multiples of n.
    sign1 = 0
    es = exp+shift
    if es < 0:
      sign1 = 1
      es = -es
    if sign1:
      shift += es%n
    else:
      shift -= es%n
    man = rshift(man, shift)
    extra = 10
    exp1 = ((exp+shift-(n-1)*prec2)//n) - extra
    rnd_shift = 0
    if flag_inverse:
        if rnd == 'u' or rnd == 'c':
            rnd_shift = 1
    else:
        if rnd == 'd' or rnd == 'f':
            rnd_shift = 1
    man = nthroot_fixed(man+rnd_shift, n, prec2, exp1)
    s = from_man_exp(man, exp1, prec, rnd)
    if flag_inverse:
        return mpf_div(fone, s, prec-extra_inverse, rnd)
    else:
        return s

def mpf_cbrt(s, prec, rnd=round_fast):
    """cubic root of a positive number"""
    return mpf_nthroot(s, 3, prec, rnd)

##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                           Mathematical constants                           #
#----------------------------------------------------------------------------#

def constant_memo(f):
    """Cache computed values of mathematical constants"""
    f.memo_prec = -1
    f.memo_val = None
    def g(prec):
        if prec == f.memo_prec:
            return f.memo_val
        if prec < f.memo_prec:
            return f.memo_val >> (f.memo_prec-prec)
        f.memo_val = f(prec)
        f.memo_prec = prec
        return f.memo_val
    g.__name__ = f.__name__
    g.__doc__ = f.__doc__
    return g

def acot(n, prec, hyperbolic):
    """Compute acot of an integer using fixed-point arithmetic. With
    hyperbolic=True, compute acoth. The standard Taylor series
    is used."""
    n = MP_BASE(n)
    s = t = (MP_ONE << prec) // n  # 1 / n
    k = 3
    while 1:
        # Repeatedly divide by k * n**2, and add
        t //= (n*n)
        term = t // k
        if not term:
            break
        # Alternate signs
        if hyperbolic or not k & 2:
            s += term
        else:
            s -= term
        k += 2
    return s

def machin(coefs, prec, hyperbolic=False):
    """Evaluate a Machin-like formula, i.e., a linear combination of
    acot(n) or acoth(n) for specific integer values of n, using fixed-
    point arithmetic. The input should be a list [(c, n), ...], giving
    c*acot[h](n) + ..."""
    extraprec = 10
    s = MP_ZERO
    for a, b in coefs:
        s += MP_BASE(a) * acot(MP_BASE(b), prec+extraprec, hyperbolic)
    return (s >> extraprec)

#----------------------------------------------------------------------------
# Pi

def agm_status(prec, step, adiff, verbose_base):
    logdiff = bitcount(adiff) * math.log(2,verbose_base)
    digits = int(prec/math.log(verbose_base,2) - logdiff)
    print "  iteration", step, ("(accuracy ~= %i base-%i digits)" % \
       (digits, verbose_base))

def pi_agm(prec, verbose=False, verbose_base=10):
    """
    Compute floor(pi * 2**prec) as a big integer using the Brent-
    Salamin algorithm based on the arithmetic-geometric mean.

    See for example Wikipedia (http://en.wikipedia.org/wiki/Brent-
    Salamin_algorithm) or "Pi and the AGM" by Jonathan and Peter
    Borwein (Wiley, 1987). The algorithm (as stated in the Wikipedia
    article) consists of setting

      a_0 = 1
      b_0 = 1/sqrt(2)
      t_0 = 1/4
      p_0 = 1

    and computing

      a_{n+1} = (a_n + b_n)/2
      b_{n+1} = sqrt(a_n * b_n)
      t_{n+1} = t_n - p_n*(a_n - a_{n+1})**2
      p_{n+1} = 2*p_n

    for n = 0, 1, 2, 3, ..., after which the approximation is given by
    pi ~= (a_n + b_n)**2 / (4*t_n). Each step roughly doubles the
    number of correct digits.
    """
    extraprec = 50
    prec += extraprec

    # Initialial values. a, b and t are fixed-point numbers
    a = MP_ONE << prec
    b = sqrt_fixed2(a >> 1, prec)
    t = a >> 2
    p = 1

    step = 1
    while 1:
        an = (a + b) >> 1
        adiff = a - an
        if verbose:
            agm_status(prec, step, adiff, verbose_base)
        # No change in a
        if p > 16 and abs(adiff) < 1000:
            break
        prod = (a * b) >> prec
        b = sqrt_fixed2(prod, prec)
        t = t - p*((adiff**2) >> prec)
        p = 2*p
        a = an
        step += 1
    if verbose:
        print "  final division"
    pi = ((((a+b)**2) >> 2) // t)
    return pi >> extraprec

# Chudnovsky's series for pi, with binary splitting
# The formulas are given in ftp://ftp.gmplib.org/pub/src/gmp-chudnovsky.c

# Constants in Chudnovsky's series
CHUD_A = MP_BASE(13591409)
CHUD_B = MP_BASE(545140134)
CHUD_C = MP_BASE(640320)
CHUD_D = MP_BASE(12)

def bs_chudnovsky(a,b,level,verbose):
    if b-a == 1:
        g = MP_BASE((6*b-5)*(2*b-1)*(6*b-1))
        p = b**3 * CHUD_C**3 // 24
        q = (-1)**b * g * (CHUD_A+CHUD_B*b)
    else:
        if verbose and level < 4:
            print "  binary splitting", a, b
        mid = (a+b)//2
        g1, p1, q1 = bs_chudnovsky(a,mid,level+1,verbose)
        g2, p2, q2 = bs_chudnovsky(mid,b,level+1,verbose)
        p = p1*p2
        g = g1*g2
        q = q1*p2 + q2*g1
    return g, p, q

def pi_chudnovsky(prec, verbose=False, verbose_base=None):
    """
    Calculate floor(pi * 2**prec) using Chudnovsky's algorithm.
    """
    # The Chudnovsky series gives 14.18 digits per term
    N = int(prec/3.3219280948/14.181647462 + 2)
    if verbose:
        print "binary splitting with N =", N
    g, p, q = bs_chudnovsky(0,N,0,verbose)
    sqrtC = sqrt_fixed(CHUD_C<<prec, prec)
    v = p*CHUD_C*sqrtC//((q+CHUD_A*p)*CHUD_D)
    return v

@constant_memo
def pi_fixed(prec):
    """
    Compute floor(pi * 2**prec) as a big integer.
    """
    #if prec < 2000:
    #    return machin([(16, 5), (-4, 239)], prec)
    #return pi_agm(prec)
    return pi_chudnovsky(prec)

def mpf_pi(prec, rnd=round_fast):
    """Compute a floating-point approximation of pi"""
    return from_man_exp(pi_fixed(prec+10), -prec-10, prec, rnd)

def mpf_degree(prec, rnd=round_fast):
    """Compute 1 degree = pi / 180."""
    return from_man_exp(pi_fixed(prec+10)//180, -prec-10, prec, rnd)

#----------------------------------------------------------------------------
# Logarithms of integers are needed for various computations involving
# logarithms, powers, radix conversion, etc
@constant_memo
def log2_fixed(prec):
    return machin([(18, 26), (-2, 4801), (8, 8749)], prec, True)

def mpf_log2(prec, rnd=round_fast):
    return from_man_exp(log2_fixed(prec+5), -prec-5, prec, rnd)

@constant_memo
def log10_fixed(prec):
    return machin([(46, 31), (34, 49), (20, 161)], prec, True)

def mpf_log10(prec, rnd=round_fast):
    return from_man_exp(log10_fixed(prec+5), -prec-5, prec, rnd)

#----------------------------------------------------------------------------
# exp(1) is computed from the Taylor series for exp, using binary splitting.

# See http://numbers.computation.free.fr/Constants/Algorithms/splitting.html
# for an explanation of this algorithm

def bspe(a,b):
    if b-a == 1:
        return MP_ONE, MP_BASE(b)
    m = (a+b)//2
    p1, q1 = bspe(a,m)
    p2, q2 = bspe(m,b)
    return p1*q2+p2, q1*q2

@constant_memo
def e_fixed(prec):
    # Slight overestimate of N needed for 1/N! < 2**(-prec)
    # This could be tightened for large N.
    N = int(1.1*prec/math.log(prec) + 20)
    p, q = bspe(0,N)
    return ((p+q)<<prec)//q

def mpf_e(prec, rnd=round_fast):
    return from_man_exp(e_fixed(prec+15), -prec-15, prec, rnd)


##############################################################################
##############################################################################

#----------------------------------------------------------------------------#
#                                                                            #
#                           Exponential function                             #
#                                                                            #
#----------------------------------------------------------------------------#

# The exponential function has a rapidly convergent Maclaurin series:
#
#     exp(x) = 1 + x + x**2/2! + x**3/3! + x**4/4! + ...
#
# The series can be summed very easily using fixed-point arithmetic.
# The convergence can be improved further, using a trick due to
# Richard P. Brent: instead of computing exp(x) directly, we choose a
# small integer r (say, r=10) and compute exp(x/2**r)**(2**r).

# The optimal value for r depends on the Python platform, the magnitude
# of x and the target precision, and has to be estimated from
# experimental timings. One test with x ~= 0.3 showed that
# r = 2.2*prec**0.42 gave a good fit to the optimal values for r for
# prec between 1 and 10000 bits, on one particular machine.

# This optimization makes the summation about twice as fast at
# low precision levels and much faster at high precision
# (roughly five times faster at 1000 decimal digits).

# If |x| is very large, we first rewrite it as t + n*log(2) with the
# integer n chosen such that |t| <= log(2), and then calculate
# exp(x) as exp(t)*(2**n), using the Maclaurin series for exp(t)
# (the multiplication by 2**n just amounts to shifting the exponent).

# Input: x * 2**prec
# Output: exp(x) * 2**(prec + r)
def exp_series(x, prec, r):
    x >>= r
    s = (MP_ONE << prec) + x
    a = x
    k = 2
    # Sum exp(x/2**r)
    while 1:
        a = ((a*x) >> prec) // k
        if not a:
            break
        s += a
        k += 1
    # Calculate s**(2**r) by repeated squaring
    while r:
        s = (s*s) >> prec
        r -= 1
    return s

def exp_series2(x, prec, r):
    x >>= r
    sign = 0
    if x < 0:
        sign = 1
        x = -x
    x2 = (x*x) >> prec
    s1 = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (k*(k-1))
        s1 += a
        k += 2
    c1 = sqrt_fixed(((s1*s1) >> prec) + (1<<prec), prec)
    if sign:
        s = c1 - s1
    else:
        s = c1 + s1
    # Calculate s**(2**r) by repeated squaring
    while r:
        s = (s*s) >> prec
        r -= 1
    return s

def mpf_exp(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if not exp:
            return fone
        if x == fninf:
            return fzero
        return x
    # Fast handling e**n. TODO: the best cutoff depends on both the
    # size of n and the precision.
    if prec > 600 and exp >= 0:
        return mpf_pow_int(mpf_e(prec+10), (-1)**sign *(man<<exp), prec, rnd)
    # extra precision needs to be similar in magnitude to log_2(|x|)
    # for the modulo reduction, plus r for the error from squaring r times
    wp = prec + max(0, bc+exp)
    if wp < 300:
        r = int(2*wp**0.4)
        if bc+exp < 0:
            r = max(1, r + bc + exp)
        wp += r + 20
        t = to_fixed(x, wp)
        # abs(x) > 1?
        if exp+bc > 1:
            lg2 = log2_fixed(wp)
            n, t = divmod(t, lg2)
        else:
            n = 0
        man = exp_series(t, wp, r)
    else:
        r = int(0.7 * wp**0.5)
        if bc+exp < 0:
            r = max(1, r + bc + exp)
        wp += r + 20
        t = to_fixed(x, wp)
        if exp+bc > 1:
            lg2 = log2_fixed(wp)
            n, t = divmod(t, lg2)
        else:
            n = 0
        man = exp_series2(t, wp, r)
    bc = wp - 2 + bctable[int(man >> (wp - 2))]
    return normalize(0, man, int(-wp+n), bc, prec, rnd)


#----------------------------------------------------------------------------#
#                                                                            #
#                                Logarithms                                  #
#                                                                            #
#----------------------------------------------------------------------------#

# The basic strategy for computing log(x) is to set r = log(x) and use
# Newton's method to solve the equation exp(r) = x.
# To obtain the Newton method one solves exp(r+h) = x, expanding in h
# which is supposed to be small; one has
# h = log(x * exp(-r)), which can be expanded in powers of
# s = (x * exp(-r) - 1) : h = s - s*s/2 + s**3/3 + ...
# The first order approximation is Newton method, the second order is
# Halley method. We use the second order approximation.
# We set the initial value r_0 to math.log(x) and then iterate
# r_{n+1} = (r_n + exp(-r_n) - 1) - (r_n + exp(-r_n) - 1)/2
# until convergence. As with square roots, we increase the working
# precision dynamically during the process so that only one full-precision
# evaluation of exp is required.

# log(x) is small for most inputs, so the r values can safely be
# computed using fixed-point arithmetic. However, when x has a very
# large or small exponent, we can improve performance through the
# normalization log(t * 2**n) = log(t) + n*log(2), choosing n such
# that 0.5 <= t <= 1 (for example).

# There are some caveats: if x is extremely close to 1, the working
# precision must be increased to maintain high relative precision in the
# output (alternatively, the series approximation for log(1+x) could
# be used in that case).

# This function performs the Newton iteration using fixed-point
# arithmetic. x is assumed to have magnitude ~= 1
def log_newton(x, prec):
    extra = 10
    # 40-bit approximation
    fx = math.log(long(x)) - 0.69314718055994529*prec
    r = MP_BASE(fx * 2.0**40)
    prevp = 40
    for p in giant_steps2(40, prec+extra):
        rb = lshift(r, p-prevp)
        # Parameters for exponential series
        r = int(2 * p**0.4)
        exp_extra = r + 10
        e = exp_series((-rb) << exp_extra, p + exp_extra, r)
        s = ((rshift(x, prec-p)*e)>>(p + exp_extra)) - (1 << p)
        s1 = -((s*s)>>(p+1))
        r = rb + s + s1
        prevp = p
    return r >> extra

def mpf_log(x, prec, rnd=round_fast):
    """Compute the natural logarithm of a positive raw mpf x."""
    sign, man, exp, bc = x
    if not man:
        if x == fzero:
            return fninf
        if x == finf:
            return finf
        return fnan
    if sign:
        raise ComplexResult("logarithm of a negative number")
    if x == fone:
        return fzero
    bc_plus_exp = bc + exp
    # Estimated precision needed for log(t) + n*log(2)
    prec2 = prec + int(math.log(1+abs(bc_plus_exp), 2)) + 15
    # Watch out for the case when x is very close to 1
    if -1 < bc_plus_exp < 2:
        near_one = mpf_abs(mpf_sub(x, fone, 53), 53)
        if near_one == 0:
            return fzero
        # estimate how close
        prec2 += -(near_one[2]) - bitcount(abs(near_one[1]))
    # Separate mantissa and exponent, calculate, join parts
    t = rshift(man, bc-prec2)
    l = log_newton(t, prec2)
    a = bc_plus_exp * log2_fixed(prec2)
    return from_man_exp(l+a, -prec2, prec, rnd)


#----------------------------------------------------------------------------#
#                                                                            #
#                          Trigonometric functions                           #
#                                                                            #
#----------------------------------------------------------------------------#

def sin_taylor(x, prec):
    x = MP_BASE(x)
    x2 = (x*x) >> prec
    s = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (k*(1-k))
        s += a
        k += 2
    return s

def cos_taylor(x, prec):
    x = MP_BASE(x)
    x2 = (x*x) >> prec
    a = c = (MP_ONE<<prec)
    k = 2
    while a:
        a = ((a * x2) >> prec) // (k*(1-k))
        c += a
        k += 2
    return c

# Input: x * 2**prec
# Output: c * 2**(prec + r), s * 2**(prec + r)
def expi_series(x, prec, r):
    x >>= r
    one = 1 << prec
    x2 = (x*x) >> prec
    s = x
    a = x
    k = 2
    while a:
        a = ((a * x2) >> prec) // (-k*(k+1))
        s += a
        k += 2
    c = sqrt_fixed(one - ((s*s)>>prec), prec)
    # Calculate (c + j*s)**(2**r) by repeated squaring
    for j in range(r):
      c, s =  (c*c-s*s) >> prec, (2*c*s ) >> prec
    return c, s

def reduce_angle(x, prec):
    """
    Let x be a nonzero, finite mpf value defining angle (measured in
    radians). Then reduce_trig(x, prec) returns (y, swaps, n) where:

      y = (man, wp) is the reduced angle as a scaled fixed-point
        number with precision wp, i.e. a floating-point number with
        exponent -wp. The mantissa is positive and has width ~equal
        to the input prec.

      swaps = (swap_cos_sin, cos_sign, sin_sign)
        Flags indicating the swaps that need to be applied
        to (cos(y), sin(y)) to obtain (cos(x), sin(x))

      n is an integer giving the original quadrant of x

    Calculation of the quadrant
    ===========================

    The integer n indices the quadrant of x. That is:

        ...
        -pi     <   x  < -pi/2     n = -2
        -pi/2   <   x  <  0        n = -1
        0       <   x  <  pi/2     n = 0
        pi/2    <   x  <  pi       n = 1
        pi      <   x  <  3*pi/2   n = 2
        3*pi/2  <   x  <  2*pi     n = 3
        2*pi    <   x  <  5*pi/2   n = 4
        ...

    Note that n does not wrap around. A quadrant index normalized to
    lie in [0, 1, 2, 3] can be found easily later on by computing
    n % 4. Keeping the extended information in n is crucial for
    interval arithmetic, as it allows one to distinguish between
    whether two points of a sine wave lie next to each other on
    a monotonic segment or are actually separated by a full
    period (or several periods).

    Note also that because is x is guaranteed to be rational, and
    all roots of the sine/cosine are irrational, all inequalities are
    strict. That is, we can always compute the correct quadrant.
    Care is required to do ensure that this is done right.

    Swaps
    =====

    The number y is a reduction of x to the first quadrant. This is
    essentially x mod pi/2. In fact, we reduce y further, to the first
    octant, by computing pi/2-x if x > pi/4.

    Due to the translation and mirror symmetries of trigonometric
    functions, this allows us to compute sin(x) or cos(x) by computing
    +/-sin(y) or +/-cos(y). The point, of course, is that if x
    is large, the Taylor series for y converges much more quickly
    than the one for x.

    """
    sign, man, exp, bc = x
    magnitude = exp + bc

    if not man:
        return (0, 0), (0, 0, 0), 0

    # Here we have abs(x) < 0.5. In this case no reduction is necessary.
    # TODO: could also handle abs(x) < 1
    if magnitude < 0:
        # Quadrant is 0 or -1
        n = -sign
        swaps = (0, 0, sign)
        fixed_exp = exp + bc - prec
        delta = fixed_exp - exp
        if delta < 0:
            man <<= (-delta)
        elif delta > 0:
            man >>= delta
        y = (man, -fixed_exp)
        return y, swaps, n

    i = 0
    while 1:
        cancellation_prec = 20 * 2**i
        wp = prec + abs(magnitude) + cancellation_prec
        pi1 = pi_fixed(wp)
        pi2 = pi1 >> 1
        pi4 = pi1 >> 2
        # Find nearest multiple
        n, y = divmod(to_fixed(x, wp), pi2)
        # Interchange cos/sin ?
        if y > pi4:
            swap_cos_sin = 1
            y = pi2 - y
        else:
            swap_cos_sin = 0
        # Now, the catch is that x might be extremely close to a
        # multiple of pi/2. This means accuracy is lost, and we may
        # even end up in the wrong quadrant, which is bad news
        # for interval arithmetic. This effect manifests by the
        # fixed-point value of y becoming small.  This is easy to check for.
        if y >> (prec + magnitude - 10):
            n = int(n)
            swaps = swap_table[swap_cos_sin^(n%2)][n%4]
            return (y>>magnitude, wp-magnitude), swaps, n
        i += 1

swap_table = ((0,0,0),(0,1,0),(0,1,1),(0,0,1)), ((1,0,0),(1,1,0),(1,1,1),(1,0,1))

def calc_cos_sin(which, y, swaps, prec, cos_rnd, sin_rnd):
    """
    Simultaneous computation of cos and sin (internal function).
    """
    y, wp = y
    swap_cos_sin, cos_sign, sin_sign = swaps

    if swap_cos_sin:
        which_compute = -which
    else:
        which_compute = which

    # XXX: assumes no swaps
    if not y:
        return fone, fzero

    # Tiny nonzero argument
    if wp > prec*2 + 30:
        y = from_man_exp(y, -wp)

        # For tiny arguments, rounding to nearest is trivial
        if cos_rnd == sin_rnd == round_nearest:
            cos, sin = fone, mpf_pos(y, prec, round_nearest)
            if swap_cos_sin:
                cos, sin = sin, cos
            if cos_sign: cos = mpf_neg(cos)
            if sin_sign: sin = mpf_neg(sin)
            return cos, sin

        # Directed rounding
        one_minus_eps = mpf_sub(fone, mpf_shift(fone, -prec-5), prec, round_down)
        y_plus_eps = mpf_add(y, mpf_shift(y, -prec-5), prec, round_up)
        y_minus_eps = mpf_sub(y, mpf_shift(y, -prec-5), prec, round_down)

        if swap_cos_sin:
            cos_rnd, sin_rnd = sin_rnd, cos_rnd
            cos_sign, sin_sign = sin_sign, cos_sign

        if cos_sign:
            cos = [mpf_neg(one_minus_eps), fnone]\
                [cos_rnd in (round_floor, round_up)]
        else:
            cos = [one_minus_eps, fone]\
                [cos_rnd in (round_ceiling, round_up)]
        if sin_sign:
            sin = [mpf_neg(y_minus_eps), mpf_neg(y_plus_eps)]\
                [sin_rnd in (round_floor, round_up)]
        else:
            sin = [y_minus_eps, y_plus_eps]\
                [sin_rnd in (round_ceiling, round_up)]

        if swap_cos_sin:
            cos, sin = sin, cos

        return cos, sin

    # Use standard Taylor series
    if prec < 600:
        if which_compute == 0:
            sin = sin_taylor(y, wp)
            # only need to evaluate one of the series
            cos = sqrt_fixed((1<<wp) - ((sin*sin)>>wp), wp)
        elif which_compute == 1:
            sin = 0
            cos = cos_taylor(y, wp)
        elif which_compute == -1:
            sin = sin_taylor(y, wp)
            cos = 0
    # Use exp(i*x) with Brent's trick 
    else:
        r = int(0.137 * prec**0.579)
        ep = r+20
        cos, sin = expi_series(y<<ep, wp+ep, r)
        cos >>= ep
        sin >>= ep

    if swap_cos_sin:
        cos, sin = sin, cos

    if cos_rnd is not round_nearest:
        # Round and set correct signs
        # XXX: this logic needs a second look
        ONE = MP_ONE << wp
        if cos_sign:
            cos += (-1)**(cos_rnd in (round_ceiling, round_down))
            cos = min(ONE, cos)
        else:
            cos += (-1)**(cos_rnd in (round_ceiling, round_up))
            cos = min(ONE, cos)
        if sin_sign:
            sin += (-1)**(sin_rnd in (round_ceiling, round_down))
            sin = min(ONE, sin)
        else:
            sin += (-1)**(sin_rnd in (round_ceiling, round_up))
            sin = min(ONE, sin)

    if which != -1:
        cos = normalize(cos_sign, cos, -wp, bitcount(cos), prec, cos_rnd)
    if which != 1:
        sin = normalize(sin_sign, sin, -wp, bitcount(sin), prec, sin_rnd)

    return cos, sin

def cos_sin(x, prec, rnd=round_fast, which=0):
    """
    Computes (cos(x), sin(x)). The parameter 'which' can disable
    evaluation of either cos or sin:

        0 -- return (cos(x), sin(x), n)
        1 -- return (cos(x), -,      n)
       -1 -- return (-,      sin(x), n)

    If only one function is wanted, this is slightly
    faster at low precision.
    """
    sign, man, exp, bc = x
    # Exact (or special) cases
    if not man:
        if exp:
            return (fnan, fnan)
        else:
            return (fone, fzero)
    y, swaps, n = reduce_angle(x, prec+10)
    return calc_cos_sin(which, y, swaps, prec, rnd, rnd)

def mpf_cos(x, prec, rnd=round_fast):
    return cos_sin(x, prec, rnd, 1)[0]

def mpf_sin(x, prec, rnd=round_fast):
    return cos_sin(x, prec, rnd, -1)[1]

def mpf_tan(x, prec, rnd=round_fast):
    c, s = cos_sin(x, prec+20)
    return mpf_div(s, c, prec, rnd)



#----------------------------------------------------------------------
# Hyperbolic functions
#

def sinh_taylor(x, prec):
    x = MP_BASE(x)
    x2 = (x*x) >> prec
    s = a = x
    k = 3
    while a:
        a = ((a * x2) >> prec) // (k*(k-1))
        s += a
        k += 2
    return s

def cosh_sinh(x, prec, rnd=round_fast, tanh=0):
    """Simultaneously compute (cosh(x), sinh(x)) for real x"""
    sign, man, exp, bc = x
    if (not man) and exp:
        if x == finf: return (finf, finf)
        if x == fninf: return (finf, fninf)
        return fnan, fnan

    if sign:
        man = -man

    high_bit = exp + bc
    prec2 = prec + 20

    if high_bit < -3:
        # Extremely close to 0, sinh(x) ~= x and cosh(x) ~= 1
        # TODO: support directed rounding
        if high_bit < -prec-2:
            return (fone, mpf_pos(x, prec, rnd))

        # Avoid cancellation when computing sinh
        # TODO: might be faster to use sinh series directly
        prec2 += (-high_bit) + 4

    # In the general case, we use
    #    cosh(x) = (exp(x) + exp(-x))/2
    #    sinh(x) = (exp(x) - exp(-x))/2
    # and note that the exponential only needs to be computed once.
    ep = mpf_exp(x, prec2)
    em = mpf_div(fone, ep, prec2)
    if tanh:
        ch = mpf_add(ep, em, prec2, rnd)
        sh = mpf_sub(ep, em, prec2, rnd)
        return mpf_div(sh, ch, prec, rnd)
    else:
        ch = mpf_shift(mpf_add(ep, em, prec, rnd), -1)
        sh = mpf_shift(mpf_sub(ep, em, prec, rnd), -1)
        return ch, sh

def mpf_cosh(x, prec, rnd=round_fast):
    """Compute cosh(x) for a real argument x"""
    return cosh_sinh(x, prec, rnd)[0]

def mpf_sinh(x, prec, rnd=round_fast):
    """Compute sinh(x) for a real argument x"""
    return cosh_sinh(x, prec, rnd)[1]

def mpf_tanh(x, prec, rnd=round_fast):
    """Compute tanh(x) for a real argument x"""
    return cosh_sinh(x, prec, rnd, tanh=1)


#----------------------------------------------------------------------
# Inverse tangent
#

"""
Near x = 0, use atan(x) = x - x**3/3 + x**5/5 - ...
Near x = 1, use atan(x) = y/x * (1 + 2/3*y + 2*4/3/5*y**2 + ...)
where y = x**2/(1+x**2).
At high precision use the Newton method  for 1.0/8 <= x < 256

TODO: see if one can improve the case x < 1.0/8;
      write a fixed-point function atan_newton
"""

def atan_taylor(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    assert not sign
    # Increase absolute precision when extremely close to 0
    diff = -(bc + exp)
    prec2 = prec
    if diff > 10:
        if 3*diff - 4 > prec:  # x**3 term vanishes; atan(x) ~x
            return from_man_exp(man, exp, prec, rnd)
        prec2 = prec + diff
    prec2 += 15  # XXX: better estimate for number of guard bits
    x = to_fixed(x, prec2)
    x2 = (x*x)>>prec2; one = 1<<prec2; s=a=x
    for n in xrange(1, 1000000):
        a = (a*x2) >> prec2
        s += a // ((-1)**n * (n+n+1))
        if -100 < a < 100:
            break
    return from_man_exp(s, -prec2, prec, rnd)

def atan_euler(x, prec, rnd=round_fast):
    prec2 = prec + 15
    x = to_fixed(x, prec2)
    one = 1<<prec2; x2 = (x*x)>>prec2; y=(x2<<prec2)//(one+x2)
    s = a = one
    for n in xrange(1, 1000000):
        a = ((a*y)>>prec2) * (2*n) // (2*n+1)
        if a < 100:
            break
        s += a
    return from_man_exp(y*s//x, -prec2, prec, rnd)

_cutoff_1 = (0, MP_FIVE, -3, 3)   # ~0.6
_cutoff_2 = (0, MP_THREE, -1, 2)   # 1.5

def mpf_atan(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if not man:
        if x == fzero: return fzero
        if x == finf: return mpf_shift(mpf_pi(prec), -1)
        if x == fninf: return mpf_neg(mpf_shift(mpf_pi(prec), -1))
        return fnan
    flag_nr = True
    if prec < 400:
        flag_nr = False
    else:
        ebc = exp + bc
    # the Newton method is used if flag_nr = True; otherwise one of the
    # two series is used.
    # This selection is based on a benchmark.
        if ebc < -2 or ebc > 8:
            flag_nr = False
        elif ebc == -2:
            if prec < 1500:
                flag_nr = False
        elif ebc <= 0:
            if prec < 600:
                flag_nr = False
        else:
            if prec < 400*ebc:
                flag_nr = False
    if not flag_nr:
        if sign:
            return mpf_neg(mpf_atan(mpf_neg(x), prec, rnd))
        if mpf_cmp(x, _cutoff_1) < 0:
            return atan_taylor(x, prec, rnd)
        if mpf_cmp(x, _cutoff_2) < 0:
            return atan_euler(x, prec, rnd)
        # For large x, use atan(x) = pi/2 - atan(1/x)
        if x[2] > 10*prec:
            pi = mpf_pi(prec, rnd)
            pihalf = mpf_shift(pi, -1)
        else:
            pi = mpf_pi(prec+4)
            pihalf = mpf_shift(pi, -1)
        t = mpf_atan(mpf_div(fone, x, prec+4), prec+4)
        return mpf_sub(pihalf, t, prec, rnd)
    # use Newton's method
    extra = 10
    extra_p = 100
    prec2 = prec + extra
    r = math.atan(to_float(x))
    r = from_float(r, 50, rnd)
    for p in giant_steps(50, prec2):
        wp = p + extra_p
        t = mpf_tan(r, wp, rnd)
        tmp1 = mpf_sub(x, t, wp, rnd)
        tmp2 = mpf_mul(t, t, wp, rnd)
        tmp2 = mpf_add(fone, tmp2, wp, rnd)
        tmp1 = mpf_div(tmp1, tmp2, wp, rnd)
        r = mpf_add(r, tmp1, wp, rnd)
    sign, man, exp, bc = r
    return normalize(sign, man, exp, bc, prec, rnd)

# TODO: cleanup the special cases
def mpf_atan2(y, x, prec, rnd=round_fast):
    xsign, xman, xexp, xbc = x
    ysign, yman, yexp, ybc = y
    if not yman:
        if y == fnan or x == fnan:
            return fnan
        if mpf_sign(x) >= 0:
            return fzero
        return mpf_pi(prec, rnd)
    if ysign:
        return mpf_neg(mpf_atan2(mpf_neg(y), x, prec, rnd))
    if not xman:
        if x == fnan:
            return fnan
        if x == finf:
            return fzero
        if x == fninf:
            return mpf_pi(prec, rnd)
        if not yman:
            return fzero
        return mpf_shift(mpf_pi(prec, rnd), -1)
    tquo = mpf_atan(mpf_div(y, x, prec+4), prec+4)
    if xsign:
        return mpf_add(mpf_pi(prec+4), tquo, prec, rnd)
    else:
        return mpf_pos(tquo, prec, rnd)

def mpf_asin(x, prec, rnd=round_fast):
    sign, man, exp, bc = x
    if bc+exp > 0 and x not in (fone, fnone):
        raise ComplexResult("asin(x) is real only for -1 <= x <= 1")
    flag_nr = True
    if prec < 1000 or exp+bc < -13:
        flag_nr = False
    else:
        ebc = exp + bc
        if ebc < -13:
            flag_nr = False
        elif ebc < -3:
            if prec < 3000:
                flag_nr = False
    if not flag_nr:
        # asin(x) = 2*atan(x/(1+sqrt(1-x**2)))
        wp = prec + 15
        a = mpf_mul(x, x)
        b = mpf_add(fone, mpf_sqrt(mpf_sub(fone, a, wp), wp), wp)
        c = mpf_div(x, b, wp)
        return mpf_shift(mpf_atan(c, prec, rnd), 1)
    # use Newton's method
    extra = 10
    extra_p = 10
    prec2 = prec + extra
    r = math.asin(to_float(x))
    r = from_float(r, 50, rnd)
    for p in giant_steps(50, prec2):
        wp = p + extra_p
        c, s = cos_sin(r, wp, rnd)
        tmp = mpf_sub(x, s, wp, rnd)
        tmp = mpf_div(tmp, c, wp, rnd)
        r = mpf_add(r, tmp, wp, rnd)
    sign, man, exp, bc = r
    return normalize(sign, man, exp, bc, prec, rnd)

def mpf_acos(x, prec, rnd=round_fast):
    # acos(x) = 2*atan(sqrt(1-x**2)/(1+x))
    sign, man, exp, bc = x
    if bc + exp > 0:
        if x not in (fone, fnone):
            raise ComplexResult("acos(x) is real only for -1 <= x <= 1")
        if x == fnone:
            return mpf_pi(prec, rnd)
    wp = prec + 15
    a = mpf_mul(x, x)
    b = mpf_sqrt(mpf_sub(fone, a, wp), wp)
    c = mpf_div(b, mpf_add(fone, x, wp), wp)
    return mpf_shift(mpf_atan(c, prec, rnd), 1)

def mpf_asinh(x, prec, rnd=round_fast):
    # asinh(x) = log(x+sqrt(x**2+1))
    wp = prec + 15
    q = mpf_sqrt(mpf_add(mpf_mul(x,x), fone, wp), wp)
    return mpf_log(mpf_add(x, q, wp), prec, rnd)

def mpf_acosh(x, prec, rnd=round_fast):
    # acosh(x) = log(x+sqrt(x**2-1))
    wp = prec + 15
    if mpf_cmp(x, fone) == -1:
        raise ComplexResult("acosh(x) is real only for x >= 1")
    q = mpf_sqrt(mpf_add(mpf_mul(x,x), fnone, wp), wp)
    return mpf_log(mpf_add(x, q, wp), prec, rnd)

def mpf_atanh(x, prec, rnd=round_fast):
    # atanh(x) = log((1+x)/(1-x))/2
    sign, man, exp, bc = x
    mag = bc + exp
    if mag > 0:
        raise ComplexResult("atanh(x) is real only for -1 < x < 1")
    wp = prec + 15
    a = mpf_add(x, fone, wp)
    b = mpf_sub(fone, x, wp)
    return mpf_shift(mpf_log(mpf_div(a, b, wp), prec, rnd), -1)


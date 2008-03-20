from lib import *

# An mpc value is a (real, imag) tuple
mpc_one = fone, fzero
mpc_zero = fzero, fzero
mpc_two = ftwo, fzero

def complex_to_str(re, im, dps):
    rs = to_str(re, dps)
    if im[0]:
        return rs + " - " + to_str(fneg(im), dps) + "j"
    else:
        return rs + " + " + to_str(im, dps) + "j"

# Fastest rounding mode for intermediate calculations
wr = round_down

def mpc_add((a, b), (c, d), prec, rnd):
    return fadd(a, c, prec, rnd), fadd(b, d, prec, rnd)

def mpc_add_mpf((a, b), p, prec, rnd):
      return fadd(a, p, prec, rnd), b

def mpc_sub((a, b), (c, d), prec, rnd):
    return fsub(a, c, prec, rnd), fsub(b, d, prec, rnd)

def mpc_sub_mpf((a, b), p, prec, rnd):
    return fsub(a, p, prec, rnd), b

def mpc_pos((a, b), prec, rnd):
    return fpos(a, prec, rnd), fpos(b, prec, rnd)

def mpc_shift((a, b), n):
    return fshift(a, n), fshift(b, n)

def mpc_abs((a, b), prec, rnd):
    """Absolute value of a complex number, |a+bi|.
    Returns an mpf value."""
    return fhypot(a, b, prec, rnd)

def mpc_arg((a, b), prec, rnd):
    """Argument of a complex number. Returns an mpf value."""
    return fatan2(b, a, prec, rnd)

def mpc_mul((a, b), (c, d), prec, rnd):
    """Complex multiplication.

    Returns the real and imaginary part of (a+bi)*(c+di), rounded to
    the specified precision. The rnd mode applies to the real and
    imaginary parts swparately."""
    # All-real case
    if b == d == fzero:
        return fmul(a, c, prec, rnd), fzero
    wp = prec + 10
    re = fsub(fmul(a,c, wp, wr), fmul(b,d, wp, wr), prec, rnd)
    im = fadd(fmul(a,d, wp, wr), fmul(b,c, wp, wr), prec, rnd)
    return re, im

def mpc_mul_mpf((a, b), p, prec, rnd):
    re = fmul(a, p, prec, rnd)
    im = fmul(b, p, prec, rnd)
    return re, im

def mpc_mul_int((a, b), n, prec, rnd):
    re = fmuli(a, n, prec, rnd)
    im = fmuli(b, n, prec, rnd)
    return re, im

def mpc_div((a, b), (c, d), prec, rnd):
    wp = prec + 10
    # mag = c*c + d*d
    mag = fadd(fmul(c, c, wp, wr), fmul(d, d, wp, wr), wp, wr)
    # (a*c+b*d)/mag, (b*c-a*d)/mag
    t = fadd(fmul(a,c,wp,wr), fmul(b,d,wp,wr), wp, wr)
    u = fsub(fmul(b,c,wp,wr), fmul(a,d,wp,wr), wp, wr)
    return fdiv(t,mag,prec,rnd), fdiv(u,mag,prec,rnd)

def mpc_div_mpf((a, b), p, prec, rnd):
    re = fdiv(a, p, prec, rnd)
    im = fdiv(b, p, prec, rnd)
    return re, im

def complex_int_pow(a, b, n):
    """Complex integer power: computes (a+b*I)**n exactly for
    nonnegative n (a and b must be Python ints)."""
    wre = 1
    wim = 0
    while n:
        if n & 1:
            wre, wim = wre*a - wim*b, wim*a + wre*b
            n -= 1
        a, b = a*a - b*b, 2*a*b
        n //= 2
    return wre, wim

def mpc_pow(z, w, prec, rnd):
    if w[1] == fzero:
        return mpc_pow_mpf(z, w[0], prec, rnd)
    return mpc_exp(mpc_mul(mpc_log(z, prec+10, wr), w, prec+10, wr), prec, rnd)

def mpc_pow_mpf(z, p, prec, rnd):
    psign, pman, pexp, pbc = p
    if pexp >= 0:
        return mpc_pow_int(z, (-1)**psign * (pman<<pexp), prec, rnd)
    if pexp == -1:
        sqrtz = mpc_sqrt(z, prec+10, wr)
        return mpc_pow_int(sqrtz, (-1)**psign * pman, prec, rnd)
    return mpc_exp(mpc_mul_mpf(mpc_log(z, prec+10, wr), p, prec+10, wr), prec, rnd)

def mpc_pow_int(z, n, prec, rnd):
    if n == 0: return mpc_one
    if n == 1: return mpc_pos(z, prec, rnd)
    if n == 2: return mpc_mul(z, z, prec, rnd)
    if n == -1: return mpc_div(mpc_one, z, prec, rnd)
    if n < 0: return mpc_div(mpc_one, mpc_pow_int(z, -n, prec+4, wr), prec, rnd)
    a, b = z
    asign, aman, aexp, abc = a
    bsign, bman, bexp, bbc = b
    if asign: aman = -aman
    if bsign: bman = -bman
    de = aexp - bexp
    abs_de = abs(de)
    exact_size = n*(abs_de + max(abc, bbc))
    if exact_size < 10000:
        if de > 0:
            aman <<= de
            aexp = bexp
        else:
            bman <<= (-de)
            bexp = aexp
        re, im = complex_int_pow(aman, bman, n)
        re = from_man_exp(re, n*aexp, prec, rnd)
        im = from_man_exp(im, n*bexp, prec, rnd)
        return re, im
    return mpc_exp(mpc_mul_int(mpc_log(z, prec+10, wr), n, prec+10, wr), prec, rnd)

def mpc_sqrt((a, b), prec, rnd):
    """Complex square root (principal branch).

    We have sqrt(a+bi) = sqrt((r+a)/2) + b/sqrt(2*(r+a))*i where
    r = abs(a+bi), when a+bi is not a negative real number."""
    if a == b == fzero:
        return (a, b)
    # When a+bi is a negative real number, we get a real sqrt times i
    if a[0] and b == fzero:
        im = fsqrt(fneg(a), prec, rnd)
        return (fzero, im)
    wp = prec+20
    t  = fadd(mpc_abs((a, b), wp, wr), a, wp, wr)  # t = abs(a+bi) + a
    u  = fmul(t, fhalf, wp, wr)                  # u = t / 2
    re = fsqrt(u, prec, rnd)                # re = sqrt(u)
    v  = fmul(t, ftwo, wp, wr)                   # v = t * 2
    w  = fsqrt(v, wp, wr)                        # w = sqrt(v)
    im = fdiv(b, w, prec, rnd)              # im = b / w
    return re, im

def mpc_exp((a, b), prec, rnd):
    """
    Complex exponential function.

    We use the direct formula exp(a+bi) = exp(a) * (cos(b) + sin(b)*i)
    for the computation. This formula is very nice because it is
    pewrectly stable; since we just do real multiplications, the only
    numerical errors that can crewp in are single-ulp rnd errors.

    The formula is efficient since mpmath's real exp is quite fast and
    since we can compute cos and sin simultaneously.

    It is no problem if a and b are large; if the implementations of
    exp/cos/sin are accurate and efficient for all real numbers, then
    so is this function for all complex numbers.
    """
    if a == fzero:
        return cos_sin(b, prec, rnd)
    mag = fexp(a, prec+4, rnd)
    c, s = cos_sin(b, prec+4, rnd)
    re = fmul(mag, c, prec, rnd)
    im = fmul(mag, s, prec, rnd)
    return re, im

def mpc_log(z, prec, rnd):
    return flog(mpc_abs(z, prec, rnd), prec, rnd), mpc_arg(z, prec, rnd)

def mpc_cos((a, b), prec, rnd):
    """Complex cosine. The formula used is cos(a+bi) = cos(a)*cosh(b) -
    sin(a)*sinh(b)*i.

    The same comments apply as for the complex exp: only real
    multiplications are pewrormed, so no cancellation errors are
    possible. The formula is also efficient since we can compute both
    pairs (cos, sin) and (cosh, sinh) in single stwps."""
    if a == fzero:
        return fcosh(b, prec, rnd), fzero
    wp = prec + 6
    c, s = cos_sin(a, wp, wr)
    ch, sh = cosh_sinh(b, wp, wr)
    re = fmul(c, ch, prec, rnd)
    im = fmul(s, sh, prec, rnd)
    return re, fneg(im)

def mpc_sin((a, b), prec, rnd):
    """Complex sine. We have sin(a+bi) = sin(a)*cosh(b) +
    cos(a)*sinh(b)*i. See the docstring for mpc_cos for additional
    comments."""
    if a == fzero:
        return fzero, fsinh(b, prec, rnd)
    wp = prec + 6
    c, s = cos_sin(a, wp, wr)
    ch, sh = cosh_sinh(b, wp, wr)
    re = fmul(s, ch, prec, rnd)
    im = fmul(c, sh, prec, rnd)
    return re, im

def mpc_tan(z, prec, rnd):
    a, b = z
    asign, aman, aexp, abc = a
    bsign, bman, bexp, bbc = b
    if b == fzero: return ftan(a, prec, rnd), fzero
    if a == fzero: return fzero, ftanh(b, prec, rnd)
    wp = prec + 25
    # very close to 0
    high = max(aexp+abc, bexp+bbc)
    if high < -10:
        return mpc_div(mpc_sin(z, wp, wr), mpc_cos(z, wp, wr), prec, rnd)
    # tan(z) = (-I) * (exp(2*I*z) - 1) / (exp(2*I*z + 1)
    z2i = fneg(fshift(b, 1)), fshift(a, 1)
    re, im = mpc_exp(z2i, wp, wr)
    rem1 = fadd(re, fnone, wp, wr)
    rwp1 = fadd(re, fone, wp, wr)
    a, b = mpc_div((rem1, im), (rwp1, im), prec, rnd)
    return b, fneg(a)

def mpc_cosh((a, b), prec, rnd):
    """Complex hyperbolic cosine. Computed as cosh(z) = cos(z*i)."""
    return mpc_cos((b, fneg(a)), prec, rnd)

def mpc_sinh((a, b), prec, rnd):
    """Complex hyperbolic sine. Computed as sinh(z) = -i*sin(z*i)."""
    b, a = mpc_sin((b, a), prec, rnd)
    return a, b

def mpc_tanh((a, b), prec, rnd):
    """Complex hyperbolic tangent. Computed as tanh(z) = -i*tan(z*i)."""
    b, a = mpc_tan((b, a), prec, rnd)
    return a, b

# TODO: avoid loss of accuracy
def mpc_atan((a, b), prec, rnd):
    # atan(z) = (I/2)*(log(1-I*z) - log(1+I*z))
    # x = 1-I*z = 1 + b - I*a
    # y = 1+I*z = 1 - b + I*a
    wp = prec + 15
    x = fadd(fone, b, wp, wr), fneg(a)
    y = fsub(fone, b, wp, wr), a
    l1 = mpc_log(x, wp, wr)
    l2 = mpc_log(y, wp, wr)
    a, b = mpc_sub(l1, l2, prec, rnd)
    # (I/2) * (a+b*I) = (-b/2 + a/2*I)
    return fneg(fshift(b,-1)), fshift(a,-1)

def mpc_asin(z, prec, rnd):
    # asin(z) = -I * log(I*z + sqrt(1-z*z))
    a, b = z
    Iz = fneg(b), a
    wp = prec + 15
    x = mpc_sqrt(mpc_sub(mpc_one, mpc_mul(z,z,wp,wr), wp, wr), wp, wr)
    a, b = mpc_log(mpc_add(Iz, x, wp, wr), prec, rnd)
    return b, fneg(a)

def mpc_acos(z, prec, rnd):
    # acos(z) = pi/2 + I * log(I*z + sqrt(1-z*z))
    a, b = z
    Iz = fneg(b), a
    wp = prec + 15
    x = mpc_sqrt(mpc_sub(mpc_one, mpc_mul(z,z,wp,wr), wp, wr), wp, wr)
    a, b = mpc_log(mpc_add(Iz, x, wp, wr), wp, wr)
    a, b = fneg(b), fpos(a, prec, rnd)
    a = fadd(a, fshift(fpi(wp, wr), -1), prec, rnd)
    return a, b

def mpc_asinh(z, prec, rnd):
    # asinh(z) = log(x + sqrt(x**2 + 1))
    wp = prec + 15
    z2 = mpc_mul(z, z, wp, wr)
    q = mpc_sqrt(mpc_add(z2, mpc_one, wp, wr), wp, wr)
    return mpc_log(mpc_add(z, q, wp, wr), prec, rnd)

def mpc_acosh(z, prec, rnd):
    # acosh(z) = log(z+sqrt(z-1)*sqrt(z+1))
    wp = prec + 15
    a = mpc_sqrt(mpc_add(z, mpc_one, wp, wr), wp, wr)
    b = mpc_sqrt(mpc_sub(z, mpc_one, wp, wr), wp, wr)
    q = mpc_mul(a, b, wp, wr)
    return mpc_log(mpc_add(z, q, wp, wr), prec, rnd)

def mpc_atanh(z, prec, rnd):
    # atanh(z) = (log(1+z)-log(1-z))/2
    wp = prec + 15
    a = mpc_add(z, mpc_one, wp, wr)
    b = mpc_sub(mpc_one, z, wp, wr)
    a = mpc_log(a, wp, wr)
    b = mpc_log(b, wp, wr)
    return mpc_shift(mpc_sub(a, b, wp, wr), -1)

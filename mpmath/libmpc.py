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
rf = round_down

def mpc_add((a, b), (c, d), prec, rnd):
    return fadd(a, c, prec, rnd), fadd(b, d, prec, rnd)

def mpc_sub((a, b), (c, d), prec, rnd):
    return fsub(a, c, prec, rnd), fsub(b, d, prec, rnd)

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
    imaginary parts separately."""
    # All-real case
    if b == d == fzero:
        return fmul(a, c, prec, rnd), fzero
    ep = prec + 10
    re = fsub(fmul(a,c, ep, rf), fmul(b,d, ep, rf), prec, rnd)
    im = fadd(fmul(a,d, ep, rf), fmul(b,c, ep, rf), prec, rnd)
    return re, im

def mpc_div((a, b), (c, d), prec, rnd):
    wp = prec + 10
    # mag = c*c + d*d
    mag = fadd(fmul(c, c, wp, rf), fmul(d, d, wp, rf), wp, rf)
    # (a*c+b*d)/mag, (b*c-a*d)/mag
    t = fadd(fmul(a,c,wp,rf), fmul(b,d,wp,rf), wp, rf)
    u = fsub(fmul(b,c,wp,rf), fmul(a,d,wp,rf), wp, rf)
    return fdiv(t,mag,prec,rnd), fdiv(u,mag,prec,rnd)

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
    ep = prec+20
    t  = fadd(mpc_abs((a, b), ep, rf), a, ep, rf)  # t = abs(a+bi) + a
    u  = fmul(t, fhalf, ep, rf)                  # u = t / 2
    re = fsqrt(u, prec, rnd)                # re = sqrt(u)
    v  = fmul(t, ftwo, ep, rf)                   # v = t * 2
    w  = fsqrt(v, ep, rf)                        # w = sqrt(v)
    im = fdiv(b, w, prec, rnd)              # im = b / w
    return re, im

def mpc_exp((a, b), prec, rnd):
    """
    Complex exponential function.

    We use the direct formula exp(a+bi) = exp(a) * (cos(b) + sin(b)*i)
    for the computation. This formula is very nice because it is
    perfectly stable; since we just do real multiplications, the only
    numerical errors that can creep in are single-ulp rnd errors.

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
    return flog(mpc_abs(z, prec, rnd), prec, rnd), \
        mpc_arg(z, prec, rnd)

def mpc_cos((a, b), prec, rnd):
    """Complex cosine. The formula used is cos(a+bi) = cos(a)*cosh(b) -
    sin(a)*sinh(b)*i.

    The same comments apply as for the complex exp: only real
    multiplications are performed, so no cancellation errors are
    possible. The formula is also efficient since we can compute both
    pairs (cos, sin) and (cosh, sinh) in single steps."""
    ep = prec + 6
    c, s = cos_sin(a, ep, rf)
    ch, sh = cosh_sinh(b, ep, rf)
    re = fmul(c, ch, prec, rnd)
    im = fmul(s, sh, prec, rnd)
    return re, fneg(im)

def mpc_sin((a, b), prec, rnd):
    """Complex sine. We have sin(a+bi) = sin(a)*cosh(b) +
    cos(a)*sinh(b)*i. See the docstring for mpc_cos for additional
    comments."""
    ep = prec + 6
    c, s = cos_sin(a, ep, rf)
    ch, sh = cosh_sinh(b, ep, rf)
    re = fmul(s, ch, prec, rnd)
    im = fmul(c, sh, prec, rnd)
    return re, im

def mpc_cosh((a, b), prec, rnd):
    """Complex hyperbolic cosine. Computed as cosh(z) = cos(z*i)."""
    return mpc_cos((b, fneg(a)), prec, rnd)

def mpc_sinh((a, b), prec, rnd):
    """Complex hyperbolic sine. Computed as sinh(z) = -i*sin(z*i)."""
    b, a = mpc_sin((b, a), prec, rnd)
    return a, b

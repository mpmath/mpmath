from lib import (
    round_down, round_up, round_floor, round_ceiling, round_nearest,
    ComplexResult,
    fnan, finf, fninf, fzero, fhalf, fone, fnone,
    fsign, flt, fle, fgt, fge, feq, fcmp, ffloor, from_int, to_int,
    to_str, prec_to_dps,
    fabs, fneg, fpos, fadd, fsub, fmul, fdiv, fshift, fpowi,
    flog, fexp, fsqrt,
    reduce_angle, calc_cos_sin)

def mpi_str(s, prec):
    sa, sb = s
    dps = prec_to_dps(prec) + 5
    return "[%s, %s]" % (to_str(sa, dps), to_str(sb, dps))

    #dps = prec_to_dps(prec)
    #m = mpi_mid(s, prec)
    #d = fshift(mpi_delta(s, 20), -1)
    #return "%s +/- %s" % (to_str(m, dps), to_str(d, 3))

def mpi_add(s, t, prec):
    sa, sb = s
    ta, tb = t
    a = fadd(sa, ta, prec, round_floor)
    b = fadd(sb, tb, prec, round_ceiling)
    if a == fnan: a = fninf
    if b == fnan: b = finf
    return a, b

def mpi_sub(s, t, prec):
    sa, sb = s
    ta, tb = t
    a = fsub(sa, tb, prec, round_floor)
    b = fsub(sb, ta, prec, round_ceiling)
    if a == fnan: a = fninf
    if b == fnan: b = finf
    return a, b

def mpi_delta(s, prec):
    sa, sb = s
    return fsub(sb, sa, prec, round_up)

def mpi_mid(s, prec):
    sa, sb = s
    return fshift(fadd(sa, sb, prec, round_nearest), -1)

def mpi_pos(s, prec):
    sa, sb = s
    a = fpos(sa, prec, round_floor)
    b = fpos(sb, prec, round_ceiling)
    return a, b

def mpi_neg(s, prec=None):
    sa, sb = s
    a = fneg(sb, prec, round_floor)
    b = fneg(sa, prec, round_ceiling)
    return a, b

def mpi_abs(s, prec):
    sa, sb = s
    sas = fsign(sa)
    sbs = fsign(sb)
    # Both points nonnegative?
    if sas >= 0:
        a = fpos(sa, prec, round_floor)
        b = fpos(sb, prec, round_ceiling)
    # Upper point nonnegative?
    elif sbs >= 0:
        a = fzero
        negsa = fneg(sa)
        if flt(negsa, sb):
            b = fpos(sb, prec, round_ceiling)
        else:
            b = fpos(negsa, prec, round_ceiling)
    # Both negative?
    else:
        a = fneg(sb, prec, round_floor)
        b = fneg(sa, prec, round_ceiling)
    return a, b

def mpi_mul(s, t, prec):
    sa, sb = s
    ta, tb = t
    sas = fsign(sa)
    sbs = fsign(sb)
    tas = fsign(ta)
    tbs = fsign(tb)
    if sas == sbs == 0:
        # Should maybe be undefined
        if ta == fninf or tb == finf:
            return fninf, finf
        return fzero, fzero
    if tas == tbs == 0:
        # Should maybe be undefined
        if sa == fninf or sb == finf:
            return fninf, finf
        return fzero, fzero
    if sas >= 0:
        # positive * positive
        if tas >= 0:
            a = fmul(sa, ta, prec, round_floor)
            b = fmul(sb, tb, prec, round_ceiling)
            if a == fnan: a = fzero
            if b == fnan: b = finf
        # positive * negative
        elif tbs <= 0:
            a = fmul(sb, ta, prec, round_floor)
            b = fmul(sa, tb, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = fzero
        # positive * both signs
        else:
            a = fmul(sb, ta, prec, round_floor)
            b = fmul(sb, tb, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = finf
    elif sbs <= 0:
        # negative * positive
        if tas >= 0:
            a = fmul(sa, tb, prec, round_floor)
            b = fmul(sb, ta, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = fzero
        # negative * negative
        elif tbs <= 0:
            a = fmul(sb, tb, prec, round_floor)
            b = fmul(sa, ta, prec, round_ceiling)
            if a == fnan: a = fzero
            if b == fnan: b = finf
        # negative * both signs
        else:
            a = fmul(sb, tb, prec, round_floor)
            b = fmul(sa, ta, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = finf
    else:
        # General case: perform all cross-multiplications and compare
        # Since the multiplications can be done exactly, we need only
        # do 4 (instead of 8: two for each rounding mode)
        cases = [fmul(sa, ta), fmul(sa, tb), fmul(sb, ta), fmul(sb, tb)]
        if fnan in cases:
            a, b = (fninf, finf)
        else:
            cases = sorted(cases, cmp=fcmp)
            a = fpos(cases[0], prec, round_floor)
            b = fpos(cases[-1], prec, round_ceiling)
    return a, b

def mpi_div(s, t, prec):
    sa, sb = s
    ta, tb = t
    sas = fsign(sa)
    sbs = fsign(sb)
    tas = fsign(ta)
    tbs = fsign(tb)
    # 0 / X
    if sas == sbs == 0:
        # 0 / <interval containing 0>
        if (tas < 0 and tbs > 0) or (tas == 0 or tbs == 0):
            return fninf, finf
        return fzero, fzero
    # Denominator contains both negative and positive numbers;
    # this should properly be a multi-interval, but the closest
    # match is the entire (extended) real line
    if tas < 0 and tbs > 0:
        return fninf, finf
    # Assume denominator to be nonnegative
    if tas < 0:
        return mpi_div(mpi_neg(s), mpi_neg(t), prec)
    # Division by zero
    # XXX: make sure all results make sense
    if tas == 0:
        # Numerator contains both signs?
        if sas < 0 and sbs > 0:
            return fninf, finf
        if tas == tbs:
            return fninf, finf
        # Numerator positive?
        if sas >= 0:
            a = fdiv(sa, tb, prec, round_floor)
            b = finf
        if sbs <= 0:
            a = fninf
            b = fdiv(sb, tb, prec, round_ceiling)
    # Division with positive denominator
    # We still have to handle nans resulting from inf/0 or inf/inf
    else:
        # Nonnegative numerator
        if sas >= 0:
            a = fdiv(sa, tb, prec, round_floor)
            b = fdiv(sb, ta, prec, round_ceiling)
            if a == fnan: a = fzero
            if b == fnan: b = finf
        # Nonpositive numerator
        elif sbs <= 0:
            a = fdiv(sa, ta, prec, round_floor)
            b = fdiv(sb, tb, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = fzero
        # Numerator contains both signs?
        else:
            a = fdiv(sa, ta, prec, round_floor)
            b = fdiv(sb, ta, prec, round_ceiling)
            if a == fnan: a = fninf
            if b == fnan: b = finf
    return a, b

def mpi_exp(s, prec):
    sa, sb = s
    # exp is monotonous
    a = fexp(sa, prec, round_floor)
    b = fexp(sb, prec, round_ceiling)
    return a, b

def mpi_log(s, prec):
    sa, sb = s
    # log is monotonous
    a = flog(sa, prec, round_floor)
    b = flog(sb, prec, round_ceiling)
    return a, b

def mpi_sqrt(s, prec):
    sa, sb = s
    # sqrt is monotonous
    a = fsqrt(sa, prec, round_floor)
    b = fsqrt(sb, prec, round_ceiling)
    return a, b

def mpi_pow_int(s, n, prec):
    sa, sb = s
    if n < 0:
        return mpi_div((fone, fone), mpi_pow_int(s, -n, prec+20), prec)
    if n == 0:
        return (fone, fone)
    if n == 1:
        return s
    # Odd -- signs are preserved
    if n & 1:
        a = fpowi(sa, n, prec, round_floor)
        b = fpowi(sb, n, prec, round_ceiling)
    # Even -- important to ensure positivity
    else:
        sas = fsign(sa)
        sbs = fsign(sb)
        # Nonnegative?
        if sas >= 0:
            a = fpowi(sa, n, prec, round_floor)
            b = fpowi(sb, n, prec, round_ceiling)
        # Nonpositive?
        elif sbs <= 0:
            a = fpowi(sb, n, prec, round_floor)
            b = fpowi(sa, n, prec, round_ceiling)
        # Mixed signs?
        else:
            a = fzero
            # max(-a,b)**n
            sa = fneg(sa)
            if fge(sa, sb):
                b = fpowi(sa, n, prec, round_ceiling)
            else:
                b = fpowi(sb, n, prec, round_ceiling)
    return a, b

def mpi_pow(s, t, prec):
    ta, tb = t
    if ta == tb and ta not in (finf, fninf):
        if ta == from_int(to_int(ta)):
            return mpi_pow_int(s, to_int(ta), prec)
        if ta == fhalf:
            return mpi_sqrt(s, prec)
    u = mpi_log(s, prec + 20)
    v = mpi_mul(u, t, prec + 20)
    return mpi_exp(v, prec)

def MIN(x, y):
    if fle(x, y):
        return x
    return y

def MAX(x, y):
    if fge(x, y):
        return x
    return y

def mpi_cos_sin(x, prec):
    a, b = x
    # Guaranteed to contain both -1 and 1
    if finf in (a, b) or fninf in (a, b):
        return (fnone, fone), (fnone, fone)
    y, yswaps, yn = reduce_angle(a, prec+20)
    z, zswaps, zn = reduce_angle(b, prec+20)
    # Guaranteed to contain both -1 and 1
    if zn - yn >= 4:
        return (fnone, fone), (fnone, fone)
    # Both points in the same quadrant -- cos and sin both strictly monotonous
    if yn == zn:
        m = yn % 4
        if m == 0:
            cb, sa = calc_cos_sin(0, y, yswaps, prec, round_ceiling, round_floor)
            ca, sb = calc_cos_sin(0, z, zswaps, prec, round_floor, round_ceiling)
        if m == 1:
            cb, sb = calc_cos_sin(0, y, yswaps, prec, round_ceiling, round_ceiling)
            ca, sa = calc_cos_sin(0, z, zswaps, prec, round_floor, round_ceiling)
        if m == 2:
            ca, sb = calc_cos_sin(0, y, yswaps, prec, round_floor, round_ceiling)
            cb, sa = calc_cos_sin(0, z, zswaps, prec, round_ceiling, round_floor)
        if m == 3:
            ca, sa = calc_cos_sin(0, y, yswaps, prec, round_floor, round_floor)
            cb, sb = calc_cos_sin(0, z, zswaps, prec, round_ceiling, round_ceiling)
        return (ca, cb), (sa, sb)
    # Intervals spanning multiple quadrants
    yn %= 4
    zn %= 4
    case = (yn, zn)
    if case == (0, 1):
        cb, sy = calc_cos_sin(0, y, yswaps, prec, round_ceiling, round_floor)
        ca, sz = calc_cos_sin(0, z, zswaps, prec, round_floor, round_floor)
        return (ca, cb), (MIN(sy, sz), fone)
    if case == (3, 0):
        cy, sa = calc_cos_sin(0, y, yswaps, prec, round_floor, round_floor)
        cz, sb = calc_cos_sin(0, z, zswaps, prec, round_floor, round_ceiling)
        return (MIN(cy, cz), fone), (sa, sb)


    raise NotImplementedError("cos/sin spanning multiple quadrants")

def mpi_cos(x, prec):
    return mpi_cos_sin(x, prec)[0]

def mpi_sin(x, prec):
    return mpi_cos_sin(x, prec)[1]

def mpi_tan(x, prec):
    cos, sin = mpi_cos_sin(x, prec+20)
    return mpi_div(sin, cos, prec)

def mpi_cot(x, prec):
    cos, sin = mpi_cos_sin(x, prec+20)
    return mpi_div(cos, sin, prec)

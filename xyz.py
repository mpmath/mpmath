import operator
import math
import sys

import mpmath

stddigits = '0123456789abcdefghijklmnopqrstuvwxyz'
stddigits_as_bytes = bytearray(stddigits.encode('ascii'))

def fpp2(f, e, p, B=10, k=None):
    cmp_low = operator.lt if f & 1 else operator.le
    cmp_high = operator.lt if f & 1 else operator.le
    ep = e - p
    R = f << max(+ep, 0) + 1
    S = 1 << max(-ep, 0) + 1
    Mminus = Mplus = 1 << max(ep, 0)
    if f == 1 << (p - 1):
        Mplus <<= 1
        R <<= 1
        S <<= 1
    k = None
    if k is None:
        k = 0
        while R + Mplus <= S:
            k -= 1
            R *= B
            Mplus *= B
            Mminus *= B
        while R + Mplus >= S*B:
            k += 1
            S *= B
    else:
        if k >= 0:
            S *= B ** k
        else:
            R *= B ** (-k)
            Mminus *= B ** (-k)
            Mplus *= B ** (-k)
    assert R + Mplus >= S
    E = k
    D = bytearray()
    low = False
    high = False
    while True:
        U, R = divmod(R, S)
        low = cmp_low(R, Mminus)
        high = cmp_high(S, R + Mplus)
        D.append(stddigits_as_bytes[U])
        if low or high:
            break
        R *= B
        Mminus *= B
        Mplus *= B
    round_up = high
    if low and high:
        round_up = 2*R >= S
        if round_up and 2*R == S:
            round_up = U & 1
    if round_up:
        assert ord('0') <= D[-1] < ord(stddigits[B - 1])
        D[-1] += 1
    return D.decode(), E


def burger_dybvig(f, prec=53, B=10, scientific=False):
    """
    Converts a positive float to its shortest unique decimal string
    using the Burger & Dybvig algorithm.
    """
    if not f or not mpmath.isfinite(f):
        return str(float(f))

    sign, man, exp, bc = f._mpf_
    man <<= prec - bc
    exp += bc

    k = math.floor(mpmath.workprec(prec + 10)(mpmath.log)(abs(f), B))
    digits, max_exponent = fpp2(man, exp, prec, 10, k)

    return format_decimal("-" if sign else "", digits, max_exponent, scientific)


def format_decimal(sign, digits, exp, scientific=False):
    """Helper to convert raw extracted digits and exponent into a standard string."""
    ndigits = len(digits)

    # Scientific Notation Needed (Very large or small numbers)
    if exp < -4 or exp > sys.float_info.dig or scientific:
        if ndigits > 1:
            return f"{sign}{digits[0]}.{digits[1:]}e{exp:+03}"
        return f"{sign}{digits}e{exp:+03}"

    # Regular Floating Point Notation
    offset = exp + 1
    if exp >= 0:
        if offset < ndigits:
            return f"{sign}{digits[:offset]}.{digits[offset:]}"
        return f"{sign}{digits}{"0" * (offset - ndigits)}.0"
    return f"{sign}0.{'0' * (-offset)}{digits}"


if __name__ == "__main__":
    import random
    import struct

    def check(x):
        mx = mpmath.mpf(x)
        if x and abs(x) < sys.float_info.min:
            return x == float(burger_dybvig(mx))
        if not x and str(x) == '-0.0':
            return True  # no signed zero
        return str(x) == burger_dybvig(mx)

    def random_subnormal():
        mantissa = random.getrandbits(52)
        x = struct.unpack('>d', mantissa.to_bytes(8))[0]
        assert abs(x) < sys.float_info.min
        return x

    # Test corner cases
    for x in [0.1, 0.3, 1.2345, 12345678.9, 0.0000456, 1.1e20, 1e20,
              0.0, float('inf'), float('nan'), float('-inf'),
              -0.25, -10.1, -1e200, -3497605547.3320312, 10.0,
              float.fromhex('1.1p10'), float.fromhex('-1.ap-13'),
              1e15, 1e16, 1e17, 5.5443656173533e-310,
              5.960464477539063e-08]:
        assert check(x)

    # Some random data
    for _ in range(10**5):
        x = random.choice([(random.random() - 0.5)*2,
                           (random.random() - 0.5)*12,
                           (random.random() - 0.5)*10**10,
                           (random.random() - 0.5)/10**7,
                           random_subnormal()])
        assert check(x)

    # All finite float16
    for n in range(1 << 16):
        x = struct.unpack('e', n.to_bytes(2))[0]
        assert check(x)

    for prec in [11, 24, 53, 113, 237]:
        mpmath.mp.prec = prec
        for _ in range(10000):
            f = random.choice([(mpmath.rand()-0.5)*2 for _ in range(10)]
                              + [(mpmath.rand()-0.5)*2*10**5 for _ in range(5)]
                              + [(mpmath.rand()-0.5)*2/10**5 for _ in range(5)]
                              + [(mpmath.rand()-0.5)*2*10**100 for _ in range(2)])
            if not f:
                continue
            s = burger_dybvig(f, prec=prec, scientific=True)
            b = mpmath.mpf(s)
            assert f == b  # round-trip
            if '.' not in s:
                continue
            # test that short repr is really minimal
            frac, *exponent = s.split('e')
            integer, frac = frac.split('.')
            exponent = 'e' + exponent[0]
            assert f == mpmath.mpf(str(integer + '.' + frac + exponent))
            frac = frac[:-random.randint(1, len(frac))]
            assert f != mpmath.mpf(str(integer + '.' + frac + exponent))

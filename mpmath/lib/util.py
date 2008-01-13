import math

# Same as standard Python float
STANDARD_PREC = 53

LOG2_10 = math.log(10,2)  # 3.3219...


def giant_steps(start, target):
    """Generate a list of precisions ranging from 'start' to 'target'
    that doubles with each step. This is used by quadratically
    convergent iterations (that is, Newton iterations), where we want
    to keep the precision at the same level as the accuracy in each
    step to minimize work.

    For example, to find a sequence of precisions to reach 1000 bits
    starting from a 53-bit estimate, giant_steps(53, 1000) gives

        [64, 126, 251, 501, 1000]

    So, if we start Newton's method with a 53-bit accurate initial
    guess, the first iteration should be carried out at 64-bit
    precision, the second at 126-bit precision, and so on.

    Note the conservative rounding (1000 to 501, etc); this is used
    guard against unit errors in the last place."""
    L = [target]
    while L[-1] > start*2:
        L = L + [L[-1]//2 + 1]
    return L[::-1]


def rshift_quick(x, n):
    """For an integer x, calculate x >> n with the fastest (floor)
    rounding. Unlike the plain Python expression (x >> n), n is
    allowed to be negative, in which case a left shift is performed."""
    if n >= 0: return x >> n
    else:      return x << (-n)


def lshift_quick(x, n):
    """For an integer x, calculate x << n. Unlike the plain Python
    expression (x << n), n is allowed to be negative, in which case a
    right shift with default (floor) rounding is performed."""
    if n >= 0: return x << n
    else:      return x >> (-n)


def make_fixed(s, prec):
    """Convert a floating-point number to a fixed-point big integer"""
    man, exp, bc = s
    offset = exp + prec
    if offset >= 0:
        return man << offset
    else:
        return man >> (-offset)



# TODO: speed up for bases 2, 4, 8, 16, ...

def bin_to_radix(x, xbits, base, bdigits):
    """
    Radix conversion for fixed-point numbers. That is, convert
    x * 2**xbits to floor(x * 10**bdigits).
    """
    return x * (base**bdigits) >> xbits

stddigits = '0123456789abcdefghijklmnopqrstuvwxyz'

def small_numeral(n, base=10, digits=stddigits):
    """
    Return the string numeral of a positive integer in an arbitrary
    base. Most efficient for small input.
    """
    if base == 10:
        return str(n)
    digs = []
    while n:
        n, digit = divmod(n, base)
        digs.append(digits[digit])
    return "".join(digs[::-1])

def numeral(n, base=10, size=0, digits=stddigits):
    """
    Represent the integer n as a string of digits in the given base.
    Recursive division is used to make this function about 3x faster
    than Python's str() for converting integers to decimal strings.

    The 'size' parameters specifies the number of digits in n; this
    number is only used to determine splitting points and need not
    be exact.
    """

    if n < 0:
        return "-" + numeral(-n, base, size, digits)

    # Fast enough to do directly
    if size < 250:
        return small_numeral(n, base, digits)

    # Divide in half
    half = (size // 2) + (size & 1)
    A, B = divmod(n, base**half)
    ad = numeral(A, base, half, digits)
    bd = numeral(B, base, half, digits).rjust(half, "0")
    return ad + bd

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

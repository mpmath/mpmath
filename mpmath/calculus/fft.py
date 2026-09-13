from .calculus import defun


def _is_power_of_two(n):
    return (n & (n - 1)) == 0

def _next_power_of_two(n):
    return 1 << (n - 1).bit_length()

def _fft_cooley_tuckey(ctx, values, inverse=False):
    """
    This function implements the Radix-2 Cooley-Tukey FFT algorithm iteratively.
    It computes the Fast Fourier Transform (or Inverse FFT) of a sequence of
    complex numbers.

    https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
    """
    n = len(values)
    if n <= 1:
        return values

    # Bit-Reversal Permutation
    transformed = [ctx.zero] * n
    num_bits = n.bit_length() - 1
    for i in range(n):
        rev = 0
        val = values[i]
        for _ in range(num_bits):
            rev <<= 1
            rev |= (i & 1)
            i >>= 1
        transformed[rev] = val

    sign = ctx.one if inverse else -ctx.one

    length = 2
    while length <= n:
        half = length // 2
        w_len = ctx.expjpi(2 * sign / length)

        for i in range(0, n, length):
            w = ctx.one
            for j in range(half):
                u = transformed[i + j]
                v = transformed[i + j + half] * w

                transformed[i + j] = u + v
                transformed[i + j + half] = u - v

                w *= w_len

        length <<= 1

    return transformed


def _fft_convolve(ctx, values_a, values_b):
    """
    Computes the convolution of two sequences using the Fast Fourier Transform (FFT).
    """
    size = _next_power_of_two(len(values_a) + len(values_b) - 1)
    padded_a = values_a + [ctx.zero] * (size - len(values_a))
    padded_b = values_b + [ctx.zero] * (size - len(values_b))

    spectrum_a = _fft_cooley_tuckey(ctx, padded_a)
    spectrum_b = _fft_cooley_tuckey(ctx, padded_b)
    spectrum_product = [a * b for a, b in zip(spectrum_a, spectrum_b)]
    result = _fft_cooley_tuckey(ctx, spectrum_product, True)
    return [value / size for value in result]

def _fft_bluestein(ctx, values, inverse=False):
    """
    This function implements Bluestein's algorithm for computing the Fast Fourier Transform (or Inverse Fast Fourier Transform) of a sequence of complex numbers of arbitrary length.

    https://en.wikipedia.org/wiki/Chirp_Z-transform
    https://edukatesengkang.com/2026/09/01/how-to-learn-bluesteins-fft-algorithm-chirp-multiplication-convolution-arbitrary-length-dfts-and-prime-size-fourier-transforms/
    """
    n = len(values)

    sign = ctx.one if inverse else -ctx.one
    chirp = [ctx.expjpi(sign * j * j / n) for j in range(n)]

    a = [values[j] * chirp[j] for j in range(n)]
    b = [ctx.zero] * (2 * n - 1)
    center = n - 1
    for j in range(n):
        b[center + j] = ctx.expjpi(-sign * j * j / n)
        b[center - j] = ctx.expjpi(-sign * j * j / n)

    convolution = _fft_convolve(ctx, a, b)
    return [convolution[n - 1 + k] * chirp[k] for k in range(n)]

@defun
def fft(ctx, values):
    r"""
    Computes the Discrete Fourier Transform (DFT) of a sequence.

    Uses the radix-2 Cooley-Tukey algorithm for power-of-two lengths and Bluestein's algorithm for all other lengths.

    **Examples**

    >>> from mpmath import mp
    >>> mp.pretty = True
    >>> mp.fft([1, 0, 0, 0])
    [1.0, (1.0 + 0.0j), 1.0, (1.0 + 0.0j)]
    >>> mp.fft([1 + 2j, 1 + 2j])
    [(2.0 + 4.0j), (0.0 + 0.0j)]
    >>> mp.fft([1, 2, 3, 4])
    [10.0, (-2.0 + 2.0j), -2.0, (-2.0 - 2.0j)]
    >>> [mp.chop(x) for x in mp.fft([1, 2, 1])]
    [4.0, (-0.5 - 0.866025403784439j), (-0.5 + 0.866025403784439j)]
    """
    n = len(values)
    if n == 0:
        return []

    converted_values = [ctx.convert(v) for v in values]
    with ctx.extraprec(10):
        if _is_power_of_two(n):
            result = _fft_cooley_tuckey(ctx, converted_values)
        else:
            result = _fft_bluestein(ctx, converted_values)
    return [+v for v in result]

@defun
def invfft(ctx, values):
    r"""
    Computes the inverse Discrete Fourier Transform (IDFT) of a sequence.

    Uses the radix-2 Cooley-Tukey algorithm for power-of-two lengths and Bluestein's algorithm for all other lengths.

    **Examples**

    >>> from mpmath import mp
    >>> mp.pretty = True
    >>> mp.invfft([1, 1, 1, 1])
    [1.0, (0.0 + 0.0j), 0.0, (0.0 + 0.0j)]
    >>> x = [1, 2, 3, 4]
    >>> mp.invfft(mp.fft(x))
    [(1.0 + 0.0j), (2.0 + 0.0j), (3.0 + 0.0j), (4.0 + 0.0j)]
    >>> mp.invfft(mp.fft([1.0 + 1.0j, 2.0 + 2.0j, 3.0 + 3.0j]))
    [(1.0 + 1.0j), (2.0 + 2.0j), (3.0 + 3.0j)]
    """
    n = len(values)
    if n == 0:
        return []

    converted_values = [ctx.convert(v) for v in values]
    with ctx.extraprec(10):
        if _is_power_of_two(n):
            result = _fft_cooley_tuckey(ctx, converted_values, True)
        else:
            result = _fft_bluestein(ctx, converted_values, True)
    return [val / n for val in result]

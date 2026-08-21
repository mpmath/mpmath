from .calculus import defun

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
        for _ in range(num_bits):
            rev <<= 1
            rev |= (i & 1)
            i >>= 1           
        transformed[rev] = values[i]

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

@defun
def fft(ctx, values):
    r"""
    Computes the Discrete Fourier Transform (DFT) of a sequence.

    Raises NotImplementedError if the input sequence length is not a power of 2.

    **Examples**

    >>> from mpmath import mp
    >>> mp.pretty = True
    >>> mp.fft([1, 0, 0, 0])
    [1.0, (1.0 + 0.0j), 1.0, (1.0 + 0.0j)]
    >>> mp.fft([1 + 2j, 1 + 2j])
    [(2.0 + 4.0j), (0.0 + 0.0j)]
    >>> mp.fft([1, 2, 3, 4])
    [10.0, (-2.0 + 2.0j), -2.0, (-2.0 - 2.0j)]
    """
    n = len(values)
    if n == 0:
        return []

    is_power_of_two = (n & (n - 1)) == 0
    if not is_power_of_two:
        raise NotImplementedError("FFT is only implemented for lengths that "
                                  f"are powers of 2, got length: {n}")

    converted_values = [ctx.convert(v) for v in values]
    with ctx.extraprec(10):
        result = _fft_cooley_tuckey(ctx, converted_values)
    return [+v for v in result]

@defun
def invfft(ctx, values):
    r"""
    Computes the inverse Discrete Fourier Transform (IDFT) of a sequence.

    Raises NotImplementedError if the input sequence length is not a power of 2.

    **Examples**

    >>> from mpmath import mp
    >>> mp.pretty = True
    >>> mp.invfft([1, 1, 1, 1])
    [1.0, (0.0 + 0.0j), 0.0, (0.0 + 0.0j)]
    >>> x = [1, 2, 3, 4]
    >>> mp.invfft(mp.fft(x))
    [(1.0 + 0.0j), (2.0 + 0.0j), (3.0 + 0.0j), (4.0 + 0.0j)]
    """
    n = len(values)
    if n == 0:
        return []

    is_power_of_two = (n & (n - 1)) == 0
    if not is_power_of_two:
        raise NotImplementedError("Inverse FFT is only implemented for lengths that "
                                  f"are powers of 2, got length: {n}")

    converted_values = [ctx.convert(v) for v in values]
    with ctx.extraprec(10):
        result = _fft_cooley_tuckey(ctx, converted_values, True)
    return [val / n for val in result]
